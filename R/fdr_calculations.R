#' @title Fit log-linear model
#' @description Fits a linear model to non-zero probabilities only.
#' @noRd
._log_lm <- function(y, x) {
    # Only keep points where hits > 0
    valid <- y > 0

    # Require at least 3 points to provide a stable slope estimate
    if (sum(valid) < 3) {
        return(NULL)
    }

    y_filt <- y[valid]
    x_filt <- x[valid]

    model <- stats::lm(log(y_filt) ~ x_filt)
    s <- summary(model)

    # Return coefficients, their standard errors, and the residual variance (MSE)
    return(list(
        slope = stats::coef(model)[2],
        intercept = stats::coef(model)[1],
        se_slope = s$coefficients[2, "Std. Error"],
        se_intercept = s$coefficients[1, "Std. Error"],
        mse = s$sigma^2
    ))
}

#' @title Prepare occupancy models from a binding profile
#' @description This internal function generates FDR estimation models for a single sample's
#'   binding profile using a permutation-based approach.
#' @param profile_scores A numeric vector of binding scores for a single sample.
#' @param profile_widths A numeric vector of widths for the genomic fragments.
#' @param iterations Integer. The number of random sampling iterations.
#' @param frag_samples A numeric vector of fragment counts to sample.
#' @param BPPARAM A BiocParallel backend instance.
#' @param seed An optional integer to ensure reproducible random number generation.
#' @noRd
._fit_occupancy_models <- function(profile_scores,
                                           profile_widths,
                                           iterations = 100000,
                                           frag_samples = c(1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 40, 80, 160),
                                           BPPARAM = BiocParallel::bpparam(),
                                           seed = NULL) {

    # Determine thresholds from the profile distribution
    probs <- c(seq(0.4, 0.8, by = 0.05), seq(0.825, 0.975, by = 0.025),0.985)
    all_quantiles <- stats::quantile(profile_scores, probs = probs)

    # Positive thresholds only
    threshold_samples <- sort(unique(all_quantiles[all_quantiles > 0.1]))
    message("   Thresholds: ", paste(sprintf("%0.2f", threshold_samples), collapse = " "))

    num_elements <- length(profile_scores)

    bpp <- BPPARAM
    bpp$RNGseed <- seed
    fdr_sim_results <- BiocParallel::bplapply(
        threshold_samples,
        BPPARAM = bpp,
        FUN = function(threshold) {
            hits <- stats::setNames(rep(0, length(frag_samples)), frag_samples)
            for (i in seq_len(iterations)) {
                for (f in frag_samples) {
                    idx <- sample.int(num_elements, size = f, replace = TRUE)
                    s <- profile_scores[idx]
                    w <- profile_widths[idx]

                    total_w <- sum(w)
                    random_mean_occupancy <- if (total_w > 0) sum(s * w) / total_w else 0

                    if (random_mean_occupancy > threshold) {
                        hits[as.character(f)] <- hits[as.character(f)] + 1
                    }
                }
            }
            return(hits)
        }
    )

    # Calculate empirical probabilities
    fdr_matrix <- do.call(rbind, fdr_sim_results) / iterations

    model_params <- lapply(seq_along(frag_samples), function(i) {
        f <- frag_samples[i]
        fdr_values <- fdr_matrix[, as.character(f)]
        m <- ._log_lm(y = fdr_values, x = threshold_samples)

        if (is.null(m)) return(NULL)

        data.frame(
            fragment_count = f,
            slope = m$slope,
            intercept = m$intercept,
            se_slope = m$se_slope,
            se_int = m$se_intercept,
            mse = m$mse
        )
    })

    model_params_df <- do.call(rbind, model_params)
    model_params_df <- stats::na.omit(model_params_df)

    # Second level of regression using Weighted Least Squares (WLS)
    # Weights are 1 / (Standard Error)

    model_params_df$log_f <- log(model_params_df$fragment_count)

    # Use a small epsilon to prevent Inf weights
    eps <- 1e-10
    w_slope <- 1 / (pmax(model_params_df$se_slope, eps))
    w_int <- 1 / (pmax(model_params_df$se_int, eps))

    # All Tier 2 models are non-linear
    slope_model <- stats::lm(slope ~ splines::ns(log_f, df = 3),
                             data = model_params_df,
                             weights = w_slope)
    intercept_model <- stats::lm(intercept ~ splines::ns(log_f, df = 3),
                                 data = model_params_df,
                                 weights = w_int)
    mse_model <- stats::lm(mse ~ splines::ns(fragment_count, df = 3),
                           data = model_params_df)

    return(list(
        models = list(
            slope_model = slope_model,
            intercept_model = intercept_model,
            mse_model = mse_model
        ),
        params_df = model_params_df,
        fdr_matrix = fdr_matrix,
        threshold_samples = threshold_samples
    ))
}

#' @title Calculate FDR for gene occupancy
#' @description uses pre-computed occupancy models to find the FDR for
#'   occupancy scores and fragment counts.
#' @param fdr_models a list containing `slope_model`, `intercept_model`, and `mse_model`.
#' @param occupancy numeric vector of mean weighted occupancy scores.
#' @param fragment_count integer vector of fragment counts.
#' @return a numeric vector of calculated FDR values, capped at 1.0.
#' @noRd
._estimate_occupancy_pvals <- function(fdr_models, occupancy, fragment_count) {
    # Initialise results vector with the default value for invalid or non-significant cases
    n_genes <- length(occupancy)
    fdr_out <- rep(1.0, n_genes)

    # Identify indices with valid data for prediction
    # Genes with 0 fragments cannot have their log probability calculated
    valid_idx <- !is.na(occupancy) & !is.na(fragment_count) & fragment_count > 0

    if (any(valid_idx)) {
        # Extract subset for calculation
        occ_v <- occupancy[valid_idx]
        frag_v <- fragment_count[valid_idx]

        # Construct prediction data frame once
        new_data <- data.frame(
            fragment_count = frag_v,
            log_f = log(frag_v)
        )

        pred_int <- stats::predict(fdr_models$intercept_model, newdata = new_data)
        pred_slope <- stats::predict(fdr_models$slope_model, newdata = new_data)
        pred_mse <- stats::predict(fdr_models$mse_model, newdata = new_data)

        # Ensure predicted MSE is non-negative to avoid issues with Jensen's correction
        pred_mse <- pmax(pred_mse, 0)

        # Calculate log-probability with Jensen's inequality correction
        log_p <- (pred_slope * occ_v) + pred_int
        res <- exp(log_p + (pred_mse / 2))

        # Clamp values and update the output vector
        res[res > 1.0] <- 1.0
        res[is.na(res)] <- 1.0
        fdr_out[valid_idx] <- res
    }

    return(fdr_out)
}



#' Calculate and add gene occupancy FDR
#'
#' @description
#' This function calculates a p-value for gene occupancy scores.  When adjusted
#' for multiple-hypothesis testing, the resulting FDR may used as a proxy for
#' determining whether a gene is actively expressed in RNA Polymerase TaDa
#' experiments. The function applies a permutation-based null model to each
#' sample's binding profile to determine empirical p-values, and returns these
#' as new columns in the input occupancy dataframe.
#'
#' These p-values are then aggregated and
#'
#' @param binding_data A `GRanges` object of binding profiles, where metadata
#'   columns represent samples. This is typically the `binding_profiles_data`
#'   element from the list returned by `load_data_genes`.
#' @param occupancy_df A data frame of gene occupancies, typically the `occupancy`
#'   element from the list returned by `load_data_genes`. It must contain
#'   sample columns and a `nfrags` column.
#' @param null_model_iterations Integer. The number of sampling iterations to build the
#'   FDR null model. Higher values are more accurate but slower. Default is 100000.
#' @param return_per_replicate_fdr Logical. Returns individual FDR scores by applying BH
#'   adjustment on each individual sample.  Use this option if not intending
#'   to apply downstream condition-level p-value aggregation via
#'   `differential_binding()`.  (Default: FALSE)
#' @param plot_diagnostics Logical.  If TRUE, will plot the Tier 2 regression fits
#'   for the p-value prediction slope, intercept and mean squared error (MSE).
#'   (Default: FALSE)
#' @param BPPARAM A `BiocParallelParam` instance specifying the parallel backend
#'   to use for computation. See `BiocParallel::bpparam()`.
#' @param seed An optional integer. If provided, it is used to set the random
#'   seed before the sampling process begins, ensuring that the FDR calculations
#'   are fully reproducible. If `NULL` (default), results may vary between runs.
#'
#' @return
#' An updated version of the input `occupancy_df` data frame, with new columns
#' appended. For each sample column (e.g., `SampleA`), a corresponding p-value
#' column (e.g., `SampleA_pval`) is added.
#'
#' @details
#' The algorithm used here is substantially revised and adapted from the
#' method described in Southall et al. (2013), Dev Cell, 26(1), 101-12.  The
#' broad principle is as follows:
#' \enumerate{
#'   \item For each sample, a null distribution of mean occupancy scores is
#'     simulated by repeatedly sampling random fragments from the genome-wide
#'     binding profile. This is done for various numbers of fragments.
#'   \item A two-tiered regression model is fitted to the simulation results. This
#'     creates a statistical model that can predict the empirical p-value for any given
#'     occupancy score and fragment count.
#'   \item This predictive model is then applied to the actual observed mean
#'     occupancy and fragment count for each gene in the `occupancy_df`.
#'   \item The final set of p-values is adjusted using BH correction to generate
#'     FDR scores for each gene.
#'   \item The FDR value for each gene in each sample is appended to the
#'     `occupancy_df` in a new column.
#' }
#'
#' Key differences from the original Southall et al. algorithm include fitting
#' natural spline models, using WLS to account for heterscedasticity of models,
#' and sampling with replacement to generate the underlying null distribution.
#'
#' This function is typically not called directly by the user, but is instead
#' invoked via `load_data_genes(calculate_occupancy_pvals = TRUE)`.
#'
#' @seealso [load_data_genes()] which uses this function internally.
#'
#' @examples
#'
#' # Prepare sample binding data (GRanges object)
#' # Here, we'll load one of the sample files included with the package
#' # (this is a TF binding profile, which would not normally be used for
#' # occupancy calculations)
#' data_dir <- system.file("extdata", package = "damidBind")
#' bgraph_file <- file.path(data_dir, "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm.gatc.2L.bedgraph.gz")
#'
#' # The internal function `import_bedgraph_as_df` is used here for convenience.
#' # The column name 'score' will be replaced with the sample name.
#' binding_df <- damidBind:::import_bedgraph_as_df(bgraph_file, colname = "Sample_1")
#'
#' binding_gr <- GenomicRanges::GRanges(
#'     seqnames = binding_df$chr,
#'     ranges = IRanges::IRanges(start = binding_df$start, end = binding_df$end),
#'     Sample_1 = binding_df$Sample_1
#' )
#'
#' # Create a mock occupancy data frame.
#' # In a real analysis, this would be generated by `calculate_occupancy()`.
#' mock_occupancy <- data.frame(
#'     name = c("geneA", "geneB", "geneC"),
#'     nfrags = c(5, 10, 2),
#'     Sample_1 = c(1.5, 0.8, -0.2),
#'     row.names = c("geneA", "geneB", "geneC")
#' )
#'
#' # Calculate pvals.
#' # We use a low null_model_iterations for speed, and a seed for reproducibility.
#' occupancy_with_pvals <- calculate_and_add_occupancy_pvals(
#'     binding_data = binding_gr,
#'     occupancy_df = mock_occupancy,
#'     null_model_iterations = 100,
#'     BPPARAM = BiocParallel::SerialParam(),
#'     seed = 42
#' )
#'
#' # View the result, which now includes a _pval column for Sample_1.
#' print(occupancy_with_pvals)
#'
#'
#' @export
calculate_and_add_occupancy_pvals <- function(binding_data,
                                  occupancy_df,
                                  null_model_iterations = 100000,
                                  return_per_replicate_fdr = FALSE,
                                  plot_diagnostics = FALSE,
                                  BPPARAM = BiocParallel::bpparam(),
                                  seed = NULL) {

    message("\nBuilding occupancy models ...")

    sample_cols <- colnames(mcols(binding_data))

    # Generate occupancy models for each sample profile
    fdr_models_list <- list()
    for (sample in sample_cols) {
        message(sprintf(" - %s", sample))
        profile_scores <- mcols(binding_data)[[sample]]
        profile_widths <- GenomicRanges::width(binding_data)

        fdr_results <- ._fit_occupancy_models(
            profile_scores = profile_scores,
            profile_widths = profile_widths,
            iterations = null_model_iterations,
            BPPARAM = BPPARAM,
            seed = seed
        )

        fdr_models_list[[sample]] <- fdr_results
    }
    message("Occupancy models prepared for all samples.\n")

    # Calculate p-values for each gene in each sample and add to occupancy table
    message("Calculating occupancy p-values ...")
    for (sample in sample_cols) {
        # Check sample exists in the occupancy table
        if (!sample %in% colnames(occupancy_df)) {
            warning(sprintf(
                "Sample '%s' from binding profiles not found in occupancy table. Skipping occupancy p-value calculation.",
                sample
            ))
            next
        } else {
            message(sprintf(" - %s", sample))
        }

        fdr_models <- fdr_models_list[[sample]]$models

        # Get occupancy and fragment counts for all genes for this sample
        occupancy_values <- occupancy_df[[sample]]
        fragment_counts <- occupancy_df$nfrags

        # Apply the p-value calculation
        p_values <- ._estimate_occupancy_pvals(
            fdr_models = fdr_models,
            occupancy = occupancy_values,
            fragment_count = fragment_counts
        )

        if (isTRUE(return_per_replicate_fdr)) {
            # Return per-replicate BH-adjusted FDR
            occupancy_df[[paste0(sample, "_FDR")]] <- stats::p.adjust(p_values, method = "BH")
        } else {
            # Return unadjusted p-values
            occupancy_df[[paste0(sample, "_pval")]] <- p_values
        }
    }
    message("FDR calculation complete.\n")

    # Call diagnostic plots if requested
    if (isTRUE(plot_diagnostics)) {
        message("Generating occupancy diagnostic plots")
        for (sample in sample_cols) {
            res <- fdr_models_list[[sample]]
            plot_occupancy_diagnostics(
                model_params_df = res$params_df,
                fdr_models = res$models,
                fdr_matrix = res$fdr_matrix,
                threshold_samples = res$threshold_samples,
                sample_name = sample
            )
        }
    }
    return(occupancy_df)
}


#' @title Filter genes by FDR within a specific condition
#' @description
#' Filters a list of genes to retain only those that meet a specified
#' False Discovery Rate (FDR) threshold. If the input is a `DamIDResults` object
#' and a combined condition-level FDR has been calculated, that value is used.
#' Otherwise, the function falls back to filtering against individual replicates.
#'
#' @param data A `DamIDResults` object or the `list` returned by `load_data_genes()`.
#' @param fdr A numeric value between 0 and 1 specifying the FDR cutoff. (Default: 0.05)
#' @param condition A character string identifying the experimental condition.
#'   This string should uniquely match the relevant sample columns (e.g., "L4" will match
#'   "L4_rep1_FDR" and "L4_rep2_FDR"). If `data` is a `DamIDResults` object, this can
#'   be either the internal identifier or the display name for the condition.
#' @param which A character string, either `"any"` or `"all"`. Only applicable
#'   when falling back to individual replicate scores. (Default: `"any"`)
#'   \itemize{
#'     \item If `"any"`, a gene is kept if it meets the `fdr` threshold in at least one
#'       replicate of the specified `condition`.
#'     \item If `"all"`, a gene is kept only if it meets the `fdr` threshold in all
#'       replicates of the specified `condition`.
#'   }
#'
#' @return A `data.frame` containing the `gene_name`, `gene_id`, `avg_occ`, and the
#'   most significant FDR value found.
#'
#' @details
#' This function is primarily used in workflows involving RNA Polymerase TaDa data, where
#' an FDR is calculated for gene occupancy to determine if a gene is actively transcribed.
#' It allows users to identify genes in a single condition that can be considered to be
#' expressed (i.e. RNA Pol occupancy is significantly greater than background).
#'
#' Note that while this is an effective proxy for gene expression, there are edge cases
#' (e.g. paused polymerase, short genes directly adjacent to an expressed gene TSS or TES)
#' where a gene may have significant occupancy but not, in fact, be transcribed.
#'
#' The function locates the relevant FDR columns in the `occupancy` table by searching for
#' column names that end with `_FDR` and also contain the `condition` string.
#'
#' @export
#' @examples
#' # Create a mock data object with an occupancy table containing FDR values,
#' # similar to the output of `load_data_genes(calculate_occupancy_pvals = TRUE)`.
#' .generate_fdr_example_results <- function() {
#'     occupancy_df <- data.frame(
#'         gene_name = c("geneA", "geneB", "geneC"),
#'         gene_id = c("FBgn01", "FBgn02", "FBgn03"),
#'         L4_rep1 = c(1.5, 0.2, 0.8),
#'         L4_rep2 = c(1.7, 0.9, 0.1),
#'         L5_rep1 = c(0.1, 0.1, 2.0),
#'         L4_rep1_FDR = c(0.01, 0.10, 0.04),
#'         L4_rep2_FDR = c(0.03, 0.02, 0.50),
#'         L5_rep1_FDR = c(0.80, 0.90, 0.01),
#'         row.names = c("geneA", "geneB", "geneC")
#'     )
#'     diff_results_base <- list(
#'         occupancy = occupancy_df,
#'         test_category = "expressed",
#'         matched_samples = list("L4" = c("L4_rep1", "L4_rep2"), "L5" = "L5_rep1")
#'     )
#'     new("DamIDResults",
#'         analysis = data.frame(row.names = rownames(occupancy_df)),
#'         upCond1 = data.frame(),
#'         upCond2 = data.frame(),
#'         cond = c("L4 mock" = "L4", "L5 mock" = "L5"),
#'         data = diff_results_base
#'     )
#' }
#' mock_data <- .generate_fdr_example_results()
#'
#' # Example 1: Get genes with FDR <= 0.05 in ANY L4 replicate.
#' # geneA (0.01, 0.03), geneB (0.02), and geneC (0.04) pass.
#' expressed_in_L4_any <- filter_genes_by_fdr(
#'     mock_data,
#'     fdr = 0.05,
#'     condition = "L4",
#'     which = "any"
#' )
#' print(expressed_in_L4_any)
#'
#' # Example 2: Get genes with FDR <= 0.05 in ALL L4 replicates.
#' # Only geneA (0.01, 0.03) passes.
#' expressed_in_L4_all <- filter_genes_by_fdr(
#'     mock_data,
#'     fdr = 0.05,
#'     condition = "L4",
#'     which = "all"
#' )
#' print(expressed_in_L4_all)
#'
#' # Example 3: Get genes with FDR <= 0.05 in any L5 replicate.
#' # geneC (0.01) and geneD (0.02) pass.
#' expressed_in_L5 <- filter_genes_by_fdr(
#'     mock_data,
#'     fdr = 0.05,
#'     condition = "L5",
#'     which = "any"
#' )
#' print(expressed_in_L5)
filter_genes_by_fdr <- function(data, fdr = 0.05, condition, which = "any") {
    # Input validation
    if (!is(data, "DamIDResults") && !(is.list(data) && "occupancy" %in% names(data))) {
        stop("'data' must be a DamIDResults object or a list from load_data_genes().")
    }
    if (!is.numeric(fdr) || fdr < 0 || fdr > 1) {
        stop("'fdr' must be a numeric value between 0 and 1.")
    }
    if (!is.character(condition) || length(condition) != 1) {
        stop("'condition' must be a single character string.")
    }
    which <- match.arg(which, c("any", "all"))

    # Initialise variables
    occupancy_df <- NULL
    analysis_df <- NULL
    internal_id <- condition

    # Resolve condition names and extract relevant tables
    if (is(data, "DamIDResults")) {
        occupancy_df <- inputData(data)$occupancy
        analysis_df <- analysisTable(data)
        cond_map <- conditionNames(data)
        match_patt <- inputData(data)$match_patterns
        if (condition %in% names(cond_map)) {
            # User provided display name
            internal_id <- unname(cond_map[condition])
        } else if (condition %in% cond_map) {
            # User provided internal name
            internal_id <- condition
        } else if (!is.null(match_patt) && condition %in% match_patt) {
            # User provided original match pattern
            idx <- which(match_patt == condition)[1]
            internal_id <- cond_map[idx]
        } else {
            # Fallback to the provided string, though likely to fail column lookup
            internal_id <- condition
        }
    } else {
        occupancy_df <- data$occupancy
        internal_id <- condition
    }

    # Ensure required annotation columns exist
    if (!all(c("gene_name", "gene_id") %in% colnames(occupancy_df))) {
        stop("The occupancy data frame must contain 'gene_name' and 'gene_id' columns.")
    }

    # Default functionality from v0.99.12: attempt to use combined condition-level FDR
    # These aggregate p-values are calculated in `differential_binding()`
    # and should be present if feeding in a DamIDResults object generated
    # post-v0.99.10.
    if (!is.null(analysis_df)) {
        combined_col <- paste0(internal_id, "_FDR")
        if (combined_col %in% colnames(analysis_df)) {
            message(sprintf("Found combined FDR for condition '%s'. Applying filter.", condition))
            idx <- which(analysis_df[[combined_col]] <= fdr)

            if (length(idx) == 0) {
                return(.empty_fdr_filter_df(fdr, condition, "combined"))
            }

            results <- analysis_df[idx, , drop = FALSE]

            # Formulate the return table consistent with package standards
            # We use the mean columns if they exist, otherwise fallback to occupancy_df
            mean_col <- paste0(internal_id, "_mean")

            # Resolve replicates via metadata for avg_occ calculation
            avg_occ <- if (mean_col %in% colnames(results)) {
                results[[mean_col]]
            } else {
                reps <- inputData(data)$matched_samples[[internal_id]]
                rowMeans(occupancy_df[rownames(results), reps, drop = FALSE], na.rm = TRUE)
            }

            result_df <- data.frame(
                gene_name = results$gene_name,
                gene_id = results$gene_id,
                avg_occ = avg_occ,
                fdr_val = results[[combined_col]],
                row.names = rownames(results)
            )
            colnames(result_df)[4] <- "combined_fdr"

            message(sprintf("Found %d genes using combined condition FDR.", nrow(result_df)))
            return(result_df)
        }
    }

    # Legacy fallback to individual replicate FDRs in the occupancy table
    # This is used if the combined FDR is missing (either data from versions v0.99.10 or older)
    # or when the user has explicity requested per-replicate FDRs on loading
    # If replicate pvals are found, and no condition aggregate is found, then currently these
    # are BH-adjusted individually and the older underpowered any/all logic is applied
    # TODO: allow on-the-fly condition aggregation of p-values?

    if (is.null(occupancy_df)) {
        stop("Occupancy data not found. Cannot perform filtering.")
    }

    # Retrieve matched samples from metadata if available
    relevant_samples <- NULL
    if (is(data, "DamIDResults")) {
        m_samples <- inputData(data)$matched_samples
        if (!is.null(m_samples) && internal_id %in% names(m_samples)) {
            relevant_samples <- m_samples[[internal_id]]
        }
    }

    # If metadata is missing (common in old objects or manual mocks), fallback to grep
    if (is.null(relevant_samples)) {
        relevant_samples <- .grep_replicate_cols(occupancy_df, internal_id)
    }

    # Look for per-replicate _FDR values or _pvals; prioritise FDRs
    fdr_cols <- paste0(relevant_samples, "_FDR")
    pval_cols <- paste0(relevant_samples, "_pval")

    if (length(fdr_cols) > 0 && all(fdr_cols %in% colnames(occupancy_df))) {
        sig_subset <- occupancy_df[, fdr_cols, drop = FALSE]
    } else if (length(pval_cols) > 0 && all(pval_cols %in% colnames(occupancy_df))) {
        message("Adjusting p-values to FDR for filtering...")
        sig_subset <- apply(occupancy_df[, pval_cols, drop = FALSE], 2, p.adjust, method = "BH")
    } else {
        warning(sprintf("No '_FDR' columns found in the data matching the condition '%s'.", condition))
        return(.empty_fdr_filter_df(fdr, condition, which))
    }

    # Apply row-wise comparison logic
    keep_indices <- if (which == "any") {
        apply(sig_subset, 1, function(r) any(r <= fdr, na.rm = TRUE))
    } else {
        apply(sig_subset, 1, function(r) !any(is.na(r)) && all(r <= fdr, na.rm = TRUE))
    }

    if (!any(keep_indices)) {
        return(.empty_fdr_filter_df(fdr, condition, which))
    }

    # Extract results and calculate summary metrics
    results <- occupancy_df[keep_indices, , drop = FALSE]
    avg_occ <- rowMeans(results[, relevant_samples, drop = FALSE], na.rm = TRUE)
    fdr_val <- apply(sig_subset[keep_indices, , drop = FALSE], 1, min, na.rm = TRUE)

    # avg_occ <- rowMeans(results[, replicate_cols, drop = FALSE], na.rm = TRUE)
    # min_fdr <- apply(results[, relevant_fdr_cols, drop = FALSE], 1, min, na.rm = TRUE)

    result_df <- data.frame(
        gene_name = results$gene_name,
        gene_id = results$gene_id,
        avg_occ = avg_occ,
        fdr_val = fdr_val,
        row.names = rownames(results)
    )

    message(sprintf("%d genes passed the FDR filter using rule '%s'.", nrow(result_df), which))
    return(result_df)
}


#' Internal helper to generate empty result dataframes
#' @noRd
.empty_fdr_filter_df <- function(fdr, cond, which) {
    message(sprintf("No genes passed the FDR <= %s filter for '%s' (%s).", fdr, cond, which))
    return(data.frame(
        gene_name = character(0),
        gene_id = character(0),
        avg_occ = numeric(0),
        fdr_val = numeric(0)
    ))
}


#' Internal helper to find sample replicate columns
#' @noRd
.grep_replicate_cols <- function(df, pattern) {
    all_cols <- colnames(df)
    exclude <- c("chr", "start", "end", "name", "gene_name", "gene_id", "nfrags")
    candidates <- setdiff(all_cols, exclude)

    # First, try to find base sample columns (excluding suffixes)
    base_samples <- candidates[!grepl("(_FDR|_pval)$", candidates)]
    hits <- base_samples[grepl(pattern, base_samples)]

    # If no base samples found (e.g. in a summary-only dataframe),
    # look for the FDR columns themselves and strip the suffix
    if (length(hits) == 0) {
        fdr_hits <- candidates[grepl(paste0(pattern, ".*_FDR$"), candidates)]
        hits <- sub("_FDR$", "", fdr_hits)
    }

    return(hits)
}

# .grep_replicate_cols <- function(df, pattern) {
#     all_cols <- colnames(df)
#     # Exclude base metadata and FDR columns
#     exclude <- c("chr", "start", "end", "name", "gene_name", "gene_id", "nfrags")
#     candidates <- setdiff(all_cols, exclude)
#     candidates <- candidates[!grepl("(_FDR|_pval)$", candidates)]
#
#     # Use fixed=FALSE to allow the user to provide a regex pattern
#     # if they are calling this on a raw list
#     hits <- candidates[grepl(pattern, candidates)]
#
#     if (length(hits) == 0) {
#         # Try fixed matching as a fallback
#         hits <- candidates[grepl(pattern, candidates, fixed = TRUE)]
#     }
#     return(hits)
# }


#' @title Filter Loci Based on Gene Expression FDR
#' @description Internal helper to filter a list of locus identifiers based on an FDR
#'   threshold applied to `_FDR` columns from the two conditions being compared.
#' @param diff_results A `DamIDResults` object from which to get occupancy data
#'   and condition information.
#' @param fdr_filter_threshold Numeric threshold for FDR filtering.
#' @param row_names_to_filter (Optional) A character vector of locus identifiers (row names)
#'   to be filtered. (Default: NULL, will use all rownames from table)
#' @return A character vector containing the subset of `row_names_to_filter`
#'   that pass the FDR threshold. If filtering is not possible (e.g., no FDR
#'   columns found), it returns the original `row_names_to_filter` vector
#'   along with a warning.
#' @noRd
filter_on_fdr <- function(diff_results, fdr_filter_threshold, row_names_to_filter=NULL) {
    analysis_tab <- analysisTable(diff_results)
    # Get internal IDs (the values in the named cond vector)
    cond_internal_ids <- as.character(conditionNames(diff_results))

    # Search for combined columns first: e.g. "NSCs_FDR"
    combined_cols <- paste0(cond_internal_ids, "_FDR")
    available_combined <- intersect(combined_cols, colnames(analysis_tab))

    if (length(available_combined) > 0) {
        # This is now the default and preferred approach as of v0.99.12
        message("Using condition FDRs from combined p-values ...")
        if (is.null(row_names_to_filter)) row_names_to_filter <- rownames(analysis_tab)

        fdr_subset <- analysis_tab[row_names_to_filter, available_combined, drop = FALSE]

        # Keep if the combined FDR for EITHER condition passes the threshold
        keep_indices <- apply(fdr_subset, 1, function(row) {
            any(row <= fdr_filter_threshold, na.rm = TRUE)
        })

        filtered_row_names <- row_names_to_filter[keep_indices]
        message(sprintf("Filtering universe: %d of %d total loci pass combined condition FDR <= %.3f.",
                        length(filtered_row_names), length(row_names_to_filter), fdr_filter_threshold))
        return(filtered_row_names)
    }

    # Legacy fallback to individual replicate logic in occupancy_df if combined columns missing
    input_data <- inputData(diff_results)
    occupancy_df <- input_data$occupancy

    if (is.null(occupancy_df)) {
        warning("fdr_filter_threshold was set, but no 'occupancy' data found. Returning all loci unfiltered.")
        return(row_names_to_filter)
    }

    all_fdr_cols <- grep("_FDR$", colnames(occupancy_df), value = TRUE)
    all_pval_cols <- grep("_pval$", colnames(occupancy_df), value = TRUE) ## NOT USED currently. TODO: apply BH correction within function?

    if (length(all_fdr_cols) == 0) {
        warning("fdr_filter_threshold was set, but no '_FDR' columns were found. Returning all loci unfiltered.")
        return(row_names_to_filter)
    }

    if (is.null(row_names_to_filter)) {
        row_names_to_filter <- rownames(analysisTable(diff_results))
    }

    relevant_fdr_cols <- character(0)
    if (!is.null(input_data$matched_samples) && is.list(input_data$matched_samples)) {
        message("Using pre-matched samples for FDR filtering.")
        matched_samples <- unlist(input_data$matched_samples)
        fdr_col_names <- paste0(matched_samples, "_FDR")
        relevant_fdr_cols <- intersect(fdr_col_names, all_fdr_cols)
    } else {
        # Fallback to old regex-based method
        warning("Falling back to regex-based FDR column matching. This may be unreliable with complex condition patterns.")
        cond_identifiers <- as.character(conditionNames(diff_results))
        pattern <- paste(cond_identifiers, collapse = "|")
        relevant_fdr_cols <- grep(pattern, all_fdr_cols, value = TRUE)
    }

    if (length(relevant_fdr_cols) == 0) {
        cond_display <- names(conditionNames(diff_results))
        warning(sprintf(
            "fdr_filter_threshold was set, but no '_FDR' columns matching conditions '%s' or '%s' were found. Returning all loci unfiltered.",
            cond_display[1], cond_display[2]
        ))
        return(row_names_to_filter)
    }

    # Align occupancy data using the provided row names
    aligned_occ_df <- occupancy_df[row_names_to_filter, , drop = FALSE]

    if (anyNA(rownames(aligned_occ_df))) {
        stop("A mismatch in row names between the data to be filtered and occupancy data was detected. Cannot perform FDR filtering.")
    }

    fdr_subset <- aligned_occ_df[, relevant_fdr_cols, drop = FALSE]

    # For each locus, check if any relevant FDR value is <= threshold
    keep_indices <- apply(fdr_subset, 1, function(row) {
        any(row <= fdr_filter_threshold, na.rm = TRUE)
    })

    filtered_row_names <- row_names_to_filter[keep_indices]

    message(sprintf(
        "Filtering universe: %d of %d total loci pass gene expression FDR <= %.3f.",
        length(filtered_row_names), length(row_names_to_filter), fdr_filter_threshold
    ))

    return(filtered_row_names)
}


#' Combine p-values using Fisher's method
#' @param p_matrix A numeric matrix of p-values (genes in rows, replicates in columns)
#' @return A numeric vector of combined p-values
#' @noRd
._combine_p_fisher <- function(p_matrix) {
    if (is.null(p_matrix) || ncol(p_matrix) == 0) return(NULL)
    if (ncol(p_matrix) == 1) return(as.numeric(p_matrix))

    message(" - Combining replicate p-values using Fisher's method")
    # Calculate Fisher statistic: -2 * sum(log(p))
    # Using pmax to avoid log(0)
    fisher_stat <- -2 * rowSums(log(pmax(p_matrix, 1e-16)), na.rm = TRUE)
    df <- 2 * rowSums(!is.na(p_matrix))

    # Combined p-value from chi-squared distribution
    p_combined <- stats::pchisq(fisher_stat, df = df, lower.tail = FALSE)
    return(p_combined)
}

#' Combine p-values using Stouffer's method
#' @param p_matrix A numeric matrix of p-values (genes in rows, replicates in columns)
#' @return A numeric vector of combined p-values
#' @noRd
._combine_p_stouffer <- function(p_matrix) {
    if (is.null(p_matrix) || ncol(p_matrix) == 0) return(NULL)
    if (ncol(p_matrix) == 1) return(as.numeric(p_matrix))

    message(" - Combining replicate p-values using Stouffer's method")
    # Ensuring p is clamped for stability
    p_matrix <- pmax(pmin(p_matrix, 1 - 1e-16), 1e-16)

    # Convert p-values to z-scores (inverse normal distribution)
    z_scores <- stats::qnorm(p_matrix, lower.tail = FALSE)

    # Sum z-scores and divide by sqrt(k)
    k <- rowSums(!is.na(z_scores))
    z_combined <- rowSums(z_scores, na.rm = TRUE) / sqrt(k)

    # Convert back to p-values
    p_combined <- stats::pnorm(z_combined, lower.tail = FALSE)
    return(p_combined)
}

