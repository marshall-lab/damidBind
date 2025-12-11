#' @title Linear Model on Log-Transformed Data
#' @description An internal helper to fit a linear model on log-transformed data.
#'   Ignores infinite values that result from log(0).
#' @param y Numeric vector. The response variable (e.g., FDR values).
#' @param x Numeric vector. The predictor variable (e.g., occupancy thresholds).
#' @return An `lm` object representing the fitted model. Returns a model with NA
#'   coefficients if fewer than two finite data points are available.
#' @noRd
._inf_log_lm <- function(y, x) {
    is_finite <- is.finite(log(y))
    if (sum(is_finite) < 2) {
        # Cannot fit a model with fewer than two data points.
        return(stats::lm(NA ~ NA))
    }
    model <- stats::lm(log(y)[is_finite] ~ x[is_finite])
    return(model)
}

#' @title Prepare FDR Models from a Binding Profile
#' @description This internal function generates FDR estimation models for a single sample's
#'   binding profile using a permutation-based approach.
#' @param profile_scores A numeric vector of binding scores for a single sample.
#' @param iterations Integer. The number of random sampling iterations.
#' @param frag_samples A numeric vector of fragment counts to sample.
#' @param threshold_samples A numeric vector of occupancy thresholds to test.
#' @param BPPARAM A BiocParallel backend instance.
#' @param seed An optional integer to ensure reproducible random number generation
#'   within the parallel jobs.
#' @return A list containing two `lm` objects: `slope_model` and `intercept_model`.
#' @details Major changes from polii.gene.call -- uses sampling with replacement.
#' @noRd
._prepare_fdr_models <- function(profile_scores,
                                 iterations = 50000,
                                 frag_samples = c(1, 2, 3, 4, 6, 8, 10, 12, 15),
                                 threshold_samples = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.5, 2.0),
                                 BPPARAM = BiocParallel::bpparam(),
                                 seed = NULL) {
    bpp <- BPPARAM
    bpp$RNGseed <- seed
    fdr_sim_results <- BiocParallel::bplapply(
        threshold_samples,
        BPPARAM = bpp,
        FUN = function(threshold) {
            # Initialise a named vector to store hit counts for this threshold
            hits <- stats::setNames(rep(0, length(frag_samples)), frag_samples)
            for (i in seq_len(iterations)) {
                # For each fragment sample size, take a random sample and check against threshold
                for (f in frag_samples) {
                    # Now sampling with replacement
                    random_mean_occupancy <- mean(sample(profile_scores, size = f, replace = TRUE))
                    if (random_mean_occupancy > threshold) {
                        hits[as.character(f)] <- hits[as.character(f)] + 1
                    }
                }
            }
            return(hits)
        }
    )

    # Calculate FDR for each combination of fragment count and threshold
    fdr_matrix <- do.call(rbind, fdr_sim_results) / iterations

    # First level of regression: log(FDR) ~ occupancy_threshold
    model_params <- lapply(seq_along(frag_samples), function(i) {
        f <- frag_samples[i]
        fdr_values <- fdr_matrix[, as.character(f)]
        model <- ._inf_log_lm(y = fdr_values, x = threshold_samples)

        if (any(is.na(stats::coef(model)))) {
            return(data.frame(fragment_count = f, slope = NA, intercept = NA))
        }
        data.frame(
            fragment_count = f,
            slope = stats::coef(model)[2],
            intercept = stats::coef(model)[1]
        )
    })
    model_params_df <- do.call(rbind, model_params)
    model_params_df <- stats::na.omit(model_params_df)

    # Second level of regression: model_parameters ~ fragment_count
    slope_model <- stats::lm(slope ~ fragment_count, data = model_params_df)
    intercept_model <- stats::lm(intercept ~ fragment_count, data = model_params_df)

    return(list(
        slope_model = slope_model,
        intercept_model = intercept_model
    ))
}


#' @title Calculate FDR for specific gene occupancy
#' @description Uses pre-computed FDR models to find the FDR for a given
#' mean occupancy and fragment count.
#' @param fdr_models A list containing `slope_model` and `intercept_model`.
#' @param occupancy Numeric. The mean weighted occupancy of the gene.
#' @param fragment_count Integer. The number of fragments covering the gene.
#' @return A single numeric value for the calculated FDR, capped at 1.0.
#' @noRd
._call_fdr_for_gene <- function(fdr_models, occupancy, fragment_count) {
    if (is.na(occupancy) || is.na(fragment_count) || fragment_count == 0) {
        return(1.0)
    }

    new_data <- data.frame(fragment_count = fragment_count)

    predicted_intercept <- stats::predict(fdr_models$intercept_model, newdata = new_data)
    predicted_slope <- stats::predict(fdr_models$slope_model, newdata = new_data)

    fdr_out <- exp(predicted_slope * occupancy + predicted_intercept)

    return(min(fdr_out, 1.0, na.rm = TRUE))
}


#' Calculate and add gene occupancy FDR
#'
#' @description
#' This function calculates a False Discovery Rate (FDR) for gene occupancy scores,
#' which is often used as a proxy for determining whether a gene is actively
#' expressed in RNA Polymerase TaDa experiments. It applies a permutation-based
#' null model to each sample's binding profile to determine empirical p-values,
#' applies the Benjamini-Hochberg (BH) adjustment, and adds the resulting FDR values
#' as new columns to the occupancy data frame.
#'
#' @param binding_data A `GRanges` object of binding profiles, where metadata
#'   columns represent samples. This is typically the `binding_profiles_data`
#'   element from the list returned by `load_data_genes`.
#' @param occupancy_df A data frame of gene occupancies, typically the `occupancy`
#'   element from the list returned by `load_data_genes`. It must contain
#'   sample columns and a `nfrags` column.
#' @param fdr_iterations Integer. The number of sampling iterations to build the
#'   FDR null model. Higher values are more accurate but slower. Default is 50000.
#'   For quick tests, a much lower value (e.g., 100) can be used.
#' @param use_legacy_pvals Logical. The original Southall et al (2013) paper reported
#'   the unadjusted p-values as 'FDR' scores. This implementation applies BH correction
#'   to return the more useful and expected corrected FDR scores for all genes.
#'   (Default: FALSE)
#' @param BPPARAM A `BiocParallelParam` instance specifying the parallel backend
#'   to use for computation. See `BiocParallel::bpparam()`. For full
#'   reproducibility, especially in examples or tests, use
#'   `BiocParallel::SerialParam()`.
#' @param seed An optional integer. If provided, it is used to set the random
#'   seed before the sampling process begins, ensuring that the FDR calculations
#'   are fully reproducible. If `NULL` (default), results will vary between runs.
#'
#' @return
#' An updated version of the input `occupancy_df` data frame, with new columns
#' appended. For each sample column (e.g., `SampleA`), a corresponding FDR
#' column (e.g., `SampleA_FDR`) is added.
#'
#' @details
#' The FDR calculation algorithm is refined and adapted from the method described in
#' Southall et al. (2013), Dev Cell, 26(1), 101-12, and re-implemented in the
#' `polii.gene.call` tool. It operates in several stages:
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
#' linear models and sampling with replacement to generate the underlying null
#' distribution.
#'
#' This function is typically not called directly by the user, but is instead
#' invoked via `load_data_genes(calculate_fdr = TRUE)`. Exporting it allows for
#' advanced use cases where FDR needs to be calculated on a pre-existing
#' occupancy table.
#'
#' @seealso [load_data_genes()] which uses this function internally.
#'
#' @examples
#' if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'     # Prepare sample binding data (GRanges object)
#'     # Here, we'll load one of the sample files included with the package
#'     # (this is a TF binding profile, which would not normally be used for
#'     # occupancy calculations)
#'     data_dir <- system.file("extdata", package = "damidBind")
#'     bgraph_file <- file.path(data_dir, "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm.gatc.2L.bedgraph.gz")
#'
#'     # The internal function `import_bedgraph_as_df` is used here for convenience.
#'     # The column name 'score' will be replaced with the sample name.
#'     binding_df <- damidBind:::import_bedgraph_as_df(bgraph_file, colname = "Sample_1")
#'
#'     binding_gr <- GenomicRanges::GRanges(
#'         seqnames = binding_df$chr,
#'         ranges = IRanges::IRanges(start = binding_df$start, end = binding_df$end),
#'         Sample_1 = binding_df$Sample_1
#'     )
#'
#'     # Create a mock occupancy data frame.
#'     # In a real analysis, this would be generated by `calculate_occupancy()`.
#'     mock_occupancy <- data.frame(
#'         name = c("geneA", "geneB", "geneC"),
#'         nfrags = c(5, 10, 2),
#'         Sample_1 = c(1.5, 0.8, -0.2),
#'         row.names = c("geneA", "geneB", "geneC")
#'     )
#'
#'     # Calculate FDR.
#'     # We use a low fdr_iterations for speed, and a seed for reproducibility.
#'     # BiocParallel::SerialParam() is used to ensure deterministic, single-core execution.
#'     occupancy_with_fdr <- calculate_and_add_fdr(
#'         binding_data = binding_gr,
#'         occupancy_df = mock_occupancy,
#'         fdr_iterations = 100,
#'         BPPARAM = BiocParallel::SerialParam(),
#'         seed = 42
#'     )
#'
#'     # View the result, which now includes an FDR column for Sample_1.
#'     print(occupancy_with_fdr)
#' }
#'
#' @export
calculate_and_add_fdr <- function(binding_data,
                                  occupancy_df,
                                  fdr_iterations = 50000,
                                  use_legacy_pvals = FALSE,
                                  BPPARAM = BiocParallel::bpparam(),
                                  seed = NULL) {

    message("Building FDR models ...")

    sample_cols <- colnames(mcols(binding_data))

    # Generate FDR models for each sample profile
    fdr_models_list <- list()
    for (sample in sample_cols) {
        message(sprintf(" - %s", sample))
        profile_scores <- mcols(binding_data)[[sample]]
        fdr_models_list[[sample]] <- ._prepare_fdr_models(
            profile_scores = profile_scores,
            iterations = fdr_iterations,
            BPPARAM = BPPARAM,
            seed = seed
        )
    }
    message("FDR models prepared for all samples.\n")

    # Calculate FDR for each gene in each sample and add to occupancy table
    message("Calculating FDR ...")
    for (sample in sample_cols) {
        # Check sample exists in the occupancy table
        if (!sample %in% colnames(occupancy_df)) {
            warning(sprintf(
                "Sample '%s' from binding profiles not found in occupancy table. Skipping FDR calculation.",
                sample
            ))
            next
        } else {
            message(sprintf(" - %s", sample))
        }

        fdr_models <- fdr_models_list[[sample]]

        # Get occupancy and fragment counts for all genes for this sample
        occupancy_values <- occupancy_df[[sample]]
        fragment_counts <- occupancy_df$nfrags

        # Use mapply to apply the FDR calculation row-wise (gene-wise)
        fdr_column_values <- mapply(
            ._call_fdr_for_gene,
            occupancy = occupancy_values,
            fragment_count = fragment_counts,
            MoreArgs = list(fdr_models = fdr_models)
        )

        # Add the new FDR column to the occupancy data frame
        fdr_col_name <- paste0(sample, "_FDR")
        true_fdr_column_values <- stats::p.adjust(fdr_column_values, method = "BH")
        final_vals <- if (isFALSE(use_legacy_pvals))
            true_fdr_column_values else fdr_column_values
        occupancy_df[[fdr_col_name]] <- final_vals
    }
    message("FDR calculation complete.")

    return(occupancy_df)
}


#' @title Filter genes by FDR within a specific condition
#' @description
#' Filters a list of genes to retain only those that meet a specified
#' False Discovery Rate (FDR) threshold. The filtering is applied to one or more
#' replicates within a single experimental condition, based on the `_FDR` columns
#' which must be present in the data (e.g., as generated by
#' `load_data_genes(calculate_fdr = TRUE)`).
#'
#' @param data A `DamIDResults` object or the `list` returned by `load_data_genes()`.
#'   The data object must contain an `occupancy` table with FDR columns (e.g., `SampleA_FDR`).
#' @param fdr A numeric value between 0 and 1 specifying the FDR cutoff. Loci with an
#'   FDR less than or equal to this value will be considered. (Default: 0.05)
#' @param condition A character string that identifies the experimental condition to filter on.
#'   This string should uniquely match the relevant sample columns (e.g., "L4" will match
#'   "L4_rep1_FDR" and "L4_rep2_FDR"). If `data` is a `DamIDResults` object, this can
#'   be either the internal identifier or the display name for the condition.
#' @param which A character string, either `"any"` (the default) or `"all"`.
#'   \itemize{
#'     \item If `"any"`, a gene is kept if it meets the `fdr` threshold in at least one
#'       replicate of the specified `condition`.
#'     \item If `"all"`, a gene is kept only if it meets the `fdr` threshold in all
#'       replicates of the specified `condition`.
#'   }
#'
#' @return A `data.frame` containing the `gene_name` and `gene_id` of the genes
#'   that passed the filter. If no genes pass, an empty `data.frame` is returned.
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
#' # similar to the output of `load_data_genes(calculate_fdr = TRUE)`.
#' .create_mock_fdr_data <- function() {
#'     occupancy_df <- data.frame(
#'         gene_name = c("geneA", "geneB", "geneC", "geneD"),
#'         gene_id = c("FBgn01", "FBgn02", "FBgn03", "FBgn04"),
#'         L4_rep1_FDR = c(0.01, 0.10, 0.04, 0.06),
#'         L4_rep2_FDR = c(0.03, 0.02, 0.50, 0.07),
#'         L5_rep1_FDR = c(0.80, 0.90, 0.01, 0.02),
#'         row.names = c("geneA", "geneB", "geneC", "geneD")
#'     )
#'     list(occupancy = occupancy_df, test_category = "expressed")
#' }
#' mock_data <- .create_mock_fdr_data()
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

    # Extract occupancy data and condition map
    occupancy_df <- if (is(data, "DamIDResults")) inputData(data)$occupancy else data$occupancy
    cond_map <- if (is(data, "DamIDResults")) conditionNames(data) else NULL

    # Centralised validation checks
    if (is.null(occupancy_df)) {
        stop("Could not find 'occupancy' data in the input object.")
    }
    if (!all(c("gene_name", "gene_id") %in% colnames(occupancy_df))) {
        stop("The 'occupancy' table must contain 'gene_name' and 'gene_id' columns.")
    }
    if (is.null(rownames(occupancy_df))) {
        stop("The 'occupancy' table is missing row names, cannot perform filtering.")
    }

    common_cols <- c("gene_id","gene_name","nfrags","name")
    sample_and_fdr_cols <- setdiff(colnames(occupancy_df),common_cols)

    # Find relevant FDR columns
    all_fdr_cols <- grep("_FDR$", sample_and_fdr_cols, value = TRUE)
    if (length(all_fdr_cols) == 0) {
        warning("No '_FDR' columns found in the data. Returning an empty data frame.")
        return(data.frame(gene_name = character(0), gene_id = character(0)))
    }
    all_sample_cols <- setdiff(sample_and_fdr_cols,all_fdr_cols)

    # Determine the internal condition identifier to search for
    internal_condition_id <- condition
    if (!is.null(cond_map) && (condition %in% names(cond_map))) {
        internal_condition_id <- cond_map[condition]
    }

    # Filter FDR columns to those belonging to the specified condition
    relevant_fdr_cols <- grep(internal_condition_id, all_fdr_cols, value = TRUE, fixed = TRUE)
    relevant_sample_cols <- grep(internal_condition_id, all_sample_cols, value = TRUE, fixed = TRUE)

    if (length(relevant_fdr_cols) == 0) {
        warning(sprintf("No FDR columns found matching the condition '%s'.", condition))
        return(data.frame(gene_name = character(0), gene_id = character(0)))
    }
    message(
        sprintf("Found %d FDR columns for condition '%s': %s",
                length(relevant_fdr_cols), condition, paste(relevant_fdr_cols, collapse = ", "))
    )

    # Perform filtering
    fdr_subset <- occupancy_df[, relevant_fdr_cols, drop = FALSE]

    filter_fun <- if (which == "any") {
        function(row) any(row <= fdr, na.rm = TRUE)
    } else {
        function(row) !any(is.na(row)) && all(row <= fdr, na.rm = TRUE)
    }

    keep_indices <- apply(fdr_subset, 1, filter_fun)

    if (sum(keep_indices, na.rm = TRUE) == 0) {
        message(sprintf("No genes passed the FDR <= %s filter for '%s' with rule '%s'.", fdr, condition, which))
        return(data.frame(gene_name = character(0), gene_id = character(0)))
    }

    # Format and return
    result_df <- occupancy_df[keep_indices, c("gene_name", "gene_id", relevant_sample_cols,relevant_fdr_cols), drop = FALSE]
    result_df$avg_occ <- apply(result_df[relevant_sample_cols],1,mean)
    result_df$min_fdr <- apply(result_df[relevant_fdr_cols],1,min)
    result_df <- result_df[c("gene_name", "gene_id","avg_occ","min_fdr")]

    message(
        sprintf("%d genes/loci passed the FDR <= %s filter for condition '%s' with rule '%s'.",
                nrow(result_df), fdr, condition, which)
    )

    return(result_df)
}


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
    input_data <- inputData(diff_results)
    occupancy_df <- input_data$occupancy

    if (is.null(occupancy_df)) {
        warning("fdr_filter_threshold was set, but no 'occupancy' data found. Returning all loci unfiltered.")
        return(row_names_to_filter)
    }

    all_fdr_cols <- grep("_FDR$", colnames(occupancy_df), value = TRUE)

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

