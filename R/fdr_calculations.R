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
#' null model to each sample's binding profile and adds the resulting FDR values
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
#' The FDR calculation algorithm is adapted from the method described in
#' Southall et al. (2013), Dev Cell, 26(1), 101-12, and implemented in the
#' `polii.gene.call` tool. It operates in several stages:
#' \enumerate{
#'   \item For each sample, a null distribution of mean occupancy scores is
#'     simulated by repeatedly sampling random fragments from the genome-wide
#'     binding profile. This is done for various numbers of fragments.
#'   \item A two-tiered regression model is fitted to the simulation results. This
#'     creates a statistical model that can predict the FDR for any given
#'     occupancy score and fragment count.
#'   \item This predictive model is then applied to the actual observed mean
#'     occupancy and fragment count for each gene in the `occupancy_df`.
#'   \item The calculated FDR value for each gene in each sample is appended to the
#'     `occupancy_df` in a new column.
#' }
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
        occupancy_df[[fdr_col_name]] <- fdr_column_values
    }
    message("FDR calculation complete.")

    return(occupancy_df)
}


