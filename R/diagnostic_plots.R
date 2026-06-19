#' @title Verify Underlying Assumptions for `limma` Analysis
#'
#' @description
#' This diagnostic function is a wrapper around the internal
#' `._plot_limma_diagnostics_internal()` function, to help assess whether the
#' assumptions of the `limma` empirical Bayes framework hold for a given dataset.
#' It generates a series of plots to check for normality of residuals, homoscedasticity,
#' and the mean-variance relationship, illustrating in particular the effect of
#' `trend` and `robust` parameters to `limma::eBayes`.
#'
#' During `limma`-based fits, the internal plot routine is called by default.  This
#' wrapper allows diagnostics to be displayed for any given log2 ratio-based `data_list`
#' object from `load_data_peaks()` or `load_data_genes()`, and the effect of moderation
#' parameters on the fit tested.
#'
#' @param data_list List. The output from `load_data_peaks` or `load_data_genes`.
#' @param cond A named character vector of length two defining the conditions for
#'   comparison, identical to the `cond` argument in `differential_binding`.
#' @param drop_samples An optional character vector of sample names or patterns to
#'   remove for this diagnostic check. Default: `NULL`.
#' @param filter_occupancy NULL or integer. See `prep_data_for_differential_analysis`.
#'   Default is `TRUE`.
#' @param filter_threshold Numeric. Threshold value for `filter_occupancy`.
#'   (default: 0)
#' @param eBayes_trend Logical. See `limma::eBayes`. Default: `TRUE`
#' @param eBayes_robust Logical. See `limma::eBayes`. Default: `TRUE`
#' @param regex Logical. Whether to use regular expressions for matching condition
#'   names. Default is `FALSE`.
#'
#' @details
#' The function first prepares the data and fits a linear model using the
#' `limma` package. It then calls an internal plotting routine to generate the
#' following checks:
#' \enumerate{
#'   \item \strong{Homoscedasticity (Residuals vs. Fitted):} A scatter plot of
#'     model residuals against fitted values. A random cloud around y=0 supports
#'     the assumption of constant variance.
#'   \item \strong{Effect of eBayes moderation:} Histograms of t-statistics
#'     before and after empirical Bayes moderation.
#'   \item \strong{Mean-variance trend (`plotSA`):} The primary diagnostic
#'     for the `eBayes` step, showing the relationship between average
#'     log2 occupancy and variance.  Points should be evenly distributed around
#'     the blue trendline; any outliers are highlighted in red.
#' }
#' The function uses the internal `prep_data_for_differential_analysis` function
#' to ensure that the data being tested is identical to that used in the main
#' differential analysis.
#'
#' @return Invisibly returns `NULL`. This function is called to
#'   generate diagnostic plots in the active graphics device.
#'
#' @examples
#' mock_genes_gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle("2L", 7),
#'     ranges = IRanges::IRanges(
#'         start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'         end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'     ),
#'     gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'     gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene")
#' )
#'
#' data_dir <- system.file("extdata", package = "damidBind")
#'
#' loaded_data <- load_data_peaks(
#'     binding_profiles_path = data_dir,
#'     peaks_path = data_dir,
#'     ensdb_genes = mock_genes_gr,
#'     quantile_norm = TRUE,
#'     plot_diagnostics = FALSE
#' )
#'
#' conditions <- c("L4 Neurons" = "L4", "L5 Neurons" = "L5")
#'
#' plot_limma_diagnostics(
#'     data_list = loaded_data,
#'     cond = conditions
#' )
#'
#' @export
plot_limma_diagnostics <- function(data_list,
                                   cond,
                                   drop_samples = NULL,
                                   filter_occupancy = TRUE,
                                   filter_threshold = 0,
                                   eBayes_trend = TRUE,
                                   eBayes_robust = TRUE,
                                   regex = FALSE) {
    if (!is.null(drop_samples)) {
        message("Temporarily dropping samples for assumption verification.")
        message("Note: occupancy data is subsetted; it is not re-calculated.")
        data_list <- ._drop_input_samples(data_list, drop_samples)
    }

    # Prepare data for analysis
    message("Preparing data using 'prep_data_for_differential_analysis'...")
    prep_results <- prep_data_for_differential_analysis(
        data_list = data_list,
        cond = cond,
        regex = regex,
        filter_occupancy = filter_occupancy,
        filter_threshold = filter_threshold
    )

    mat <- as.matrix(prep_results$mat)
    factors <- prep_results$factors
    cond_internal <- prep_results$cond_internal

    if (nrow(mat) < 10) {
        warning(
            "Fewer than 10 loci available for analysis. ",
            "Diagnostic plots may not be meaningful."
        )
    }

    # Fit the linear model
    message("Fitting linear model with 'lmFit'...")
    group <- factor(factors$condition, levels = cond_internal)
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- cond_internal
    fit <- limma::lmFit(mat, design)

    # Define contrasts and run empirical Bayes moderation
    message("Running empirical Bayes moderation with 'eBayes'...")
    contrast_str <- sprintf("%s-%s", cond_internal[1], cond_internal[2])
    contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    fit2 <- limma::eBayes(fit2, trend = eBayes_trend, robust = eBayes_robust)

    # Call the internal plotting function with the computed fit objects
    ._plot_limma_diagnostics_internal(fit = fit, fit2 = fit2, mat = mat)

    invisible(NULL)
}


#' @title Internal Routine for Plotting Limma Diagnostics
#' @description
#' This internal function generates a panel of diagnostic plots from pre-computed
#' limma fit objects.
#'
#' @param fit An MArrayLM object from `limma::lmFit`.
#' @param fit2 An MArrayLM object from `limma::eBayes`.
#' @param mat The numeric matrix of expression/occupancy values that was used to
#'   generate the `fit` object. Required for calculating residuals.
#'
#' @return Invisibly returns NULL.
#'
#' @noRd
._plot_limma_diagnostics_internal <- function(fit, fit2, mat) {

    # Homoscedasticity (Residuals vs. Fitted)
    fitted_vals <- stats::fitted(fit)
    resids <- stats::residuals(fit, y = mat)
    resids_vs_fitted_df <- data.frame(
        residual = as.vector(resids),
        fitted = as.vector(fitted_vals)
    )

    p_rvf <- ggplot2::ggplot(resids_vs_fitted_df, ggplot2::aes(x = .data$fitted, y = .data$residual)) +
        ggplot2::geom_point(alpha = 0.05, shape = 20) +
        ggplot2::geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
        ggplot2::labs(
            title = "Residuals vs Fitted",
            x = "Fitted values (average occupancy)",
            y = "Residuals"
        ) +
        ggplot2::theme_bw()

    # Effect of eBayes moderation (histograms of t-statistics)
    raw_t <- fit2$coefficients[, 1] / (fit2$stdev.unscaled[, 1] * fit2$sigma)
    moderated_t <- fit2$t[, 1]

    t_df <- data.frame(
        value = c(raw_t, moderated_t),
        Type = factor(rep(c("Unmoderated", "Moderated"), each = length(raw_t)),
                      levels = c("Unmoderated", "Moderated"))
    )

    p_thist <- ggplot2::ggplot(t_df, ggplot2::aes(x = .data$value)) +
        ggplot2::geom_histogram(bins = 75, fill = "grey70", color = "black") +
        ggplot2::facet_wrap(~Type, scales = "free_x") +
        ggplot2::labs(
            title = "eBayes Moderation",
            x = "t-statistic",
            y = "Frequency"
        ) +
        ggplot2::theme_bw()

    # Mean-variance trend (plotSA)
    p_sa <- ggPlotSA(fit2)

    # Combine and display plots
    message("  Displaying diagnostic plots...")
    final_plot <- (p_rvf | p_sa) / p_thist +
        patchwork::plot_annotation(title = "Diagnostic plots for limma eBayes")

    print(final_plot)

    invisible(NULL)
}


#' @title Create a ggplot2 Mean-Variance (SA) Plot
#'
#' @description
#' This function generates a `ggplot2` version of the mean-variance `plotSA` from a
#' `limma` `MArrayLM` object after running `eBayes`. It visualises the relationship
#' between the average log-expression and the square root of the residual
#' standard deviation (`sigma`).
#'
#' @param fit An `MArrayLM` object, typically the output from `eBayes`.
#' @param main Character. The main title for the plot.
#' @param xlab Character.
#' @param ylab Character.
#' @param point_size Numeric.
#' @param point_alpha Numeric.
#' @param outlier_col Character. Colour for points flagged as potential outliers.
#' @param base_col Character. Colour for points not flagged as outliers.
#' @param trend_col Character. Colour for the fitted trend line.
#' @param subtitle Character. Optional subtitle for the plot.
#'
#' @return A `ggplot` object.
#' @noRd
ggPlotSA <- function(
        fit,
        main = "Mean-variance trend (SA plot)",
        xlab = "Average log2 occupancy",
        ylab = expression(sqrt(sigma)), # Use expression for nice sigma symbol
        point_size = 0.5,
        point_alpha = 0.2,
        outlier_col = "red",
        base_col = "black",
        trend_col = "blue"
) {

    # Input validation
    if (!is(fit, "MArrayLM")) {
        stop("Input 'fit' must be an MArrayLM object from limma.")
    }
    if (is.null(fit$Amean) || is.null(fit$sigma) || is.null(fit$s2.prior)) {
        stop("Input 'fit' object is missing required components (Amean, sigma, s2.prior). Did you run eBayes?")
    }

    plot_df <- data.frame(
        Amean = fit$Amean,
        sqrt_sigma = sqrt(fit$sigma)
    )

    # Identify outliers (mimics plotSA logic for robust eBayes)
    plot_df$color <- base_col
    plot_df$is_outlier <- "Normal"

    if (length(fit$df.prior) > 1L) {
        df2 <- max(fit$df.prior, na.rm = TRUE)
        s2_ratio <- fit$sigma^2 / fit$s2.prior

        # Calculate two-sided p-values from the F-distribution
        p_down <- pf(s2_ratio, df1 = fit$df.residual, df2 = df2, lower.tail = TRUE)
        p_up <- pf(s2_ratio, df1 = fit$df.residual, df2 = df2, lower.tail = FALSE)
        p_two_sided <- 2 * pmin(p_down, p_up)

        FDR <- p.adjust(p_two_sided, method = "BH")

        # Apply the same 0.5 FDR cutoff as the original plotSA
        plot_df$color[FDR <= 0.5 & !is.na(FDR)] <- outlier_col
        plot_df$is_outlier[FDR <= 0.5 & !is.na(FDR)] <- "Outlier"
    }

    p <- ggplot(plot_df, aes(x = .data$Amean, y = .data$sqrt_sigma)) +
        geom_point(
            aes(color = .data$is_outlier),
            size = point_size,
            alpha = point_alpha
        ) +
        scale_colour_manual(
            name = "Variance",
            values = c("Normal" = base_col, "Outlier" = outlier_col)
        ) +
        labs(
            title = main,
            x = xlab,
            y = ylab
        ) +
        theme_bw() +
        theme(legend.position = if (length(fit$df.prior) > 1L) "topright" else "none")

    # Add trend Line
    if (length(fit$s2.prior) == 1L) {
        # trend=FALSE: add a horizontal line for the common prior
        # The s2.prior is the prior variance, so sqrt(sigma) is sqrt(sqrt())
        p <- p + geom_hline(yintercept = sqrt(sqrt(fit$s2.prior)), color = trend_col, linewidth = 1)

    } else {
        # trend=TRUE: Add a line plot for the mean-dependent prior
        # Create a data frame for the trend line to ensure correct ordering
        trend_df <- data.frame(
            Amean = fit$Amean,
            s2_prior_sqrt_sqrt = sqrt(sqrt(fit$s2.prior))
        )

        # Remove NAs and order by Amean to draw the line correctly
        trend_df <- trend_df[complete.cases(trend_df), ]
        trend_df <- trend_df[order(trend_df$Amean), ]

        p <- p + geom_line(data = trend_df, aes(x = .data$Amean, y = .data$s2_prior_sqrt_sqrt), color = trend_col, linewidth = 1)
    }

    return(p)
}


#' Display diagnostic plots for input data
#'
#' @description This function creates and displays diagnostic plots
#' (PCA and correlation heatmap) for both occupancy and raw binding data. It is
#' called by `load_data_peaks` and `load_data_genes`.
#'
#' @param loaded_data A list object, the output of `load_data_peaks`
#'   or `load_data_genes`.
#' @param drop_samples An optional character vector of sample names or patterns to
#'   remove for this diagnostic check. When used, the occupancy data is subsetted,
#'   not recalculated, providing an approximation of the effect of dropping samples.
#'   Default: `NULL`.
#'
#' @return Returns the input `loaded_data` object invisibly
#'
#' @examples
#' # Mock ensdb data to avoid network access
#' mock_genes_gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle("2L", 7),
#'     ranges = IRanges::IRanges(
#'         start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'         end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'     ),
#'     gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'     gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene")
#' )
#' data_dir <- system.file("extdata", package = "damidBind")
#'
#' # Load the example package data
#' loaded_data <- load_data_peaks(
#'     binding_profiles_path = data_dir,
#'     peaks_path = data_dir,
#'     ensdb_genes = mock_genes_gr,
#'     plot_diagnostics = FALSE # don't call the function here ...
#' )
#'
#' # Plot diagnostics
#' plot_input_diagnostics(loaded_data) # ... so that we can call it explicity :/
#'
#' @export
plot_input_diagnostics <- function(loaded_data, drop_samples = NULL) {

    if (!is.null(drop_samples)) {
        message("Temporarily dropping samples for diagnostic plotting.")
        message("Note: occupancy data is subsetted from the full analysis; it is not re-calculated.")
        loaded_data <- ._drop_input_samples(loaded_data, drop_samples)
    }

    message("Generating diagnostic plots...")

    # Process occupancy data
    occupancy_df <- loaded_data$occupancy
    if (!is.null(occupancy_df)) {
        base_cols <- c("chr", "start", "end", "name", "gene_name", "gene_id", "nfrags")
        sample_cols <- setdiff(colnames(occupancy_df), base_cols)
        sample_cols <- sample_cols[!grepl("(_FDR|_pval)$", sample_cols)]

        if (length(sample_cols) > 1) {
            occ_matrix <- as.matrix(occupancy_df[, sample_cols, drop = FALSE])
            colnames(occ_matrix) <- extract_unique_sample_ids(colnames(occ_matrix))

            occ_matrix <- stats::na.omit(occ_matrix)

            if (nrow(occ_matrix) > length(sample_cols)) {
                pca_occ <- .make_pca_plot(occ_matrix, title = "PCA of Samples")
                heatmap_occ <- .make_corr_heatmap(occ_matrix, title = "Sample Correlation")

                heatmap_grob_occ <- grid::grid.grabExpr(ComplexHeatmap::draw(heatmap_occ))

                combined_plot_occ <- pca_occ + heatmap_grob_occ +
                    patchwork::plot_annotation(title = "Diagnostic plots for locus occupancy data")

                print(combined_plot_occ)
            } else {
                message(" - Skipped diagnostics for occupancy data: not enough features/rows for analysis.")
            }
        } else {
            message(" - Skipped diagnostics for occupancy data: not enough sample columns.")
        }
    }

    # Process raw genomic data
    binding_profiles_gr <- loaded_data$binding_profiles_data
    if (!is.null(binding_profiles_gr)) {
        raw_matrix <- as.matrix(mcols(binding_profiles_gr))
        raw_matrix <- stats::na.omit(raw_matrix)

        # Check for at least two samples before attempting diagnostics
        if (ncol(raw_matrix) > 1) {
            colnames(raw_matrix) <- extract_unique_sample_ids(colnames(raw_matrix))

            max_rows <- 100000
            was_sampled <- FALSE
            if (nrow(raw_matrix) > max_rows) {
                message(sprintf(" - Raw binding data has > %d rows, sampling %d rows for diagnostics.", max_rows, max_rows))
                raw_matrix <- raw_matrix[sample(nrow(raw_matrix), max_rows), , drop = FALSE]
                was_sampled <- TRUE
            }

            if (nrow(raw_matrix) > ncol(raw_matrix)) {
                pca_raw <- .make_pca_plot(raw_matrix, title = "PCA of Samples")
                heatmap_raw <- .make_corr_heatmap(raw_matrix, title = "Sample Correlation")
                heatmap_grob_raw <- grid::grid.grabExpr(ComplexHeatmap::draw(heatmap_raw))

                plot_title <- "Diagnostic plots for raw genomic binding data"
                if (was_sampled) {
                    plot_title <- paste(plot_title, "(Sampled)")
                }

                combined_plot_raw <- pca_raw + heatmap_grob_raw +
                    patchwork::plot_annotation(title = plot_title)

                print(combined_plot_raw)
            } else {
                message(" - Skipped diagnostics for raw binding data: not enough features/rows.")
            }
        } else {
            message(" - Skipped diagnostics for raw binding data: not enough sample columns.")
        }
    }

    return(invisible(loaded_data))
}


#' PCA plot for data diagnostics
#' @param data_matrix A numeric matrix with features in rows and samples in columns
#' @param title The title for the plot.
#' @return A ggplot object.
#' @noRd
.make_pca_plot <- function(data_matrix, title = "PCA Plot") {

    # Identify features with zero variance and remove
    variances <- apply(data_matrix, 1, stats::var, na.rm = TRUE)
    non_constant_rows <- variances > 1e-10

    num_filtered_rows <- sum(!non_constant_rows)

    if (sum(non_constant_rows) < 2) {
        message(sprintf(" - Skipped PCA '%s': Fewer than 2 features with non-zero variance.", title))
        return(ggplot() + theme_void() + labs(title = title, subtitle = "Skipped: Insufficient data variance"))
    }

    if (num_filtered_rows>0) {
        message(sprintf("  %d rows with zero variance were filtered.",num_filtered_rows))
    }
    filtered_matrix <- data_matrix[non_constant_rows, , drop = FALSE]

    pca_res <- stats::prcomp(t(filtered_matrix), scale. = TRUE, center = TRUE)
    pca_data <- as.data.frame(pca_res$x)

    var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

    p <- ggplot(pca_data, aes(x = .data$PC1, y = .data$PC2)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(aes(label = rownames(pca_data)), point.padding = 0.2, box.padding = 0.2) +
        labs(
            title = title,
            x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
            y = sprintf("PC2 (%.1f%% variance)", var_explained[2])
        ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        coord_fixed()

    return(p)
}


#' Correlation heatmap for diagnostics
#' @param data_matrix A numeric matrix with features in rows and samples in columns
#' @param title The title for the plot.
#' @return A Heatmap object.
#' @noRd
.make_corr_heatmap <- function(data_matrix, title = "Correlation Heatmap") {
    cor_matrix <- stats::cor(data_matrix, use = "pairwise.complete.obs")

    # Define a robust color ramp
    min_cor <- min(cor_matrix, na.rm = TRUE)
    max_cor <- max(cor_matrix, na.rm = TRUE)

    if (isTRUE(all.equal(min_cor, max_cor))) {
        # Handle case where all correlations are identical
        breaks <- c(min_cor - 0.01, min_cor, min_cor + 0.01)
    } else {
        breaks <- c(min_cor, (min_cor + max_cor) / 2, max_cor)
    }

    col_fun <- circlize::colorRamp2(breaks, c("deepskyblue4", "white", "firebrick3"))

    ht <- ComplexHeatmap::Heatmap(
        cor_matrix,
        name = "Corr.",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        column_title = title,
        heatmap_legend_param = list(title = "Correlation"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 5.5))
        }
    )
    return(ht)
}

#' Extract unique sample names from complex labels
#'
#' This function takes a vector of complex sample labels and iteratively
#' constructs a simplified, unique name for each. It identifies all blocks of
#' text that differ across the sample set and progressively adds them to a base
#' name until the combination of the base name and a replicate identifier is
#' unique for every sample.
#'
#' @param sample_names A character vector of sample labels.
#' @param delimiter A regular expression used as a delimiter to split labels
#'        into blocks. (Default: `[-_\\.]`)
#' @param replicate_pattern A regular expression used to identify the
#'        replicate block. (Default: `^(n|N|r|rep|replicate|sample)\\d+`)
#'
#' @return A vector of simplified, unique names. If a unique name cannot be
#'   formed or essential information is missing for a sample, the original
#'   label for that sample is returned as a fallback.
#'
#' @export
#'
#' @examples
#' labels <- c(
#'   "RNAPII_elav-GSE77860-n1-SRR3164378-2017-vs-Dam.scaled.kde-norm",
#'   "RNAPII_elav-GSE77860-n2-SRR3164379-2017-vs-Dam.scaled.kde-norm",
#'   "RNAPII_elav-GSE77860-n4-SRR3164380-2017-vs-Dam.scaled.kde-norm",
#'   "RNAPII_Wor-GSE77860-n1-SRR3164346-2017-vs-Dam.scaled.kde-norm",
#'   "RNAPII_Wor-GSE77860-n2-SRR3164347-2017-vs-Dam.scaled.kde-norm",
#'   "RNAPII_Wor-GSE77860-sample1-SRR2038537-2017-vs-Dam.scaled.kde-norm"
#' )
#' extract_unique_sample_ids(labels)
#'
extract_unique_sample_ids <- function(sample_names,
                                      delimiter = "[-_\\.]",
                                      replicate_pattern = "^(n|N|r|rep|replicate|sample)\\d+") {

    # Input validation
    if (!is.character(sample_names) || length(sample_names) < 2) {
        warning("Input 'sample_names' must be a character vector with at least two elements. Returning original names.")
        return(sample_names)
    }

    # Tokenise and create matrix
    tokenised_list <- strsplit(sample_names, split = delimiter, perl = TRUE)
    max_len <- max(vapply(tokenised_list, length, FUN.VALUE = integer(1)))
    padded_list <- lapply(tokenised_list, function(x) c(x, rep(NA, max_len - length(x))))
    token_matrix <- do.call(rbind, padded_list)

    if (ncol(token_matrix) == 0) {
        warning("Could not tokenise sample names. Returning original names.")
        return(sample_names)
    }

    # Identify replicate information
    replicate_info <- lapply(tokenised_list, function(tokens) {
        match_indices <- grep(replicate_pattern, tokens, perl = TRUE, value = FALSE)
        if (length(match_indices) > 0) {
            # Return the first match's index and value
            list(index = match_indices[1], value = tokens[match_indices[1]])
        } else {
            list(index = NA_integer_, value = NA_character_)
        }
    })
    replicate_values <- vapply(replicate_info, `[[`, "value", FUN.VALUE = character(1))
    replicate_col_indices <- vapply(replicate_info, `[[`, "index", FUN.VALUE = integer(1))

    # Identify all differing columns (excluding replicates)
    num_unique_per_col <- apply(token_matrix, 2, function(col) length(unique(na.omit(col))))
    all_diff_indices <- which(num_unique_per_col > 1)

    unique_replicate_cols <- unique(na.omit(replicate_col_indices))
    # These are the columns we will iterate through to build the base name
    iterative_diff_indices <- setdiff(all_diff_indices, unique_replicate_cols)

    # Handle no differing cols
    if (length(iterative_diff_indices) == 0) {
        # If replicates exist, use them. Otherwise, fall back.
        final_names <- ifelse(
            !is.na(replicate_values),
            replicate_values,
            sample_names
        )
        if (length(unique(final_names)) < length(sample_names)) {
            warning("Could not generate unique names from replicates alone. Returning original names where collisions occurred.")
            return(make.unique(final_names))
        }
        return(final_names)
    }


    # Iteratively build the unique id basenames
    base_name_components <- list(token_matrix[, iterative_diff_indices[1], drop = FALSE])

    for (i in seq_along(iterative_diff_indices)) {

        # Construct current base name from components
        current_base_name <- apply(do.call(cbind, base_name_components), 1, paste, collapse = "_")

        # Combine with replicate for uniqueness check
        temp_full_names <- paste(current_base_name, replicate_values, sep = "_")

        # Check for uniqueness
        if(length(unique(temp_full_names[!is.na(replicate_values)])) == length(temp_full_names[!is.na(replicate_values)])) {
            break # Names are unique
        }

        if ((i + 1) <= length(iterative_diff_indices)) {
            next_col_idx <- iterative_diff_indices[i + 1]
            base_name_components[[length(base_name_components) + 1]] <- token_matrix[, next_col_idx, drop = FALSE]
        }
    }

    # Final basename assembly
    final_base_names <- apply(do.call(cbind, base_name_components), 1, paste, collapse = "_")

    # Combine the final base name with the replicate value
    adjusted_names <- paste(final_base_names, replicate_values, sep = "_")

    # If any sample is missing a replicate, the name will end in "_NA".
    # We fall back to the original name for any sample that could not be fully processed.
    missing_info_indices <- is.na(replicate_values)
    adjusted_names[missing_info_indices] <- sample_names[missing_info_indices]

    # Final check: if the process still resulted in duplicates, use make.unique as a last resort
    if(length(unique(adjusted_names)) < length(sample_names)) {
        warning("Iterative process failed to produce fully unique names. Applying `make.unique` as a fallback.")
        adjusted_names <- make.unique(adjusted_names, sep = ".")
    }

    return(adjusted_names)
}


#' @title Plot occupancy model diagnostics
#'
#' @description
#' Generates diagnostic plots for the non-linear regression models used in occupancy
#' estimation to assess model fit and the appropriateness of the WLS approach.
#'
#' @param model_params_df Dataframe. The parameters extracted from first-level
#'   simulations.
#' @param fdr_models List. The WLS models for slope, intercept, and MSE.
#' @param sample_name Character. The name of the sample being visualised.
#'
#' @return Invisibly returns NULL.
#'
#' @noRd
plot_occupancy_diagnostics <- function(model_params_df, fdr_models, fdr_matrix=NULL, threshold_samples=NULL, sample_name) {

    # Create a sequence of fragment counts for smooth trend lines
    # Using length.out = 200 provides a smoother curve for the spline fit
    frag_seq <- seq(
        min(model_params_df$fragment_count),
        max(model_params_df$fragment_count),
        length.out = 200
    )

    # Prepare the prediction data frame
    trend_df <- data.frame(fragment_count = frag_seq)
    trend_df$log_f <- log(trend_df$fragment_count)

    # Predict values across the sequence using the models
    trend_df$slope_pred <- stats::predict(fdr_models$slope_model, newdata = trend_df)
    trend_df$int_pred <- stats::predict(fdr_models$intercept_model, newdata = trend_df)
    trend_df$mse_pred <- stats::predict(fdr_models$mse_model, newdata = trend_df)

    # Slope fit plot
    # Size aesthetic is 1/SE to show WLS weighting influence
    p_slope <- ggplot2::ggplot(model_params_df, ggplot2::aes(x = .data$fragment_count, y = .data$slope)) +
        ggplot2::geom_point(ggplot2::aes(size = 1 / .data$se_slope), alpha = 0.5, colour = "#333399") +
        ggplot2::geom_line(data = trend_df, ggplot2::aes(y = .data$slope_pred), colour = "#CC0000", linewidth = 0.8, alpha = 0.6) +
        ggplot2::labs(
            title = "Slope Regression",
            #subtitle = "Natural spline fit (log-scale)",
            x = "Fragment Count",
            y = "Log-Linear Slope"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")

    # Intercept fit plot
    p_int <- ggplot2::ggplot(model_params_df, ggplot2::aes(x = .data$fragment_count, y = .data$intercept)) +
        ggplot2::geom_point(ggplot2::aes(size = 1 / .data$se_int), alpha = 0.5, colour = "#226622") +
        ggplot2::geom_line(data = trend_df, ggplot2::aes(y = .data$int_pred), colour = "#CC0000", linewidth = 0.8, alpha = 0.6) +
        ggplot2::labs(
            title = "Intercept Regression",
            #subtitle = "Natural spline fit (log-scale)",
            x = "Fragment Count",
            y = "Log-Linear Intercept"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")

    # MSE fit plot
    # Visualises the variance component used for the Jensen's correction factor
    p_mse <- ggplot2::ggplot(model_params_df, ggplot2::aes(x = .data$fragment_count, y = .data$mse)) +
        ggplot2::geom_point(alpha = 0.5, colour = "#996633") +
        ggplot2::geom_line(data = trend_df, ggplot2::aes(y = .data$mse_pred), colour = "#CC0000", linewidth = 0.8, alpha = 0.6) +
        ggplot2::labs(
            title = "MSE Regression",
            #subtitle = "Natural spline fit (log-scale)",
            x = "Fragment Count",
            y = "Residual Variance"
        ) +
        ggplot2::theme_bw()

    # Heteroscedasticity plot
    p_wls <- ggplot2::ggplot(model_params_df, ggplot2::aes(x = .data$fragment_count, y = .data$se_slope)) +
        ggplot2::geom_point(colour = "#663399", alpha = 0.6) +
        ggplot2::geom_smooth(
            method = "lm",
            formula = y ~ splines::ns(x, df = 3),
            colour = "#663399",
            linetype = "dashed",
            linewidth = 0.5,
            alpha = 0.4,
            se = FALSE
        ) +
        ggplot2::labs(
            title = "Heteroscedasticity",
            #subtitle = "Estimate precision by fragment count",
            x = "Fragment Count",
            y = "SE (Slope)"
        ) +
        ggplot2::theme_bw()

    final_plot <- (p_slope + p_int + p_mse) +
        patchwork::plot_annotation(
            title = "Occupancy model diagnostics",
            subtitle = sample_name,
            theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
        )

    print(final_plot)
    invisible(NULL)
}

#' Plot mean-variance relationship for CATaDa results
#'
#' @description
#' This diagnostic function assesses whether CATaDa occupancy signal follows
#' the theoretical mean-variance relationship assumed by Negative Binomial
#' models. It is only applicable to objects generated via the
#' `differential_accessibility()` function.
#'
#' @param diff_results A DamIDResults object.
#'
#' @return A ggplot object.
#'
#' @export
plot_catada_mean_variance <- function(diff_results) {
    # Validate input class
    if (!is(diff_results, "DamIDResults")) {
        stop("Input 'diff_results' must be a DamIDResults S4 object.")
    }

    # Extract metadata
    cond_map <- conditionNames(diff_results)
    data_list <- inputData(diff_results)
    category <- data_list$test_category

    # Only allow execution for CATaDa data (not valid for log2 ratio data)
    if (is.null(category) || category != "accessible") {
        stop(sprintf(
            "This plot is only relevant for CATaDa accessibility data (test_category: 'accessible'). Current category is '%s'.",
            if (is.null(category)) "NULL" else category
        ))
    }

    # Extract coordinates and columns
    occ_df <- data_list$occupancy
    matched_list <- data_list$matched_samples
    analysis_loci <- rownames(analysisTable(diff_results))

    if (is.null(occ_df)) {
        stop("The results object does not contain the required occupancy data.")
    }

    # Restrict to the set of loci used in the final analysis
    occ_df <- occ_df[analysis_loci, , drop = FALSE]

    # Map display names to internal sample list keys
    cond_display <- names(cond_map)
    cond_internal <- as.character(cond_map)

    # Determine empirical mean and variance per condition
    stats_list <- lapply(seq_along(cond_display), function(i) {
        display_name <- cond_display[i]
        internal_id <- cond_internal[i]

        # Get replicates identified during the differential prep stage
        samples <- matched_list[[internal_id]]

        if (length(samples) < 2) {
            warning(sprintf("Condition '%s' has fewer than 2 replicates; skipping variance calculation", display_name))
            return(NULL)
        }

        # Calculate mean and variance
        cond_mat <- as.matrix(occ_df[, samples, drop = FALSE])
        data.frame(
            mean_val = rowMeans(cond_mat, na.rm = TRUE),
            var_val = apply(cond_mat, 1, stats::var, na.rm = TRUE),
            condition = display_name
        )
    })

    plot_df <- do.call(rbind, stats_list)

    # Keep condition display names
    plot_df$condition <- factor(plot_df$condition, levels = cond_display)

    # Filter non-positive values to prevent log transformation errors
    plot_df <- plot_df[plot_df$mean_val > 0 & plot_df$var_val > 0, ]

    if (nrow(plot_df) == 0) {
        stop("No data with positive mean and variance found for plotting.")
    }

    # Create theoretical reference lines
    # Var = mean + alpha * mean^2
    mean_seq <- exp(seq(log(min(plot_df$mean_val)), log(max(plot_df$mean_val)), length.out = 100))
    theo_df <- do.call(rbind, lapply(c(0.01, 0.1, 0.5), function(a) {
        data.frame(
            mean_val = mean_seq,
            var_val = mean_seq + a * (mean_seq^2),
            type = paste0("NB (alpha=", a, ")")
        )
    }))

    # Add the Poisson line (var = mean)
    poisson_df <- data.frame(
        mean_val = mean_seq,
        var_val = mean_seq,
        type = "Poisson"
    )
    theo_df <- rbind(theo_df, poisson_df)

    # Plot
    ggplot2::ggplot(plot_df, aes(x = .data$mean_val, y = .data$var_val)) +
        ggplot2::geom_point(alpha = 0.2, size = 0.5, colour = "grey30", shape = 20) +
        ggplot2::geom_line(data = theo_df, aes(colour = .data$type), linewidth = 1) +
        ggplot2::scale_x_log10(labels = scales::label_log()) +
        ggplot2::scale_y_log10(labels = scales::label_log()) +
        ggplot2::facet_wrap(~condition) +
        ggplot2::labs(
            title = "CATaDa mean-variance relationship",
            subtitle = "Empirical variance vs parametric model assumptions",
            x = "Mean Intensity (log scale)",
            y = "Variance (log scale)",
            colour = "Model assumption"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
        )
}

#' @title Internal helper to fit a robust fragment-length bias model
#' @description Fits a GAM to GATC fragment signals against log10-transformed widths.
#' @param widths Numeric vector of fragment widths.
#' @param signals Numeric vector of fragment signals.
#' @return An mgcv::gam object.
#' @noRd
._fit_bias_gam <- function(widths, signals) {
    # Filter for valid data
    valid <- is.finite(signals) & widths > 0

    training_df <- data.frame(
        signal = signals[valid],
        width = widths[valid]
    )

    # Fit GAM with cubic regression spline.
    model <- mgcv::gam(signal ~ s(log10(width), bs = "cr"), data = training_df)
    return(model)
}


#' @title Internal helper to fit a robust fragment-length bias model using medians
#' @description Bins fragments by length and fits a GAM to the median signal per bin.
#'   Handles empty bins by mapping midpoints to observed bin indices.
#' @param widths Numeric vector of fragment widths.
#' @param signals Numeric vector of fragment signals.
#' @param n_bins Integer. Number of bins for the median calculation.
#' @return An mgcv::gam object representing the technical floor.
#' @noRd
._fit_bias_median_gam <- function(widths, signals, n_bins = 100) {
    # Ensure we only work with finite signal and positive widths
    valid <- is.finite(signals) & widths > 0
    if (sum(valid) < n_bins) {
        # Fallback if data is extremely sparse
        n_bins <- max(5, floor(sum(valid) / 10))
    }

    df <- data.frame(w = log10(widths[valid]), s = signals[valid])

    # Define the binning breaks and midpoints
    w_min <- min(df$w)
    w_max <- max(df$w)

    # Add a tiny epsilon to handle the max value correctly in cut()
    breaks <- seq(w_min, w_max, length.out = n_bins + 1)
    midpoints <- breaks[-1] - (diff(breaks) / 2)

    # Assign fragments to bin indices (labels=FALSE returns 1, 2, ..., n_bins)
    df$bin_idx <- cut(df$w, breaks = breaks, labels = FALSE, include.lowest = TRUE)

    # Calculate median signal per bin index
    # aggregate() will drop empty bins, returning only bin_idxs that have data
    bin_stats <- stats::aggregate(s ~ bin_idx, data = df, FUN = stats::median)

    # Correctly map the corresponding midpoints to the non-empty bins
    bin_stats$w_mid <- midpoints[bin_stats$bin_idx]

    # Fit GAM to the binned medians
    model <- mgcv::gam(s ~ s(w_mid, bs = "cr"), data = bin_stats)
    return(model)
}


#' Test weighting discrepancy against systematic length bias
#'
#' @description
#' Tests whether the discrepancy between weighted and simple means is driven by
#' a systematic fragment-length score bias or by underlying biological/stochastic
#' variance. It models any length-based bias by fitting a Generalized
#' Additive Model (GAM) to all non-peak fragments, predicts the expected bias
#' discrepancy per region, and compares it to the observed discrepancy.
#'
#' @param diff_results A `DamIDResults` object
#' @param sample_name Character string of the sample to analyse
#' @param plot Logical. Whether to print the summary plots (Default: TRUE)
#'
#' @return A list containing the discrepancy data, the summary statistics,
#'   and a two-panel patchwork plot for verification and results.
#' @export
test_weighting_vs_bias_artifact <- function(diff_results, sample_name, plot=TRUE) {

    if (!is(diff_results, "DamIDResults")) {
        stop("Input 'diff_results' must be a DamIDResults S4 object.")
    }

    data_list <- inputData(diff_results)
    binding_gr <- data_list$binding_profiles_data
    occ_df <- data_list$occupancy
    category <- data_list$test_category

    # Restrict to loci that were actually tested for differential analysis
    analysis_loci <- rownames(analysisTable(diff_results))
    foreground_occ_df <- occ_df[analysis_loci, , drop = FALSE]

    # Convert test region IDs back to GRanges for overlap analysis
    matches <- stringr::str_match(rownames(foreground_occ_df), "^(.*?):(\\d+)-(\\d+)")
    foreground_gr <- GenomicRanges::GRanges(
        seqnames = matches[, 2],
        ranges = IRanges::IRanges(as.integer(matches[, 3]), as.integer(matches[, 4])),
        name = rownames(foreground_occ_df)
    )

    message("Defining signal-free background...")
    signal_gr <- NULL

    if (category == "bound") {
        # Try to find per-sample peaks
        peaks_list <- data_list$peaks
        # Match using simplified IDs to account for normalisation suffixes
        sample_id_clean <- extract_unique_sample_ids(sample_name)
        peak_match_idx <- grep(sample_id_clean, names(peaks_list))

        if (length(peak_match_idx) > 0) {
            message(sprintf(" - Excluding fragments overlapping peaks from: %s", names(peaks_list)[peak_match_idx[1]]))
            signal_gr <- peaks_list[[peak_match_idx[1]]]
        } else {
            message(" - No specific peaks found for sample; falling back to reduced union peaks (pr)")
            signal_gr <- data_list$pr
        }
    } else if (category == "expressed") {
        # Find genes with FDR < 0.05 in this specific sample to exclude
        fdr_col <- paste0(sample_name, "_FDR")
        if (!fdr_col %in% colnames(occ_df)) {
            # Fallback to pval if FDR is missing (unadjusted)
            fdr_col <- paste0(sample_name, "_pval")
        }

        if (fdr_col %in% colnames(occ_df)) {
            sig_expr_indices <- which(occ_df[[fdr_col]] < 0.05)
            if (length(sig_expr_indices) > 0) {
                message(sprintf(" - Excluding fragments overlapping %d expressed genes (FDR < 0.05)", length(sig_expr_indices)))
                # Extract coordinates from names (chr:start-end)
                expr_names <- rownames(occ_df)[sig_expr_indices]
                expr_matches <- stringr::str_match(expr_names, "^(.*?):(\\d+)-(\\d+)")
                signal_gr <- GenomicRanges::GRanges(
                    seqnames = expr_matches[, 2],
                    ranges = IRanges::IRanges(as.integer(expr_matches[, 3]), as.integer(expr_matches[, 4]))
                )
            }
        }
    }

    # Identify background fragments (no overlap with signal regions)
    all_widths <- GenomicRanges::width(binding_gr)
    all_signals <- S4Vectors::mcols(binding_gr)[[sample_name]]

    bg_mask <- rep(TRUE, length(binding_gr))
    if (!is.null(signal_gr)) {
        signal_overlaps <- GenomicRanges::findOverlaps(binding_gr, signal_gr)
        bg_mask[S4Vectors::queryHits(signal_overlaps)] <- FALSE
    }

    if (sum(bg_mask) < 100) {
        warning("Insufficient background fragments found; using all fragments for bias model.")
        bg_mask <- rep(TRUE, length(binding_gr))
    }

    # Model the purified technical bias using the background fragments
    message("Fitting purified bias models...")
    bias_model_null <- ._fit_bias_gam(all_widths[bg_mask], all_signals[bg_mask])

    # Global median illustrates any zero bias on plot (not used otherwise)
    bias_model_median <- ._fit_bias_median_gam(all_widths, all_signals)

    # Predict expected bias for all fragments
    expected_bias_signals <- stats::predict(bias_model_null,
                                            newdata = data.frame(width = all_widths))

    # Identify overlapping fragments within tested foreground regions
    overlaps <- GenomicRanges::findOverlaps(foreground_gr, binding_gr)
    q_hits <- S4Vectors::queryHits(overlaps)
    s_hits <- S4Vectors::subjectHits(overlaps)

    # Partition variance per genomic region
    message("Calculating observed vs. predicted bias discrepancy per region...")
    w_list <- split(all_widths[s_hits], q_hits)
    s_list <- split(all_signals[s_hits], q_hits)
    b_list <- split(expected_bias_signals[s_hits], q_hits)

    eval_results <- lapply(names(w_list), function(i) {
        w <- w_list[[i]]
        s <- s_list[[i]]
        b <- b_list[[i]]

        if (length(w) < 2 || any(is.na(s))) return(NULL)

        width_cv <- stats::sd(w) / mean(w)
        D_obs <- abs(stats::weighted.mean(s, w) - mean(s))
        D_bias <- abs(stats::weighted.mean(b, w) - mean(b))

        data.frame(width_cv = width_cv, D_obs = D_obs, D_bias = D_bias)
    })

    eval_df <- do.call(rbind, eval_results)

    # Global Summary Statistics
    obs_median  <- stats::median(eval_df$D_obs)
    bias_median <- stats::median(eval_df$D_bias)

    wilcox_res <- stats::wilcox.test(
        eval_df$D_obs,
        eval_df$D_bias,
        paired = TRUE,
        alternative = "greater"
    )

    boot_ratios <- vapply(seq_len(1000), function(i) {
        idx <- sample(nrow(eval_df), replace = TRUE)
        stats::median(eval_df$D_obs[idx]) / max(stats::median(eval_df$D_bias[idx]), 1e-6)
    }, numeric(1))

    final_ratio <- stats::median(boot_ratios)
    ratio_ci    <- stats::quantile(boot_ratios, probs = c(0.025, 0.975))

    summary_stats <- list(
        median_obs_discrepancy = obs_median,
        median_bias_induced    = bias_median,
        ratio                  = final_ratio,
        ratio_ci               = ratio_ci,
        p_value                = wilcox_res$p.value
    )

    message("\n# Weighted bias test (signal-free baseline)")
    message(sprintf("Median Observed Discrepancy:     %0.4f", obs_median))
    message(sprintf("Median Bias-Induced Discrepancy: %0.4f", bias_median))
    message(sprintf("The weighting corrects %0.1fx more variance than can be explained by any systematic length bias.", final_ratio))
    message(sprintf("Wilcoxon p-value (D_obs > D_bias): %e", wilcox_res$p.value))
    message(sprintf("95%% CI for Adjustment Ratio: [%.2f - %.2f]", ratio_ci[1], ratio_ci[2]))

    # Visualisation
    plot_verify_df <- data.frame(width = all_widths, signal = all_signals)
    plot_verify_df <- plot_verify_df[is.finite(plot_verify_df$signal), ]

    fit_range <- seq(min(log10(all_widths[all_widths > 0])), max(log10(all_widths)), length.out = 300)
    fit_line_df <- data.frame(
        width = 10^fit_range,
        null_sig = stats::predict(bias_model_null, newdata = data.frame(width = 10^fit_range)),
        floor_sig = stats::predict(bias_model_median, newdata = data.frame(w_mid = fit_range))
    )

    p1 <- ggplot2::ggplot(plot_verify_df, ggplot2::aes(x = .data$width, y = .data$signal)) +
        ggplot2::geom_point(alpha = 0.01, size = 0.2, colour = "grey30") +
        ggplot2::geom_line(data = fit_line_df, ggplot2::aes(y = .data$null_sig, colour = "Bias model"), linewidth = 1) +
        ggplot2::geom_line(data = fit_line_df, ggplot2::aes(y = .data$floor_sig, colour = "Median"), linewidth = 1, linetype = "dashed") +
        ggplot2::scale_x_log10() +
        ggplot2::scale_colour_manual(values = c("Bias model" = "firebrick", "Median" = "royalblue")) +
        ggplot2::labs(
            title = "DamID-seq fragment distribution by length",
            subtitle = sprintf("Sample: %s", sample_name),
            x = "Fragment length",
            y = "Log2 ratio signal",
            colour = "GAM Model"
        ) +
        ggplot2::theme_bw()+
        ggplot2::theme(
            legend.position = "inside",
            legend.position.inside = c(0.05, 0.95),
            legend.justification = c(0, 1)
        )

    plot_discrep_df <- data.frame(
        cv = rep(eval_df$width_cv, 2),
        val = c(eval_df$D_obs, eval_df$D_bias),
        type = factor(rep(c("Observed Discrepancy", "Bias-Only Prediction"), each = nrow(eval_df)),
                      levels = c("Observed Discrepancy", "Bias-Only Prediction"))
    )

    p2 <- ggplot2::ggplot(plot_discrep_df, ggplot2::aes(x = .data$cv, y = .data$val, colour = .data$type)) +
        ggplot2::geom_point(alpha = 0.1, size = 1, shape = 16) +
        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.2) +
        ggplot2::scale_colour_manual(
            values = c(
                "Observed Discrepancy" = "black",
                "Bias-Only Prediction" = "firebrick"
            ),
            breaks = c(
                "Observed Discrepancy",
                "Bias-Only Prediction"
            ),
            labels = c(
                "Observed",
                "Bias-only predicted"
            )
        ) +
        ggplot2::labs(
            title = "Effect of fragment-length-weighted mean",
            subtitle = "Observed vs bias-only discrepancies from the simple mean",
            x = "Fragment width heterogeneity (CV)",
            y = "|Weighted mean - Simple mean|",
            colour = "Source"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "inside",
            legend.position.inside = c(0.05, 0.95),
            legend.justification = c(0, 1)
        )

    if (isTRUE(plot)) {
        combined_plot <- (p1 + p2) +
            patchwork::plot_annotation(tag_levels = "A") &
            theme(
                plot.tag = element_text(face = "bold", size = 14)
            )

        print(combined_plot)
    }

    return(invisible(list(
        data = eval_df,
        stats = summary_stats,
        plot = if (isTRUE(plot)) combined_plot else NULL
    )))
}


#' Test weighting discrepancy against systematic length bias (batch analysis)
#'
#' @description
#' A wrapper for \code{\link{test_weighting_vs_bias_artifact}} that iterates through
#' all samples present in a \code{DamIDResults} object and returns a summary
#' table of the results.
#'
#' @param diff_results A \code{DamIDResults} object.
#' @param plot_all Logical. Whether to print summary plots for every sample in the batch.
#'   (Default: FALSE)
#'
#' @return A \code{data.frame} summarising the bias test for each sample,
#'   including the median discrepancies, adjustment ratios, and individual p-values.
#' @export
test_weighting_vs_bias_batch <- function(diff_results, plot_all = FALSE) {

    if (!is(diff_results, "DamIDResults")) {
        stop("'diff_results' must be a DamIDResults object.")
    }

    # Extract sample names
    binding_gr <- inputData(diff_results)$binding_profiles_data
    sample_names <- colnames(S4Vectors::mcols(binding_gr))

    if (length(sample_names) == 0) {
        stop("No sample columns found in the results object metadata.")
    }

    message(sprintf("Running batch weighting-bias test for %d samples...", length(sample_names)))

    results_list <- lapply(sample_names, function(sn) {
        res <- test_weighting_vs_bias_artifact(diff_results, sample_name = sn, plot = plot_all)

        data.frame(
            sample = sn,
            median_obs = res$stats$median_obs_discrepancy,
            median_bias = res$stats$median_bias_induced,
            ratio = res$stats$ratio,
            p_value = res$stats$p_value,
            ci_lower = res$stats$ratio_ci[1],
            ci_upper = res$stats$ratio_ci[2]
        )
    })

    summary_df <- do.call(rbind, results_list)
    rownames(summary_df) <- NULL

    # Calculate aggregate stats
    ratios <- summary_df$ratio
    n <- length(ratios)

    if (n > 1) {
        log_ratios <- log(ratios)
        mean_log <- mean(log_ratios)
        se_log <- stats::sd(log_ratios) / sqrt(n)

        # Two-tailed critical t-value
        t_crit <- stats::qt(0.975, df = n - 1)

        # Geometric mean and CI
        geom_mean <- exp(mean_log)
        ci_lower  <- exp(mean_log - (t_crit * se_log))
        ci_upper  <- exp(mean_log + (t_crit * se_log))

        message("\n## Batch Summary (Log-Normal Model)")
        message(sprintf("Samples tested: %d", n))
        message(sprintf("Geometric Mean Adjustment Ratio: %0.2fx", geom_mean))
        message(sprintf("95%% Population CI: [%.2f - %.2f]", ci_lower, ci_upper))
    }

    return(summary_df)
}



