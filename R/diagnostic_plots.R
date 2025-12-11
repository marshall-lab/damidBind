#' @title Verify Underlying Assumptions for `limma` Analysis
#'
#' @description
#' This diagnostic function is a wrapper around the internal
#' `._plot_limma_diagnostics_internal()` funtion, to help assess whether the
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
            title = "Homoscedasticity",
            x = "Fitted Values (Average Occupancy)",
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
        xlab = "Average log-expression",
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
        # trend=FALSE -- add a horizontal line for the common prior
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

    # Process pccupancy data
    occupancy_df <- loaded_data$occupancy
    if (!is.null(occupancy_df)) {
        base_cols <- c("chr", "start", "end", "name", "gene_name", "gene_id", "nfrags")
        sample_cols <- setdiff(colnames(occupancy_df), base_cols)
        sample_cols <- sample_cols[!grepl("_FDR$", sample_cols)]

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

        colnames(raw_matrix) <- extract_unique_sample_ids(colnames(raw_matrix))

        max_rows <- 100000
        was_sampled <- FALSE
        if (nrow(raw_matrix) > max_rows) {
            message(sprintf(" - Raw binding data has > %d rows, sampling %d rows for diagnostics.", max_rows, max_rows))
            raw_matrix <- raw_matrix[sample(nrow(raw_matrix), max_rows), ]
            was_sampled <- TRUE
        }

        if (ncol(raw_matrix) > 1 && nrow(raw_matrix) > ncol(raw_matrix)) {
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
            message(" - Skipped diagnostics for raw binding data: not enough features/rows or samples.")
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
