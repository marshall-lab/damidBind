#' Differential binding/expression analysis (`limma`)
#'
#' Setup and differential analysis for occupancy/binding experiments
#' using `limma`. Accepts output from `load_data_peaks` or `load_data_genes`,
#' prepares an experiment matrix, fits linear models, and returns DE loci.
#'
#' @param data_list List. Output from load_data_peaks or load_data_genes.
#' @param cond A named or unnamed character vector of length two. The values are
#'   strings or regular expressions used to identify samples for each condition.
#'   If the vector is named, the names are used as user-friendly display names
#'   for the conditions in plots and outputs. If unnamed, the match strings are
#'   used as display names. The order determines the contrast, e.g., `cond[1]` vs `cond[2]`.
#' @param regex Logical. If `TRUE`, the strings in `cond` are treated as
#'   regular expressions for matching sample names. If `FALSE` (the default),
#'   fixed string matching is used.
#' @param fdr Numeric. FDR threshold for significance (default 0.05).
#' @param eBayes_trend Logical. If `TRUE`, the analysis will account for data
#'   heteroscedasticity, which is common in DamID-seq data. (default: TRUE)
#' @param eBayes_robust Logical. If `TRUE`, the fitted trend should be robust
#'   to outliers.  Only useful when `eBayes_trend = TRUE`.  Recommended for DamID-seq
#'   data. (default: TRUE)
#' @param plot_diagnostics Logical.  If `TRUE` (default for interactive sessions),
#'   plots limma diagnostics to assess eBayes moderation, using `plot_limma_diagnostics`.
#' @param filter_occupancy NULL or integer. Filters out any locus with
#'   occupancy > `filter_threshold` in fewer than this number of samples of any single condition
#'   when set.  If set to TRUE, defaults to the minimum length of the two conditions.
#'    If FALSE or NULL, no filtering is applied. (default: TRUE)
#' @param filter_threshold Numeric.  `filter_occupancy` uses this value for thresholding
#'   the input data.  (default: 0)
#' @param filter_positive_enrichment Logical. If `TRUE` (default), regions
#'   are only considered significantly enriched if the mean score in the
#'   enriched condition is greater than zero. For example, for a region to be
#'   in `upCond1`, its logFC must be positive and its mean score in condition 1
#'   must be > 0. Set to `FALSE` to include all statistically significant changes.
#' @return A `DamIDResults` object containing the results. Access slots using
#'   accessors:
#'   \item{enrichedCond1()}{data.frame of regions enriched in condition 1}
#'   \item{enrichedCond2()}{data.frame of regions enriched in condition 2}
#'   \item{analysisTable()}{data.frame of full results for all tested regions}
#'   \item{conditionNames()}{A named character vector mapping display names to internal condition names}
#'   \item{inputData()}{The original `data_list` input}
#'
#' @examples
#' # Create a mock GRanges object for gene annotations
#' # This object, based on the package's unit tests, avoids network access.
#' mock_genes_gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle("2L", 7),
#'     ranges = IRanges::IRanges(
#'         start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'         end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'     ),
#'     strand = S4Vectors::Rle(GenomicRanges::strand(c("+", "-", "+", "+", "-", "-", "+"))),
#'     gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'     gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene")
#' )
#'
#' # Get path to sample data files included with the package
#' data_dir <- system.file("extdata", package = "damidBind")
#'
#' # Load data
#' loaded_data <- load_data_peaks(
#'     binding_profiles_path = data_dir,
#'     peaks_path = data_dir,
#'     ensdb_genes = mock_genes_gr,
#'     quantile_norm = TRUE
#' )
#'
#' # Run differential binding analysis
#' diff_results <- differential_binding(
#'     loaded_data,
#'     cond = c("L4 Neurons" = "L4", "L5 Neurons" = "L5")
#' )
#'
#' # View the results summary
#' diff_results
#'
#' @export
differential_binding <- function(
        data_list,
        cond,
        regex = FALSE,
        fdr = 0.05,
        eBayes_trend = TRUE,
        eBayes_robust = TRUE,
        plot_diagnostics = interactive(),
        filter_occupancy = TRUE,
        filter_threshold = 0,
        filter_positive_enrichment = TRUE) {

    # Prep data for analysis
    prep_results <- prep_data_for_differential_analysis(
        data_list = data_list,
        cond = cond,
        regex = regex,
        filter_occupancy = filter_occupancy,
        filter_threshold = filter_threshold
    )



    mat <- prep_results$mat
    factors <- prep_results$factors
    cond_internal <- prep_results$cond_internal
    cond_display <- prep_results$cond_display
    cond_matches <- prep_results$cond_matches
    occupancy_df <- prep_results$occupancy_df
    data_list <- prep_results$data_list # Get the updated data_list

    # Ensure 'mat' is a numeric matrix for limma
    mat <- as.matrix(mat)

    # Ensure the factor levels are explicitly set to the desired contrast order
    group <- factor(factors$condition, levels = cond_internal)
    design <- model.matrix(~ 0 + group)
    colnames(design) <- cond_internal # Use the internal condition names for design matrix

    fit <- lmFit(mat, design)

    # Make the contrast string: ensures cond[1] vs cond[2]
    contrast_str <- sprintf("%s-%s", cond_internal[1], cond_internal[2])
    message(sprintf("limma contrasts: %s", contrast_str))
    contrast_mat <- makeContrasts(
        contrasts = contrast_str,
        levels = design
    )

    fit2 <- contrasts.fit(fit, contrast_mat)

    message(sprintf("Applying eBayes moderation with `trend = %s, robust = %s`",as.character(eBayes_trend), as.character(eBayes_robust)))
    fit2 <- eBayes(fit2, trend = eBayes_trend, robust = eBayes_robust)

    result_table <- topTable(fit2, number = Inf, sort.by = "B", adjust.method = "BH")

    # Handle case where topTable could return no rows (e.g., empty input).
    if (nrow(result_table) == 0) {
        message("limma::topTable returned no results. Returning empty DamIDResults object.")
        mapping_cond <- setNames(cond_matches, cond_display)
        return(new("DamIDResults",
                   upCond1 = result_table[0, ],
                   upCond2 = result_table[0, ],
                   analysis = result_table,
                   cond = mapping_cond,
                   data = data_list
        ))
    }

    # Condition means
    mean1_name <- sprintf("%s_mean", cond_internal[1])
    mean2_name <- sprintf("%s_mean", cond_internal[2])

    mat1_samples <- rownames(factors[factors$condition == cond_internal[1], , drop = FALSE])
    mat2_samples <- rownames(factors[factors$condition == cond_internal[2], , drop = FALSE])

    # Reorder 'mat' by the rownames of result_table before calculating means
    reordered_mat <- mat[rownames(result_table), , drop = FALSE]

    result_table[[mean1_name]] <- rowMeans(reordered_mat[, mat1_samples, drop = FALSE], na.rm = TRUE)
    result_table[[mean2_name]] <- rowMeans(reordered_mat[, mat2_samples, drop = FALSE], na.rm = TRUE)

    # Gene annotation
    if ("gene_name" %in% colnames(occupancy_df) && "gene_id" %in% colnames(occupancy_df)) {
        gene_name_annot <- occupancy_df[rownames(result_table), "gene_name"]
        gene_id_annot <- occupancy_df[rownames(result_table), "gene_id"]
        result_table[, "gene_name"] <- gene_name_annot
        result_table[, "gene_id"] <- gene_id_annot
    } else {
        result_table[, "gene_name"] <- NA_character_
        result_table[, "gene_id"] <- NA_character_
    }

    result_table$minuslogp <- -log10(result_table$adj.P.Val) # Adjusted P values

    # Up/down regulation at FDR
    upCond1_all <- result_table[result_table$adj.P.Val <= fdr & result_table$logFC > 0, ]
    upCond2_all <- result_table[result_table$adj.P.Val <= fdr & result_table$logFC < 0, ]

    # Filter positive enrichment
    if (isTRUE(filter_positive_enrichment)) {
        message("\nFiltering for positive enrichment (mean score > 0 in the enriched condition).")
        upCond1 <- upCond1_all[upCond1_all[[mean1_name]] > 0, ]
        upCond2 <- upCond2_all[upCond2_all[[mean2_name]] > 0, ]

        # Report discarded results
        n_discarded1 <- nrow(upCond1_all) - nrow(upCond1)
        n_discarded2 <- nrow(upCond2_all) - nrow(upCond2)
        if (n_discarded1 > 0) {
            message(sprintf(" - Discarded %d regions enriched in '%s' due to negative mean scores.", n_discarded1, cond_display[1]))
        }
        if (n_discarded2 > 0) {
            message(sprintf(" - Discarded %d regions enriched in '%s' due to negative mean scores.", n_discarded2, cond_display[2]))
        }
    } else {
        upCond1 <- upCond1_all
        upCond2 <- upCond2_all
    }

    # Prepare mapping for output: Display Name -> Match Pattern
    mapping_cond <- setNames(cond_matches, cond_display)

    # User-friendly output summary and top genes
    ._report_results(cond_display[1], upCond1)
    ._report_results(cond_display[2], upCond2)

    # Return a formal S4 object
    results_object <- new("DamIDResults",
        upCond1 = upCond1,
        upCond2 = upCond2,
        analysis = result_table,
        cond = mapping_cond,
        data = data_list
    )

    if (isTRUE(plot_diagnostics)) {
        ._plot_limma_diagnostics_internal(fit = fit, fit2 = fit2, mat = mat)
    }

    return(results_object)
}
