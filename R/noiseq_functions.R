#' Differential accessibility analysis for CATaDa (`NOISeq` based)
#'
#' Setup and differential analysis for CATaDa chromatin accessibility experiments
#' using `NOISeq`. Accepts output from `load_data_peaks`, prepares a count matrix,
#' performs `NOISeq` analysis, and returns differentially-accessible loci.
#'
#' @param data_list List. Output from load_data_peaks.
#' @param cond A named or unnamed character vector of length two. The values are
#'   strings or regular expressions used to identify samples for each condition.
#'   If the vector is named, the names are used as user-friendly display names
#'   for the conditions in plots and outputs. If unnamed, the match strings are
#'   used as display names. The order determines the contrast, e.g., `cond[1]` vs `cond[2]`.
#' @param regex Logical. If `TRUE`, the strings in `cond` are treated as
#'   regular expressions for matching sample names. If `FALSE` (the default),
#'   fixed string matching is used.
#' @param q Numeric. Q-value threshold for NOISeq significance (default 0.8).
#' @param norm Normalisation method passed to NOISeq.  Defaults to "n" (no normalisation), but "uqua"
#'   (upper quantile) or "tmm" (trimmed mean of M) are options if needed
#' @return A `DamIDResults` object containing the results. Access slots using accessors (e.g., `analysisTable(results)`). The object includes:
#'   \item{upCond1}{data.frame of regions enriched in condition 1}
#'   \item{upCond2}{data.frame of regions enriched in condition 2}
#'   \item{analysis}{data.frame of full results for all tested regions}
#'   \item{cond}{A named character vector mapping display names to internal condition names}
#'   \item{data}{The original `data_list` input}
#'
#' @examples
#' # NOTE: This example uses mock counts data, as the package's sample
#' # data is in log2-ratio format.
#'
#' # Create a mock data_list with plausible count data
#' mock_occupancy_counts <- data.frame(
#'     name = c("peak1", "peak2", "peak3"),
#'     gene_name = c("GeneA", "GeneB", "GeneC"),
#'     gene_id = c("ID_A", "ID_B", "ID_C"),
#'     GroupA_rep1 = c(100, 20, 50), GroupA_rep2 = c(110, 25, 45),
#'     GroupB_rep1 = c(10, 200, 55), GroupB_rep2 = c(15, 220, 60),
#'     row.names = c("peak1", "peak2", "peak3")
#' )
#'
#' mock_data_list <- list(
#'     occupancy = mock_occupancy_counts,
#'     test_category = "accessible"
#' )
#'
#' # Run differential accessibility analysis
#' diff_access_results <- differential_accessibility(
#'     mock_data_list,
#'     cond = c("Group A" = "GroupA", "Group B" = "GroupB")
#' )
#'
#' # View the results summary
#' diff_access_results
#'
#' @export
differential_accessibility <- function(
        data_list,
        cond,
        regex = FALSE,
        norm = "n",
        q = 0.8) {
    # Prep data for analysis
    prep_results <- prep_data_for_differential_analysis(
        data_list = data_list,
        cond = cond,
        regex = regex,
        filter_occupancy = NULL # no negatives in count data
    )

    mat <- prep_results$mat
    factors <- prep_results$factors
    cond_internal <- prep_results$cond_internal
    cond_display <- prep_results$cond_display
    cond_matches <- prep_results$cond_matches
    occupancy_df <- prep_results$occupancy_df
    data_list <- prep_results$data_list

    # Ensure 'mat' is a numeric matrix
    mat <- as.matrix(mat)

    # Prepare NOISeq input
    noiseq_data <- readData(data = mat, factors = factors)

    # Using capture.output() to remove the automatic warning to user NOIseqBIO with biological replicates
    # NOIseqBIO uses assumptions relevant to RNAseq data, which are not appropriate
    # for CATaDa, so should not be used here.  The usage below is correct for CATaDa data.
    capture.output(
        {
            noiseq_res <- noiseq(
                noiseq_data,
                conditions = cond_internal,
                factor = "condition",
                norm = norm,
                replicates = "biological"
            )
        },
        file = NULL
    )

    noiseq_df <- as.data.frame(noiseq_res@results[[1]])

    # Condition means
    mean1_name <- sprintf("%s_mean", cond_internal[1])
    mean2_name <- sprintf("%s_mean", cond_internal[2])

    mat1_samples <- rownames(factors[factors$condition == cond_internal[1], , drop = FALSE])
    mat2_samples <- rownames(factors[factors$condition == cond_internal[2], , drop = FALSE])

    noiseq_df[[mean1_name]] <- rowMeans(mat[rownames(noiseq_df), mat1_samples, drop = FALSE], na.rm = TRUE)
    noiseq_df[[mean2_name]] <- rowMeans(mat[rownames(noiseq_df), mat2_samples, drop = FALSE], na.rm = TRUE)

    # Log2 ratio and -log10(p) for consistency with limma output
    noiseq_df$logFC <- noiseq_df$M
    noiseq_df$minuslogp <- -log10(1 - noiseq_df$prob)

    # Gene annotation
    if ("gene_name" %in% colnames(occupancy_df) && "gene_id" %in% colnames(occupancy_df)) {
        idx <- match(rownames(noiseq_df), rownames(occupancy_df))
        noiseq_df$gene_name <- occupancy_df$gene_name[idx]
        noiseq_df$gene_id <- occupancy_df$gene_id[idx]
    } else {
        noiseq_df$gene_name <- NA_character_
        noiseq_df$gene_id <- NA_character_
    }

    # Identify significant DE loci using degenes
    # Using capture.output to suppress redundant messaging
    capture.output(
        {
            sig_up_res <- degenes(noiseq_res, q = q, M = "up")
            sig_down_res <- degenes(noiseq_res, q = q, M = "down")
        },
        file = NULL
    )

    upCond1 <- noiseq_df[rownames(noiseq_df) %in% rownames(sig_up_res), , drop = FALSE]
    upCond2 <- noiseq_df[rownames(noiseq_df) %in% rownames(sig_down_res), , drop = FALSE]

    # Prepare mapping for output: Display Name -> Match Pattern
    mapping_cond <- setNames(cond_matches, cond_display)

    # User-friendly output summary and top genes
    ._report_results(cond_display[1], upCond1)
    ._report_results(cond_display[2], upCond2)

    # adjust test_category to reflect accessibility
    data_list$test_category <- "accessible"

    # Return a formal S4 object
    new("DamIDResults",
        upCond1    = upCond1,
        upCond2    = upCond2,
        analysis   = noiseq_df,
        cond       = mapping_cond,
        data       = data_list
    )
}
