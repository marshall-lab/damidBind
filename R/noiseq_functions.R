#' Differential Accessibility Analysis for CATaDa (NOISeq based)
#'
#' Setup and differential analysis for CATaDa chromatin accessibility experiments using NOISeq.
#' Accepts output from load_data_peaks, prepares count matrix,
#' performs NOISeq analysis, and returns differentially-accessible loci.
#'
#' @param data_list List. Output from load_data_peaks.
#' @param cond Vector (character). Two strings identifying the two conditions to compare.
#'   The order matters: `cond[1]` is used as Condition 1, `cond[2]` as Condition 2.
#' @param cond_names Vector (character, optional). Custom display names for
#'   the two conditions in outputs/plots. Order maps to `cond`.
#' @param q Numeric. Q-value threshold for NOISeq significance (default 0.8).
#' @param norm Normalisation method passed to NOISeq.  Defaults to "n" (no normalisation), but "uqua"
#'   (upper quantile) or "tmm" (trimmed mean of M) are options if needed
#' @return A `DamIDResults` object containing the results. Access slots using the `@`
#'   accessor (e.g., `results@analysis`). The object includes:
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
#'     gene_names = c("GeneA", "GeneB", "GeneC"),
#'     gene_ids = c("ID_A", "ID_B", "ID_C"),
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
#'     cond = c("GroupA", "GroupB")
#' )
#'
#' # View the results summary
#' diff_access_results
#'
#' @export
differential_accessibility <- function(
    data_list,
    cond,
    cond_names = NULL,
    norm = "n",
    q = 0.8) {
    # Prep data for analysis
    prep_results <- prep_data_for_differential_analysis(
        data_list = data_list,
        cond = cond,
        cond_names = cond_names
    )

    mat <- prep_results$mat
    factors <- prep_results$factors
    cond_internal <- prep_results$cond_internal
    cond_display <- prep_results$cond_display
    occupancy_df <- prep_results$occupancy_df

    # Ensure 'mat' is a numeric matrix for limma
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
    mean1_name <- sprintf("%s_mean", gsub("-", ".", cond_internal[1]))
    mean2_name <- sprintf("%s_mean", gsub("-", ".", cond_internal[2]))

    mat1_samples <- rownames(factors[factors$condition == cond_internal[1], , drop = FALSE])
    mat2_samples <- rownames(factors[factors$condition == cond_internal[2], , drop = FALSE])

    noiseq_df[[mean1_name]] <- rowMeans(mat[, mat1_samples, drop = FALSE])
    noiseq_df[[mean2_name]] <- rowMeans(mat[, mat2_samples, drop = FALSE])

    # Log2 ratio and -log10(p) for consistency with limma output
    noiseq_df$logFC <- noiseq_df$M
    noiseq_df$minuslogp <- -log10(1 - noiseq_df$prob)

    # Gene annotation
    if ("gene_names" %in% colnames(occupancy_df) && "gene_ids" %in% colnames(occupancy_df)) {
        idx <- match(rownames(noiseq_df), rownames(occupancy_df))
        noiseq_df$gene_names <- occupancy_df$gene_names[idx]
        noiseq_df$gene_ids <- occupancy_df$gene_ids[idx]
    } else {
        noiseq_df$gene_names <- NA_character_
        noiseq_df$gene_ids <- NA_character_
    }

    # Identify significant DE loci using degenes
    sig_up_res <- degenes(noiseq_res, q = q, M = "up")
    sig_down_res <- degenes(noiseq_res, q = q, M = "down")

    upCond1 <- noiseq_df[rownames(noiseq_df) %in% rownames(sig_up_res), , drop = FALSE]
    upCond2 <- noiseq_df[rownames(noiseq_df) %in% rownames(sig_down_res), , drop = FALSE]

    # Prepare mapping for output
    mapping_cond <- setNames(cond_internal, cond_display)

    # User-friendly output summary and top genes (uses display names)
    report_results <- function(ctype_name, loci) {
        num <- nrow(loci)
        num.adj <- sum(loci$gene_names != "" & !is.na(loci$gene_names)) # Only count non-empty, non-NA gene names
        message(sprintf("\n%d loci enriched in %s", num, ctype_name))
        if (num.adj > 0) {
            top_genes <- loci$gene_names[loci$gene_names != "" & !is.na(loci$gene_names)]
            message(sprintf("Highest-ranked genes:\n%s", paste0(top_genes[seq_len(min(10, length(top_genes)))], collapse = ", ")))
        }
    }
    report_results(cond_display[1], upCond1)
    report_results(cond_display[2], upCond2)

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
