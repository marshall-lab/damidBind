#' Differential Binding/Expression Analysis (limma)
#'
#' Setup and differential analysis for occupancy/binding experiments
#' using limma. Accepts output from `load_data.peaks` or `load_data.genes`,
#' prepares experiment matrix, fits linear models, and returns DE loci.
#'
#' @param data_list List. Output from load_data.peaks or load_data.genes.
#' @param cond Vector (character). Two strings identifying the two conditions to compare.
#'   The order matters: `cond[1]` is used as Condition 1, `cond[2]` as Condition 2.
#' @param cond_names Vector (character, optional). Custom display names for
#'   the two conditions in outputs/plots. Order maps to `cond`.
#' @param fdr Numeric. FDR threshold for significance (default 0.05).
#' @return A `DamIDResults` object containing the results. Access slots using
#'   the `@` accessor (e.g., `results@analysis`). The object includes:
#'   \item{upCond1}{data.frame of regions enriched in condition 1}
#'   \item{upCond2}{data.frame of regions enriched in condition 2}
#'   \item{analysis}{data.frame of full results for all tested regions}
#'   \item{cond}{A named character vector mapping display names to internal condition names}
#'   \item{data}{The original `data_list` input}
#' @export
differential_binding <- function(
    data_list,
    cond,
    cond_names = NULL,
    fdr = 0.05
) {

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
  fit2 <- eBayes(fit2)
  result_table <- topTable(fit2, number = Inf, sort.by = "B", adjust.method = "BH")

  # Handle case where topTable could return no rows (e.g., empty input).
  if (nrow(result_table) == 0) {
    message("limma::topTable returned no results. Returning empty results list.")
    mapping_cond <- setNames(cond_internal, cond_display)
    return(list(
      upCond1 = result_table[0, ],
      upCond2 = result_table[0, ],
      analysis = result_table,
      cond = mapping_cond,
      data_list = data_list
    ))
  }

  # Condition means
  mean1_name <- sprintf("%s_mean", gsub("-", ".", cond_internal[1]))
  mean2_name <- sprintf("%s_mean", gsub("-", ".", cond_internal[2]))

  mat1_samples <- rownames(factors[factors$condition == cond_internal[1], , drop = FALSE])
  mat2_samples <- rownames(factors[factors$condition == cond_internal[2], , drop = FALSE])

  # Reorder 'mat' by the rownames of result_table before calculating means
  reordered_mat <- mat[rownames(result_table), , drop = FALSE]

  result_table[[mean1_name]] <- rowMeans(reordered_mat[, mat1_samples, drop = FALSE])
  result_table[[mean2_name]] <- rowMeans(reordered_mat[, mat2_samples, drop = FALSE])

  # Gene annotation
  occupancy_df <- data_list$occupancy
  if ("gene_names" %in% colnames(occupancy_df) && "gene_ids" %in% colnames(occupancy_df)) {
    gene_names_annot <- occupancy_df[rownames(result_table), "gene_names"]
    gene_ids_annot   <- occupancy_df[rownames(result_table), "gene_ids"]
    result_table[,"gene_names"] <- gene_names_annot
    result_table[,"gene_ids"]   <- gene_ids_annot
  } else {
    result_table[,"gene_names"] <- NA_character_
    result_table[,"gene_ids"]   <- NA_character_
  }

  result_table$minuslogp <- -log10(result_table$adj.P.Val) # Adjusted P values

  # Up/down regulation at FDR
  upCond1 <- result_table[result_table$adj.P.Val <= fdr & result_table$logFC > 0,]
  upCond2 <- result_table[result_table$adj.P.Val <= fdr & result_table$logFC < 0,]

  # Prepare mapping for output
  mapping_cond <- setNames(cond_internal, cond_display)

  # Report based on custom names
  report_results <- function(ctype_name, loci) {
    num <- nrow(loci)
    num.adj <- sum(loci$gene_names != "" & !is.na(loci$gene_names))
    message(sprintf("\n%d loci enriched in %s", num, ctype_name))
    if (num.adj > 0) {
      top_genes <- loci$gene_names[loci$gene_names != "" & !is.na(loci$gene_names)]
      message(sprintf("Highest-ranked genes:\n%s", paste0(top_genes[seq_len(min(10, length(top_genes)))], collapse = ", ")))
    }
  }
  report_results(cond_display[1], upCond1) # Use display names
  report_results(cond_display[2], upCond2) # Use display names

  # Return a formal S4 object
  new("DamIDResults",
      upCond1 = upCond1,
      upCond2 = upCond2,
      analysis = result_table,
      cond = mapping_cond,
      data = data_list
  )
}
