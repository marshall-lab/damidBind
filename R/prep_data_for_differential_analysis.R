#' Prepare Data for Differential Analysis
#'
#' Internal function to extract and prepare data, identify conditions, and set up
#' factors for downstream differential analysis (limma or NOISeq).
#'
#' @param data_list List. Output from load_data_peaks or load_data_genes.
#' @param cond Vector (character). Two strings, each uniquely identifying one
#'   of the two conditions to compare. The order matters; cond[1] will be the
#'   first condition, cond[2] the second.
#' @param cond_names Vector (character, optional). Custom display names for the
#'   two conditions (defaults to `cond` values if not provided or NA).
#'   Order maps to `cond`.
#' @return A list containing:
#'   \itemize{
#'     \item `mat`: The matrix of occupancy/counts with sample columns only.
#'     \item `factors`: A data frame with `condition` factors for each sample.
#'     \item `cond_internal`: The actual identified condition strings from
#'           sample names, in the order they map to `cond`.
#'     \item `cond_display`: The display names for the conditions, in the order
#'           they map to `cond`.
#'     \item `occupancy_df`: The original occupancy data frame.
#'   }
#' @importFrom stringr str_detect
#' @importFrom stats setNames
#' @noRd
prep_data_for_differential_analysis <- function(data_list, cond, cond_names = NULL) {
    if (!is.list(data_list)) {
        stop("'data_list' must be the full output of load_data_peaks or load_data_genes.")
    }

    occupancy_df <- data_list$occupancy
    if (is.null(occupancy_df)) {
        stop("Could not find 'occupancy' component in input data_list.")
    }

    # Identify sample columns from the occupancy data frame
    base_cols <- c("chr", "start", "end", "name", "gene_name", "gene_id", "nfrags")
    all_sample_cols <- setdiff(colnames(occupancy_df), base_cols)

    # Explicitly remove any FDR columns from the list of potential sample columns
    all_sample_cols <- all_sample_cols[!grepl("_FDR$", all_sample_cols)]

    if (length(all_sample_cols) < 2) {
        stop("At least two sample columns are required for differential analysis.")
    }

    # Validate 'cond' input
    if (is.null(cond) || length(cond) != 2 || !is.character(cond)) {
        stop("`cond` must be a character vector of exactly two strings, each uniquely identifying one condition.")
    }

    # Create logical masks for each condition
    cond1_matches <- stringr::str_detect(all_sample_cols, fixed(cond[1]))
    cond2_matches <- stringr::str_detect(all_sample_cols, fixed(cond[2]))

    # Check for overlapping samples (a sample matching both conditions)
    if (any(cond1_matches & cond2_matches)) {
        overlapping_samples <- all_sample_cols[cond1_matches & cond2_matches]
        stop(sprintf(
            "Conditions '%s' and '%s' overlap in sample names. The following samples belong to both: %s. Please ensure `cond` values uniquely identify samples.",
            cond[1], cond[2], paste(overlapping_samples, collapse = ", ")
        ))
    }

    # Check for unassigned samples (samples matching neither condition)
    unassigned_mask <- !(cond1_matches | cond2_matches)
    if (any(unassigned_mask)) {
        unassigned_samples <- all_sample_cols[unassigned_mask]
        warning(sprintf("The following samples could not be assigned to either condition '%s' or '%s':\n    %s\n\n", cond[1], cond[2], paste(unassigned_samples, collapse = "\n    ")))
    }

    # Get sample names based on masks
    samples_cond1 <- all_sample_cols[cond1_matches]
    samples_cond2 <- all_sample_cols[cond2_matches]
    sel_cols <- c(samples_cond1, samples_cond2)

    if (length(sel_cols) < 2) {
        stop(sprintf(
            "Fewer than two sample columns matched for conditions '%s' and '%s'. Please check your 'cond' values against sample names. Found %d matched samples.",
            cond[1], cond[2], length(sel_cols)
        ))
    }

    # Prepare the matrix for analysis
    mat <- occupancy_df[, sel_cols, drop = FALSE]
    rownames(mat) <- if ("name" %in% names(occupancy_df)) make.unique(occupancy_df$name) else rownames(occupancy_df)

    # Create factors data frame
    factors <- data.frame(
        condition = factor(
            rep(cond, times = c(length(samples_cond1), length(samples_cond2))),
            levels = unique(cond)
        ),
        row.names = sel_cols
    )

    # Process cond_names for display
    if (is.null(cond_names)) {
        cond_display <- cond # Use cond values as display names if not provided
    } else if (length(cond_names) != 2 || !is.character(cond_names)) {
        warning("`cond_names` must be a character vector of exactly two strings. Defaulting to `cond` values for display names.")
        cond_display <- cond
    } else {
        # Ensure any NA in cond_names defaults to the corresponding cond
        cond_display <- ifelse(is.na(cond_names), cond, cond_names)
    }

    # Ensure the internal condition names match the order of cond
    cond_internal_sorted <- levels(factors$condition)

    message("Differential analysis setup:\n")
    message(sprintf("Condition 1: '%s' (display as '%s')", cond_internal_sorted[1], cond_display[1]))
    samples_cond1_feedback <- colnames(mat)[factors$condition == cond_internal_sorted[1]]
    message(sprintf("  Found %d replicates:\n    %s\n", length(samples_cond1_feedback), paste(samples_cond1_feedback, collapse = "\n    ")))

    message(sprintf("Condition 2: '%s' (display as '%s')", cond_internal_sorted[2], cond_display[2]))
    samples_cond2_feedback <- colnames(mat)[factors$condition == cond_internal_sorted[2]]
    message(sprintf("  Found %d replicates:\n    %s\n", length(samples_cond2_feedback), paste(samples_cond2_feedback, collapse = "\n    ")))

    # Check for sufficient replicates after identification
    if (length(samples_cond1_feedback) < 1 || length(samples_cond2_feedback) < 1) {
        stop("Each condition must have at least one sample/replicate. Please check 'cond' values or data.")
    }

    return(list(
        mat = mat,
        factors = factors,
        cond_internal = cond_internal_sorted,
        cond_display = cond_display,
        occupancy_df = occupancy_df,
        test_category = data_list$test_category
    ))
}


#' @title Report differential analysis results summary
#' @description An internal helper to report analysis results on console output
#'
#' @param ctype_name The condition name to display.
#' @param loci Differentially-enriched loci for the condition.
#'
#' @return (Returns NULL invisibly; not used)
#' @noRd
._report_results <- function(ctype_name, loci) {
    num <- nrow(loci)
    num.adj <- sum(loci$gene_name != "" & !is.na(loci$gene_name))
    message(sprintf("\n%d loci enriched in %s", num, ctype_name))
    if (num.adj > 0) {
        top_genes <- loci$gene_name[loci$gene_name != "" & !is.na(loci$gene_name)]
        message(sprintf("Highest-ranked genes:\n%s", paste0(top_genes[seq_len(min(10, length(top_genes)))], collapse = ", ")))
    }
    invisible(NULL)
}
