#' Prepare Data for Differential Analysis
#'
#' Internal function to extract and prepare data, identify conditions, and set up
#' factors for downstream differential analysis (`limma` or `NOISeq`).
#'
#' @param data_list List. Output from load_data_peaks or load_data_genes.
#' @param cond A named or unnamed character vector of length two. The values are
#'   strings or regular expressions used to identify samples for each condition.
#'   If the vector is named, the names are used as user-friendly display names
#'   for the conditions in plots and outputs. If unnamed, the match strings are
#'   used as display names. The order determines the contrast, e.g., `cond[1]` vs `cond[2]`.
#' @param regex logical. If `TRUE`, the strings in `cond` are treated as
#'   regular expressions for matching sample names. Default is `FALSE` (fixed string matching).
#' @param filter_occupancy NULL or integer. Filters out any locus with
#'   occupancy > `filter_threshold` in fewer than this number of samples of any single condition
#'   when set.  If set to TRUE, defaults to the minimum length of the two conditions.
#'    If FALSE or NULL, no filtering is applied. (default: TRUE)
#' @param filter_threshold Numeric.  `filter_occupancy` uses this value for thresholding
#'   the input data.  (default: 0)
#' @return A list containing:
#'   \itemize{
#'     \item `mat`: The matrix of occupancy/counts with sample columns only.
#'     \item `factors`: A data frame with `condition` factors for each sample.
#'     \item `cond_internal`: The santisied, internal names for the conditions, used in the model matrix.
#'     \item `cond_display`: The final display names for the conditions.
#'     \item `cond_matches`: The original match strings/patterns for the conditions.
#'     \item `occupancy_df`: The original, potentially filtered, occupancy data frame.
#'     \item `data_list`: The input `data_list`, updated with matched sample info.
#'   }
#' @importFrom stats setNames
#' @noRd
prep_data_for_differential_analysis <- function(data_list, cond, regex = FALSE, filter_occupancy = TRUE, filter_threshold = 0) {
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
    if (is.null(cond) || length(cond) != 2 || !is.character(cond) || any(duplicated(cond))) {
        stop("`cond` must be a character vector of two unique strings.")
    }

    # If 'cond' is unnamed, use its values as names for display
    if (is.null(names(cond)) || any(names(cond) == "")) {
        message("Input 'cond' vector is unnamed. Using condition strings as display names.")
        names(cond) <- cond
    }

    cond_match_strings <- unname(cond)
    cond_display_names <- names(cond)

    # Create safe names for downstream differential analysis matrix
    cond_internal_safe <- make.names(cond_display_names, unique = TRUE)

    if (!identical(cond_display_names, cond_internal_safe)) {
        message("Condition display names were sanitized for internal data:")
        for (i in seq_along(cond_display_names)) {
            if (cond_display_names[i] != cond_internal_safe[i]) {
                message(sprintf("  '%s' -> '%s'", cond_display_names[i], cond_internal_safe[i]))
            }
        }
    }

    # Create logical masks for each condition
    match_fun <- if (isTRUE(regex)) {
        message("Using regex matching for conditions.")
        function(string, pattern) grepl(pattern, string)
    } else {
        function(string, pattern) grepl(pattern, string, fixed = TRUE)
    }
    cond1_matches <- match_fun(all_sample_cols, cond_match_strings[1])
    cond2_matches <- match_fun(all_sample_cols, cond_match_strings[2])


    # Check for overlapping samples (a sample matching both conditions)
    if (any(cond1_matches & cond2_matches)) {
        overlapping_samples <- all_sample_cols[cond1_matches & cond2_matches]
        stop(sprintf(
            "Conditions '%s' and '%s' overlap in sample names. The following samples belong to both: %s. Please ensure `cond` values uniquely identify samples.",
            cond_match_strings[1], cond_match_strings[2], paste(overlapping_samples, collapse = ", ")
        ))
    }

    # Check for unassigned samples (samples matching neither condition)
    unassigned_mask <- !(cond1_matches | cond2_matches)
    if (any(unassigned_mask)) {
        unassigned_samples <- all_sample_cols[unassigned_mask]
        warning(sprintf("The following samples could not be assigned to either condition '%s' or '%s':\n    %s\n\n", cond_match_strings[1], cond_match_strings[2], paste(unassigned_samples, collapse = "\n    ")))
    }

    # Get sample names based on masks
    samples_cond1 <- all_sample_cols[cond1_matches]
    samples_cond2 <- all_sample_cols[cond2_matches]
    sel_cols <- c(samples_cond1, samples_cond2)

    if (length(sel_cols) < 2) {
        stop(sprintf(
            "Fewer than two sample columns matched for conditions '%s' and '%s'. Please check your 'cond' values against sample names. Found %d matched samples.",
            cond_match_strings[1], cond_match_strings[2], length(sel_cols)
        ))
    }

    # Prepare the matrix for analysis
    mat <- occupancy_df[, sel_cols, drop = FALSE]
    rownames(mat) <- if ("name" %in% names(occupancy_df)) make.unique(occupancy_df$name) else rownames(occupancy_df)

    # Store matched samples in data_list for FDR filtering later
    data_list$matched_samples <- list(
        stats::setNames(list(samples_cond1), cond_display_names[1]),
        stats::setNames(list(samples_cond2), cond_display_names[2])
    )

    # Filter loci with negative occupancy if requested
    if (isTRUE(filter_occupancy)) filter_occupancy <- min(length(cond1_matches),length(cond2_matches))
    if (!is.null(filter_occupancy) && is.numeric(filter_occupancy) && filter_occupancy > 0) {
        initial_loci_count <- nrow(mat)

        filter_occupancy <- as.integer(filter_occupancy)
        if (filter_occupancy <= 0) {
            stop("'filter_occupancy' must be a positive integer.")
        }

        minimum_above_threshold <- min(
            filter_occupancy,
            length(samples_cond1),
            length(samples_cond2)
        )

        message(sprintf("Applying filter: Loci must have occupancy > %g in at least %d samples of at least one condition.", filter_threshold, minimum_above_threshold))

        keep_cond1 <- if (length(samples_cond1) >= 1) {
            rowSums(mat[, samples_cond1, drop = FALSE] > filter_threshold, na.rm=TRUE) >= minimum_above_threshold
        } else {
            rep(FALSE, nrow(mat))
        }

        keep_cond2 <- if (length(samples_cond2) >= 1) {
            rowSums(mat[, samples_cond2, drop = FALSE] > filter_threshold, na.rm=TRUE) >= minimum_above_threshold
        } else {
            rep(FALSE, nrow(mat))
        }

        rows_to_keep <- keep_cond1 | keep_cond2

        mat <- mat[rows_to_keep, , drop = FALSE]
        occupancy_df <- occupancy_df[rows_to_keep, , drop = FALSE]
        data_list$occupancy <- occupancy_df # Update the occupancy in the data_list

        filtered_count <- initial_loci_count - nrow(mat)
        if (filtered_count > 0) {
            message(sprintf(
                "  Filtered out %d loci. %d loci remain for analysis.",
                filtered_count,
                nrow(mat)
            ))
        } else {
            message("  No loci were filtered out.")
        }

        if (nrow(mat) == 0) {
            stop("No loci remained after filtering. The data may have very low occupancy or the filter may be too stringent.")
        }
    }

    # Create factors data frame
    factors <- data.frame(
        condition = factor(
            rep(cond_internal_safe, times = c(length(samples_cond1), length(samples_cond2))),
            levels = cond_internal_safe
        ),
        row.names = sel_cols
    )

    # Order factor levels
    cond_internal_ordered <- levels(factors$condition)

    message("Differential analysis setup:\n")
    message(sprintf("Condition 1 Display Name: '%s' (Internal: '%s', Match Pattern: '%s')",
                    cond_display_names[1], cond_internal_safe[1], cond_match_strings[1]))
    samples_cond1_feedback <- colnames(mat)[factors$condition == cond_internal_safe[1]]
    message(sprintf("  Found %d replicates:\n    %s\n", length(samples_cond1_feedback), paste(samples_cond1_feedback, collapse = "\n    ")))

    message(sprintf("Condition 2 Display Name: '%s' (Internal: '%s', Match Pattern: '%s')",
                    cond_display_names[2], cond_internal_safe[2], cond_match_strings[2]))
    samples_cond2_feedback <- colnames(mat)[factors$condition == cond_internal_safe[2]]
    message(sprintf("  Found %d replicates:\n    %s\n", length(samples_cond2_feedback), paste(samples_cond2_feedback, collapse = "\n    ")))

    # Check for sufficient replicates after identification
    if (length(samples_cond1_feedback) < 1 || length(samples_cond2_feedback) < 1) {
        stop("Each condition must have at least one sample/replicate. Please check 'cond' values or data.")
    }

    return(list(
        mat = mat,
        factors = factors,
        cond_internal = cond_internal_ordered,
        cond_display = cond_display_names,
        cond_matches = cond_match_strings,
        occupancy_df = occupancy_df,
        data_list = data_list
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



#' @title Drop specified input samples from data objects
#' @description An internal helper to remove samples based on user-provided names
#'   or patterns. It operates on the list structure created by `load_data_*`
#'   functions and can filter binding profiles, peak lists, and occupancy tables.
#'
#' @param data_list A list that may contain `binding_profiles_data` (a GRanges
#'   object), `peaks` (a list of GRanges), and/or `occupancy` (a data.frame).
#' @param drop_samples A character vector of sample names or patterns to remove.
#'   The pattern is matched against the sample names of binding profiles, the names
#'   of peak lists, and the sample columns of the occupancy table.
#'
#' @return An updated version of `data_list` with specified samples removed. If
#'   `drop_samples` is NULL or empty, it returns the original list.
#' @noRd
._drop_input_samples <- function(data_list, drop_samples) {
    # If no samples to drop, return immediately
    if (is.null(drop_samples) || length(drop_samples) == 0) {
        return(data_list)
    }

    # Identify sample names from all potential sources within the data_list
    binding_profile_names <- if ("binding_profiles_data" %in% names(data_list) && length(S4Vectors::mcols(data_list$binding_profiles_data)) > 0) {
        colnames(S4Vectors::mcols(data_list$binding_profiles_data))
    } else {
        character(0)
    }

    has_peaks <- "peaks" %in% names(data_list) && is.list(data_list$peaks) && length(data_list$peaks) > 0
    peak_names <- if (has_peaks) names(data_list$peaks) else NULL

    has_occupancy <- "occupancy" %in% names(data_list) && is.data.frame(data_list$occupancy)
    occupancy_sample_names <- if (has_occupancy) {
        base_cols <- c("chr", "start", "end", "name", "gene_name", "gene_id", "nfrags")
        occ_cols <- colnames(data_list$occupancy)
        occ_cols <- occ_cols[!grepl("_FDR$", occ_cols)] # Exclude FDR columns
        setdiff(occ_cols, base_cols)
    } else {
        character(0)
    }

    # For each pattern, find all matching samples across the sources
    to_drop_binding <- character()
    to_drop_peaks <- character()
    to_drop_occupancy <- character()

    for (pattern in drop_samples) {
        matched_binding <- binding_profile_names[grepl(pattern, binding_profile_names)]
        matched_peaks <- if (has_peaks) peak_names[grepl(pattern, peak_names)] else character(0)
        matched_occupancy <- if (has_occupancy) occupancy_sample_names[grepl(pattern, occupancy_sample_names)] else character(0)

        # Issue a warning if a pattern matches nothing at all
        if (length(matched_binding) == 0 && length(matched_peaks) == 0 && length(matched_occupancy) == 0) {
            warning(sprintf("The pattern '%s' in 'drop_samples' did not match any binding profiles, peak lists, or occupancy table columns.", pattern), call. = FALSE)
        }

        to_drop_binding <- c(to_drop_binding, matched_binding)
        if (has_peaks) to_drop_peaks <- c(to_drop_peaks, matched_peaks)
        if (has_occupancy) to_drop_occupancy <- c(to_drop_occupancy, matched_occupancy)
    }

    unique_binding_to_drop <- unique(to_drop_binding)
    unique_peaks_to_drop <- unique(to_drop_peaks)
    unique_occupancy_to_drop <- unique(to_drop_occupancy)

    all_dropped_samples <- unique(c(unique_binding_to_drop, unique_peaks_to_drop, unique_occupancy_to_drop))

    if (length(all_dropped_samples) > 0) {
        message("The following matched samples will be excluded:")
        message(paste0(" - ", all_dropped_samples, collapse = "\n"))

        # Drop from binding profiles GRanges metadata columns
        if (length(binding_profile_names) > 0 && length(unique_binding_to_drop) > 0) {
            cols_to_keep <- setdiff(binding_profile_names, unique_binding_to_drop)
            if (length(cols_to_keep) > 0) {
                S4Vectors::mcols(data_list$binding_profiles_data) <- S4Vectors::mcols(data_list$binding_profiles_data)[, cols_to_keep, drop = FALSE]
            } else {
                S4Vectors::mcols(data_list$binding_profiles_data) <- S4Vectors::DataFrame(row.names = seq_along(data_list$binding_profiles_data))
                warning("All binding profile samples have been dropped.", call. = FALSE)
            }
        }

        # Drop from peaks list
        if (has_peaks && length(unique_peaks_to_drop) > 0) {
            peaks_to_keep <- setdiff(peak_names, unique_peaks_to_drop)
            if (length(peaks_to_keep) > 0) {
                data_list$peaks <- data_list$peaks[peaks_to_keep]
            } else {
                data_list$peaks <- list()
                warning("All peak lists have been dropped.", call. = FALSE)
            }
        }

        # Drop from occupancy dataframe
        if (has_occupancy && length(unique_occupancy_to_drop) > 0) {
            fdr_cols_to_drop <- paste0(unique_occupancy_to_drop, "_FDR")
            all_cols_to_drop <- c(unique_occupancy_to_drop, fdr_cols_to_drop)

            cols_to_keep <- setdiff(colnames(data_list$occupancy), all_cols_to_drop)
            data_list$occupancy <- data_list$occupancy[, cols_to_keep, drop = FALSE]
        }
    } else {
        message("No samples matched the patterns in 'drop_samples'. No samples were dropped.")
    }

    return(data_list)
}
