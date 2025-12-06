#' @title Check and coerce config inputs
#'
#' @description
#' This helper function validates and processes user-provided configuration inputs.
#' It ensures the input is a list and, if so, merges it with a set of default
#' values. It also handles non-list inputs gracefully, treating specific values
#' (e.g., `FALSE`, `0`) as `NULL` (meaning "no configuration"), and warning the
#' user if an unexpected non-list value is provided, while still merging valid
#' list inputs with defaults.
#'
#' @param default_list A `list` representing the default configuration keys and
#'   their values. This will be the base for merging or the fallback if
#'   `config_input` is invalid.
#' @param config_input The user-provided configuration input. This is expected
#'   to be a `list`, but can also be `NULL`, `FALSE`, or `0`.
#'   \itemize{
#'     \item If `config_input` is a `list`, it is merged with `default_list` using
#'       `modifyList`, and the merged list is returned.
#'     \item If `config_input` is `NULL`, `FALSE`, or `0`, the function returns `NULL`,
#'       indicating that no configuration should be applied (e.g., to turn off a feature).
#'     \item If `config_input` is any other non-list single value, a message is printed
#'       informing the user that the input was not a list, and the `default_list`
#'       is returned, allowing the plot to proceed with defaults.
#'   }
#'
#' @return A `list` (the merged configuration or `default_list` if `config_input`
#'   was a non-list invalid type) or `NULL` (if `config_input` was `NULL`, `FALSE`, or `0`).
#' @noRd
check_list_input <- function(default_list, config_input) {
    if (is.list(config_input)) {
        return(modifyList(default_list, config_input))
    } else if (is.null(config_input) || isFALSE(config_input) || identical(config_input, 0)) {
        return(NULL)
    } else {
        message(
            "Input value was ", dQuote(as.character(config_input)),
            ", which is not a list. Using default settings for this option."
        )
        return(default_list)
    }
}

#' Sample data points based on local isolation
#'
#' @description
#' An issue with labelling points on dense plots (e.g., volcano plots) is that
#' high point density prevents clear labelling, even with tools like `ggrepel`.
#' This function addresses this by retaining isolated points while sampling from
#' points in higher-density regions. It takes a dataframe with Cartesian coordinates
#' and returns a logical vector identifying which points to select for labelling.
#' The result is a less cluttered plot where labels are present even in crowded areas,
#' providing a better representation of the underlying data.
#'
#' @details
#' Algorithm in detail:
#' 1.  If `scale = TRUE`, the coordinate data is centred and scaled.
#' 2.  An approximate k-nearest neighbour (KNN) search for all points is conducted
#'     using the HNSW algorithm.
#' 3.  A priority score is calculated for each point, defined as the median distance
#'     to its `k` nearest neighbours, where higher scores signify greater isolation.
#' 4.  The function iterates through the sorted list of points in descending order:
#'     a. If a point has not yet been processed, it is marked to be 'kept'.
#'     b. All neighbours of this point within the specified exclusion radius `r`
#'        are then marked as 'processed' and will be ignored in subsequent iterations.
#' 5.  A logical vector is returned, where `TRUE` corresponds to a point
#'     that should be kept for labelling.
#'
#' @param df A dataframe containing the point coordinates.
#' @param x_col A character string with the name of the column containing x-coordinates.
#' @param y_col A character string with the name of the column containing y-coordinates.
#' @param scale A logical value. If `TRUE`, the coordinate data is centred and scaled
#'   (using `scale(center=TRUE, scale=TRUE)`) before distance calculations.
#'   Defaults: `TRUE`.
#' @param r The exclusion radius. This can be a positive numeric value or the
#'   string "auto". If "auto", the radius is calculated as the median distance to
#'   the `k_for_r`-th nearest neighbour across all points. A smaller `r` will
#'   result in more points being kept. Note: The interpretation of `r` depends
#'   on whether `scale` is `TRUE`.
#' @param k_priority An integer for calculating the isolation priority score.
#'   Must be less than or equal to `k_search`. Default: 50.
#' @param k_search The maximum number of neighbours to find in the initial KNN search.
#'   This value must be greater than or equal to both `k` and `k_for_r`.
#' @param k_for_r An integer specifying which neighbour to use for the 'auto' `r`
#'   calculation. Default: 5.
#'
#' @return A logical vector of length `nrow(df)`. `TRUE` indicates the point
#'   at that index should be kept for labelling.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(ggrepel)
#'
#' # Generate sample data with a dense cluster
#' set.seed(42)
#' n_points <- 1000
#' cluster_data <- data.frame(
#'   x = rnorm(n_points, mean = 5, sd = 1),
#'   y = rnorm(n_points, mean = 5, sd = 1),
#'   label = paste("Point", 1:n_points)
#' )
#'
#' # Use the function to get a logical vector for filtering
#' kept_labels <- sample_labels_by_isolation(
#'   df = cluster_data,
#'   x_col = "x",
#'   y_col = "y",
#'   scale = FALSE,
#'   r = "auto",
#'   k = 5,
#'   k_search = 100,
#'   k_for_r = 10
#' )
#'
#' # Create the label dataframe for ggplot
#' label_df <- cluster_data[kept_labels, ]
#'
#' # Plot the results
#' ggplot(cluster_data, aes(x = x, y = y)) +
#'   geom_point(colour = "grey70", alpha = 0.7) +
#'   geom_point(data = label_df, colour = "firebrick") +
#'   geom_text_repel(
#'     data = label_df,
#'     aes(label = label),
#'     min.segment.length = 0,
#'     box.padding = 0.25,
#'     max.overlaps = Inf
#'   ) +
#'   coord_fixed() +
#'   labs(
#'     title = "Sampled Labels",
#'     subtitle = paste(sum(kept_labels), "of", nrow(cluster_data), "points labelled"),
#'     caption = "Red points are selected for labelling."
#'   ) +
#'   theme_bw()
sample_labels_by_isolation <- function(df, x_col, y_col, r, k_priority = 30, scale = TRUE, k_search = 30, k_for_r = 5) {

    # Input validation
    if (!is.data.frame(df)) stop("'df' must be a data frame.", call. = FALSE)
    if (!all(c(x_col, y_col) %in% names(df))) {
        stop("'", x_col, "' and '", y_col, "' must be valid column names in 'df'.", call. = FALSE)
    }
    if (!is.logical(scale) || length(scale) != 1) {
        stop("'scale' must be a logical value (TRUE or FALSE).", call. = FALSE)
    }

    k_priority <- as.integer(k_priority)
    k_search <- as.integer(k_search)
    k_for_r <- as.integer(k_for_r)

    if (is.na(k_priority) || k_priority <= 0 ) {
        stop("'k_priority' must be a positive integer.", call. = FALSE)
    }
    if (is.na(k_search) || k_search <= 0 ) {
        stop("'k_search' must be a positive integer.", call. = FALSE)
    }
    if (is.na(k_for_r) || k_for_r <= 0 ) {
        stop("'k_for_r' must be a positive integer.", call. = FALSE)
    }

    # param validation
    n_points <- nrow(df)

    # Sanity check
    if (2*k_search > n_points) {
        return(rep(TRUE, n_points))
        message("  Label sampling: too few points to effectively sample; all will be labelled.")
    }
    if (k_priority > k_search) {
        k_priority <- k_search
        message("  Label sampling: 'k_priority' was too large, is set to 'k_search'.")
    }
    if (k_for_r > k_priority) {
        k_for_r <- k_priority
        message("  Label sampling: 'k_for_r' was too large, is set to 'k'.")
    }


    coords <- as.matrix(df[, c(x_col, y_col)])

    # Coordinate scaling
    if (isTRUE(scale)) {
        coords <- scale(coords, center = TRUE, scale = TRUE)
        message("  Label sampling: coordinates have been centred and scaled.")
    }

    # KNN search
    nn_data <- RcppHNSW::hnsw_knn(X = coords, k = k_search)

    # Exclusion radius (r)
    if (identical(r, "auto")) {
        if (k_for_r > k_search) {
            stop("'k_for_r' cannot be greater than 'k_search'.", call. = FALSE)
        }
        # Use the k_for_r-th neighbour's distance for auto-calculation
        distances_to_kth_nn <- nn_data$dist[, k_for_r]
        r_calculated <- median(distances_to_kth_nn)
        message(sprintf(
            "  Label sampling: auto-calculated radius 'r' (median distance to %d-th nearest neighbour): %0.3f",
            k_for_r, r_calculated)
        )
        r <- r_calculated
    } else if (!is.numeric(r) || r <= 0) {
        stop("'r' must be a positive numeric value or the string 'auto'.", call. = FALSE)
    }

    # Priority rank score (median distance to first 'k' neighbours)
    k_priority <- min(k_priority, ncol(nn_data$dist))
    priority_scores <- apply(nn_data$dist[, seq_len(k_priority), drop = FALSE], 1, median)

    # Get the order of points to process: from most isolated to least.
    sorted_indices <- order(priority_scores, decreasing = TRUE)

    processed <- rep(FALSE, n_points)
    kept <- rep(FALSE, n_points)

    for (i in sorted_indices) {
        if (processed[i]) {
            next
        }

        kept[i] <- TRUE
        processed[i] <- TRUE

        neighbour_indices <- nn_data$idx[i, ]
        neighbour_distances <- nn_data$dist[i, ]

        indices_to_discard <- neighbour_indices[neighbour_distances < r & neighbour_distances > 0] # Exclude self-distance if 0
        processed[indices_to_discard] <- TRUE
    }

    return(kept)
}


#' @title Manage volcano plot configurations
#' @description Merges user-provided configurations with defaults and returns a
#'   single list of all settings.
#' @param diff_results A `DamIDResults` object.
#' @param plot_config,label_config,highlight_config,save User-provided config lists.
#' @return A list containing the final, merged configurations (`plot`, `label`,
#'   `highlight`, `save`).
#' @noRd
.manage_volcano_configs <- function(diff_results, plot_config, label_config, highlight_config, label_display, save) {
    cond_display <- names(conditionNames(diff_results))
    test_category <- inputData(diff_results)$test_category

    plot_defaults <- list(
        title = sprintf("Differentially %s loci", test_category),
        xlab = bquote(log[2] * "FC (" * .(cond_display[1]) * " / " * .(cond_display[2]) * ")"),
        ylab = "",
        ystat = ifelse(test_category == "accessible", "minuslogp", "B"),
        base_size = 18,
        sig_colour = "orange",
        nonsig_colour = rgb(0.4, 0.4, 0.4),
        sig_alpha = 0.4, sig_size = 1,
        nonsig_alpha = 0.1, nonsig_size = 1
    )
    final_plot_config <- check_list_input(plot_defaults, plot_config)
    if (final_plot_config$ylab == "") {
        final_plot_config$ylab <- if (final_plot_config$ystat == "minuslogp") "-log(p)" else final_plot_config$ystat
    }

    label_defaults <- list(
        genes = NULL, label_size = 3, clean_names = FALSE,
        names_clean = "snoRNA|snRNA|^CR|tRNA|RNA", names_clean_extra = NULL,
        max_overlaps = 10
    )
    final_label_config <- check_list_input(label_defaults, label_config)

    highlight_defaults <- list(
        alpha = 1, size = 2, label = TRUE, colour = NULL, label_size = 4,
        max_overlaps = 10,
        legend = TRUE,
        legend_inside = TRUE,
        legend_inside_pos = 'r',
        legend_position_override = NULL,
        legend_justification_override = NULL,
        sig_labels_only = FALSE,
        label_fill = FALSE, text_col = FALSE,
        text_luminosity = 0
    )
    final_highlight_config <- check_list_input(highlight_defaults, highlight_config)

    label_display_defaults <- list(
        k_priority = 30,
        r = 0.2,
        k_for_r = 10,
        k_search = 30,
        scale = TRUE
    )
    final_label_display <- check_list_input(label_display_defaults, label_display)

    save_defaults <- list(
        filename = "damidBind_volcano_plot", format = "pdf",
        width = 5, height = 4
    )
    final_save_config <- check_list_input(save_defaults, save)

    list(
        plot = final_plot_config,
        label = final_label_config,
        highlight = final_highlight_config,
        label_display = final_label_display,
        save = final_save_config
    )
}


#' @title Find indices of elements within comma-separated label strings
#' @description Finds which `labels` (character vector of comma-separated strings)
#' contain any of the given `elements`. It performs a whole-word match.
#' @param labels Character vector of labels, e.g., `c("geneA,geneB", "geneC")`.
#' @param elements Character vector of elements to search for, e.g., `c("geneA", "geneD")`.
#' @return Integer vector of indices in `labels` that contain at least one of the `elements`.
#' @noRd
.find_element_indices <- function(labels, elements) {
    if (is.null(elements) || length(elements) == 0) {
        return(integer(0))
    }
    escaped_elements <- gsub("([[:punct:]])", "\\\\\\1", elements)
    patt <- stringr::str_c(escaped_elements, collapse = "|")
    which(stringr::str_detect(labels, sprintf("(^|,)(%s)(,|$)", patt)))
}


#' @title Prepare data for highlight and label layers
#' @description Generates dataframes for ggplot highlight and label layers.
#' @return A list with `highlight_df` (for geom_point), `label_df` (for ggrepel),
#'   and `highlight_colours` (named vector).
#' @noRd
.prepare_highlight_and_label_data <- function(plot_df, gene_labels_all, highlight, configs) {
    highlight_dfs <- list()
    label_dfs <- list()
    highlighted_ids <- character(0)
    highlight_config <- configs$highlight
    label_config <- configs$label
    sig_indices <- which(plot_df$sig)

    highlight_colours <- NULL
    highlight_text_colours <- NULL

    if (!is.null(highlight) && length(highlight) > 0) {
        num_groups <- length(highlight)
        pal <- if (is.null(highlight_config$colour) || length(highlight_config$colour) < num_groups) {
            scales::hue_pal(l = 50)(num_groups)
        } else {
            highlight_config$colour
        }
        names(pal) <- names(highlight)
        highlight_colours <- pal

        if (is.numeric(highlight_config$text_luminosity) && highlight_config$text_luminosity > 0) {
            luminosity_reduction <- highlight_config$text_luminosity / 100
            highlight_text_colours <- colorspace::darken(highlight_colours, amount = luminosity_reduction)
            names(highlight_text_colours) <- names(highlight_colours)
        } else {
            highlight_text_colours <- highlight_colours
        }

        for (i in seq_along(highlight)) {
            group_name <- names(highlight)[i]
            indices <- .find_element_indices(gene_labels_all, highlight[[i]])

            indices_remove_labels <- NULL
            if (isTRUE(highlight_config$sig_labels_only)) {
                indices_remove_labels <- setdiff(indices, sig_indices)
            }

            if (length(indices) > 0) {
                group_df <- plot_df
                group_df$highlight_group_name <- group_name
                split_labels <- strsplit(gene_labels_all, split = ",")
                group_df$label_to_display <- vapply(split_labels, function(vec) {
                    matched <- intersect(trimws(vec), highlight[[i]])
                    stringr::str_c(if (length(matched) == 0) vec else matched, collapse = ",")
                }, FUN.VALUE = character(1L))

                # Remove unwanted label
                group_df$label_to_display[indices_remove_labels] <- ""

                # Now filter on highlight indices
                group_df <- group_df[indices, , drop = FALSE]
                highlight_dfs[[group_name]] <- group_df
                if (isTRUE(highlight_config$label)) {
                    label_dfs[[group_name]] <- group_df
                    highlighted_ids <- c(highlighted_ids, group_df$id)
                }
            }
        }
    }

    # Process general labels, avoid overlap with highlight labels
    if (!is.null(label_config)) {
        candidate_indices <- if (!is.null(label_config$genes)) {
            intersect(sig_indices, .find_element_indices(gene_labels_all, label_config$genes))
        } else {
            sig_indices
        }

        final_indices <- setdiff(candidate_indices, which(plot_df$id %in% unique(highlighted_ids)))
        if (length(final_indices) > 0) {
            general_label_df <- plot_df[final_indices, , drop = FALSE]
            general_label_df$label_to_display <- gene_labels_all[final_indices]
            general_label_df$highlight_group_name <- NA_character_

            if (isTRUE(label_config$clean_names)) {
                clean_regex <- stringr::str_c(label_config$names_clean, label_config$names_clean_extra, sep = "|")
                keep_indices <- which(!stringr::str_detect(general_label_df$label_to_display, clean_regex))
                general_label_df <- general_label_df[keep_indices, , drop = FALSE]
            }
            if (nrow(general_label_df) > 0) label_dfs[["general"]] <- general_label_df
        }
    }

    # Combine dataframes and return
    final_highlight_df <- if (length(highlight_dfs) > 0) do.call(rbind, highlight_dfs) else NULL
    final_label_df <- if (length(label_dfs) > 0) do.call(rbind, label_dfs) else NULL

    list(
        highlight_df = final_highlight_df,
        label_df = final_label_df,
        highlight_colours = highlight_colours,
        highlight_text_colours = highlight_text_colours
    )
}


#' Volcano plot of differentially bound/expressed loci
#'
#' @description
#' Creates a volcano plot from the results of a differential analysis. The plot
#' shows the log-fold change against a measure of statistical significance.
#' The function offers extensive customisation for point appearance, gene
#' labelling, and highlighting specific groups of loci.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param fdr_filter_threshold Numeric or NULL. If a value (e.g., 0.05) is provided,
#'   the volcano plot will only include loci that have an FDR value less than or
#'   equal to this threshold in at least one replicate of the two conditions
#'   being plotted. This requires that the data was loaded using
#'   `load_data_genes` with `calculate_fdr = TRUE`, which generates the
#'   necessary `_FDR` columns. If `NULL` (default), no FDR-based filtering is
#'   performed.
#' @param plot_config List. Names to override plot details (title, axes, size,
#'   colours, etc); see details.
#'   \itemize{
#'     \item \code{title}, \code{xlab}, \code{ylab} (character)
#'     \item \code{ystat} (character): The column name from `analysisTable(diff_results)`
#'       to use for the y-axis (e.g., "minuslogp" or "B"). Default is "B".
#'     \item \code{base_size} (integer): ggplot theme base font size.
#'     \item \code{sig_colour}, \code{nonsig_colour} (colours)
#'     \item \code{sig_alpha}, \code{sig_size}: alpha and size for significant points.
#'     \item \code{nonsig_alpha}, \code{nonsig_size}: alpha and size for non-significant points.
#'   }
#' @param label_config List. Fine-grained label controls; if missing or `NULL`,
#'   no labels are added (see details).
#'   \itemize{
#'     \item \code{genes}: character vector to restrict labels to a subset
#'       (default: label all significant).
#'     \item \code{label_size}: label size (numeric).
#'     \item \code{clean_names}: logical; if `TRUE`, applies regex filtering to labels.
#'     \item \code{names_clean}, \code{names_clean_extra}: regex to exclude from labels
#'       when \code{clean_names} is `TRUE`.
#'     \item \code{max_overlaps}: integer; maximum ggrepel overlaps. (default: 10)
#'   }
#' @param highlight List. A simple list where each element is a character vector
#'   of genes/loci to highlight. Each element of this list will correspond to a
#'   separate highlight group. If `NULL`, no highlight overlays are drawn.
#' @param highlight_config List. Additional highlight configuration options,
#'   applied consistently across all highlight groups. If missing or `NULL`,
#'   defaults are used.
#'   \itemize{
#'     \item \code{alpha}: Numeric; transparency for highlight points (default: 1).
#'     \item \code{size}: Numeric; size for highlight points (default: 2).
#'     \item \code{label}: Logical; if `TRUE`, labels are added for all highlight groups (default: `FALSE`).
#'     \item \code{colour}: A list of colours, where each element corresponds to a highlight
#'       group in the `highlight` list. If not specified or not enough colours are
#'       provided, a default hue palette is used.
#'     \item \code{label_size}: Numeric; label size (default: 4).
#'     \item \code{max_overlaps}: Integer; maximum ggrepel overlaps for highlight labels (default: 10).
#'     \item \code{sig_labels_only}: Logical; whether to only label significant loci in the set
#'     \item \code{legend}: Logical; whether to draw a plot legend for the highlight groups (default: TRUE).
#'     \item \code{legend_inside}: Logical; whether to draw the plot legend for
#'       the highlight groups inside the plot (default: TRUE).
#'     \item \code{legend_inside_pos}: String, either 'r' (right) or 'l' (left).
#'       Presets for internal legend position in the bottom right or left corners
#'       of the plot.  (default: 'r')
#'     \item \code{legend_position_override}: Numeric. Manual override for internal
#'       legend positioning when not set to the default, NULL.
#'     \item \code{legend_justification_override}: Numeric. Manual override for internal
#'       legend justification when not set to the default, NULL.
#'     \item \code{label_fill}: logical; if `TRUE`, uses `geom_label_repel`, else `geom_text_repel` (default: FALSE)
#'     \item \code{text_col}: logical; if `TRUE`, text is coloured as per points, else black (default: FALSE)
#'     \item \code{text_luminosity}: Numeric (0-100).  When using `text_col`, setting
#'       a non-zero value will darken the luminosity of the highlight colour on text
#'       labels for increased contrast.  0 = no change; 100 = black.  (default: 0)
#'   }
#' @param label_display List. Additional label display options for sampling dense
#'   labels in all groups.  Uses KNN-based sampling to optimise display when not
#'   NULL.
#'   \itemize{
#'     \item \code{scale}: Logical; if \code{TRUE}, labelled coordinate data are centred
#'       and scaled (using \code{scale(center = TRUE, scale = TRUE)}) before
#'       sampling.  Note: this does not affect plotted values. (default: \code{TRUE}).
#'     \item \code{r}: Numeric or \code{"auto"}; the sampling exclusion radius. If \code{"auto"},
#'        \code{r} is set to the median distance to the \code{k_for_r}-th
#'        nearest neighbour across all points. A smaller \code{r} keeps more
#'        points. (Default: 0.2).
#'     \item \code{k_search}: Integer; maximum number of neighbours to find in the initial KNN search. Must be greater than or equal to both \code{k} and \code{k_for_r} (default: 30).
#'     \item \code{k_priority}: Integer; number of neighbours used to infer the
#'       isolation priority score. Must be less than or equal to \code{k_search} (default: 30).
#'     \item \code{k_for_r}: Integer; which neighbour to use for the
#'       \code{"auto"} \code{r} calculation (default: 30).
#'   }
#' @param save List or `NULL`. Controls saving the plot to a file.
#'   If `NULL`, `FALSE`, or `0`, the plot is not saved.
#'   If a `list`, it specifies saving parameters:
#'   \itemize{
#'     \item \code{filename} (character): The path and base name for the output file
#'       (e.g., "my_volcano_plot"). If not specified, a default is used.
#'     \item \code{format} (character): File format ("pdf", "svg", or "png").
#'       Default is "pdf".
#'     \item \code{width} (numeric): Width of the plot in inches. Default is 5.
#'     \item \code{height} (numeric): Height of the plot in inches. Default is 4.
#'   }
#'
#' @return A `ggplot` object
#'
#' @examples
#' # Helper function to create a sample DamIDResults object
#' .generate_example_results <- function() {
#'     mock_genes_gr <- GenomicRanges::GRanges(
#'         seqnames = S4Vectors::Rle("2L", 7),
#'         ranges = IRanges::IRanges(
#'             start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'             end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'         ),
#'         gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'         gene_name = c("ap", "dpr1", "side", "mav", "geneE", "geneF", "LargeTestGene")
#'     )
#'     data_dir <- system.file("extdata", package = "damidBind")
#'     loaded_data <- load_data_peaks(
#'         binding_profiles_path = data_dir,
#'         peaks_path = data_dir,
#'         ensdb_genes = mock_genes_gr,
#'         quantile_norm = TRUE
#'     )
#'     diff_results <- differential_binding(
#'         loaded_data,
#'         cond = c("L4 Neurons" = "L4",
#'                  "L5 Neurons" = "L5")
#'     )
#'     return(diff_results)
#' }
#' diff_results <- .generate_example_results()
#'
#' # Generate a default volcano plot
#' plot_volcano(diff_results)
#'
#' @export
plot_volcano <- function(
        diff_results,
        fdr_filter_threshold = NULL,
        plot_config = list(),
        label_config = list(),
        highlight = NULL,
        highlight_config = list(),
        label_display = list(),
        save = NULL) {
    stopifnot(is(diff_results, "DamIDResults"))

    # Prepare data and configurations
    analysis_table <- analysisTable(diff_results)

    # Filter loci by gene expression FDR, if requested
    if (!is.null(fdr_filter_threshold) && is.numeric(fdr_filter_threshold)) {
        kept_loci <- filter_on_fdr(diff_results, fdr_filter_threshold)
        analysis_table <- analysis_table[kept_loci, , drop = FALSE]
        if (nrow(analysis_table) == 0) {
            warning("After FDR filtering, no loci remain to be plotted. The plot will be empty.")
        }
    }

    configs <- .manage_volcano_configs(diff_results, plot_config, label_config,
                                       highlight_config, label_display, save)
    plot_cfg <- configs$plot
    if (!plot_cfg$ystat %in% names(analysis_table)) {
        stop(sprintf(
            "ystat ('%s') is not a valid column name. Valid columns are: %s",
            plot_cfg$ystat, stringr::str_c(names(analysis_table), collapse = ", ")
        ))
    }

    plot_df <- analysis_table
    plot_df$id <- rownames(plot_df)
    plot_df$sig <- plot_df$id %in% c(rownames(enrichedCond1(diff_results)), rownames(enrichedCond2(diff_results)))
    gene_labels_all <- if ("gene_name" %in% names(plot_df)) plot_df$gene_name else plot_df$id
    layer_data <- .prepare_highlight_and_label_data(plot_df, gene_labels_all, highlight, configs)

    if (!is.null(layer_data$label_df) && nrow(layer_data$label_df) > 0 && !is.null(configs$label_display)) {
        kept_mask <- sample_labels_by_isolation(
            df = layer_data$label_df,
            x_col = "logFC",
            y_col = plot_cfg$ystat,
            r = configs$label_display$r,
            k_priority = configs$label_display$k_priority,
            k_for_r = configs$label_display$k_for_r,
            k_search = configs$label_display$k_search
        )
        layer_data$label_df <- layer_data$label_df[kept_mask, ]
    }

    # Build plot layers
    p <- ggplot(plot_df, aes(x = .data$logFC, y = .data[[plot_cfg$ystat]])) +
        geom_point(data = ~ subset(., !sig), colour = plot_cfg$nonsig_colour, alpha = plot_cfg$nonsig_alpha, size = plot_cfg$nonsig_size, shape = 20) +
        geom_point(data = ~ subset(., sig), colour = plot_cfg$sig_colour, alpha = plot_cfg$sig_alpha, size = plot_cfg$sig_size, shape = 20)

    if (!is.null(layer_data$highlight_df)) {
        p <- p +
            geom_point(data = layer_data$highlight_df, aes(colour = .data$highlight_group_name), alpha = configs$highlight$alpha, size = configs$highlight$size, shape = 20) +
            scale_colour_manual(
                name = NULL,
                values = layer_data$highlight_colours,
                guide = if (isTRUE(configs$highlight$legend)) "legend" else "none"
            )
    }

    if (
        !is.null(layer_data$highlight_df) &&
        isTRUE(configs$highlight$legend) &&
        isTRUE(configs$highlight$legend_inside)
        ) {
        # This awkward kludge exists because of issues with ggnewscale
        p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3), position = "inside"))
    }

    # Add labels
    if (!is.null(layer_data$label_df) && nrow(layer_data$label_df) > 0) {
        label_size <- if (!is.null(configs$label)) configs$label$label_size else configs$highlight$label_size
        max_overlaps <- if (!is.null(configs$label)) configs$label$max_overlaps else configs$highlight$max_overlaps

        if (isTRUE(configs$highlight$label_fill)) {
            p <- p + ggrepel::geom_label_repel(data = layer_data$label_df, aes(label = .data$label_to_display, fill = .data$highlight_group_name), size = label_size, max.overlaps = max_overlaps, min.segment.length = 0, box.padding = 0.1, point.padding = 0.1, color = "black") +
                ggplot2::scale_fill_manual(name = NULL, values = layer_data$highlight_colours, guide = "none")
        } else if (isTRUE(configs$highlight$text_col)) {
            # Using a separate colour scale for text to allow darker text
            p <- p +
                ggnewscale::new_scale_colour() +
                ggrepel::geom_text_repel(
                    data = layer_data$label_df,
                    aes(label = .data$label_to_display, colour = .data$highlight_group_name),
                    size = label_size, max.overlaps = max_overlaps,
                    min.segment.length = 0, box.padding = 0.1, point.padding = 0.1
                ) +
                scale_colour_manual(
                    name = NULL,
                    values = layer_data$highlight_text_colours,
                    guide = "none"
                )
        } else {
            p <- p + ggrepel::geom_text_repel(data = layer_data$label_df, aes(label = .data$label_to_display), size = label_size, max.overlaps = max_overlaps, min.segment.length = 0, box.padding = 0.1, point.padding = 0.1)
        }
    }

    # Finalise labels and legend
    p <- p + labs(title = plot_cfg$title, x = plot_cfg$xlab, y = plot_cfg$ylab) + theme_bw(base_size = plot_cfg$base_size)
    if (!is.null(layer_data$highlight_df)) {
        if (isTRUE(configs$highlight$legend)) {
            if (isTRUE(configs$highlight$legend_inside)) {

                if (configs$highlight$legend_inside_pos == 'r') {
                    legend.position <- c(0.97, 0.05)
                    legend.justification <- c(1, 0)
                } else if (configs$highlight$legend_inside_pos == 'l') {
                    legend.position <- c(0.03, 0.05)
                    legend.justification <- c(0, 0)
                }

                if (!is.null(configs$highlight$legend_position_override)) {
                    legend.position <- configs$highlight$legend_position_override
                }
                if (!is.null(configs$highlight$legend_justification_override)) {
                    legend.justification <- configs$highlight$legend_justification_override
                }

                p <- p +
                    theme(
                        legend.position.inside = legend.position,
                        legend.justification = legend.justification,
                        legend.background = ggplot2::element_rect(fill = alpha("white", 0.7), colour = "grey20", linewidth = 0.3),
                        legend.key = ggplot2::element_rect(fill = "transparent", colour = NA)
                    )
            } else {
                p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))
            }
        }
    }

    # Save plot if requested
    if (!is.null(configs$save)) {
        sc <- configs$save
        fn <- sprintf("%s.%s", sc$filename, sc$format)
        tryCatch(
            {
                ggsave(fn, p, width = sc$width, height = sc$height, units = "in", device = sc$format)
                message(sprintf("Volcano plot saved to: %s", fn))
            },
            error = function(e) warning(sprintf("Failed to save volcano plot: %s", conditionMessage(e)))
        )
    }
    p
}


