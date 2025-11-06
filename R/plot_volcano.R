#' @title Check and Coerce List Input for Configuration
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


#' @title Manage Volcano Plot Configurations
#' @description Merges user-provided configurations with defaults and returns a
#'   single list of all settings.
#' @param diff_results A `DamIDResults` object.
#' @param plot_config,label_config,highlight_config,save User-provided config lists.
#' @return A list containing the final, merged configurations (`plot`, `label`,
#'   `highlight`, `save`).
#' @noRd
.manage_volcano_configs <- function(diff_results, plot_config, label_config, highlight_config, save) {
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
        max_overlaps = 20
    )
    final_label_config <- check_list_input(label_defaults, label_config)

    highlight_defaults <- list(
        alpha = 1, size = 2, label = TRUE, colour = NULL, label_size = 4,
        max_overlaps = 10, legend = TRUE, legend_inside = TRUE,
        legend_pos = list(
            legend.position = c(0.97, 0.05),
            legend.justification = c(1, 0)
        ),
        label_fill = TRUE, text_col = FALSE
    )
    final_highlight_config <- check_list_input(highlight_defaults, highlight_config)

    save_defaults <- list(
        filename = "damidBind_volcano_plot", format = "pdf",
        width = 5, height = 4
    )
    final_save_config <- check_list_input(save_defaults, save)

    list(
        plot = final_plot_config,
        label = final_label_config,
        highlight = final_highlight_config,
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


#' @title Prepare Data for Highlight and Label Layers
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

    # Process highlight groups
    if (!is.null(highlight) && length(highlight) > 0) {
        num_groups <- length(highlight)
        pal <- if (is.null(highlight_config$colour) || length(highlight_config$colour) < num_groups) {
            scales::hue_pal(l = 50)(num_groups)
        } else {
            highlight_config$colour
        }
        names(pal) <- names(highlight)

        for (i in seq_along(highlight)) {
            group_name <- names(highlight)[i]
            indices <- .find_element_indices(gene_labels_all, highlight[[i]])

            if (length(indices) > 0) {
                group_df <- plot_df[indices, , drop = FALSE]
                group_df$highlight_group_name <- group_name
                split_labels <- strsplit(gene_labels_all[indices], split = ",")
                group_df$label_to_display <- vapply(split_labels, function(vec) {
                    matched <- intersect(trimws(vec), highlight[[i]])
                    stringr::str_c(if (length(matched) == 0) vec else matched, collapse = ",")
                }, FUN.VALUE = character(1L))

                highlight_dfs[[group_name]] <- group_df
                if (isTRUE(highlight_config$label)) {
                    label_dfs[[group_name]] <- group_df
                    highlighted_ids <- c(highlighted_ids, group_df$id)
                }
            }
        }
        highlight_colours <- pal
    } else {
        highlight_colours <- NULL
    }

    # Process general labels, avoid overlap with highlight labels
    if (!is.null(label_config)) {
        sig_indices <- which(plot_df$sig)
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

    list(highlight_df = final_highlight_df, label_df = final_label_df, highlight_colours = highlight_colours)
}


#' Volcano Plot of Differentially Bound/Expressed Loci
#'
#' @description
#' Creates a volcano plot from the results of a differential analysis. The plot
#' shows the log-fold change against a measure of statistical significance.
#' The function offers extensive customisation for point appearance, gene
#' labelling, and highlighting specific groups of loci.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
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
#'     \item \code{genes}: character vector to restrict labels to a subset (default: label all significant).
#'     \item \code{label_size}: label size (numeric).
#'     \item \code{clean_names}: logical; if `TRUE`, applies regex filtering to labels.
#'     \item \code{names_clean}, \code{names_clean_extra}: regex to exclude from labels
#'       when \code{clean_names} is `TRUE`.
#'     \item \code{max_overlaps}: integer; maximum ggrepel overlaps.
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
#'     \item \code{legend}: Logical; whether to draw a plot legend for the highlight groups (default: TRUE).
#'     \item \code{legend_inside}: Logical; whether to draw the plot legend for the highlight groups inside the plot (default: TRUE).
#'     \item \code{legend_pos}: list; when legend_inside is TRUE, internal position for the legend box (default: list(legend.position = c(0.97, 0.05), legend.justification = c(1, 0)) ).
#'     \item \code{label_fill}: logical; if `TRUE`, uses `geom_label_repel`, else `geom_text_repel` (default: FALSE)
#'     \item \code{text_col}: logical; if `TRUE`, text is coloured as per points, else black (default: FALSE)
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
#'         cond = c("L4", "L5"),
#'         cond_names = c("L4 Neurons", "L5 Neurons")
#'     )
#'     return(diff_results)
#' }
#' diff_results <- .generate_example_results()
#'
#' # Generate a default volcano plot
#' plot_volcano(diff_results)
#'
#' # Generate a plot with a highlighted gene group, but no other labels
#' L4_genes_to_highlight <- c("ap", "dpr1", "side", "mav")
#' plot_volcano(
#'     diff_results,
#'     label_config = NULL,
#'     highlight = list("Key L4 Genes" = L4_genes_to_highlight)
#' )
#'
#' @export
plot_volcano <- function(
    diff_results,
    plot_config = list(),
    label_config = list(),
    highlight = NULL,
    highlight_config = list(),
    save = NULL) {
    stopifnot(is(diff_results, "DamIDResults"))

    # Prepare data and configurations
    analysis_table <- analysisTable(diff_results)
    configs <- .manage_volcano_configs(diff_results, plot_config, label_config, highlight_config, save)
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
    gene_labels_all <- if ("gene_names" %in% names(plot_df)) plot_df$gene_names else plot_df$id
    layer_data <- .prepare_highlight_and_label_data(plot_df, gene_labels_all, highlight, configs)

    # Build plot layers
    p <- ggplot(plot_df, aes(x = .data$logFC, y = .data[[plot_cfg$ystat]])) +
        geom_point(data = ~ subset(., !sig), colour = plot_cfg$nonsig_colour, alpha = plot_cfg$nonsig_alpha, size = plot_cfg$nonsig_size, shape = 20) +
        geom_point(data = ~ subset(., sig), colour = plot_cfg$sig_colour, alpha = plot_cfg$sig_alpha, size = plot_cfg$sig_size, shape = 20)

    if (!is.null(layer_data$highlight_df)) {
        p <- p + geom_point(data = layer_data$highlight_df, aes(colour = .data$highlight_group_name), alpha = configs$highlight$alpha, size = configs$highlight$size, shape = 20) +
            scale_colour_manual(name = NULL, values = layer_data$highlight_colours, guide = if (isTRUE(configs$highlight$legend)) "legend" else "none")
    }

    if (!is.null(layer_data$label_df)) {
        label_size <- if (!is.null(configs$label)) configs$label$label_size else configs$highlight$label_size
        max_overlaps <- if (!is.null(configs$label)) configs$label$max_overlaps else configs$highlight$max_overlaps

        if (isTRUE(configs$highlight$label_fill)) {
            p <- p + ggrepel::geom_label_repel(data = layer_data$label_df, aes(label = .data$label_to_display, fill = .data$highlight_group_name), size = label_size, max.overlaps = max_overlaps, min.segment.length = 0, box.padding = 0.1, point.padding = 0.1, color = "black")+
            ggplot2::scale_fill_manual(name = NULL, values = layer_data$highlight_colours, guide = "none")
        } else if (isTRUE(configs$highlight$label_col)) {
            p <- p + ggrepel::geom_text_repel(data = layer_data$label_df, aes(label = .data$label_to_display, colour = .data$highlight_group_name), size = label_size, max.overlaps = max_overlaps, min.segment.length = 0, box.padding = 0.1, point.padding = 0.1)
        } else {
            p <- p + ggrepel::geom_text_repel(data = layer_data$label_df, aes(label = .data$label_to_display), size = label_size, max.overlaps = max_overlaps, min.segment.length = 0, box.padding = 0.1, point.padding = 0.1)
        }
    }

    # Finalise labels and legend
    p <- p + labs(title = plot_cfg$title, x = plot_cfg$xlab, y = plot_cfg$ylab) + theme_bw(base_size = plot_cfg$base_size)
    if (!is.null(layer_data$highlight_df)) {
        if (isTRUE(configs$highlight$legend)) {
            if (isTRUE(configs$highlight$legend_inside)) {
                p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3), position = "inside"))
                p <- p +
                    theme(
                        legend.position.inside = configs$highlight$legend_pos$legend.position,
                        legend.justification = configs$highlight$legend_pos$legend.justification,
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
