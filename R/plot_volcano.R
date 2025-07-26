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

#' Volcano Plot of Differentially Bound/Expressed Loci
#'
#' Plots a volcano plot with optional layered highlights and overlays,
#' given the results from `differential_binding()` or `differential_analysis()`.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param plot_config List. Names to override plot details (title, axes, size,
#'   colours, etc); see details.
#'   \itemize{
#'     \item \code{title}, \code{xlab}, \code{ylab} (character)
#'     \item \code{ystat} (character): The column name from `diff_results@analysis`
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
#' @return A `ggplot` object, invisibly.
#'
#' @examples
#' # ---- Helper function to create a sample DamIDResults object ----
#' .generate_example_results <- function() {
#'   mock_genes_gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle("2L", 7),
#'     ranges = IRanges::IRanges(
#'       start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'       end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'     ),
#'     gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'     gene_name = c("ap", "dpr1", "side", "mav", "geneE", "geneF", "LargeTestGene")
#'   )
#'   data_dir <- system.file("extdata", package = "damidBind")
#'   loaded_data <- load_data_peaks(
#'     binding_profiles_path = data_dir,
#'     peaks_path = data_dir,
#'     ensdb_genes = mock_genes_gr,
#'     quantile_norm = TRUE
#'   )
#'   diff_results <- differential_binding(
#'      loaded_data,
#'      cond = c("L4", "L5"),
#'      cond_names = c("L4 Neurons", "L5 Neurons")
#'   )
#'   return(diff_results)
#' }
#' diff_results <- .generate_example_results()
#' # ---- End of helper section ----
#'
#' # Generate a default volcano plot
#' plot_volcano(diff_results)
#'
#' # Generate a plot with a highlighted gene group, but no other labels
#' L4_genes_to_highlight <- c("ap", "dpr1", "side", "mav")
#' plot_volcano(
#'   diff_results,
#'   label_config = NULL,
#'   highlight = list("Key L4 Genes" = L4_genes_to_highlight)
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

  analysis_table <- diff_results@analysis
  upCond1 <- rownames(diff_results@upCond1)
  upCond2 <- rownames(diff_results@upCond2)
  all_sig <- unique(c(upCond1, upCond2))
  plot_data <- analysis_table
  gene_labels_all <- if ("gene_names" %in% names(plot_data)) plot_data$gene_names else rownames(plot_data)

  # Get custom condition names from diff_results@cond
  cond_display <- names(diff_results@cond)
  test_category <- diff_results@data$test_category

  # Set up base plotting config as before, e.g.:
  plot_defaults <- list(
    title = sprintf("Differentially %s loci", test_category),
    xlab = bquote(log[2] * "FC (" * .(cond_display[1]) * " / " * .(cond_display[2]) * ")"),
    ylab = "",
    ystat = ifelse(test_category == "accessible", "minuslogp", "B"),
    base_size = 18,
    sig_colour = "orange",
    nonsig_colour = rgb(0.4, 0.4, 0.4),
    sig_alpha = 0.4,
    sig_size = 1,
    nonsig_alpha = 0.1,
    nonsig_size = 1
  )

  plot_config <- check_list_input(plot_defaults, plot_config)
  if (is.null(plot_config)) {
    # should not occur!
    plot_config <- plot_defaults
  }

  label_defaults <- list(
    genes = NULL,
    label_size = 3,
    clean_names = FALSE,
    names_clean = "snoRNA|snRNA|^CR|tRNA|RNA",
    names_clean_extra = NULL,
    max_overlaps = 20
  )
  label_config <- check_list_input(label_defaults, label_config)

  highlight_defaults <- list(
    alpha = 1,
    size = 2,
    label = TRUE,
    colour = NULL,
    label_size = 4,
    max_overlaps = 10
  )
  highlight_config <- check_list_input(highlight_defaults, highlight_config)


  save_defaults <- list(
    filename = "damidBind_volcano_plot",
    format = "pdf",
    width = 5,
    height = 4
  )
  save_config <- check_list_input(save_defaults, save)


  # Plot dataframe
  plot_df <- as.data.frame(plot_data, stringsAsFactors = FALSE)
  plot_df$id <- rownames(plot_data)
  plot_df$sig <- rownames(plot_data) %in% all_sig

  # Helper functions
  clean_names <- function(x, clean, clean_extra = NULL) {
    if (is.null(x)) {
      return(character(0))
    }
    regex <- if (is.null(clean_extra)) clean else paste0(clean, "|", clean_extra)
    x[!grepl(regex, x, ignore.case = FALSE)]
  }

  indices_to_clean <- function(x, clean, clean_extra = NULL) {
    if (is.null(x)) {
      return(integer(0))
    }
    regex <- if (is.null(clean_extra)) clean else paste0(clean, "|", clean_extra)
    regex <- gsub("\\|\\|", "|", x = regex) # Clean up potential double pipes
    grep(regex, x, ignore.case = FALSE)
  }

  find_element_indices <- function(labels, elements) {
    if (is.null(elements) || length(elements) == 0) {
      return(integer())
    }
    # Escape special regex characters in the element names to prevent unintended regex patterns
    escaped_elements <- gsub("([[:punct:]])", "\\\\\\1", elements)
    patt <- paste(escaped_elements, collapse = "|")
    which(
      grepl(sprintf("(^|,)(%s)(,|$)", patt), labels, perl = TRUE, ignore.case = TRUE)
    )
  }

  # Store IDs of genes already labelled by highlight groups
  highlight_labelled_gene_ids <- character(0)

  # Highlight overlays
  highlight_layers <- list()
  highlight_labels <- list()

  if (!is.null(highlight) && length(highlight) > 0) {
    # Generate default colours if not provided or insufficient
    if (is.null(highlight_config$colour) || length(highlight_config$colour) < length(highlight)) {
      highlight_default_pal <- hue_pal(l = 50)(length(highlight))
    } else {
      highlight_default_pal <- highlight_config$colour
    }


    # Create names for highlight groups for legend and messaging
    highlight_names <- names(highlight)
    if (is.null(highlight_names)) {
      highlight_names <- paste0("Highlight Group ", seq_along(highlight))
    }

    # Ensure highlight_config$colour is a list if it contains multiple colours
    if (!is.list(highlight_config$colour) && length(highlight_config$colour) > 1) {
      highlight_config$colour <- as.list(highlight_config$colour)
    }

    for (i in seq_along(highlight)) {
      group_elements <- highlight[[i]]
      indices <- find_element_indices(gene_labels_all, group_elements)

      if (length(indices) > 0) {
        layer_df <- plot_df[indices, , drop = FALSE]
        layer_df$highlight_group_name <- highlight_names[i]

        # Assign the specific labels to the layer_df for ggrepel
        # layer_df$label_for_repel <- gene_labels_all[indices]

        # Split each locus label string by commas, trim whitespace
        split_labels <- strsplit(gene_labels_all[indices], split = ",")

        # For each split label vector, intersect with the highlight group genes (group_elements)
        # and paste back together (could be more than one gene)
        matched_labels <- vapply(
          split_labels,
          function(lab_vec) {
            matched <- intersect(trimws(lab_vec), group_elements)
            if (length(matched) == 0) {
              # fallback: you could keep full label or NA; here, use full label
              paste(lab_vec, collapse = ",")
            } else {
              paste(matched, collapse = ",")
            }
          },
          FUN.VALUE = character(1L)
        )
        layer_df$label_for_repel <- matched_labels

        # Determine colour for this group
        group_colour <- if (!is.null(highlight_config$colour) && i <= length(highlight_config$colour)) {
          if (is.list(highlight_config$colour)) highlight_config$colour[[i]] else highlight_config$colour[i]
        } else {
          highlight_default_pal[i]
        }

        highlight_layers[[i]] <- geom_point(
          data = layer_df,
          aes(x = logFC, y = .data[[plot_config$ystat]], colour = highlight_group_name),
          alpha = highlight_config$alpha,
          size = highlight_config$size
        )

        if (isTRUE(highlight_config$label)) {
          message(sprintf("Highlight group '%s' will label: %s", highlight_names[i], paste(gene_labels_all[indices], collapse = ", ")))
          highlight_labels[[i]] <- geom_text_repel(
            data = layer_df,
            aes(label = label_for_repel),
            size = highlight_config$label_size,
            box.padding = 0.2,
            point.padding = 0.3,
            max.time = 5,
            max.iter = 10000,
            min.segment.length = 0.1,
            force = 0.5,
            force_pull = 3,
            max.overlaps = highlight_config$max_overlaps
          )
          # Store IDs of genes labelled by this highlight group
          highlight_labelled_gene_ids <- unique(c(highlight_labelled_gene_ids, rownames(layer_df)))
        }
      }
    }
  }

  # Y-stat logic
  if (plot_config$ystat %in% names(plot_df)) {
    if (plot_config$ylab == "") {
      plot_config$ylab <- if (plot_config$ystat == "minuslogp") "-log(p)" else plot_config$ystat
    }
  } else {
    stop(sprintf("ystat ('%s') is not a valid column in plot data. Valid columns are: %s", plot_config$ystat, paste(names(plot_df), collapse = ", ") ))
  }

  # Volcano plot construction
  p <- ggplot(plot_df, aes(x = logFC, y = .data[[plot_config$ystat]])) +
    theme_bw(base_size = plot_config$base_size) +
    geom_point(
      data = subset(plot_df, !sig),
      colour = plot_config$nonsig_colour,
      alpha = plot_config$nonsig_alpha,
      size = plot_config$nonsig_size,
      shape = 16
    ) +
    geom_point(
      data = subset(plot_df, sig),
      colour = plot_config$sig_colour,
      alpha = plot_config$sig_alpha,
      size = plot_config$sig_size,
      shape = 16
    ) +
    labs(
      title = plot_config$title,
      x = plot_config$xlab,
      y = plot_config$ylab
    )

  # Add highlight layers
  if (length(highlight_layers) > 0) {
    for (hl_layer in highlight_layers) p <- p + hl_layer

    # Add highlight colour scales for the legend
    all_highlight_colours <- character(length(highlight_names))
    names(all_highlight_colours) <- highlight_names
    for (idx in seq_along(highlight_names)) {
      # Use the colour determined in the loop above for consistency
      group_colour <- if (!is.null(highlight_config$colour) && idx <= length(highlight_config$colour)) {
        if (is.list(highlight_config$colour)) highlight_config$colour[[idx]] else highlight_config$colour[idx]
      } else {
        highlight_default_pal[idx]
      }
      all_highlight_colours[idx] <- group_colour
    }
    p <- p + scale_colour_manual(name = NULL, values = all_highlight_colours)

    # Adjust legend appearance for highlights (e.g., make points larger)
    p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))
  }

  # Add highlight labels
  if (length(highlight_labels) > 0) {
    for (hl_label in highlight_labels) p <- p + hl_label
  }

  # Labels for significant points (general)
  if (!is.null(label_config)) {
    candidate_general_label_indices <-
      if (!is.null(label_config$genes)) {
        restricted_genes_indices <- find_element_indices(gene_labels_all, label_config$genes)
        intersect(which(plot_df$sig), restricted_genes_indices)
      } else {
        which(plot_df$sig)
      }

    # Get the IDs of genes for potential general labelling
    candidate_general_labels_ids <- plot_df$id[candidate_general_label_indices]

    # Filter out any IDs that were already labelled by highlights
    filtered_general_labels_ids <- setdiff(candidate_general_labels_ids, highlight_labelled_gene_ids)

    # Convert back to indices in plot_df
    final_general_label_indices <- which(plot_df$id %in% filtered_general_labels_ids)

    label_df_general <- plot_df[final_general_label_indices, , drop = FALSE]
    label_df_general$label_to_display <- gene_labels_all[final_general_label_indices]

    # Apply cleaning filter
    if (isTRUE(label_config$clean_names) && nrow(label_df_general) > 0) {
      indices_to_remove <- indices_to_clean(
        label_df_general$label_to_display,
        label_config$names_clean,
        label_config$names_clean_extra
      )
      if (length(indices_to_remove) > 0) { # Explicitly use length
        label_df_general <- label_df_general[seq_len(nrow(label_df_general))[-indices_to_remove], , drop = FALSE]
      }
    }

    if (nrow(label_df_general) > 0) {
      p <- p +
        geom_text_repel(
          data = label_df_general,
          aes(label = label_to_display),
          size = label_config$label_size,
          box.padding = 0.2,
          point.padding = 0.3,
          max.time = 5,
          max.iter = 10000,
          min.segment.length = 0.1,
          force = 0.5,
          force_pull = 3,
          max.overlaps = label_config$max_overlaps
        )
    }
  }

  # Save to file if requested
  if (!is.null(save_config)) {
    tryCatch(
      {
        full_filename <- paste0(save_config$filename, ".", save_config$format)
        ggsave(
          filename = full_filename,
          plot = p,
          width = save_config$width,
          height = save_config$height,
          units = "in",
          device = save_config$format
        )
      },
      error = function(e) {
        message("An error occurred while saving the plot file(s): ", e$message)
      }
    )
  }

  invisible(p)
}
