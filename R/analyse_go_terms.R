#' Perform Gene Ontology (GO) Enrichment Analysis for Differentially Bound/Expressed Regions
#'
#' This function performs Gene Ontology (GO) enrichment analysis using `clusterProfiler`
#' for either the upCond1ulated or upCond2ulated regions/genes identified by
#' `differential_binding()` or `differential_accessibility()`. It automatically extracts the relevant Flybase IDs (FBgnIDs)
#' and the background universe from the input `DamIDResults` object.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param direction Character string. Specifies which set of genes to analyse, either using condition names,
#'   "cond1" or "cond2", or "all" (for all significantly enriched genes from either direction).
#'   Default is "cond1".
#' @param org_db An OrgDb object specifying the organism's annotation database.
#'   For Drosophila, use `org.Dm.eg.db::org.Dm.eg.db`.
#' @param ontology Character string. The GO ontology to use: "BP" (Biological Process),
#'   "MF" (Molecular Function), or "CC" (Cellular Component). Default is "BP".
#' @param pvalue_cutoff Numeric. Adjusted p-value cutoff for significance. Default: 0.05.
#' @param qvalue_cutoff Numeric. Q-value cutoff for significance. Default: 0.2.
#' @param plot_title Character string. Title for the generated dot plot.
#' @param show_category Integer. Number of top enriched GO categories to display in the plot. Default: 12.
#' @param label_format_width Integer. Max character length for GO term labels on the plot. Default: 30.
#' @param save List or `NULL`. Controls saving the plot to a file (dot plot).
#'   If `NULL`, `FALSE`, or `0`, the plot is not saved.
#'   If a `list`, it specifies saving parameters:
#'   \itemize{
#'     \item \code{filename} (character): The path and base name for the output file. If not specified, the default name "damidBind_GSEA_dotplot" is used.
#'     \item \code{format} (character): File format ("pdf", "svg", or "png"). Default is "pdf".
#'     \item \code{width} (numeric): Width of the plot in inches. Default is 6.
#'     \item \code{height} (numeric): Height of the plot in inches. Default is 6.
#'   }
#' @param save_results_path Character string or NULL. If a path is provided (e.g., "go_results.csv"), the
#'   enrichment results table will be saved to this CSV file.
#' @param maxGSSize Integer. Maximum size of gene sets to consider. Default: 1000.
#' @param minGSSize Integer. Minimum size of gene sets to consider. Default: 10.
#'
#' @return A list containing:
#'   \item{enrich_go_object}{`enrichResult` object from `clusterProfiler`.}
#'   \item{results_table}{Data frame of enrichment results.}
#'   \item{dot_plot}{`ggplot` object of the dot plot.}
#'   NULL if no significant enrichment is found or if input validation fails.
#'
#' @details
#' This function assumes that the `analysis` slot in the `diff_results`
#' object contains a `gene_ids` column.
#' If this column is not present, or cannot be processed, the function will return NULL.
#'
#' The function includes an internal helper `clean_gene_symbols` which filters common ambiguous
#' gene symbols (snoRNA, snRNA, tRNA) that may not be useful for GO enrichment.
#'
#' @export
analyse_go_enrichment <- function(
    diff_results,
    direction = "cond1",
    org_db = org.Dm.eg.db::org.Dm.eg.db,
    ontology = "BP",
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.2,
    plot_title = NULL,
    show_category = 12,
    label_format_width = 30,
    save = NULL,
    save_results_path = NULL,
    maxGSSize = 1000,
    minGSSize = 10) {
  # Input validation
  stopifnot(is(diff_results, "DamIDResults"))

  if (!inherits(org_db, "OrgDb")) {
    stop("'org_db' must be a valid OrgDb object (e.g., org.Dm.eg.db::org.Dm.eg.db).")
  }

  analysis_df <- diff_results@analysis

  cond <- diff_results@cond
  cond_display_names <- names(cond)
  cond1_display_name <- cond_display_names[1]
  cond2_display_name <- cond_display_names[2]

  # Determine which gene list to use based on 'direction'
  selected_display_name <- ""
  if (direction %in% c("cond1", cond[1], cond1_display_name)) {
    selected_rownames <- rownames(diff_results@upCond1)
    selected_display_name <- cond1_display_name
  } else if (direction %in% c("cond2", cond[2], cond2_display_name)) {
    selected_rownames <- rownames(diff_results@upCond2)
    selected_display_name <- cond2_display_name
  } else if (direction == "all") {
    selected_rownames <- unique(c(rownames(diff_results@upCond1), rownames(diff_results@upCond2)))
    selected_display_name <- "All Significant"
  } else {
    stop("Invalid 'direction' specified. Must be 'cond1' (or condition 1 name), 'cond2' (or condition 2 name), or 'all'.")
  }

  if (length(selected_rownames) == 0) {
    message(sprintf("No significant genes found for direction '%s' (representing %s). Returning NULL.", direction, selected_display_name))
    return(NULL)
  }

  # Ensure the selected_rownames are present in the analysis_df
  selected_rownames <- selected_rownames[selected_rownames %in% rownames(analysis_df)]
  if (length(selected_rownames) == 0) {
    message(sprintf("Selected genes for direction '%s' (representing %s) not found in analysis_data. Returning NULL.", direction, selected_display_name))
    return(NULL)
  }

  query_df <- analysis_df[selected_rownames, , drop = FALSE]

  # Extract query gene list (FBgnIDs), handling potentially comma-separated values
  if (!"gene_ids" %in% colnames(query_df)) {
    stop("Column 'gene_ids' not found in analysis_data. Ensure it is included.")
  }
  gene_list_fbgnid <- unique(unlist(strsplit(query_df$gene_ids, ",")))
  gene_list_fbgnid <- gene_list_fbgnid[nchar(gene_list_fbgnid) > 0]

  if (length(gene_list_fbgnid) == 0) {
    message("No valid Flybase IDs found in the query gene list. Returning NULL.")
    return(NULL)
  }

  # Extract universe gene list (all FBgnIDs in the input data)
  if (!"gene_ids" %in% colnames(analysis_df)) {
    stop("Column 'gene_ids' not found in analysis_data for universe. Ensure it is included.")
  }
  universe_fbgnid <- unique(unlist(strsplit(analysis_df$gene_ids, ",")))
  universe_fbgnid <- universe_fbgnid[nchar(universe_fbgnid) > 0]

  if (length(universe_fbgnid) < length(gene_list_fbgnid) || length(universe_fbgnid) < 11) {
    message("Universe of Flybase IDs is too small or invalid for enrichment analysis. Returning NULL.")
    return(NULL)
  }

  # Map FBgnID to SYMBOL for clusterProfiler for both query and universe
  gene_universe_map <- tryCatch(
    {
      bitr(
        geneID = unique(c(gene_list_fbgnid, universe_fbgnid)),
        fromType = "FLYBASE",
        toType = "SYMBOL",
        OrgDb = org_db
      )
    },
    error = function(e) {
      stop("Failed to map Flybase IDs to symbols using bitr: ", e$message)
    }
  )

  # Get symbols corresponding to the query FBgnIDs
  query_symbols_all <- gene_universe_map$SYMBOL[match(gene_list_fbgnid, gene_universe_map$FLYBASE)]
  query_symbols_all <- na.omit(query_symbols_all)

  # Clean the gene symbols
  cleaned_query_symbols <- clean_gene_symbols(query_symbols_all)

  if (length(cleaned_query_symbols) == 0) {
    message("No valid gene symbols remaining after filtering and mapping. Returning NULL.")
    return(NULL)
  }

  # Get symbols corresponding to the universe FBgnIDs
  universe_symbols <- gene_universe_map$SYMBOL[match(universe_fbgnid, gene_universe_map$FLYBASE)]
  universe_symbols <- na.omit(universe_symbols)
  universe_symbols <- unique(universe_symbols)

  if (length(universe_symbols) < length(cleaned_query_symbols)) {
    message("Mapped universe symbols count is less than query symbols. Adjust universe or check mapping. Returning NULL.")
    return(NULL)
  }

  # Perform GO enrichment
  ego <- tryCatch(
    {
      enrichGO(
        gene          = cleaned_query_symbols,
        OrgDb         = org_db,
        universe      = universe_symbols,
        keyType       = "SYMBOL",
        ont           = ontology,
        pAdjustMethod = "BH",
        maxGSSize = maxGSSize,
        minGSSize = minGSSize,
        pvalueCutoff  = pvalue_cutoff,
        qvalueCutoff  = qvalue_cutoff,
        readable      = TRUE # Retrieves gene symbols in result instead of IDs
      )
    },
    error = function(e) {
      message("Error during GO enrichment: ", e$message)
      return(NULL)
    })

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No significant GO terms found. Returning NULL.")
    return(NULL)
  }

  # Convert to dataframe for plotting
  # ego <- mutate(ego, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  ego.df <- ego %>% as.data.frame()
  ego.df$Description <- ego.df$Description %>% gsub("regulation", "reg.", .)
  ego.df$GeneRatio <- sapply(strsplit(ego.df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), USE.NAMES = FALSE)


  # Generate dot plot
  if (is.null(plot_title)) {
    plot_title <- sprintf("GO Enrichment for %s (%s)", selected_display_name, ontology)
  }

  if (nrow(ego.df) > show_category) {
    plot_df <- ego.df %>% head(n = show_category)
  } else {
    plot_df <- ego.df
  }

  max_x_value <- max(plot_df$GeneRatio, na.rm = TRUE)

  dplot <- ggplot(
    plot_df,
    aes(GeneRatio, fct_reorder(Description, GeneRatio))
  ) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
      colours = c("#f7ca64", "#46bac2", "#7e62a3"),
      name = bquote(italic(p)[adj]),
      trans = "log2",
      guide = guide_colorbar(reverse = TRUE, order = 1),
      labels = scientific_format(digits = 1)
    ) +
    scale_size_continuous(
      range = c(1, 10),
      labels = number_format(accuracy = 1),
      breaks = c(
        (min(ego.df$Count) + 1),
        round(mean(c(min(ego.df$Count), max(ego.df$Count)))),
        max(ego.df$Count)
      )
    ) +
    labs(x = "Gene Ratio", y = NULL, title = plot_title) + # Use labs for clarity
    theme_bw(14) +
    scale_x_continuous(
      limits = c(0, max_x_value),
      expand = expansion(mult = c(0.02, 0.1), add = 0),
      breaks = function(x) {
        f <- as.numeric(as.character(sprintf("%0.1f", x[1])))
        t <- as.numeric(as.character(sprintf("%0.1f", x[2])))
        return(seq(f, t, by = ((t - f) / 2)))
      }
    ) +
    scale_y_discrete(
      labels = function(y) str_wrap(y, width = label_format_width),
      expand = expansion(add = c(0.5, 1))
    ) +
    theme(
      axis.text.y = element_text(lineheight = 0.8),
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, vjust = 1)
    )

  # Save results table if path provided
  if (!is.null(save_results_path)) {
    tryCatch(
      {
        write.csv(ego.df, file = save_results_path, row.names = FALSE)
        message("GO enrichment results table saved to: ", save_results_path)
      },
      error = function(e) {
        warning("Failed to save GO results table to CSV: ", e$message)
      }
    )
  }

  # Plot saving
  save_defaults <- list(
    filename = "damidBind_GSEA_dotplot",
    format = "pdf",
    width = 6,
    height = 6
  )
  save_config <- check_list_input(save_defaults, save)

  if (!is.null(save_config)) {
    tryCatch(
      {
        full_filename <- paste0(save_config$filename, ".", save_config$format)
        ggsave(
          filename = full_filename,
          plot = dplot,
          width = save_config$width,
          height = save_config$height,
          units = "in",
          device = save_config$format
        )
        message("GO dot plot saved to: ", full_filename)
      },
      error = function(e) {
        message("An error occurred while saving the GO dot plot: ", e$message)
      }
    )
  } else {
    # Only print the plot if not saving to a file
    print(dplot)
  }

  invisible(list(
    enrich_go_object = ego,
    results_table = ego.df,
    dot_plot = dplot
  ))
}


#' GO Enrichment Analysis (US spelling alias)
#'
#' This is a direct alias for `analyse_go_enrichment()`; it is provided for users who prefer US English spelling.
#' See `?analyse_go_enrichment` for documentation.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param direction Character string. Specifies which set of genes to analyse, either using condition names, "cond1" or "cond2", or "all" (for all significantly enriched genes from either direction).  Default is "cond1".
#' @param org_db An OrgDb object specifying the organism's annotation database.
#'   For Drosophila, use `org.Dm.eg.db::org.Dm.eg.db`.
#' @param ontology Character string. The GO ontology to use: "BP" (Biological Process),
#'   "MF" (Molecular Function), or "CC" (Cellular Component). Default is "BP".
#' @param pvalue_cutoff Numeric. Adjusted p-value cutoff for significance. Default: 0.05.
#' @param qvalue_cutoff Numeric. Q-value cutoff for significance. Default: 0.2.
#' @param plot_title Character string. Title for the generated dot plot.
#' @param show_category Integer. Number of top enriched GO categories to display in the plot. Default: 12.
#' @param label_format_width Integer. Max character length for GO term labels on the plot. Default: 30.
#' @param save List or `NULL`. Controls saving the plot to a file (dot plot).
#'   If `NULL`, `FALSE`, or `0`, the plot is not saved.
#'   If a `list`, it specifies saving parameters:
#'   - \code{filename} (character): The path and base name for the output file. If not specified, the default name "damidBind_GSEA_dotplot" is used.
#'   - \code{format} (character): File format ("pdf", "svg", or "png").
#'     Default is "pdf".
#'   - \code{width} (numeric): Width of the plot in inches. Default is 6
#'   - \code{height} (numeric): Height of the plot in inches. Default is 6.
#' @param save_results_path Character string or NULL. If a path is provided (e.g., "go_results.csv"), the
#'   enrichment results table will be saved to this CSV file.
#' @param maxGSSize Integer. Maximum size of gene sets to consider. Default: 1000.
#' @param minGSSize Integer. Minimum size of gene sets to consider. Default: 10.
#'
#' @aliases analyze_go_enrichment
#' @export
analyze_go_enrichment <- analyse_go_enrichment


# Internal helper function (not exported)
# This function helps filter out gene symbols that are typically not useful for GO enrichment
clean_gene_symbols <- function(x, clean_regex = "snoRNA|snRNA|tRNA", clean_extra = NULL) {
  if (is.null(x) || length(x) == 0) {
    return(character(0))
  }
  if (!is.null(clean_extra)) {
    clean_regex <- paste(clean_regex, clean_extra, sep = "|")
  }
  keep <- !grepl(clean_regex, as.character(x), ignore.case = TRUE)
  return(x[keep])
}
