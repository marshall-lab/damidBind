#' @title Filter Gene Symbols
#' @description Internal helper to filter out common ambiguous or uninformative gene
#'   symbols (e.g., snoRNA, snRNA) from a list before GO analysis.
#' @param x A character vector of gene symbols to be filtered.
#' @param clean_regex A character string containing the base regular expression for symbols to remove.
#' @param clean_extra An optional character string containing an additional regex pattern to be
#'   appended to `clean_regex`.
#' @return A character vector containing the gene symbols that did not match the filter patterns.
#' @noRd
._clean_gene_symbols <- function(x, clean_regex = "snoRNA|snRNA|tRNA", clean_extra = NULL) {
  if (is.null(x) || length(x) == 0) {
    return(character(0))
  }
  if (!is.null(clean_extra)) {
    clean_regex <- paste(clean_regex, clean_extra, sep = "|")
  }
  keep <- !grepl(clean_regex, as.character(x), ignore.case = TRUE)
  return(x[keep])
}

#' @title Prepare Gene ID Lists for GO Analysis
#' @description Internal helper to select significant genes based on `direction`,
#' and extract query and universe FBgnID lists.
#' @param diff_results A `DamIDResults` object.
#' @param direction The user-specified direction ("cond1", "cond2", etc.).
#' @return A list with `query_fbgnids`, `universe_fbgnids`, and `selected_display_name`,
#' or NULL if no valid genes are found.
#' @noRd
._prepare_go_gene_ids <- function(diff_results, direction) {
  analysis_df <- analysisTable(diff_results)
  cond <- conditionNames(diff_results)
  cond_display_names <- names(cond)
  cond1_display_name <- cond_display_names[1]
  cond2_display_name <- cond_display_names[2]

  # Determine which gene list to use
  if (direction %in% c("cond1", cond[1], cond1_display_name)) {
    selected_rownames <- rownames(enrichedCond1(diff_results))
    selected_display_name <- cond1_display_name
  } else if (direction %in% c("cond2", cond[2], cond2_display_name)) {
    selected_rownames <- rownames(enrichedCond2(diff_results))
    selected_display_name <- cond2_display_name
  } else if (direction == "all") {
    selected_rownames <- c(rownames(enrichedCond1(diff_results)), rownames(enrichedCond2(diff_results)))
    selected_display_name <- "All Significant"
  } else {
    stop("Invalid 'direction'. Must be 'cond1' (or its name), 'cond2' (or its name), or 'all'.")
  }

  if (length(selected_rownames) == 0) {
    message(sprintf("No significant genes for direction '%s' (%s). Returning NULL.", direction, selected_display_name))
    return(NULL)
  }

  # Ensure gene_ids column exists
  if (!"gene_ids" %in% colnames(analysis_df)) {
    stop("Column 'gene_ids' not found in analysis data. Cannot perform GO enrichment.")
  }

  # Extract query gene list (FBgnIDs), handling comma-separated values
  query_df <- analysis_df[rownames(analysis_df) %in% unique(selected_rownames), , drop = FALSE]
  query_fbgnids <- unique(unlist(strsplit(query_df$gene_ids, ",")))
  query_fbgnids <- query_fbgnids[nchar(query_fbgnids) > 0]

  if (length(query_fbgnids) == 0) {
    message("No valid Flybase IDs in the query gene list. Returning NULL.")
    return(NULL)
  }

  # Extract universe gene list
  universe_fbgnids <- unique(unlist(strsplit(analysis_df$gene_ids, ",")))
  universe_fbgnids <- universe_fbgnids[nchar(universe_fbgnids) > 0]

  if (length(universe_fbgnids) < length(query_fbgnids) || length(universe_fbgnids) < 11) {
    message("Universe of Flybase IDs is too small. Returning NULL.")
    return(NULL)
  }

  list(
    query_fbgnids = query_fbgnids,
    universe_fbgnids = universe_fbgnids,
    selected_display_name = selected_display_name
  )
}

#' @title Map and Clean Gene Symbols
#' @description Internal helper to map FBgnIDs to symbols using bitr and clean them.
#' @param query_fbgnids,universe_fbgnids Character vectors of Flybase IDs.
#' @param org_db An OrgDb object.
#' @param clean_gene_symbols Logical, whether to remove snoRNAs/tRNAs prior to enrichment analysis
#' @return A list with `query_symbols` and `universe_symbols`, or NULL on failure.
#' @noRd
._map_and_clean_go_symbols <- function(query_fbgnids, universe_fbgnids, org_db, clean_gene_symbols) {
  # Map FBgnID to SYMBOL for both query and universe
  gene_universe_map <- tryCatch({
    clusterProfiler::bitr(
      geneID = unique(c(query_fbgnids, universe_fbgnids)),
      fromType = "FLYBASE", toType = "SYMBOL", OrgDb = org_db
    )
  }, error = function(e) {
    stop("Failed to map Flybase IDs to symbols using bitr: ", conditionMessage(e))
  })

  # Get symbols for query genes and clean them
  query_symbols_all <- gene_universe_map$SYMBOL[match(query_fbgnids, gene_universe_map$FLYBASE)]
  cleaned_query_symbols <- if (isTRUE(clean_gene_symbols))
    ._clean_gene_symbols(stats::na.omit(query_symbols_all))
  else stats::na.omit(query_symbols_all)

  if (length(cleaned_query_symbols) == 0) {
    message("No valid gene symbols remained after mapping and filtering. Returning NULL.")
    return(NULL)
  }

  # Get symbols for universe genes
  universe_symbols <- gene_universe_map$SYMBOL[match(universe_fbgnids, gene_universe_map$FLYBASE)]
  universe_symbols <- unique(stats::na.omit(universe_symbols))

  if (length(universe_symbols) < length(cleaned_query_symbols)) {
    message("Mapped universe symbol count is less than query symbols. Returning NULL.")
    return(NULL)
  }

  list(query_symbols = cleaned_query_symbols, universe_symbols = universe_symbols)
}

#' @title Run GO Enrichment Analysis
#' @description Internal helper to execute the `enrichGO` call.
#' @param ... Arguments passed to `clusterProfiler::enrichGO`.
#' @return An `enrichResult` object, or NULL on failure.
#' @noRd
._run_go_enrichment <- function(gene, universe, org_db, ontology, pvalue_cutoff,
                                qvalue_cutoff, maxGSSize, minGSSize) {
  ego <- tryCatch({
    clusterProfiler::enrichGO(
      gene = gene, OrgDb = org_db, universe = universe,
      keyType = "SYMBOL", ont = ontology, pAdjustMethod = "BH",
      maxGSSize = maxGSSize, minGSSize = minGSSize,
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
      readable = TRUE
    )
  }, error = function(e) {
    warning("GO enrichment failed: ", conditionMessage(e))
    return(NULL)
  })
  return(ego)
}


#' @title Create GO Dot Plot
#' @description Internal helper to generate a ggplot object from processed GO results.
#' @param results_df A processed data.frame of GO enrichment results.
#' @param plot_title,show_category,label_format_width Plotting parameters.
#' @return A `ggplot` object.
#' @importFrom scales scientific_format number_format
#' @noRd
._create_go_dotplot <- function(results_df, plot_title, show_category, label_format_width) {
  plot_df <- if (nrow(results_df) > show_category) {
    head(results_df, n = show_category)
  } else {
    results_df
  }
  max_x_value <- max(plot_df$GeneRatio, na.rm = TRUE)

  ggplot(
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
      labels = scales::scientific_format(digits = 1)
    ) +
    scale_size_continuous(
      range = c(1, 10),
      labels = scales::number_format(accuracy = 1),
      breaks = tryCatch({
        c(
          (min(results_df$Count) + 1),
          round(mean(c(min(results_df$Count), max(results_df$Count)))),
          max(results_df$Count)
        )
      }, error = function(e) waiver()) # Handle cases with few points
    ) +
    labs(x = "Gene Ratio", y = NULL, title = plot_title) +
    theme_bw(18) +
    scale_x_continuous(
      limits = c(0, max_x_value),
      expand = expansion(mult = c(0.02, 0.1), add = 0),
      breaks = function(x) {
        f <- as.numeric(as.character(sprintf("%0.1f", x[1])))
        t <- as.numeric(as.character(sprintf("%0.1f", x[2])))
        # handle cases where f and t are the same
        if (f == t) return(f)
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
}

#' @title Save GO Results Table
#' @description Internal helper to save enrichment results to a CSV file.
#' @param results_df The data.frame of GO results.
#' @param save_results_path The file path for the CSV.
#' @return Invisibly returns NULL.
#' @noRd
._save_go_table <- function(results_df, save_results_path) {
  if (!is.null(save_results_path)) {
    tryCatch({
      write.csv(results_df, file = save_results_path, row.names = FALSE)
      message("GO enrichment results table saved to: ", save_results_path)
    }, error = function(e) {
      warning("Failed to save GO results table to CSV: ", conditionMessage(e))
    })
  }
  invisible(NULL)
}

#' @title Save or Print GO Plot
#' @description Internal helper to save a ggplot object or print it.
#' @param dplot The ggplot object to process.
#' @param save_config A list of save parameters, or NULL to print.
#' @return Invisibly returns NULL.
#' @noRd
._save_or_print_go_plot <- function(dplot, save_config) {
  if (!is.null(save_config)) {
    tryCatch({
      full_filename <- paste0(save_config$filename, ".", save_config$format)
      ggsave(
        filename = full_filename, plot = dplot,
        width = save_config$width, height = save_config$height,
        units = "in", device = save_config$format
      )
      message("GO dot plot saved to: ", full_filename)
    }, error = function(e) {
      message("Failed to save the GO dot plot: ", conditionMessage(e))
    })
  } else {
    print(dplot)
  }
  invisible(NULL)
}


#' Perform Gene Ontology (GO) Enrichment Analysis for Differentially Bound/Expressed Regions
#'
#' This function performs Gene Ontology (GO) enrichment analysis using `clusterProfiler`
#' for either the up-regulated or down-regulated regions/genes identified by
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
#' @param clean_gene_symbols Logical. Removes snoRNAs and tRNAs (common sources of accidental bias
#'   between different NGS methods) from the gene lists prior to enrichment analysis. Default: TRUE.
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
#' The function includes an internal helper `._clean_gene_symbols` which filters common ambiguous
#' gene symbols (snoRNA, snRNA, tRNA) that may not be useful for GO enrichment.
#'
#' @examples
#' # This example requires the 'org.Dm.eg.db' package
#' if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
#'     # Helper function to create a sample DamIDResults object
#'     .generate_example_results <- function() {
#'         # Define a mock gene object. Note: Real, mappable FlyBase IDs are
#'         # used for the 'gene_id' column to ensure the example runs.
#'         mock_genes_gr <- GenomicRanges::GRanges(
#'             seqnames = S4Vectors::Rle("2L", 7),
#'             ranges = IRanges::IRanges(
#'                 start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'                 end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'             ),
#'             gene_id = c(
#'                 "FBgn0034439", "FBgn0031267", "FBgn0051138", "FBgn0031265",
#'                 "FBgn0004655", "FBgn0000251", "FBgn0000252"
#'             ),
#'             gene_name = c("ap", "dpr1", "side", "dpr2", "eg", "bi", "br")
#'         )
#'         data_dir <- system.file("extdata", package = "damidBind")
#'         loaded_data <- load_data_peaks(
#'             binding_profiles_path = data_dir,
#'             peaks_path = data_dir,
#'             ensdb_genes = mock_genes_gr,
#'             quantile_norm = TRUE
#'         )
#'         diff_results <- differential_binding(
#'             loaded_data,
#'             cond = c("L4", "L5"),
#'             cond_names = c("L4 Neurons", "L5 Neurons")
#'         )
#'         return(diff_results)
#'     }
#'     diff_results <- .generate_example_results()
#'
#'     # Run GO Enrichment for genes enriched in the first condition ("L4")
#'     # Note: with tiny sample data, this may not find significant terms.
#'     go_results <- analyse_go_enrichment(
#'         diff_results,
#'         direction = "L4",
#'         org_db = org.Dm.eg.db::org.Dm.eg.db
#'     )
#'
#'     # Print the results table if any enrichment was found
#'     if (!is.null(go_results)) {
#'         print(go_results$results_table)
#'     }
#' }
#'
#' @aliases analyze_go_enrichment
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
    minGSSize = 10,
    clean_gene_symbols = TRUE) {

  # Input validation
  stopifnot(is(diff_results, "DamIDResults"))
  if (!inherits(org_db, "OrgDb")) {
    stop("'org_db' must be a valid OrgDb object (e.g., org.Dm.eg.db::org.Dm.eg.db).")
  }

  # Prepare gene ID lists
  gene_lists <- ._prepare_go_gene_ids(diff_results, direction)
  if (is.null(gene_lists)) { return(NULL) }

  # Map and clean gene symbols
  gene_symbols <- ._map_and_clean_go_symbols(
    gene_lists$query_fbgnids, gene_lists$universe_fbgnids, org_db, clean_gene_symbols
  )
  if (is.null(gene_symbols)) { return(NULL) }

  # Perform GO enrichment analysis
  ego <- ._run_go_enrichment(
    gene = gene_symbols$query_symbols, universe = gene_symbols$universe_symbols,
    org_db = org_db, ontology = ontology, pvalue_cutoff = pvalue_cutoff,
    qvalue_cutoff = qvalue_cutoff, maxGSSize = maxGSSize, minGSSize = minGSSize
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No significant GO terms found. Returning NULL.")
    return(NULL)
  }

  ego_df <- as.data.frame(ego)
  ego_df$Description <- gsub("regulation", "reg.", ego_df$Description)
  ego_df$GeneRatio <- vapply(strsplit(ego_df$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  }, FUN.VALUE = numeric(1))

  # Plot
  if (is.null(plot_title)) {
    plot_title <- sprintf(
      "GO Enrichment for %s (%s)",
      gene_lists$selected_display_name, ontology
    )
  }
  dplot <- ._create_go_dotplot(ego_df, plot_title, show_category, label_format_width)

  # Save results table and plot (if requested)
  ._save_go_table(ego_df, save_results_path)

  save_defaults <- list(filename = "damidBind_GSEA_dotplot", format = "pdf", width = 6, height = 6)
  save_config <- check_list_input(save_defaults, save) # Assumes check_list_input is available
  ._save_or_print_go_plot(dplot, save_config)

  # Return
  list(
    dot_plot = dplot,
    enrich_go_object = ego,
    results_table = ego_df
  )
}

#' @export
analyze_go_enrichment <- analyse_go_enrichment
