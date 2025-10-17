#' Check for Shiny/IGV Dependencies helper function
#' @description This internal function stops execution if the required packages
#' for the interactive browser are not installed.
#' @return `NULL` invisibly, or stops with an error.
#' @noRd
.check_igv_dependencies <- function() {
  required <- c("igvShiny", "shiny", "DT")
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package '%s' is required for browse_igv_regions().\n  Install with: BiocManager::install('%s')",
        pkg, pkg
      ), call. = FALSE)
    }
  }
}

#' Prepare Data for IGV Shiny App helper function
#' @description This internal function extracts and transforms data from the
#' `DamIDResults` object into the formats needed by the Shiny app, including
#' the main region table, peak data, and sample information.
#' @param diff_results A `DamIDResults` object.
#' @param samples Character vector of specific samples to include.
#' @return A list containing prepared data frames and vectors.
#' @noRd
.prepare_igv_data <- function(diff_results, samples) {
  # Get custom condition names
  cond_display_names <- names(conditionNames(diff_results))

  # Process peaks data if it exists
  binding_profiles_data <- inputData(diff_results)$binding_profiles_data
  peaks_bed <- NULL
  if ("pr" %in% names(inputData(diff_results))) {
    peaks <- inputData(diff_results)$pr
    peaks_bed <- as.data.frame(peaks) %>%
      select(c("seqnames", "start", "end")) %>%
      rename(chr = "seqnames")
  }

  # Determine which sample columns to use
  all_sample_cols <- setdiff(names(binding_profiles_data), c("chr", "start", "end"))
  use_samples <- if (is.null(samples)) all_sample_cols else intersect(as.character(samples), all_sample_cols)
  if (length(use_samples) == 0) stop("None of the requested samples are present in the data.")

  # Start with the full analysis table, preserving its correct row names
  occ_tab <- analysisTable(diff_results)

  # Use the original row names as the definitive region identifier
  occ_tab$region_name <- rownames(occ_tab)

  # Parse coordinates from the region name.
  base_region_name <- sub("\\.\\d+$", "", occ_tab$region_name)
  matches <- str_match(base_region_name, "^(.*?):(\\d+)-(\\d+)")

  # Add coordinate columns, which will be NA if parsing fails for any reason
  occ_tab$chr <- matches[, 2]
  occ_tab$start <- as.integer(matches[, 3])
  occ_tab$end <- as.integer(matches[, 4])

  # Get row names of significantly regulated regions
  up_rows <- rownames(enrichedCond1(diff_results))
  down_rows <- rownames(enrichedCond2(diff_results))

  # Filter `occ_tab` to create the final table for the interactive display.
  # Subsetting is done on the row names, which are guaranteed to match.
  region_tab <- occ_tab[rownames(occ_tab) %in% c(up_rows, down_rows), , drop = FALSE]

  if (nrow(region_tab) == 0) {
    stop("No differentially bound regions found to display.")
  }

  # Assign the 'enriched' status based on which significant list the region belongs to
  region_tab$enriched <- ifelse(
    rownames(region_tab) %in% up_rows,
    cond_display_names[1],
    cond_display_names[2]
  )

  list(
    region_tab = region_tab,
    peaks_bed = peaks_bed,
    binding_profiles_data = binding_profiles_data,
    use_samples = use_samples,
    cond_display_names = cond_display_names
  )
}


#' Configure IGV and Shiny options helper function
#' @description This internal function sets up configuration lists for `igvShiny`,
#' Shiny connection options, and calculates scales for track visualisation.
#' @param prepped_data The list returned by `.prepare_igv_data()`.
#' @param final_genome Character string for the IGV genome.
#' @param host,port Connection options for Shiny.
#' @return A list containing `igv_options`, `shiny_opts`, and `track_scales`.
#' @noRd
.configure_igv_shiny <- function(prepped_data, final_genome, host, port) {
  # Set genome and location
  igvshiny_options <- parseAndValidateGenomeSpec(genomeName = final_genome, initialLocus = "all")

  # Calculate min/max for bedgraph tracks to ensure consistent scaling
  bp_min <- floor(min(prepped_data$binding_profiles_data[, prepped_data$use_samples]))
  bp_max <- ceiling(max(prepped_data$binding_profiles_data[, prepped_data$use_samples]))
  enrich_max <- ceiling(max(abs(prepped_data$region_tab$logFC)))
  track_scales <- list(bp_min = bp_min, bp_max = bp_max, enrich_max = enrich_max)

  list(
    igv_options = igvshiny_options,
    track_scales = track_scales
  )
}


#' Build the IGV Shiny App UI helper function
#' @description This internal function constructs the `fluidPage` UI for the
#' interactive browser, including custom JavaScript for SVG saving.
#' @param diff_results A `DamIDResults` object to extract the title category.
#' @return A Shiny UI object.
#' @noRd
.build_igv_shiny_ui <- function(diff_results) {
  # Custom JavaScript to handle saving the IGV panel as an SVG file
  svg_save_js <- tags$head(tags$script(HTML("
      Shiny.addCustomMessageHandler('igvShiny:saveSVG', function(message) {
        const elementID = message.elementID;
        if (!elementID) {
          console.error('damidBind: Message did not contain elementID for saving SVG.');
          return;
        }
        const igvDiv = document.getElementById(elementID);
        if (!igvDiv || !igvDiv.igvBrowser) {
          console.error('damidBind: IGV browser instance not found for ID: ' + elementID);
          return;
        }
        const igvBrowser = igvDiv.igvBrowser;
        const now = new Date();
        const timestamp = now.getFullYear().toString() +
                          (now.getMonth() + 1).toString().padStart(2, '0') +
                          now.getDate().toString().padStart(2, '0') + '_' +
                          now.getHours().toString().padStart(2, '0') +
                          now.getMinutes().toString().padStart(2, '0') +
                          now.getSeconds().toString().padStart(2, '0');
        const filename = `igv_snapshot_${timestamp}.svg`;
        igvBrowser.saveSVGtoFile({ filename: filename });
        console.log(`damidBind: SVG snapshot saved as ${filename}`);
      });
    ")))

  fluidPage(
    svg_save_js,
    titlePanel(sprintf("damidBind: Differentially-%s regions", inputData(diff_results)$test_category)),
    sidebarLayout(
      sidebarPanel(
        width = 6,
        h4("Differentially-bound regions"),
        DT::DTOutput("region_table", width = "100%"),
        hr(),
        p("Click a row to navigate IGV to that region."),
        hr(),
        actionButton("saveSVGButton", "Save as SVG", icon = icon("camera"))
      ),
      mainPanel(width = 6, igvShinyOutput("igv"))
    )
  )
}

#' Load Tracks into IGV helper function
#' @description This internal function is called by the server to load all
#' BED, bedGraph, and annotation tracks into the IGV browser.
#' @param session The Shiny server session object.
#' @param prepped_data,track_scales,colours Data and configuration objects.
#' @return `NULL` invisibly
#' @noRd
.load_igv_tracks <- function(session, prepped_data, track_scales, colours) {
  message("IGV browser initialized. Loading tracks...")
  preptbl <- function(tbl) { # A small helper to ensure 'chr' column exists
    if ("seqnames" %in% names(tbl)) tbl <- rename(tbl, chr = "seqnames")
    tbl$chr <- as.character(tbl$chr)
    return(tbl)
  }
  # Add binding peaks track if present
  if (!is.null(prepped_data$peaks_bed)) {
    loadBedTrack(session, "igv",
                 tbl = preptbl(prepped_data$peaks_bed), trackHeight = 50,
                 trackName = "Binding peaks", color = "darkgreen"
    )
    message(" - Added 'Binding peaks' track")
  }
  # Add sample quantitative tracks
  for (sample in prepped_data$use_samples) {
    bprof <- prepped_data$binding_profiles_data[, c("chr", "start", "end", sample)]
    loadBedGraphTrack(session, "igv",
                      tbl = preptbl(bprof), autoscale = FALSE,
                      min = track_scales$bp_min, max = track_scales$bp_max,
                      trackHeight = 65, trackName = sample, color = "#6666cc"
    )
    message(sprintf(" - Added sample track: %s", sample))
  }
  # Add differentially-enriched region tracks
  cond1_name <- prepped_data$cond_display_names[1]
  cond2_name <- prepped_data$cond_display_names[2]
  upCond1_df <- prepped_data$region_tab[prepped_data$region_tab$enriched == cond1_name, ]
  upCond2_df <- prepped_data$region_tab[prepped_data$region_tab$enriched == cond2_name, ]

  if (nrow(upCond1_df) > 0) {
    loadBedGraphTrack(session, "igv",
                      tbl = preptbl(upCond1_df[, c("chr", "start", "end", "logFC")]), autoscale = FALSE,
                      min = 0, max = track_scales$enrich_max,
                      trackName = sprintf("Enriched (%s)", cond1_name), color = colours$cond1
    )
    message(sprintf(" - Added 'Enriched in %s' track", cond1_name))
  }
  if (nrow(upCond2_df) > 0) {
    loadBedGraphTrack(session, "igv",
                      tbl = preptbl(upCond2_df[, c("chr", "start", "end", "logFC")]), autoscale = FALSE,
                      min = -track_scales$enrich_max, max = 0,
                      trackName = sprintf("Enriched (%s)", cond2_name), color = colours$cond2
    )
    message(sprintf(" - Added 'Enriched in %s' track", cond2_name))
  }
}

#' Define the IGV Shiny App server helper function
#' @description This internal function creates the server function for the Shiny
#' app, defining reactive outputs and event observers.
#' @param prepped_data,shiny_configs,colours,padding_width Data and config objects.
#' @return A Shiny server function.
#' @noRd
.define_igv_shiny_server <- function(prepped_data, shiny_configs, colours, padding_width) {
  server_func <- function(input, output, session) {
    output$igv <- renderIgvShiny({ igvShiny(shiny_configs$igv_options) })

    observeEvent(input$igvReady,
                 { .load_igv_tracks(session, prepped_data, shiny_configs$track_scales, colours) },
                 once = TRUE
    )

    output$region_table <- DT::renderDT({
      show_cols <- intersect(c("region_name", "logFC", "enriched", "gene_names", "gene_ids"), names(prepped_data$region_tab))
      DT::datatable(
        prepped_data$region_tab[, show_cols, drop = FALSE],
        rownames = FALSE, selection = "single",
        options = list(pageLength = 12, order = list(list(which(show_cols == "logFC")-1, "desc")), scrollX = TRUE),
        width = "100%"
      ) %>% DT::formatRound(columns = "logFC", digits = 2)
    })

    observeEvent(input$region_table_rows_selected, {
      if (length(input$region_table_rows_selected)) {
        sel_row <- prepped_data$region_tab[input$region_table_rows_selected[1], ]
        region_string <- sprintf("%s:%d-%d", sel_row$chr, max(sel_row$start - padding_width, 1), sel_row$end + padding_width)
        showGenomicRegion(session, "igv", region_string)
      }
    })

    observeEvent(input$saveSVGButton, {
      message("Saving IGV view as SVG...")
      session$sendCustomMessage(type = "igvShiny:saveSVG", message = list(elementID = "igv"))
    })
  }
  return(server_func)
}


#' Interactive IGV visualisation (Shiny + igvShiny) of differential regions
#'
#' Launches a Shiny app with an embedded IGV browser and an interactive table listing
#' differentially-bound regions (from `differential_binding()` or `differential_accessibility()` results).
#' Clicking on a region in the table will pan IGV to that locus. Sample coverage and region tracks are loaded as quantitative/annotation tracks.
#' A dedicated "Save as SVG" button is provided to export the current IGV view.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param samples Optional character vector of sample names to display (default: all in dataset).
#' @param colour_cond1,colour_cond2  Colours for differentially enriched region tracks.
#' @param use_genome IGV genome name (inferred from peak annotations if not given).
#' @param padding_width Width to pad browser viewbox on either side of the peak.
#' @param host Hostname for the server location (defaults to localhost).
#' @param port Port for connection (if NULL (default) the port is assigned by Shiny).
#' @return Invisibly returns the Shiny app object created by `shinyApp()`.
#'
#' @examples
#' \donttest{
#' # This example launches an interactive Shiny app and is not run by
#' # automated checks. It requires an internet connection for IGV.
#'
#' .generate_example_results <- function() {
#'     mock_genes_gr <- GenomicRanges::GRanges(
#'         seqnames = S4Vectors::Rle("2L", 7),
#'         ranges = IRanges::IRanges(
#'             start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'             end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'         ),
#'         gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'         gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene")
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
#' # Launch the interactive browser (requires network access; uncomment to run)
#' # browse_igv_regions(diff_results)
#' }
#'
#' @importFrom shiny actionButton HTML tags icon h4 hr p
#' @export
browse_igv_regions <- function(
    diff_results,
    samples = NULL,
    colour_cond1 = "#ff6600",
    colour_cond2 = "#2288dd",
    use_genome = NULL,
    padding_width = 20000,
    host = "localhost",
    port = NULL) {

  # Validate inputs and dependencies
  stopifnot(is(diff_results, "DamIDResults"))
  .check_igv_dependencies()

  # Prepare data structures
  prepped_data <- .prepare_igv_data(diff_results, samples)

  # Determine IGV genome name from peak object, or use argument, or fallback to Drosophila dm6 by default
  final_genome <- use_genome
  if (is.null(final_genome)) {
    if ("pr" %in% names(inputData(diff_results))) {
      final_genome <- as.character(genome(inputData(diff_results)$pr))[1]
    }
    if (is.na(final_genome) || length(final_genome) == 0 || is.null(final_genome)) {
      final_genome <- "dm6" # Fallback genome
    }
  }

  shiny_configs <- .configure_igv_shiny(prepped_data, final_genome, host, port)

  # Shiny options
  shiny_connection_opts <- list(host = host)
  if (!is.null(port) && is.numeric(port)) {
    shiny_connection_opts[["port"]] <- as.integer(port)
  }

  # Build UI and Server
  app_ui <- .build_igv_shiny_ui(diff_results)
  app_server <- .define_igv_shiny_server(
    prepped_data,
    shiny_configs,
    colours = list(cond1 = colour_cond1, cond2 = colour_cond2),
    padding_width = padding_width
  )

  # Create and run the Shiny app
  igvShinyApp <- shinyApp(
    ui = app_ui,
    server = app_server,
    options = list(shiny_connection_opts)
  )

  invisible(runApp(igvShinyApp, launch.browser = TRUE))
}

