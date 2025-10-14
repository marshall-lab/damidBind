#' Interactive IGV visualisation (Shiny + igvShiny) of differential regions
#'
#' Launches a Shiny app with an embedded IGV browser and an interactive table listing
#' differentially-bound regions (from `differential_binding()` or `differential_acccessibility()` results).
#' Clicking on a region in the table will pan IGV to that locus. Sample coverage and region tracks are loaded as quantitative/annotation tracks.
#' A dedicated "Save as SVG" button is provided to export the current IGV view.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param samples Optional character vector of sample names to display (default: all in dataset).
#' @param colour_cond1,colour_cond2  Colours for differentially enriched region tracks.
#' @param use_genome IGV genome name (inferred from peak annotations if not given).
#' @param padding_width Width to pad browser viewbox on either side of the peak
#' @param host Hostname for the server location (defaults to localhost)
#' @param port Port for connection (if NULL (default) the port is assigned by Shiny)
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
#' @importFrom shiny actionButton HTML tags
#' @export
browse_igv_regions <- function(
    diff_results,
    samples = NULL,
    colour_cond1 = "#ff6600",
    colour_cond2 = "#2288dd",
    use_genome = NULL,
    padding_width = 20000,
    host = "localhost",
    port = NULL,
    default_svg_scale = 2) {
  stopifnot(is(diff_results, "DamIDResults"))

  required <- c("igvShiny", "shiny", "DT")
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for browse_igv_regions().\n  Install with: BiocManager::install('%s')", pkg, pkg), call. = FALSE)
    }
  }

  preptbl <- function(tbl) {
    tbl$chr <- as.character(tbl$chr)
    return(tbl)
  }

  # Get the custom condition names
  cond_display_names <- names(conditionNames(diff_results))
  cond1_display_name <- cond_display_names[1]
  cond2_display_name <- cond_display_names[2]

  # Check for peaks and process if present
  binding_profiles_data <- inputData(diff_results)$binding_profiles_data
  peaks_incl <- FALSE
  if ("pr" %in% names(inputData(diff_results))) {
    peaks <- inputData(diff_results)$pr
    peaks_bed <- peaks %>%
      as.data.frame() %>%
      select(c("seqnames", "start", "end")) %>%
      rename(chr = seqnames)
    peaks_incl <- TRUE
  }

  # Sample columns: default is all, or user selection
  all_sample_cols <- setdiff(names(binding_profiles_data), c("chr", "start", "end"))
  if (is.null(samples)) {
    use_samples <- all_sample_cols
  } else {
    use_samples <- intersect(as.character(samples), all_sample_cols)
    if (length(use_samples) == 0) stop("None of the requested samples present in data.")
  }

  # Prepare data for table/region finding
  occupancy <- analysisTable(diff_results)

  # Extract region locations from rownames
  rn <- rownames(occupancy)
  matches <- str_match(rn, "^(.*?):(\\d+)-(\\d+)")
  genloc <- data.frame(
    region_name = matches[, 1],
    chr = matches[, 2],
    start = as.integer(matches[, 3]),
    end = as.integer(matches[, 4]),
    stringsAsFactors = FALSE
  )
  occ_tab <- cbind(genloc, occupancy)

  # Differentially-regulated regions
  up_idx <- match(rownames(enrichedCond1(diff_results)), occ_tab$region_name)
  down_idx <- match(rownames(enrichedCond2(diff_results)), occ_tab$region_name)

  up_idx <- up_idx[!is.na(up_idx)]
  down_idx <- down_idx[!is.na(down_idx)]

  # Table of differentially-bound regions to display
  region_tab <- occ_tab[c(up_idx, down_idx), , drop = FALSE]
  if (nrow(region_tab) == 0) stop("No differentially bound regions found to display.")

  region_tab$enriched <- c(rep(cond1_display_name, length(up_idx)), rep(cond2_display_name, length(down_idx)))
  region_tab <- region_tab

  # Determine IGV genome name from peak object, or use argument, or fallback to Drosophila dm6 by default
  if (is.null(use_genome)) {
    if (peaks_incl) use_genome <- as.character(genome(peaks))[1]
    if (is.na(use_genome) || length(use_genome) == 0 || is.null(use_genome)) use_genome <- "dm6"
  }

  # Shiny options
  igvshiny_options <- parseAndValidateGenomeSpec(genomeName = use_genome, initialLocus = "all")
  shiny_connection_opts <- list(host = host)
  if (!is.null(port) && is.numeric(port)) {
    shiny_connection_opts[["port"]] <- as.integer(port)
  }

  # Get min/maxes of profiles and logFC to display all at a set scale without clipping
  bp_min <- floor(min(binding_profiles_data[, use_samples]))
  bp_max <- ceiling(max(binding_profiles_data[, use_samples]))
  enrich_max <- ceiling(max(abs(region_tab$logFC)))

  # JavaScript for triggering SVG save
  js_code <- HTML("
      Shiny.addCustomMessageHandler('save-igv-svg', function(message) {
        let browser = document.getElementById(message.igvId).igvBrowser;
        if (browser) {
          browser.saveSVG();
        } else {
          console.error('Could not find IGV browser instance with ID: ' + message.igvId);
        }
      });
    ")

  # Build Shiny app
  igvShinyApp <- shinyApp(
    ui = fluidPage(
      tags$head(tags$script(HTML(
        "
            Shiny.addCustomMessageHandler('igvShiny:saveSVG', function(message) {
                const elementID = message.elementID;
                if (!elementID) {
                    console.error('damidBind: Message did not contain an elementID for saving SVG.');
                    return;
                }
                const igvDiv = document.getElementById(elementID);
                if (!igvDiv || !igvDiv.igvBrowser) {
                    console.error('damidBind: IGV browser instance not found for ID: ' + elementID);
                    return;
                }
                const igvBrowser = igvDiv.igvBrowser;

                // Generate a filename with a timestamp to prevent overwriting
                const now = new Date();
                const timestamp = now.getFullYear().toString() +
                                  (now.getMonth() + 1).toString().padStart(2, '0') +
                                  now.getDate().toString().padStart(2, '0') + '_' +
                                  now.getHours().toString().padStart(2, '0') +
                                  now.getMinutes().toString().padStart(2, '0') +
                                  now.getSeconds().toString().padStart(2, '0');
                const filename = 'igv_snapshot_' + timestamp + '.svg';

                // Call the igv.js API method
                igvBrowser.saveSVGtoFile({ filename: filename });
                console.log('damidBind: SVG snapshot saved as ' + filename);
            });
            "
      ))),
      titlePanel(sprintf("damidBind: Differentially-%s regions", inputData(diff_results)$test_category)),

      sidebarLayout(
        sidebarPanel(
          width = 6,
          h4("Differentially-bound regions"),
          DTOutput("region_table", width = "100%"),
          hr(),
          p("Click a row to navigate IGV to that region."),
          hr(),
          actionButton("saveSVGButton", "Save as SVG", icon = icon("camera")),
        ),
        mainPanel(
          width = 6,
          igvShinyOutput("igv")
        )
      )
    ),
    server = function(input, output, session) {
      # Set up IGV browser
      output$igv <- renderIgvShiny({
        igvShiny(igvshiny_options)
      })

      # Add tracks once ready
      observeEvent(input$igvReady,
                   {
                     message("IGV browser initialized.  Please note that the browser will block further interactive R analysis until stopped.")

                     # Add peaks track if present
                     if (peaks_incl) {
                       loadBedTrack(
                         session,
                         "igv",
                         tbl = preptbl(peaks_bed),
                         trackHeight = 50,
                         trackName = "Binding peaks",
                         color = "darkgreen"
                       )
                       message(" - Added 'Binding peaks' track")
                     }

                     # Add sample quantitative tracks
                     for (sample in use_samples) {
                       bprof <- binding_profiles_data[, c("chr", "start", "end", sample)]
                       loadBedGraphTrack(
                         session,
                         "igv",
                         tbl = preptbl(bprof),
                         autoscale = FALSE,
                         min = bp_min,
                         max = bp_max,
                         trackHeight = 65,
                         trackName = sample,
                         color = "#6666cc"
                       )
                       message(sprintf(" - Added sample track: %s", sample))
                     }

                     regcol_sel <- c("chr", "start", "end", "logFC")
                     # differentially-bound, significant cond1 regions track
                     if (length(up_idx) > 0) {
                       upCond1_df <- region_tab[region_tab$enriched == cond1_display_name, regcol_sel, drop = FALSE] # Use custom name
                       loadBedGraphTrack(
                         session,
                         "igv",
                         tbl = preptbl(upCond1_df),
                         autoscale = FALSE,
                         min = 0,
                         max = enrich_max,
                         trackName = sprintf("Enriched (%s)", cond1_display_name), # Use custom name
                         color = colour_cond1
                       )
                       message(sprintf(" - Added 'Enriched (%s) regions track", cond1_display_name))
                     }

                     # differentially-bound, significant cond2 regions track
                     if (length(down_idx) > 0) {
                       upCond2_df <- region_tab[region_tab$enriched == cond2_display_name, regcol_sel, drop = FALSE] # Use custom name
                       loadBedGraphTrack(
                         session,
                         "igv",
                         tbl = preptbl(upCond2_df),
                         autoscale = FALSE,
                         min = -enrich_max,
                         max = 0,
                         trackName = sprintf("Enriched (%s)", cond2_display_name), # Use custom name
                         color = colour_cond2
                       )
                       message(sprintf(" - Added 'Enriched (%s) regions track", cond2_display_name))
                     }
                   },
                   once = TRUE
      )

      # Region datatable for display
      output$region_table <- renderDT({
        show_cols <- intersect(
          c("region_name", "logFC", "enriched", "gene_names", "gene_ids"),
          names(region_tab)
        )
        datatable(
          region_tab[, show_cols, drop = FALSE],
          rownames = FALSE,
          selection = "single",
          options = list(pageLength = 12, order = list(list(2, "desc")), scrollX = TRUE),
          width = "100%"
        ) %>% formatRound(columns = 2, digits = 2)
      })

      # Pan IGV when a region is selected in table
      observeEvent(input$region_table_rows_selected, {
        if (length(input$region_table_rows_selected)) {
          sel_row <- input$region_table_rows_selected[1]
          sel_reg <- region_tab[sel_row, ]
          region_string <- sprintf("%s:%d-%d", sel_reg$chr, max(sel_reg$start - padding_width, 1), sel_reg$end + padding_width)
          showGenomicRegion(session, "igv", region_string)
        }
      })

      observeEvent(input$saveSVGButton, {
        message("Saving as SVG...")
        session$sendCustomMessage(
          type = "igvShiny:saveSVG",
          message = list(elementID = "igv")
        )
      })

    },
    options = list(shiny_connection_opts)
  )
  invisible(runApp(igvShinyApp, launch.browser = TRUE))
}
