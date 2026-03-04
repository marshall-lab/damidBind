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
#' @param average_tracks Logical. If `TRUE`, displays the average signal for each
#'   condition instead of individual replicate tracks. (Default: `FALSE`)
#' @export
#' @return A list containing prepared data frames and vectors.
#' @noRd
.prepare_igv_data <- function(diff_results, samples, average_tracks = FALSE) {
    # Get custom condition names
    cond_map <- conditionNames(diff_results)
    cond_display_names <- names(cond_map)
    cond_internal_ids <- as.character(cond_map)

    # Process binding profiles and peaks data if it exists
    binding_profiles_data <- inputData(diff_results)$binding_profiles_data

    # igvShiny does not place nicely with GRanges, so we shift everything over to
    # dataframes from here onwards
    drop_gr_cols <- c("width","strand")
    binding_profiles_data <- as.data.frame(binding_profiles_data) %>%
        subset(select = setdiff(names(.), drop_gr_cols)) %>%
        rename(chr = "seqnames")

    peaks_bed <- NULL
    if ("pr" %in% names(inputData(diff_results))) {
        peaks <- inputData(diff_results)$pr
        peaks_bed <- as.data.frame(peaks) %>%
            select(c("seqnames", "start", "end")) %>%
            rename(chr = "seqnames")
    }

    # Identify all available sample columns
    all_raw_cols <- setdiff(names(binding_profiles_data), c("chr", "start", "end"))

    if (isTRUE(average_tracks)) {
        message("Averaging replicates for IGV display tracks...")
        matched <- inputData(diff_results)$matched_samples
        avg_sample_names <- character(0)

        for (i in seq_along(matched)) {
            internal_id <- names(matched)[i]
            # These are the potentially hyphenated names from metadata
            replicates_metadata <- matched[[i]]

            # Apply make.names to match the sanitized dot format R uses for data frames
            replicates_sanitised <- make.names(replicates_metadata)

            # Map back to display name for the track label
            display_idx <- which(cond_internal_ids == internal_id)
            display_name <- if (length(display_idx) > 0) cond_display_names[display_idx] else internal_id

            avg_name <- sprintf("%s (Avg)", display_name)

            # Find which sanitised names exist in our data frame columns
            valid_reps <- intersect(replicates_sanitised, all_raw_cols)

            if (length(valid_reps) > 0) {
                # Calculate row means using the matched column names
                binding_profiles_data[[avg_name]] <- rowMeans(
                    binding_profiles_data[, valid_reps, drop = FALSE],
                    na.rm = TRUE
                )
                avg_sample_names <- c(avg_sample_names, avg_name)
            } else {
                # Diagnostic warning if the match still fails
                warning(sprintf("Could not find columns for condition: %s", internal_id))
            }
        }
        use_samples <- avg_sample_names
    } else {
        # Default behaviour for individual replicates
        use_samples <- if (is.null(samples)) all_raw_cols else intersect(as.character(samples), all_raw_cols)
    }

    if (length(use_samples) == 0) {
        stop("None of the requested samples or averaged conditions are present in the data.")
    }

    # Prepare navigation table
    occ_tab <- analysisTable(diff_results)
    occ_tab$region_name <- rownames(occ_tab)

    # Parse coordinates for IGV
    base_region_name <- sub("\\.\\d+$", "", occ_tab$region_name)
    matches <- str_match(base_region_name, "^(.*?):(\\d+)-(\\d+)")
    occ_tab$chr <- matches[, 2]
    occ_tab$start <- as.integer(matches[, 3])
    occ_tab$end <- as.integer(matches[, 4])

    # Filter to significant rows
    up_rows <- rownames(enrichedCond1(diff_results))
    down_rows <- rownames(enrichedCond2(diff_results))
    region_tab <- occ_tab[rownames(occ_tab) %in% c(up_rows, down_rows), , drop = FALSE]

    if (nrow(region_tab) == 0) {
        stop("No differentially bound regions found to display.")
    }

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
#' @param prepped_data,track_scales,colours,trackheight Data and configuration objects.
#' @param peakCol,trackCol Colours for tracks.
#' @return `NULL` invisibly
#' @noRd
.load_igv_tracks <- function(session, prepped_data, track_scales, colours, trackheight, peakCol, trackCol) {
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
                     trackName = "Binding peaks", color = peakCol
        )
        message(" - Added 'Binding peaks' track")
    }

    # Use the display names prepared upstream in browse_igv_regions
    sample_display_names <- prepped_data$sample_display_names

    # Add sample quantitative tracks
    for (sample in prepped_data$use_samples) {
        display_name <- sample_display_names[[sample]]
        bprof <- prepped_data$binding_profiles_data[, c("chr", "start", "end", sample)]
        names(bprof)[names(bprof) == sample] <- display_name

        loadBedGraphTrack(session, "igv",
                          tbl = preptbl(bprof), autoscale = FALSE,
                          min = track_scales$bp_min, max = track_scales$bp_max,
                          trackHeight = trackheight, trackName = display_name, color = trackCol
        )
        message(sprintf(" - Added track: %s", display_name))
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
#' @param prepped_data,shiny_configs,colours,padding_width,trackheight Data and config objects.
#' @param peakCol,trackCol Colours for tracks.
#' @return A Shiny server function.
#' @noRd
.define_igv_shiny_server <- function(prepped_data, shiny_configs, colours, padding_width, trackheight, peakCol, trackCol) {
    server_func <- function(input, output, session) {
        output$igv <- renderIgvShiny({
            igvShiny(shiny_configs$igv_options)
        })

        observeEvent(input$igvReady,
                     {
                         .load_igv_tracks(
                             session,
                             prepped_data,
                             shiny_configs$track_scales,
                             colours,
                             trackheight,
                             peakCol,
                             trackCol
                         )
                     },
                     once = TRUE
        )

        output$region_table <- DT::renderDT({
            show_cols <- intersect(c("region_name", "logFC", "enriched", "gene_name", "gene_id"), names(prepped_data$region_tab))
            DT::datatable(
                prepped_data$region_tab[, show_cols, drop = FALSE],
                rownames = FALSE, selection = "single",
                options = list(pageLength = 12, order = list(list(which(show_cols == "logFC") - 1, "desc")), scrollX = TRUE),
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
#' Alternatively, if `export_data_archive` is provided, the function will bypass the browser and
#' save all tracks (including averaged tracks) as a zip archive of bedGraph and BED files.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param samples Optional character vector of sample names to display (default: all in dataset).
#' @param use_unique_ids Logical.  When `TRUE` (default), simplified unique sample names will be
#'   displayed.  Set as `FALSE` to use the full sample file names from loading.
#' @param average_tracks Logical. If `TRUE`, displays the average signal for each
#'   condition instead of individual replicate tracks. (Default: `FALSE`)
#' @param export_data_archive Character or NULL. If a filename is provided (e.g. "my_data"),
#'   the function exports all IGV tracks as a zip archive ("my_data.zip") and exits
#'   without launching the browser. (Default: NULL)
#' @param colour_cond1,colour_cond2  Colours for differentially-enriched region tracks.
#' @param use_genome IGV genome name (inferred from peak annotations if not provided).
#' @param padding_width Width to pad browser viewbox on either side of the peak (Default: 20000)
#' @param trackHeight Height of bedGraph tracks (Default: 65)
#' @param peakColour Colour for significant peaks track (Default: "darkgreen")
#' @param trackColour Colour for bedGraph tracks (Default: "#6666ff")
#' @param host Hostname for the server location (Default: "localhost").
#' @param port Port for connection (if NULL (default) the port is assigned by Shiny).
#'
#' @return Invisibly returns the Shiny app object if launching the browser, or the path
#'   to the zip archive if exporting data.
#'
#' @examples
#' \dontrun{
#' if (isTRUE(curl::has_internet()) && interactive()) {
#' # This example launches an interactive Shiny app in a web browser
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
#'         cond = c("L4 Neurons" = "L4",
#'                  "L5 Neurons" = "L5")
#'     )
#'     return(diff_results)
#' }
#' diff_results <- .generate_example_results()
#'
#' # Launch the interactive browser
#' browse_igv_regions(diff_results)
#'
#' # Export data instead of launching browser
#' browse_igv_regions(diff_results, export_data_archive = "L4_vs_L5_tracks")
#' }
#' }
#'
#' @export
browse_igv_regions <- function(
        diff_results,
        samples = NULL,
        use_unique_ids = TRUE,
        average_tracks = FALSE,
        export_data_archive = NULL,
        colour_cond1 = "#ff6600",
        colour_cond2 = "#2288dd",
        use_genome = NULL,
        padding_width = 20000,
        trackHeight = 65,
        peakColour = "darkgreen",
        trackColour = "#6666ff",
        host = "localhost",
        port = NULL) {

    # Validate inputs and dependencies
    stopifnot(is(diff_results, "DamIDResults"))

    # Dependencies only needed for browser mode
    if (is.null(export_data_archive)) {
        .check_igv_dependencies()
    }

    # Prepare data structures
    prepped_data <- .prepare_igv_data(diff_results, samples, average_tracks = average_tracks)

    # Determine unique sample display names for tracks
    sample_names <- prepped_data$use_samples
    sample_display_names <- if (isTRUE(use_unique_ids) && !isTRUE(average_tracks)) {
        extract_unique_sample_ids(sample_names)
    } else {
        sample_names
    }
    prepped_data$sample_display_names <- stats::setNames(sample_display_names, sample_names)

    # Handle data export if requested
    if (!is.null(export_data_archive)) {
        archive_path <- ._export_igv_data_to_zip(prepped_data, export_data_archive)
        return(invisible(archive_path))
    }

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
        padding_width = padding_width,
        trackheight = trackHeight,
        peakCol = peakColour,
        trackCol = trackColour
    )

    # Create and run the Shiny app
    igvShinyApp <- shinyApp(
        ui = app_ui,
        server = app_server,
        options = list(shiny_connection_opts)
    )

    invisible(runApp(igvShinyApp, launch.browser = TRUE))
}


#' Export processed IGV tracks to a zip archive
#' @description Internal helper that writes dataframes to disk as BED and bedGraph
#'   files before zipping them.
#' @param prepped_data List. The result of `.prepare_igv_data`.
#' @param filename Character. The base name or path of the zip archive.
#' @return Character. Path to the created zip file.
#' @noRd
._export_igv_data_to_zip <- function(prepped_data, filename) {
    message("Exporting IGV tracks to zip archive...")

    # Check if system zip is available
    if (Sys.which("zip") == "" && .Platform$OS.type == "unix") {
        stop("The 'zip' utility is not found on your system PATH. Please install zip or use a different environment.")
    }

    # Ensure filename has .zip extension
    if (!grepl("\\.zip$", filename, ignore.case = TRUE)) {
        filename <- paste0(filename, ".zip")
    }

    # Resolve the absolute destination path immediately
    # This prevents path doubling when we later change working directories
    wd_orig <- getwd()
    zip_path <- if (fs::is_absolute_path(filename)) {
        filename
    } else {
        file.path(wd_orig, filename)
    }

    # Verify that the destination directory exists
    dest_dir <- dirname(zip_path)
    if (!dir.exists(dest_dir)) {
        stop(sprintf("Destination directory does not exist: %s", dest_dir))
    }

    # Create a staging directory in temp
    temp_dir <- file.path(tempdir(), paste0("damidBind_export_", round(as.numeric(Sys.time()))))
    if (!dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Failed to create temporary staging directory for export.")
    }

    # Helper to write bedGraph with 0-based start
    write_bg <- function(df, val_col, out_name) {
        out_df <- df[, c("chr", "start", "end", val_col)]
        # Convert 1-based start to 0-based for bedGraph format
        out_df$start <- out_df$start - 1

        write.table(
            out_df,
            file = out_name,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    }

    exported_files <- character(0)

    # Export binding peaks
    if (!is.null(prepped_data$peaks_bed)) {
        peak_file <- file.path(temp_dir, "unified_binding_peaks.bed")
        peak_out <- prepped_data$peaks_bed
        peak_out$start <- peak_out$start - 1
        write.table(peak_out, file = peak_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        exported_files <- c(exported_files, "unified_binding_peaks.bed")
    }

    # Export profile tracks
    for (sample in names(prepped_data$sample_display_names)) {
        disp_name <- prepped_data$sample_display_names[[sample]]
        # Sanitise label for filename
        safe_name <- make.names(disp_name)
        bg_filename <- paste0("track_signal_", safe_name, ".bedGraph")
        bg_full_path <- file.path(temp_dir, bg_filename)

        write_bg(prepped_data$binding_profiles_data, sample, bg_full_path)
        exported_files <- c(exported_files, bg_filename)
    }

    # Export enrichment tracks
    cond1_name <- prepped_data$cond_display_names[1]
    cond2_name <- prepped_data$cond_display_names[2]
    upCond1_df <- prepped_data$region_tab[prepped_data$region_tab$enriched == cond1_name, ]
    upCond2_df <- prepped_data$region_tab[prepped_data$region_tab$enriched == cond2_name, ]

    if (nrow(upCond1_df) > 0) {
        enr1_filename <- paste0("enrichment_", make.names(cond1_name), ".bedGraph")
        write_bg(upCond1_df, "logFC", file.path(temp_dir, enr1_filename))
        exported_files <- c(exported_files, enr1_filename)
    }

    if (nrow(upCond2_df) > 0) {
        enr2_filename <- paste0("enrichment_", make.names(cond2_name), ".bedGraph")
        write_bg(upCond2_df, "logFC", file.path(temp_dir, enr2_filename))
        exported_files <- c(exported_files, enr2_filename)
    }

    # Run zip command from within the temp directory
    # This ensures paths inside the zip archive are flat
    on.exit(setwd(wd_orig))
    setwd(temp_dir)

    # Capture success status
    zip_status <- utils::zip(zip_path, files = exported_files)

    # Revert to original directory before checking existence or reporting success
    setwd(wd_orig)

    if (zip_status != 0) {
        stop(sprintf("Zip command failed with exit code %d. Check if the output path is writable.", zip_status))
    }

    if (!file.exists(zip_path)) {
        stop("Zip archive was not found after export command. Check system permissions.")
    }

    message(sprintf("Successfully exported %d tracks to %s", length(exported_files), zip_path))
    return(zip_path)
}
