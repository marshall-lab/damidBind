#' Draw Proportional Venn Diagram for differential binding analysis
#'
#' Generates a two-set Venn/proportional diagram summarising the results of the differential binding analysis.
#' The set union represents significant binding peaks that fail to show significant differences in occupancy;
#' the exclusive regions of each set represent regions with enriched differential binding in that condition.
#' Note that regions can be bound in both conditions, and still show differential occupancy.
#'
#' @param diff_results A `DamIDResults` object, as returned by
#'   `differential_binding()` or `differential_accessibility()`.
#' @param title Plot title to use.
#' @param subtitle Subtitle to use (default is empty).
#' @param set_labels Character vector of length 2. Names for the two sets/circles (defaults to the analysis condition names).
#' @param filename Character. Path at which to save the diagram, if not NULL.
#' @param font Font name to use (default is "sans")
#' @param format Character. Output plot format, "pdf" or "svg" (default "pdf").
#' @param region_colours Character vector of length 2 or 3. Fill colours for each set region (default: c("#FFA500", "#2288DD", "#CCCCCC")).
#' @return The function is called to generating a plot. It invisibly returns `NULL`.
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
#' # Generate the Venn diagram
#' plot_venn(diff_results)
#'
#' @export
plot_venn <- function(
    diff_results,
    title = "Enriched binding at loci",
    subtitle = "",
    set_labels = NULL,
    filename = NULL,
    font = "sans",
    format = c("pdf", "svg"),
    region_colours = c("#FFA500", "#2288DD", "#CCCCCC")) {
    # Argument and field checks
    stopifnot(is(diff_results, "DamIDResults"))

    if (is.null(set_labels)) {
        set_labels <- names(conditionNames(diff_results))
    }

    upCond1 <- rownames(enrichedCond1(diff_results))
    upCond2 <- rownames(enrichedCond2(diff_results))
    all_ids <- rownames(analysisTable(diff_results))
    # Defensive: remove NAs or empty
    upCond1 <- upCond1[!is.na(upCond1) & nchar(upCond1) > 0]
    upCond2 <- upCond2[!is.na(upCond2) & nchar(upCond2) > 0]
    all_ids <- all_ids[!is.na(all_ids) & nchar(all_ids) > 0]

    # Check minimum viable dataset
    if (length(all_ids) < 2) {
        stop("Not enough loci for Venn diagram (need at least 2 in 'all').")
    }

    # Non-significant set
    nonsig <- setdiff(all_ids, union(upCond1, upCond2))
    Cond1_full <- union(upCond1, nonsig)
    Cond2_full <- union(upCond2, nonsig)

    # Warn if one group is empty
    if (length(upCond1) == 0) message(sprintf("Note: No loci were differentially enriched in condition 1 ('%s')", set_labels[1]))
    if (length(upCond2) == 0) message(sprintf("Note: No loci were differentially enriched in condition 2 ('%s')", set_labels[2]))

    biovenn_params <- list(
        list_x = Cond1_full,
        list_y = Cond2_full,
        list_z = NULL,
        xtitle = set_labels[1],
        ytitle = set_labels[2],
        ztitle = NULL,
        title = title,
        subtitle = subtitle,
        x_c = region_colours[1],
        y_c = region_colours[2],
        z_c = region_colours[3],
        t_f = font,
        st_f = font,
        xt_f = font,
        yt_f = font,
        zt_f = font,
        nr_f = font,
        t_fb = 1,
        st_fb = 1,
        xt_fb = 1,
        yt_fb = 1,
        zt_fb = 1,
        nr_fb = 1,
        t_s = 1,
        st_s = 1,
        width = 600,
        height = 600
    )

    # File handling/format
    if (!(is.null(filename))) {
        format <- match.arg(format)
        if (!grepl(paste0("\\.", format, "$"), filename, ignore.case = TRUE)) {
            filename <- paste0(tools::file_path_sans_ext(filename), ".", format)
        }

        biovenn_params$output <- format
        biovenn_params$filename <- filename
    }

    # Default colours: orange, blue, neutral
    reg_col <- rep(region_colours, length.out = 3)

    # Call BioVenn (discarding the verbose set numbers messaging from this package with capture.output())
    capture.output(
        {
           venn <- do.call(BioVenn::draw.venn, biovenn_params)
        },
        file = NULL
    )
    return(invisible(NULL))
}
