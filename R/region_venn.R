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
#' @return Invisibly, the output of `BioVenn::draw.venn()`.
#' @export
plot_venn <- function(
    diff_results,
    title = "Enriched binding at loci",
    subtitle = "",
    set_labels = NULL,
    filename = NULL,
    font = "sans",
    format = c("pdf", "svg"),
    region_colours = c("#FFA500", "#2288DD", "#CCCCCC")
    ) {
  # Argument and field checks
  stopifnot(is(diff_results, "DamIDResults"))

  if (is.null(set_labels))
    set_labels <- names(diff_results@cond)

  upCond1 <- rownames(diff_results@upCond1)
  upCond2 <- rownames(diff_results@upCond2)
  all_ids <- rownames(diff_results@analysis)
  # Defensive: remove NAs or empty
  upCond1 <- upCond1[!is.na(upCond1) & nchar(upCond1) > 0]
  upCond2 <- upCond2[!is.na(upCond2) & nchar(upCond2) > 0]
  all_ids <- all_ids[!is.na(all_ids) & nchar(all_ids) > 0]

  # Check minimum viable dataset
  if (length(all_ids) < 2)
    stop("Not enough loci for Venn diagram (need at least 2 in 'all').")

  # Non-significant set
  nonsig <- setdiff(all_ids, union(upCond1, upCond2))
  Cond1_full <- union(upCond1, nonsig)
  Cond2_full <- union(upCond2, nonsig)

  # Warn if one group is empty
  if (length(upCond1) == 0) warning("No loci present in upCond1; circle 1 will not have a non-overlap region.")
  if (length(upCond2) == 0) warning("No loci present in upCond2; circle 2 will not have a non-overlap region.")

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
    nr_f = font
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

  # Call BioVenn (discarding the verbose set numbers messaging from this package)
  venn <- suppressMessages(do.call(BioVenn::draw.venn, biovenn_params))
  return(invisible(NULL))
}
