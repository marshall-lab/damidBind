#' The DamIDResults Class
#'
#' @description
#' An S4 class to store the results of a differential analysis, as generated
#' by \code{\link{differential_binding}} or \code{\link{differential_accessibility}}. It contains
#' the full analysis table, subsets of significantly changed regions,
#' and associated metadata.
#'
#' @slot analysis A `data.frame` containing the full differential analysis table from
#'   `limma` or `NOISeq`.
#' @slot upCond1 A `data.frame` of regions significantly enriched in the first condition.
#' @slot upCond2 A `data.frame` of regions significantly enriched in the second condition.
#' @slot cond A named `character` vector mapping user-friendly display names (the names)
#'   to the internal condition identifiers (the values) used in the analysis.
#' @slot data A `list` containing the initial input data used for the analysis,
#'   including the occupancy `data.frame` and other metadata.
#'
#' @section Accessor Methods:
#' The following accessor functions are available for a \code{DamIDResults} object.
#' \itemize{
#'   \item \code{\link{analysisTable}(object)}: Returns the full differential analysis table (a \code{data.frame}).
#'   \item \code{\link{enrichedCond1}(object)}: Returns a \code{data.frame} of regions significantly enriched in the first condition.
#'   \item \code{\link{enrichedCond2}(object)}: Returns a \code{data.frame} of regions significantly enriched in the second condition.
#'   \item \code{\link{conditionNames}(object)}: Returns a named \code{character} vector mapping display names to internal identifiers.
#'   \item \code{\link{inputData}(object)}: Returns a \code{list} containing the original input data used for the analysis.
#'   \item \code{\link{expressed}(object, condition, fdr = 0.05, which = "any")}: Returns a \code{data.frame} of genes considered expressed in `condition`, based on an FDR threshold of significantly enriched occupancy.  Only available for analyses with FDR calculations, generated via \code{load_data_genes(calculate_fdr = TRUE)}.
#' }
#'
#' @section Generic Methods:
#' The generic \code{plot()} function is also S4-enabled for this class.
#' Calling \code{plot(object)} is equivalent to calling
#' \code{plot_volcano(diff_results = object)}.
#'
#' @seealso
#' For more powerful and specific plotting functions,
#' see \code{\link{plot_volcano}}, \code{\link{plot_venn}},
#' and \code{\link{analyse_go_enrichment}}.
#'
#' To explore the differential analysis results in an interactive IGV browser
#' window, see \code{\link{browse_igv_regions}}
#'
#' @name DamIDResults-class
#' @aliases DamIDResults
#' @docType class
#' @keywords classes
#' @exportClass DamIDResults
#'
#' @examples
#' # Helper function to create a sample DamIDResults object for examples
#' .generate_example_results <- function() {
#'     analysis_df <- data.frame(
#'         logFC = c(2, -2, 0.1), P.Value = c(0.01, 0.01, 0.9), B = c(4, 3, -1),
#'         gene_name = c("GeneA", "GeneB", "GeneC"),
#'         row.names = c("chr1:1-100", "chr1:101-200", "chr1:201-300")
#'     )
#'     new("DamIDResults",
#'         analysis = analysis_df,
#'         upCond1 = analysis_df[1, , drop = FALSE],
#'         upCond2 = analysis_df[2, , drop = FALSE],
#'         cond = c("Condition 1" = "C1", "Condition 2" = "C2"),
#'         data = list(test_category = "bound")
#'     )
#' }
#' mock_results <- .generate_example_results()
#'
#' # Show the object summary
#' mock_results
#'
#' # Access different parts of the object
#' analysisTable(mock_results)
#' enrichedCond1(mock_results)
#' conditionNames(mock_results)
NULL

#' @rdname DamIDResults-class
setClass("DamIDResults",
         slots = c(
             analysis = "data.frame",
             upCond1 = "data.frame",
             upCond2 = "data.frame",
             cond = "character",
             data = "list"
         )
)
setMethod(
    "show", "DamIDResults",
    function(object) {
        cat("An object of class 'DamIDResults'\n")
        cond_display <- names(object@cond)
        cat(sprintf("Differentially %s regions\n", object@data$test_category))
        cat(sprintf("Comparison: '%s' vs '%s'\n", cond_display[1], cond_display[2]))
        cat(sprintf("- %d regions enriched in %s\n", nrow(object@upCond1), cond_display[1]))
        cat(sprintf("- %d regions enriched in %s\n", nrow(object@upCond2), cond_display[2]))
        cat(sprintf("- %d total regions tested\n", nrow(object@analysis)))
        invisible(object)
    }
)
setMethod(
    "plot", signature(x = "DamIDResults", y = "missing"),
    function(x, y, ...) {
        # This method simply calls the main volcano plot function by default.
        plot_volcano(diff_results = x, ...)
    }
)
