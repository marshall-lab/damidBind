#' @title Access the differential binding analysis results
#' @description This function returns the full differential analysis table from the
#'   \code{DamIDResults} object.
#' @param object A \code{DamIDResults} object.
#' @return A \code{data.frame} with the full analysis results.
#' @seealso \code{\link{DamIDResults-class}} for an overview of the class and all its methods.
#' @aliases analysisTable,DamIDResults-method
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
#' analysisTable(mock_results)
#' @export
setGeneric("analysisTable", function(object) standardGeneric("analysisTable"))
setMethod("analysisTable", "DamIDResults", function(object) object@analysis)

#' @title Access Condition 1 enriched regions
#' @description This function returns the subset of regions significantly enriched in the first condition.
#' @param object A \code{DamIDResults} object.
#' @return A \code{data.frame}.
#' @seealso \code{\link{DamIDResults-class}}
#' @aliases enrichedCond1,DamIDResults-method
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
#' enrichedCond1(mock_results)
#' @export
setGeneric("enrichedCond1", function(object) standardGeneric("enrichedCond1"))
setMethod("enrichedCond1", "DamIDResults", function(object) object@upCond1)

#' @title Access Condition 2 enriched regions
#' @description This function returns the subset of regions significantly enriched in the second condition.
#' @param object A \code{DamIDResults} object.
#' @return A \code{data.frame}.
#' @seealso \code{\link{DamIDResults-class}}
#' @aliases enrichedCond2,DamIDResults-method
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
#' enrichedCond2(mock_results)
#' @export
setGeneric("enrichedCond2", function(object) standardGeneric("enrichedCond2"))
setMethod("enrichedCond2", "DamIDResults", function(object) object@upCond2)

#' @title Access condition name mapping
#' @description This function returns the mapping of user-friendly display names to internal condition identifiers.
#' @param object A \code{DamIDResults} object.
#' @return A named \code{character} vector.
#' @seealso \code{\link{DamIDResults-class}}
#' @aliases conditionNames,DamIDResults-method
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
#' conditionNames(mock_results)
#' @export
setGeneric("conditionNames", function(object) standardGeneric("conditionNames"))
setMethod("conditionNames", "DamIDResults", function(object) object@cond)

#' @title Access original input data and metadata
#' @description This function returns the original list of input data used to generate the results.
#' @param object A \code{DamIDResults} object.
#' @return A \code{list}.
#' @seealso \code{\link{DamIDResults-class}}
#' @aliases inputData,DamIDResults-method
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
#' inputData(mock_results)
#' @export
setGeneric("inputData", function(object) standardGeneric("inputData"))
setMethod("inputData", "DamIDResults", function(object) object@data)

#' @title Get expressed genes/loci by FDR
#' @description A method to filter genes or loci that are considered 'expressed'
#'   in a specific condition, based on a False Discovery Rate (FDR) threshold.
#'   The method is a wrapper around the \code{\link{filter_genes_by_fdr}} function.
#'
#' @param object A \code{DamIDResults} object. This object must have been generated
#'   from data loaded with \code{load_data_genes(calculate_fdr = TRUE)} for the
#'   underlying FDR columns to be present.
#' @param condition A character string identifying the experimental condition to filter.
#'   This can be the internal identifier or the user-friendly display name.
#' @param fdr Numeric. The FDR cutoff. Defaults to 0.05.
#' @param which Character string, either 'any' or 'all'. Determines whether a gene
#'   must pass the FDR threshold in any or all replicates of the condition.
#'   Defaults to 'any'.
#'
#' @return A \code{data.frame} containing the \code{gene_name} and \code{gene_id} of
#'   genes that pass the filter.
#'
#' @seealso \code{\link{filter_genes_by_fdr}}, \code{\link{DamIDResults-class}},
#'   \code{\link{load_data_genes}}
#' @aliases expressed,DamIDResults-method
#'
#' @examples
#' # Helper function to create a sample DamIDResults object with FDR data
#' .generate_fdr_example_results <- function() {
#'     occupancy_df <- data.frame(
#'         gene_name = c("geneA", "geneB", "geneC"),
#'         gene_id = c("FBgn01", "FBgn02", "FBgn03"),
#'         L4_rep1_FDR = c(0.01, 0.10, 0.04),
#'         L4_rep2_FDR = c(0.03, 0.02, 0.50),
#'         L5_rep1_FDR = c(0.80, 0.90, 0.01),
#'         row.names = c("geneA", "geneB", "geneC")
#'     )
#'     diff_results_base <- list(occupancy = occupancy_df, test_category = "expressed")
#'     new("DamIDResults",
#'         analysis = data.frame(row.names = rownames(occupancy_df)),
#'         upCond1 = data.frame(), upCond2 = data.frame(),
#'         cond = c("L4 Neurons" = "L4", "L5 Neurons" = "L5"),
#'         data = diff_results_base
#'     )
#' }
#' mock_fdr_results <- .generate_fdr_example_results()
#'
#' # Get genes expressed in L4 neurons (FDR <= 0.05 in any replicate)
#' expressed(mock_fdr_results, condition = "L4 Neurons")
#'
#' # Get genes expressed in L5 neurons with a stricter fdr
#' expressed(mock_fdr_results, condition = "L5", fdr = 0.02)
#'
#' @export
setGeneric("expressed", function(object, condition, fdr = 0.05, which = "any")
    standardGeneric("expressed"))

setMethod("expressed", "DamIDResults",
          function(object, condition, fdr = 0.05, which = "any") {
              filter_genes_by_fdr(data = object,
                                  condition = condition,
                                  fdr = fdr,
                                  which = which)
          }
)
