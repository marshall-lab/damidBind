library(testthat)
library(damidBind)

context("DamIDResults S4 Methods")

#' Helper to create a DamIDResults object with FDR data for methods testing
.create_mock_fdr_damidresults <- function() {
    occupancy_df <- data.frame(
        gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE"),
        gene_id = c("FBgn01", "FBgn02", "FBgn03", "FBgn04", "FBgn05"),
        L4_rep1 = 1:5, L4_rep2 = 1:5, L5_rep1 = 1:5,
        L4_rep1_FDR = c(0.01, 0.10, 0.06, 0.50, 0.04),
        L4_rep2_FDR = c(0.03, 0.02, 0.07, 0.60, NA),
        L5_rep1_FDR = c(0.80, 0.90, 0.70, 0.01, 0.85),
        row.names = c("geneA", "geneB", "geneC", "geneD", "geneE")
    )

    # Create the analysis table subset
    analysis_df <- data.frame(
        gene_name = occupancy_df$gene_name,
        gene_id = occupancy_df$gene_id,
        row.names = rownames(occupancy_df)
    )

    new("DamIDResults",
        analysis = analysis_df,
        upCond1 = analysis_df[1, , drop = FALSE],
        upCond2 = analysis_df[4, , drop = FALSE],
        cond = c("L4 Neurons" = "L4", "L5 Neurons" = "L5"),
        data = list(
            occupancy = occupancy_df,
            test_category = "expressed",
            matched_samples = list("L4" = c("L4_rep1", "L4_rep2"), "L5" = "L5_rep1")
        )
    )
}


test_that("expressed() method correctly filters genes", {
    mock_results <- .create_mock_fdr_damidresults()

    # Test using the display name with default arguments
    res_default <- expressed(mock_results, condition = "L4 Neurons")
    expect_s3_class(res_default, "data.frame")
    expect_equal(nrow(res_default), 3)
    expect_equal(sort(res_default$gene_name), c("geneA", "geneB", "geneE"))

    # Test with which = "all"
    res_all <- expressed(mock_results, condition = "L4", which = "all")
    expect_equal(nrow(res_all), 1)
    expect_equal(res_all$gene_name, "geneA")

    # Test with a custom FDR threshold
    res_fdr_high <- expressed(mock_results, condition = "L4", fdr = 0.1)
    expect_equal(nrow(res_fdr_high), 4)

    # Test filtering on the second condition
    res_cond2 <- expressed(mock_results, condition = "L5 Neurons")
    expect_equal(nrow(res_cond2), 1)
    expect_equal(res_cond2$gene_name, "geneD")
})

test_that("Existing accessor methods complete without error", {
    mock_results <- .create_mock_fdr_damidresults()
    expect_no_error(analysisTable(mock_results))
    expect_no_error(enrichedCond1(mock_results))
    expect_no_error(enrichedCond2(mock_results))
    expect_no_error(conditionNames(mock_results))
    expect_no_error(inputData(mock_results))
})
