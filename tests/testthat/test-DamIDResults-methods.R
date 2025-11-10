library(testthat)
library(damidBind)

context("DamIDResults S4 Methods")

# Helper function to create mock DamIDResults with FDR data
.create_mock_fdr_damidresults <- function() {
    occupancy_df <- data.frame(
        gene_name = c("geneA", "geneB", "geneC", "geneD"),
        gene_id = c("FBgn01", "FBgn02", "FBgn03", "FBgn04"),
        L4_rep1_FDR = c(0.01, 0.10, 0.04, 0.06),
        L4_rep2_FDR = c(0.03, 0.02, 0.50, 0.07),
        L5_rep1_FDR = c(0.80, 0.90, 0.01, 0.02),
        row.names = c("geneA", "geneB", "geneC", "geneD")
    )
    diff_results_base <- list(occupancy = occupancy_df, test_category = "expressed")
    new("DamIDResults",
        analysis = data.frame(row.names = rownames(occupancy_df)),
        upCond1 = data.frame(), upCond2 = data.frame(),
        cond = c("L4 Neurons" = "L4", "L5 Neurons" = "L5"),
        data = diff_results_base
    )
}

test_that("expressed() method correctly filters genes", {
    mock_results <- .create_mock_fdr_damidresults()

    # Test using the display name with default arguments (fdr=0.05, which="any")
    res_default <- expressed(mock_results, condition = "L4 Neurons")
    expect_s3_class(res_default, "data.frame")
    expect_equal(nrow(res_default), 3)
    expect_equal(sort(res_default$gene_name), c("geneA", "geneB", "geneC"))

    # Test with which = "all"
    res_all <- expressed(mock_results, condition = "L4", which = "all")
    expect_s3_class(res_all, "data.frame")
    expect_equal(nrow(res_all), 1)
    expect_equal(res_all$gene_name, "geneA")

    # Test with a custom FDR threshold
    res_fdr_high <- expressed(mock_results, condition = "L4", fdr = 0.1)
    expect_equal(nrow(res_fdr_high), 4)

    # Test filtering on the second condition
    res_cond2 <- expressed(mock_results, condition = "L5 Neurons")
    expect_equal(nrow(res_cond2), 2)
    expect_equal(sort(res_cond2$gene_name), c("geneC", "geneD"))
})

# Regression test for existing accessors
test_that("Existing accessor methods complete without error", {
    mock_results <- .create_mock_fdr_damidresults()
    expect_no_error(analysisTable(mock_results))
    expect_no_error(enrichedCond1(mock_results))
    expect_no_error(enrichedCond2(mock_results))
    expect_no_error(conditionNames(mock_results))
    expect_no_error(inputData(mock_results))
})
