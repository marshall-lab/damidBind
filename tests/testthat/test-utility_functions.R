library(testthat)
library(damidBind)

context("Utility Functions: filter_genes_by_fdr")

# Revised helper to respect parameters
.create_mock_data_for_filter_tests <- function(include_fdr_cols = TRUE) {
    occupancy_df <- data.frame(
        gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE"),
        gene_id = c("FBgn01", "FBgn02", "FBgn03", "FBgn04", "FBgn05"),
        L4_rep1 = 1:5, L4_rep2 = 1:5, L5_rep1 = 1:5,
        row.names = c("geneA", "geneB", "geneC", "geneD", "geneE")
    )

    if (include_fdr_cols) {
        occupancy_df$L4_rep1_FDR = c(0.01, 0.10, 0.06, 0.50, 0.04)
        occupancy_df$L4_rep2_FDR = c(0.03, 0.02, 0.07, 0.60, NA)
        occupancy_df$L5_rep1_FDR = c(0.80, 0.90, 0.70, 0.01, 0.85)
    }

    list(
        occupancy = occupancy_df,
        test_category = "expressed",
        matched_samples = list("L4" = c("L4_rep1", "L4_rep2"), "L5" = "L5_rep1")
    )
}



test_that("filter_genes_by_fdr works correctly with list input", {
    mock_data <- .create_mock_data_for_filter_tests()

    # Test 'any' logic
    res_any <- filter_genes_by_fdr(mock_data, fdr = 0.05, condition = "L4", which = "any")
    expect_s3_class(res_any, "data.frame")
    expect_equal(nrow(res_any), 3)
    expect_equal(sort(res_any$gene_id), sort(c("FBgn01", "FBgn02", "FBgn05")))

    # Test 'all' logic
    res_all <- filter_genes_by_fdr(mock_data, fdr = 0.05, condition = "L4", which = "all")
    expect_s3_class(res_all, "data.frame")
    expect_equal(nrow(res_all), 1)
    expect_equal(res_all$gene_id, "FBgn01")

    # Test filtering on the second condition
    res_l5 <- filter_genes_by_fdr(mock_data, fdr = 0.05, condition = "L5", which = "any")
    expect_s3_class(res_l5, "data.frame")
    expect_equal(nrow(res_l5), 1)
    expect_equal(res_l5$gene_id, "FBgn04")

    # Test case where no genes pass the filter
    res_none <- filter_genes_by_fdr(mock_data, fdr = 0.001, condition = "L4", which = "any")
    expect_s3_class(res_none, "data.frame")
    expect_equal(nrow(res_none), 0)
    expect_named(res_none, c("gene_name", "gene_id", "avg_occ", "fdr_val"))

})


test_that("filter_genes_by_fdr works correctly with DamIDResults object input", {
    mock_list_data <- .create_mock_data_for_filter_tests()

    # Create a simple DamIDResults object from the mock list data
    mock_damid_results <- new("DamIDResults",
                              analysis = data.frame(gene_name = "dummy", gene_id = "dummy"),
                              upCond1 = data.frame(),
                              upCond2 = data.frame(),
                              cond = c("L4 Neurons" = "L4", "L5 Neurons" = "L5"),
                              data = mock_list_data
    )

    # Test using the display name for the condition
    res_any <- filter_genes_by_fdr(mock_damid_results, fdr = 0.05, condition = "L4 Neurons", which = "any")
    expect_equal(nrow(res_any), 3)
    expect_equal(sort(res_any$gene_id), sort(c("FBgn01", "FBgn02", "FBgn05")))

    # Test using the internal identifier
    res_all <- filter_genes_by_fdr(mock_damid_results, fdr = 0.05, condition = "L4", which = "all")
    expect_equal(nrow(res_all), 1)
    expect_equal(res_all$gene_id, "FBgn01")
})


test_that("filter_genes_by_fdr handles errors and warnings correctly", {
    mock_data <- .create_mock_data_for_filter_tests()

    # Data with no FDR columns
    data_no_fdr <- .create_mock_data_for_filter_tests(include_fdr_cols = FALSE)
    expect_warning(
        res <- filter_genes_by_fdr(data_no_fdr, fdr = 0.05, "L4"),
        regexp = "No '_FDR' columns found in the data"
    )
    expect_equal(nrow(res), 0)

    # Condition not found
    expect_warning(
        res <- filter_genes_by_fdr(mock_data, fdr = 0.05, "NonExistentCond"),
        regexp = "No '_FDR' columns found in the data matching the condition 'NonExistentCond'"
    )
    expect_equal(nrow(res), 0)
})
