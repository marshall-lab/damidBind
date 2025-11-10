# FILE: tests/testthat/test-noiseq_functions.R

library(testthat)
library(damidBind)
library(NOISeq)
library(dplyr)

context("Differential Analysis: differential_accessibility")

# Helper to create a dummy data_list for CATaDa counts
make_dummy_catada_data_list <- function() {
    occupancy_df <- data.frame(
        name = c("peak1", "peak2", "peak3", "peak4", "peak5", "peak6"),
        chr = c("chrA", "chrA", "chrB", "chrB", "chrC", "chrC"),
        start = c(10, 110, 20, 120, 30, 130),
        end = c(100, 200, 110, 210, 120, 220),
        gene_name = c("GeneX", "GeneY", "", "GeneZ", "GeneA", "GeneB"),
        gene_id = c("ID_X", "ID_Y", "ID_none", "ID_Z", "ID_A", "ID_B"),
        nfrags = c(5, 8, 3, 6, 7, 9),
        CondA_rep1 = c(100, 50, 10, 20, 150, 30), # High in A
        CondA_rep2 = c(110, 55, 12, 22, 160, 33),
        CondB_rep1 = c(20, 200, 50, 40, 30, 100), # High in B
        CondB_rep2 = c(25, 210, 52, 45, 35, 110),
        stringsAsFactors = FALSE
    )
    rownames(occupancy_df) <- occupancy_df$name

    list(
        occupancy = occupancy_df,
        test_category = "bound" # Start with bound, function will change it to accessible
    )
}


test_that("differential_accessibility returns expected structure and values", {
    dl_catada <- make_dummy_catada_data_list()
    cond <- c("CondA", "CondB")
    cond_names <- c("Condition_Alpha", "Condition_Beta")

    res <- differential_accessibility(dl_catada, cond, cond_names, q = 0.8)

    # Check overall output structure
    expect_s4_class(res, "DamIDResults")
    expect_s3_class(analysisTable(res), "data.frame")
    expect_true(all(rownames(analysisTable(res)) == rownames(dl_catada$occupancy)))
    expect_true(all(c("logFC", "minuslogp", "gene_name", "gene_id") %in% colnames(analysisTable(res))))

    # Check gene annotations are transferred
    expect_equal(analysisTable(res)["peak1", "gene_name"], "GeneX")
    expect_equal(analysisTable(res)["peak3", "gene_name"], "")

    expect_equal(nrow(enrichedCond1(res)), 2)
    expect_true(all(c("peak1", "peak5") %in% rownames(enrichedCond1(res))))
    expect_equal(nrow(enrichedCond2(res)), 4)
    expect_true(all(c("peak2", "peak3", "peak4", "peak6") %in% rownames(enrichedCond2(res))))

    # Check calculated condition means
    expect_equal(analysisTable(res)["peak1", "CondA_mean"], mean(c(100, 110)))
    expect_equal(analysisTable(res)["peak1", "CondB_mean"], mean(c(20, 25)))

    expect_equal(analysisTable(res)["peak1", "logFC"], log2(105 / 22.5))
    # In the real run, strong evidence leads to prob=1, which correctly becomes Inf.
    expect_equal(analysisTable(res)["peak1", "minuslogp"], Inf)

    # Check 'cond' mapping
    expect_equal(conditionNames(res), c("Condition_Alpha" = "CondA", "Condition_Beta" = "CondB"))

    expect_equal(inputData(res)$test_category, "accessible")
})


test_that("differential_accessibility handles cases with no significant results", {
    dl_catada <- make_dummy_catada_data_list()
    cond <- c("CondA", "CondB")

    # Modify the input data so that there are no real differences between the conditions.
    dl_catada$occupancy$CondB_rep1 <- dl_catada$occupancy$CondA_rep1 + 1
    dl_catada$occupancy$CondB_rep2 <- dl_catada$occupancy$CondA_rep2 + 1

    # With no real difference, `noiseq` should not find anything significant.
    res <- differential_accessibility(dl_catada, cond, q = 0.8)

    expect_equal(nrow(enrichedCond1(res)), 0)
    expect_equal(nrow(enrichedCond2(res)), 0)
})
