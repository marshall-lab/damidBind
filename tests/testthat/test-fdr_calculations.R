library(testthat)
library(damidBind)
library(GenomicRanges)
library(IRanges)

context("FDR Calculations")

# Helper to prepare minimal data for FDR tests.
# This uses a real data file to ensure the test is realistic.
prepare_fdr_test_data <- function() {
    # Load one of the sample files included with the package
    data_dir <- system.file("extdata", package = "damidBind")
    bgraph_file <- file.path(data_dir, "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm.gatc.2L.bedgraph.gz")

    # Use the internal import function
    binding_df <- damidBind:::import_bedgraph_as_df(bgraph_file, colname = "L4_r1")

    # The function expects a GRanges object for binding data.
    binding_gr <- GRanges(
        seqnames = binding_df$chr,
        ranges = IRanges(start = binding_df$start, end = binding_df$end),
        L4_r1 = binding_df$L4_r1
    )

    # A mock occupancy table. The function will calculate FDR for `L4_r1`.
    occupancy_table <- data.frame(
        name = c("geneA", "geneB", "geneC", "geneD"),
        nfrags = c(5, 10, 2, 20),
        L4_r1 = c(1.8, 0.9, -0.1, 2.5),
        row.names = c("geneA", "geneB", "geneC", "geneD")
    )

    list(binding_gr = binding_gr, occupancy_table = occupancy_table)
}


test_that("FDR calculation is reproducible with a seed and not without", {
    # This test requires BiocParallel
    skip_if_not_installed("BiocParallel")
    test_data <- prepare_fdr_test_data()

    # Run 1: with seed 42
    res1 <- calculate_and_add_fdr(
        binding_data = test_data$binding_gr,
        occupancy_df = test_data$occupancy_table,
        fdr_iterations = 500, # Low iterations for test speed
        seed = 42,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Run 2: with the same seed 42
    res2 <- calculate_and_add_fdr(
        binding_data = test_data$binding_gr,
        occupancy_df = test_data$occupancy_table,
        fdr_iterations = 500,
        seed = 42,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Run 3: with a different seed
    res3 <- calculate_and_add_fdr(
        binding_data = test_data$binding_gr,
        occupancy_df = test_data$occupancy_table,
        fdr_iterations = 500,
        seed = 101,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Run 4: with no seed (should not be identical to Run 5)
    res4 <- calculate_and_add_fdr(
        binding_data = test_data$binding_gr,
        occupancy_df = test_data$occupancy_table,
        fdr_iterations = 500,
        seed = NULL,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Run 5: with no seed again
    res5 <- calculate_and_add_fdr(
        binding_data = test_data$binding_gr,
        occupancy_df = test_data$occupancy_table,
        fdr_iterations = 500,
        seed = NULL,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Reproducibility with seed
    expect_identical(res1, res2, info = "Results with the same seed should be identical.")

    # Seed effectiveness
    expect_false(identical(res1, res3), info = "Results with different seeds should not be identical.")

    # Non-reproducibility without seed
    expect_false(identical(res4, res5), info = "Results without a seed should generally not be identical.")
})


test_that("FDR calculation produces correct output structure and values", {
    # This test requires BiocParallel
    skip_if_not_installed("BiocParallel")
    test_data <- prepare_fdr_test_data()

    res <- calculate_and_add_fdr(
        binding_data = test_data$binding_gr,
        occupancy_df = test_data$occupancy_table,
        fdr_iterations = 50,
        seed = 123, # Use a fixed seed for predictable output
        BPPARAM = BiocParallel::SerialParam()
    )

    # Check output structure
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), nrow(test_data$occupancy_table))
    expect_true("L4_r1_FDR" %in% colnames(res))

    # Check value constraints
    fdr_values <- res$L4_r1_FDR
    expect_type(fdr_values, "double")
    expect_true(all(fdr_values >= 0 & fdr_values <= 1), info = "FDR values should be between 0 and 1.")
})
