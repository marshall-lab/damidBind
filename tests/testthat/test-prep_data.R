library(testthat)
library(damidBind)

context("prep_data_for_differential_analysis")

# Create a dummy data_list for testing
make_dummy_data_list <- function() {
    occupancy_df <- data.frame(
        name = c("chr1:1-100", "chr1:200-300", "chr2:50-150", "chr1:400-500"),
        chr = c("chr1", "chr1", "chr2", "chr1"),
        start = c(1, 200, 50, 400),
        end = c(100, 300, 150, 500),
        gene_names = c("GeneA", "GeneB", "GeneC", "GeneD"),
        gene_ids = c("FBgn001", "FBgn002", "FBgn003", "FBgn004"),
        nfrags = c(5, 8, 3, 6),
        L4_rep1 = c(10, 20, 15, 25),
        L4_rep2 = c(12, 21, 14, 26),
        L5_rep1 = c(5, 10, 7, 12),
        L5_rep2 = c(6, 11, 8, 13),
        stringsAsFactors = FALSE
    )
    list(occupancy = occupancy_df, test_category = "bound")
}

test_that("prep_data_for_differential_analysis prepares data correctly for valid input", {
    dl <- make_dummy_data_list()
    cond <- c("L4", "L5")
    cond_names <- c("L4 neurons", "L5 neurons")

    res <- damidBind:::prep_data_for_differential_analysis(dl, cond, cond_names)

    # Check output structure and types
    expect_type(res, "list")
    expect_s3_class(res$mat, "data.frame")
    expect_s3_class(res$factors, "data.frame")
    expect_type(res$cond_internal, "character")
    expect_type(res$cond_display, "character")
    expect_equal(res$occupancy_df, dl$occupancy)
    expect_equal(res$test_category, "bound")

    # Check 'mat' content
    expect_equal(colnames(res$mat), c("L4_rep1", "L4_rep2", "L5_rep1", "L5_rep2"))
    expect_equal(rownames(res$mat), dl$occupancy$name)

    # Check 'factors' content
    expect_equal(rownames(res$factors), colnames(res$mat))
    expect_s3_class(res$factors$condition, "factor")
    expect_equal(levels(res$factors$condition), cond)
    expect_equal(as.character(res$factors$condition), c("L4", "L4", "L5", "L5"))

    # Check condition names
    expect_equal(res$cond_internal, cond)
    expect_equal(res$cond_display, cond_names)
})

test_that("prep_data_for_differential_analysis handles default cond_names", {
    dl <- make_dummy_data_list()
    cond <- c("L4", "L5")

    res <- damidBind:::prep_data_for_differential_analysis(dl, cond, cond_names = NULL)
    expect_equal(res$cond_display, cond) # Should default to internal cond names
})

test_that("prep_data_for_differential_analysis throws error for invalid 'cond' length", {
    dl <- make_dummy_data_list()
    expect_error(damidBind:::prep_data_for_differential_analysis(dl, "L4"), "`cond` must be a character vector of exactly two strings")
    expect_error(damidBind:::prep_data_for_differential_analysis(dl, c("L4", "L5", "L6")), "`cond` must be a character vector of exactly two strings")
})

test_that("prep_data_for_differential_analysis throws error for overlapping samples", {
    dl_overlap <- make_dummy_data_list()
    # Create a sample name that will be matched by both conditions "L4" and "L5"
    colnames(dl_overlap$occupancy)[colnames(dl_overlap$occupancy) == "L4_rep1"] <- "L4_and_L5_rep1"

    expect_error(
        damidBind:::prep_data_for_differential_analysis(dl_overlap, cond = c("L4", "L5")),
        "Conditions 'L4' and 'L5' overlap"
    )
})

test_that("prep_data_for_differential_analysis throws error if samples not found", {
    dl <- make_dummy_data_list()
    cond_missing <- c("NonExistent1", "NonExistent2")

    expect_error(
        expect_warning(
            damidBind:::prep_data_for_differential_analysis(dl, cond_missing),
            "The following samples could not be assigned to either condition"
        ),
        "Fewer than two sample columns matched"
    )
})




test_that("prep_data_for_differential_analysis throws error if not all samples assigned", {
    dl_unassigned <- make_dummy_data_list()
    # Add a sample that doesn't match either cond
    dl_unassigned$occupancy$Unassigned_Sample <- c(100, 101, 102, 103)
    expect_warning(
        damidBind:::prep_data_for_differential_analysis(dl_unassigned, cond = c("L4", "L5")),
        "The following samples could not be assigned to either condition"
    )
})

test_that("Function displays appropriate messages", {
    dl <- make_dummy_data_list()
    cond <- c("L4", "L5")
    cond_names <- c("L4 neurons", "L5 neurons")

    expect_message(
        damidBind:::prep_data_for_differential_analysis(dl, cond, cond_names),
        "Differential analysis setup"
    )
    expect_message(
        damidBind:::prep_data_for_differential_analysis(dl, cond, cond_names),
        "Found 2 replicates"
    )
})
