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
        gene_name = c("GeneA", "GeneB", "GeneC", "GeneD"),
        gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004"),
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
    # Define new named vector format for conditions
    cond_display_names <- c("L4 neurons", "L5 neurons")
    cond_patterns <- c("L4", "L5")
    cond_vec <- setNames(cond_patterns, cond_display_names)

    res <- damidBind:::prep_data_for_differential_analysis(dl, cond = cond_vec)

    # Check output structure and types
    expect_type(res, "list")
    expect_s3_class(res$mat, "data.frame")
    expect_s3_class(res$factors, "data.frame")
    expect_type(res$cond_internal, "character")
    expect_type(res$cond_display, "character")
    expect_equal(res$occupancy_df, dl$occupancy)
    expect_equal(res$data_list$test_category, "bound")

    # Check 'mat' content
    expect_equal(colnames(res$mat), c("L4_rep1", "L4_rep2", "L5_rep1", "L5_rep2"))
    expect_equal(rownames(res$mat), dl$occupancy$name)

    # Check 'factors' content
    expect_equal(rownames(res$factors), colnames(res$mat))
    expect_s3_class(res$factors$condition, "factor")
    internal_cond_names <- make.names(cond_display_names)
    expect_equal(levels(res$factors$condition), internal_cond_names)
    expect_equal(as.character(res$factors$condition), rep(internal_cond_names, each = 2))

    # Check condition names
    expect_equal(res$cond_internal, internal_cond_names)
    expect_equal(res$cond_display, cond_display_names)
    expect_equal(res$cond_matches, cond_patterns)
})

test_that("prep_data_for_differential_analysis handles default display names", {
    dl <- make_dummy_data_list()
    cond_unnamed <- c("L4", "L5")

    res <- damidBind:::prep_data_for_differential_analysis(dl, cond = cond_unnamed)
    # When unnamed, display names should default to the patterns themselves
    expect_equal(res$cond_display, cond_unnamed)
})

test_that("prep_data_for_differential_analysis throws error for invalid 'cond' length", {
    dl <- make_dummy_data_list()
    expect_error(damidBind:::prep_data_for_differential_analysis(dl, "L4"), "`cond` must be a character vector of two unique strings.")
    expect_error(damidBind:::prep_data_for_differential_analysis(dl, c("L4", "L5", "L6")), "`cond` must be a character vector of two unique strings.")
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


test_that("prep_data_for_differential_analysis warns if not all samples assigned", {
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
    cond_vec <- c("L4 neurons" = "L4", "L5 neurons" = "L5")

    expect_message(
        damidBind:::prep_data_for_differential_analysis(dl, cond = cond_vec),
        "Differential analysis setup"
    )
    expect_message(
        damidBind:::prep_data_for_differential_analysis(dl, cond = cond_vec),
        "Found 2 replicates"
    )
})

test_that("filter_occupancy works as expected", {
    # Corrected data.frame creation
    occupancy_df_filter_test <- data.frame(
        name = paste0("locus_", 1:7),
        chr = "chr1",
        start = seq(1, 601, by = 100),
        end = seq(100, 700, by = 100),
        L4_rep1 = c(10, 10,  0, 10,  0, 10, 0),
        L4_rep2 = c(10, 10, 10,  0,  0, 10, 0),
        L5_rep1 = c(10,  0, 10, 10, 10, 10, 0),
        L5_rep2 = c(10, 10, 10,  0,  0, 10, 0),
        stringsAsFactors = FALSE
    )

    dl_filter <- list(occupancy = occupancy_df_filter_test, test_category = "test")
    cond_vec <- c("L4", "L5")

    # Test with default filter (TRUE, which means minimum of 2)
    res_filtered <- expect_message(
        damidBind:::prep_data_for_differential_analysis(
            dl_filter,
            cond_vec,
            filter_occupancy = TRUE
        ),
        regexp = "Filtered out 3 loci. 4 loci remain for analysis."
    )

    # Check that the matrix and occupancy_df were correctly filtered
    expect_equal(nrow(res_filtered$mat), 4)
    expect_equal(nrow(res_filtered$occupancy_df), 4)

    # Check that the correct rows were kept
    expected_rows_kept <- c("locus_1", "locus_2", "locus_3", "locus_6")
    expect_equal(res_filtered$occupancy_df$name, expected_rows_kept)


    # Test when filter is explicitly off (NULL or FALSE)
    res_unfiltered <- damidBind:::prep_data_for_differential_analysis(
            dl_filter,
            cond_vec,
            filter_occupancy = NULL
        )

    # Check that nothing was filtered
    expect_equal(nrow(res_unfiltered$mat), 7)
    expect_equal(nrow(res_unfiltered$occupancy_df), 7)
    expect_equal(res_unfiltered$occupancy_df$name, occupancy_df_filter_test$name)


    # Test edge case where a condition has only one sample
    dl_single_rep <- dl_filter
    dl_single_rep$occupancy <- dl_single_rep$occupancy[, !colnames(dl_single_rep$occupancy) %in% "L4_rep2"]

    res_single_rep <- expect_message(
        damidBind:::prep_data_for_differential_analysis(
            dl_single_rep,
            cond_vec,
            filter_occupancy = TRUE
        ),
        regexp = "Filtered out 1 loci. 6 loci remain for analysis."
    )

    expected_rows_kept_single_rep <- c("locus_1", "locus_2", "locus_3", "locus_4", "locus_5", "locus_6")
    expect_equal(nrow(res_single_rep$mat), 6)
    expect_equal(res_single_rep$occupancy_df$name, expected_rows_kept_single_rep)
})
