library(testthat)
library(damidBind)

context("differential_binding")

# Helper to create a minimal, complete data_list for differential analysis testing
make_minimal_data_list <- function() {
    occupancy_df <- data.frame(
        name = c("loc1", "loc2", "loc3", "loc4"),
        chr = c("chrA", "chrA", "chrB", "chrB"),
        start = c(1, 101, 1, 101),
        end = c(100, 200, 100, 200),
        gene_name = c("GeneA", "GeneB", "", "GeneC"),
        gene_id = c("ID_A", "ID_B", "ID_none", "ID_C"),
        CondA_rep1 = c(10, 50, 10, 20),
        CondA_rep2 = c(12, 52, 11, 22),
        CondB_rep1 = c(5, 10, 20, 40),
        CondB_rep2 = c(6, 12, 21, 41),
        stringsAsFactors = FALSE
    )
    rownames(occupancy_df) <- occupancy_df$name
    list(occupancy = occupancy_df, test_category = "bound")
}

test_that("differential_binding processes data and returns correct structure", {
    dl <- make_minimal_data_list()
    # Define the new named vector for conditions
    cond_vec <- c("CondA_display" = "CondA", "CondB_display" = "CondB")

    # Test against the real, deterministic output of limma
    res <- suppressMessages(differential_binding(dl, cond = cond_vec, fdr = 0.05))

    # Check overall object and slot structure
    expect_s4_class(res, "DamIDResults")
    expect_s3_class(analysisTable(res), "data.frame")
    expect_true(all(c("upCond1", "upCond2", "analysis", "cond", "data") %in% slotNames(res)))

    # Check that the analysis table contains all expected loci, regardless of order
    expected_loci <- c("loc1", "loc2", "loc3", "loc4")
    expect_equal(sort(rownames(analysisTable(res))), sort(expected_loci))
    expect_equal(nrow(analysisTable(res)), 4)

    # Check essential columns are present and correctly derived
    # Constructing internal mean names based on new logic (sanitised display names)
    internal_mean_names <- paste0(make.names(names(cond_vec)), "_mean")
    expect_true(all(c("logFC", "adj.P.Val", "minuslogp", "gene_name", "gene_id", internal_mean_names) %in% colnames(analysisTable(res))))
    # Check that minuslogp is correctly calculated from the adj.P.Val in the results of this specific run
    expect_equal(analysisTable(res)$minuslogp, -log10(analysisTable(res)$adj.P.Val))

    # Check that significant sets contain the correct members, ignoring order
    # In this dataset, all 4 loci should be significant at FDR < 0.05
    expect_equal(nrow(enrichedCond1(res)), 2)
    expect_equal(sort(rownames(enrichedCond1(res))), sort(c("loc1", "loc2")))
    expect_equal(nrow(enrichedCond2(res)), 2)
    expect_equal(sort(rownames(enrichedCond2(res))), sort(c("loc3", "loc4")))

    # Check gene annotations and condition means
    expect_equal(analysisTable(res)["loc1", "gene_name"], "GeneA")
    expect_equal(analysisTable(res)["loc3", "gene_name"], "")
    expect_equal(analysisTable(res)["loc2", internal_mean_names[1]], mean(c(50, 52)))
    expect_equal(analysisTable(res)["loc4", internal_mean_names[2]], mean(c(40, 41)))

    # Check cond mapping
    expect_equal(conditionNames(res), setNames(c("CondA", "CondB"), c("CondA_display", "CondB_display")))
})


test_that("differential_binding handles cases with no significant results", {
    # Create a data list where there is no real difference between conditions
    dl_no_diff <- make_minimal_data_list()
    dl_no_diff$occupancy$CondB_rep1 <- dl_no_diff$occupancy$CondA_rep1
    dl_no_diff$occupancy$CondB_rep2 <- dl_no_diff$occupancy$CondA_rep2

    cond_vec <- c("CondA" = "CondA", "CondB" = "CondB")

    # Run the function. The real limma pipeline should find no significant differences.
    res <- suppressMessages(differential_binding(dl_no_diff, cond = cond_vec, fdr = 0.05))

    # With no real difference, upCond1 and upCond2 should be empty data frames.
    expect_equal(nrow(enrichedCond1(res)), 0)
    expect_equal(nrow(enrichedCond2(res)), 0)
})


test_that("differential_binding can optionally filter for positive enrichment", {
    # Mock occupancy data:
    # 'neg_sig': significant difference, both means negative
    # 'pos_sig': significant difference, both means positive (standard case)
    # 'non_sig': not significant
    occupancy_df <- data.frame(
        name = c("neg_sig", "pos_sig", "non_sig"),
        gene_name = c("GeneNeg", "GenePos", "GeneNon"),
        gene_id = c("ID_Neg", "ID_Pos", "ID_Non"),
        CondA_rep1 = c(-0.9, 3.1, 1.0),
        CondA_rep2 = c(-1.1, 2.9, 1.1),
        CondB_rep1 = c(-2.9, 1.1, 1.0),
        CondB_rep2 = c(-3.1, 0.9, 0.9),
        stringsAsFactors = FALSE
    )
    rownames(occupancy_df) <- occupancy_df$name
    dl <- list(occupancy = occupancy_df, test_category = "bound")
    cond_vec <- c("CondA_name" = "CondA", "CondB_name" = "CondB")
    internal_mean_names <- paste0(make.names(names(cond_vec)), "_mean")

    # Test condition without filter
    # Both 'neg_sig' and 'pos_sig' should be in upCond1 (logFC > 0).
    res_unfiltered <- suppressMessages(
        differential_binding(dl, cond = cond_vec, fdr = 0.1,
                             filter_positive_enrichment = FALSE,
                             filter_negative_occupancy = NULL) # Disable pre-filtering
    )

    expect_s4_class(res_unfiltered, "DamIDResults")

    # Check that both significant loci are found when not filtering.
    expect_true("neg_sig" %in% rownames(enrichedCond1(res_unfiltered)))
    expect_true("pos_sig" %in% rownames(enrichedCond1(res_unfiltered)))
    expect_equal(nrow(enrichedCond1(res_unfiltered)), 2)
    expect_equal(nrow(enrichedCond2(res_unfiltered)), 0)

    # Test with the filter (default)
    res_filtered <- suppressMessages(
        differential_binding(dl, cond = cond_vec, fdr = 0.1,
                             filter_positive_enrichment = TRUE,
                             filter_negative_occupancy = NULL) # Disable pre-filtering
    )

    expect_s4_class(res_filtered, "DamIDResults")
    # 'neg_sig' should now be filtered out from the final "enriched" set.
    expect_false("neg_sig" %in% rownames(enrichedCond1(res_filtered)))
    # 'pos_sig' should remain.
    expect_true("pos_sig" %in% rownames(enrichedCond1(res_filtered)))
    expect_equal(nrow(enrichedCond1(res_filtered)), 1)
    expect_equal(nrow(enrichedCond2(res_filtered)), 0)

    # Verify 'neg_sig' is still present in the full analysis table.
    expect_true("neg_sig" %in% rownames(analysisTable(res_filtered)))
    # Double-check that the means were negative.
    expect_lt(analysisTable(res_filtered)["neg_sig", internal_mean_names[1]], 0)
    expect_lt(analysisTable(res_filtered)["neg_sig", internal_mean_names[2]], 0)
})
