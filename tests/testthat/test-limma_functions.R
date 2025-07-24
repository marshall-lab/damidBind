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
    gene_names = c("GeneA", "GeneB", "", "GeneC"),
    gene_ids = c("ID_A", "ID_B", "ID_none", "ID_C"),
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
  cond <- c("CondA", "CondB")
  cond_names <- c("CondA_display", "CondB_display")

  # Test against the real, deterministic output of limma
  res <- suppressMessages(differential_binding(dl, cond, cond_names, fdr = 0.05))

  # Check overall object and slot structure
  expect_s4_class(res, "DamIDResults")
  expect_s3_class(res@analysis, "data.frame")
  expect_true(all(c("upCond1", "upCond2", "analysis", "cond", "data") %in% slotNames(res)))

  # Check that the analysis table contains all expected loci, regardless of order
  expected_loci <- c("loc1", "loc2", "loc3", "loc4")
  expect_equal(sort(rownames(res@analysis)), sort(expected_loci))
  expect_equal(nrow(res@analysis), 4)

  #   Check essential columns are present and correctly derived
  expect_true(all(c("logFC", "adj.P.Val", "minuslogp", "gene_names", "gene_ids", "CondA_mean", "CondB_mean") %in% colnames(res@analysis)))
  # Check that minuslogp is correctly calculated from the adj.P.Val in the results of this specific run
  expect_equal(res@analysis$minuslogp, -log10(res@analysis$adj.P.Val))

  #  Check that significant sets contain the correct members, ignoring order
  # In this dataset, all 4 loci should be significant at FDR < 0.05
  expect_equal(nrow(res@upCond1), 2)
  expect_equal(sort(rownames(res@upCond1)), sort(c("loc1", "loc2")))
  expect_equal(nrow(res@upCond2), 2)
  expect_equal(sort(rownames(res@upCond2)), sort(c("loc3", "loc4")))

  # Check gene annotations and condition means
  expect_equal(res@analysis["loc1", "gene_names"], "GeneA")
  expect_equal(res@analysis["loc3", "gene_names"], "")
  expect_equal(res@analysis["loc2", "CondA_mean"], mean(c(50, 52)))
  expect_equal(res@analysis["loc4", "CondB_mean"], mean(c(40, 41)))

  # Check cond mapping
  expect_equal(res@cond, c("CondA_display" = "CondA", "CondB_display" = "CondB"))
})


test_that("differential_binding handles cases with no significant results", {
  # Create a data list where there is no real difference between conditions
  dl_no_diff <- make_minimal_data_list()
  dl_no_diff$occupancy$CondB_rep1 <- dl_no_diff$occupancy$CondA_rep1
  dl_no_diff$occupancy$CondB_rep2 <- dl_no_diff$occupancy$CondA_rep2

  cond <- c("CondA", "CondB")

  # Run the function. The real limma pipeline should find no significant differences.
  res <- suppressMessages(differential_binding(dl_no_diff, cond, fdr = 0.05))

  # With no real difference, upCond1 and upCond2 should be empty data frames.
  expect_equal(nrow(res@upCond1), 0)
  expect_equal(nrow(res@upCond2), 0)
})
