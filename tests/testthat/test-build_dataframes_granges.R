library(testthat)
library(GenomicRanges)
library(damidBind)

context("build_dataframes_from_granges")

test_that("build_dataframes_from_granges merges GRanges correctly", {
  # Create dummy GRanges objects
  gr1 <- GRanges("chr1", IRanges(c(1, 10, 20), width = 5), score = c(0.1, 0.5, 0.9))
  gr2 <- GRanges("chr1", IRanges(c(1, 10, 20), width = 5), score = c(0.2, 0.6, 1.0))
  gr3 <- GRanges("chr1", IRanges(c(5, 10, 25), width = 5), score = c(0.3, 0.7, 1.1)) # Partial overlap

  gr_list <- list(sampleA = gr1, sampleB = gr2)
  result <- damidBind:::build_dataframes_from_granges(gr_list) # Use ::: for internal functions

  # Expected output structure
  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result), c("chr", "start", "end", "sampleA", "sampleB"))
  expect_equal(nrow(result), 3) # Should have 3 common regions

  # Expected values for merged regions
  expect_equal(result$sampleA, c(0.1, 0.5, 0.9)) # Values from gr1
  expect_equal(result$sampleB, c(0.2, 0.6, 1.0)) # Values from gr2

  # Check sorting
  expect_equal(result$start, c(1, 10, 20))
})

test_that("build_dataframes_from_granges handles non-overlapping regions gracefully", {
  gr_a <- GRanges("chr1", IRanges(1, width = 10), score = 100)
  gr_b <- GRanges("chr1", IRanges(20, width = 10), score = 200)
  gr_list_no_overlap <- list(A = gr_a, B = gr_b)

  result <- damidBind:::build_dataframes_from_granges(gr_list_no_overlap)
  expect_equal(nrow(result), 0) # No common regions should result in 0 rows
})

test_that("build_dataframes_from_granges handles GRanges with missing numeric mcols", {
  gr_no_score <- GRanges("chr1", IRanges(1, width = 10)) # No score MCol
  gr_list_invalid <- list(sampleA = gr_no_score)
  expect_error(
        damidBind:::build_dataframes_from_granges(gr_list_invalid),
        "has no numeric metadata column|invalid subscript type 'list'|subscript out of bounds"
  )
})

test_that("build_dataframes_from_granges handles GRanges with multiple numeric mcols", {
  gr_multi_score <- GRanges("chr1", IRanges(1, width = 10), score1 = 1, score2 = 2)
  gr_list_invalid <- list(sampleA = gr_multi_score)
  expect_error(damidBind:::build_dataframes_from_granges(gr_list_invalid), "has multiple numeric metadata columns")
})

test_that("build_dataframes_from_granges handles empty input list", {
  expect_error(damidBind:::build_dataframes_from_granges(list()), "Empty GRanges list supplied")
})

test_that("build_dataframes_from_granges maintains original order if no merge needed", {
  gr1 <- GRanges("chr1", IRanges(c(1, 10, 20), width = 5), score = c(0.1, 0.5, 0.9))
  gr_list <- list(sampleA = gr1)
  result <- data.frame(
    chr = as.character(seqnames(gr1)),
      start = start(gr1),
      end = end(gr1),
      sampleA = mcols(gr1)$score,
      stringsAsFactors = FALSE
    )
  expect_equal(damidBind:::build_dataframes_from_granges(gr_list), result, ignore_attr = "row.names")
})
