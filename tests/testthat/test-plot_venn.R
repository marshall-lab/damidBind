library(testthat)
library(damidBind)

context("Visualization: plot_venn")

# Dummy diff_results object for testing
make_dummy_diff_results_for_venn <- function() {
  analysis_table <- data.frame(
    logFC = c(2.5, -3.0, 0.5, 1.8, -1.2, 0.1, 0.01),
    adj.P.Val = c(0.001, 0.005, 0.1, 0.002, 0.08, 0.5, 0.2), # L1,L2,L4 sig; L5 non-sig; L3,L6 non-sig
    row.names = c("Locus1", "Locus2", "Locus3", "Locus4", "Locus5", "Locus6", "Locus7")
  )

  upCond1 <- analysis_table[c("Locus1", "Locus4"), ] # Sig for Cond1
  upCond2 <- analysis_table[c("Locus2"), ] # Sig for Cond2 (note, not Locus5)

  new("DamIDResults",
    upCond1 = upCond1,
    upCond2 = upCond2,
    analysis = analysis_table,
    cond = c(CondA = "Treatment_A", CondB = "Treatment_B"),
    data = list(test_category = "bound")
  )
}

test_that("plot_venn correctly prepares data and calls BioVenn::draw.venn", {
  diff_res <- make_dummy_diff_results_for_venn()

  # Loci: Locus1,2,3,4,5,6,7
  # upCond1 (CondA): Locus1, Locus4
  # upCond2 (CondB): Locus2
  # All unique sig: Locus1, Locus2, Locus4
  # Non-sig: Locus3, Locus5, Locus6, Locus7

  # Cond1_full: union(upCond1, nonsig) = {Locus1, Locus4, Locus3, Locus5, Locus6, Locus7}
  # Cond2_full: union(upCond2, nonsig) = {Locus2, Locus3, Locus5, Locus6, Locus7}

  # Define a mock for BioVenn::draw.venn
  # This mock will store the parameters it was called with.
  mock_draw_venn_args <- NULL
  mock_draw_venn <- function(list_x, list_y, ..., filename = NULL, output = NULL) {
    mock_draw_venn_args <<- list(
      list_x = list_x,
      list_y = list_y,
      filename = filename,
      output = output,
      other_args = list(...) # Capture other arguments
    )
    message("BioVenn::draw.venn mocked successfully!") # Confirm mock was hit
    invisible(NULL) # Match invisible return of original
  }

  local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

  # Define specific set labels for the test
  custom_set_labels <- c("Set A Features", "Set B Features")

  # Call plot_venn
  expect_no_error({
    plot_venn(
      diff_res,
      title = "My Venn Plot",
      subtitle = "Test Subtitle",
      set_labels = custom_set_labels,
      filename = "test_venn.pdf",
      format = "pdf"
    )
  })

  # Verify mock was called
  expect_true(!is.null(mock_draw_venn_args))

  # Verify list_x and list_y contents based on diff_res's logic
  # As per your function's logic:
  # Cond1_full are all ids *not only* in Cond2.
  # Cond2_full are all ids *not only* in Cond1.
  # This matches the BioVenn logic of:
  # A_only = upCond1_only
  # B_only = upCond2_only
  # AB_overlap = non_sig

  # Expected non-significant loci:
  nonsig_expected <- setdiff(rownames(diff_res@analysis), union(rownames(diff_res@upCond1), rownames(diff_res@upCond2)))
  expect_equal(sort(mock_draw_venn_args$list_x), sort(union(rownames(diff_res@upCond1), nonsig_expected)))
  expect_equal(sort(mock_draw_venn_args$list_y), sort(union(rownames(diff_res@upCond2), nonsig_expected)))

  # Verify other parameters
  expect_equal(mock_draw_venn_args$other_args$xtitle, custom_set_labels[1])
  expect_equal(mock_draw_venn_args$other_args$ytitle, custom_set_labels[2])
  expect_equal(mock_draw_venn_args$other_args$title, "My Venn Plot")
  expect_equal(mock_draw_venn_args$other_args$subtitle, "Test Subtitle")
  expect_equal(mock_draw_venn_args$filename, "test_venn.pdf")
  expect_equal(mock_draw_venn_args$output, "pdf")
})

test_that("plot_venn handles cases with no significant regions", {
  diff_res_no_sig <- make_dummy_diff_results_for_venn()
  # Override with no significant results
  diff_res_no_sig@upCond1 <- diff_res_no_sig@upCond1[FALSE, ]
  diff_res_no_sig@upCond2 <- diff_res_no_sig@upCond2[FALSE, ]

  mock_draw_venn_args <- NULL # Reset mock args
  mock_draw_venn <- function(list_x, list_y, ...) {
    mock_draw_venn_args <<- list(list_x = list_x, list_y = list_y)
    invisible(NULL)
  }
  local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

  # Use evaluate_promise for robustly testing multiple warnings
  res <- evaluate_promise(plot_venn(diff_res_no_sig))
  expect_length(res$warnings, 2)
  expect_match(res$warnings[1], "No loci present in upCond1", fixed = TRUE)
  expect_match(res$warnings[2], "No loci present in upCond2", fixed = TRUE)

  expect_true(!is.null(mock_draw_venn_args))
  # All loci should be in both 'full' sets, as they're all non-significant
  expect_equal(sort(mock_draw_venn_args$list_x), sort(rownames(diff_res_no_sig@analysis)))
  expect_equal(sort(mock_draw_venn_args$list_y), sort(rownames(diff_res_no_sig@analysis)))
})


test_that("plot_venn warns if one condition has no loci", {
  diff_res_partial_sig <- make_dummy_diff_results_for_venn()
  diff_res_partial_sig@upCond1 <- diff_res_partial_sig@upCond1[FALSE, ] # No Cond1 sig

  mock_draw_venn <- function(list_x, list_y, ...) invisible(NULL)
  local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

  expect_warning(plot_venn(diff_res_partial_sig), "No loci present in upCond1")
})

test_that("plot_venn uses default set_labels if not provided", {
  diff_res <- make_dummy_diff_results_for_venn()
  mock_draw_venn_args <- NULL
  mock_draw_venn <- function(list_x, list_y, xtitle, ytitle, ...) {
    mock_draw_venn_args <<- list(xtitle = xtitle, ytitle = ytitle)
    invisible(NULL)
  }
  local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

  expect_no_error(plot_venn(diff_res, set_labels = NULL))
  expect_true(!is.null(mock_draw_venn_args))
  expect_equal(mock_draw_venn_args$xtitle, names(diff_res@cond)[1])
  expect_equal(mock_draw_venn_args$ytitle, names(diff_res@cond)[2])
})
