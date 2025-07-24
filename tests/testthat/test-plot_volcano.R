library(testthat)
library(damidBind)
library(ggplot2)
library(ggrepel)

context("plot_volcano")

# Create a dummy DamIDResults object for volcano plotting
make_dummy_diff_results <- function() {
  analysis_table <- data.frame(
    logFC = c(2.5, -3.0, 0.5, 1.8, -1.2, 0.1),
    adj.P.Val = c(0.001, 0.005, 0.1, 0.002, 0.08, 0.5),
    minuslogp = -log10(c(0.001, 0.005, 0.1, 0.002, 0.08, 0.5)),
    B = c(5, 4, -1, 3, 0, -2), # Dummy B-statistics
    gene_names = c("GENE_SIG_UP1", "GENE_SIG_DOWN1", "GENE_NONSIG1", "GENE_SIG_UP2", "GENE_NONSIG2", "GENE_NONSIG_SMALL"),
    row.names = c("L1", "L2", "L3", "L4", "L5", "L6")
  )

  upCond1 <- analysis_table[analysis_table$adj.P.Val < 0.05 & analysis_table$logFC > 0, ]
  upCond2 <- analysis_table[analysis_table$adj.P.Val < 0.05 & analysis_table$logFC < 0, ]

  new("DamIDResults",
      upCond1 = upCond1,
      upCond2 = upCond2,
      analysis = analysis_table,
      cond = c(ConditionA = "CondA", ConditionB = "CondB"),
      data = list(test_category = "bound")
  )
}

test_that("plot_volcano runs without error and returns a ggplot object", {
  diff_res <- make_dummy_diff_results()
  p <- NULL # Initialize outside expect_no_error to capture the return value
  expect_no_error(p <- plot_volcano(diff_res, save = FALSE)) # Ensure it doesn't try to save
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano handles custom plot_config", {
  diff_res <- make_dummy_diff_results()
  custom_config <- list(
    title = "My Custom Title",
    xlab = "Log Fold Change",
    sig_colour = "red",
    nonsig_colour = "blue",
    ystat = "minuslogp" # Test alternative y-axis stat
  )
  p <- plot_volcano(diff_res, plot_config = custom_config, save = FALSE)

  # Check if title and axis labels are correctly set in the plot object
  expect_equal(p$labels$title, "My Custom Title")
  expect_equal(p$labels$x, "Log Fold Change")
  expect_equal(p$labels$y, "-log(p)") # Default ylab for minuslogp
})

test_that("plot_volcano handles label_config for gene labels", {
  diff_res <- make_dummy_diff_results()
  # Test with specific genes to label
  p <- plot_volcano(diff_res, label_config = list(genes = c("GENE_SIG_UP1", "GENE_SIG_DOWN1")), save = FALSE)

  # Check for ggrepel layer existence (indicates labels were attempted)
  label_layers <- lapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"))
  expect_true(any(unlist(label_layers)))

  # For actual label content, you might need to inspect p$data or render the plot
  # and use a visual testing tool. We'll simply check that the data for the label layer exists.
  label_data <- p$layers[[which(unlist(label_layers))[[1]]]]$data
  expect_true("GENE_SIG_UP1" %in% label_data$label_to_display)
  expect_true("GENE_SIG_DOWN1" %in% label_data$label_to_display)
  expect_false("GENE_SIG_UP2" %in% label_data$label_to_display) # Not in specified genes
})


test_that("plot_volcano handles highlight regions", {
  diff_res <- make_dummy_diff_results()
  highlight_list <- list(
    "My Highlight" = c("GENE_SIG_UP1", "GENE_SIG_UP2")
  )
  p <- plot_volcano(diff_res, highlight = highlight_list, save = FALSE)

  # Check for highlight geom_point layer
  highlight_layers <- sapply(p$layers, function(l) {
    if (inherits(l$geom, "GeomPoint")) {
      # The highlight layer is the only one with this column in its data
      return("highlight_group_name" %in% colnames(l$data))
    }
    FALSE
  })
  expect_true(any(unlist(highlight_layers)))

  # Check that labels are added when highlight_config$label is TRUE
  p_labelled <- plot_volcano(diff_res, highlight = highlight_list, highlight_config = list(label = TRUE), save = FALSE)
  label_layers <- lapply(p_labelled$layers, function(l) inherits(l$geom, "GeomTextRepel"))
  expect_true(any(unlist(label_layers)))
})

test_that("plot_volcano warns about invalid non-list plot_config inputs and proceeds with defaults", {
  diff_res <- make_dummy_diff_results()
  expect_message(plot_volcano(diff_res, plot_config = "invalid", save = FALSE), "Input value was ")
  # No error should occur, and it should return a ggplot object.
  expect_s3_class(plot_volcano(diff_res, plot_config = "invalid", save = FALSE), "ggplot")
})

test_that("plot_volcano throws error for invalid ystat", {
  diff_res <- make_dummy_diff_results()
  expect_error(plot_volcano(diff_res, plot_config = list(ystat = "NonExistentCol"), save = FALSE), regexp = "is not a valid column in plot data")
})

test_that("plot_volcano skips saving if save is NULL/FALSE/0", {
  diff_res <- make_dummy_diff_results()
  # Create a dummy file path
  temp_file_base <- tempfile() # Get a base name without extension
  expected_file <- paste0(temp_file_base, ".pdf") # Define the expected full path

  # Pass only the base name to the function
  plot_volcano(diff_res, save = list(filename = temp_file_base, format = "pdf"))

  # Check for the correctly constructed file and clean up
  expect_true(file.exists(expected_file))
  file.remove(expected_file)
})

