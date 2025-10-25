library(testthat)
library(damidBind)

context("Quantile Normalisation")

test_that("quantile_normalisation handles basic numeric matrices", {
    mat2 <- matrix(c(1, 3, 2, 4, 6, 5, 7, 9, 8), nrow = 3, byrow = TRUE)
    # Sorted columns:
    # col1: 1, 4, 7
    # col2: 3, 6, 9
    # col3: 2, 5, 8
    # Quantile means (row means of sorted cols):
    # row1: mean(1,3,2) = 2
    # row2: mean(4,6,5) = 5
    # row3: mean(7,9,8) = 8
    # Final matrix has each value replaced by the quantile mean corresponding to its rank within its original column.
    # Since each original column is just a permutation of (1,2,3), (4,5,6), (7,8,9), ranks will map to the same means.
    expected_norm2 <- matrix(rep(c(2, 5, 8), 3), nrow = 3)

    # Run the normalisation
    actual_norm <- quantile_normalisation(mat2)

    # Use expect_equal for numeric comparison with tolerance
    expect_equal(actual_norm, expected_norm2, tolerance = 1e-6)
    expect_equal(dim(actual_norm), dim(mat2))
    expect_true(is.numeric(actual_norm))
})


test_that("quantile_normalisation handles edge cases (single row/col)", {
    mat_single_row <- matrix(c(1, 5, 2), nrow = 1)
    expect_equal(quantile_normalisation(mat_single_row), mat_single_row)

    mat_single_col <- matrix(c(1, 5, 2), ncol = 1)
    expect_equal(quantile_normalisation(mat_single_col), mat_single_col)
})

test_that("quantile_normalisation throws errors for invalid inputs", {
    # Test non-matrix input
    expect_error(quantile_normalisation(c(1, 2, 3)), "Input 'x' must be a matrix.")
    # Test non-numeric input
    expect_error(quantile_normalisation(matrix(c("a", 1), nrow = 1)), "Input 'x' must be a numeric matrix.")
    # Test NA values
    expect_error(quantile_normalisation(matrix(c(1, NA, 3, 4), nrow = 2)), "Input 'x' contains missing values", fixed = TRUE)
    # Test infinite values
    expect_error(quantile_normalisation(matrix(c(1, Inf, 3, 4), nrow = 2)), "Input 'x' contains non-finite values", fixed = TRUE)
    # Test NaN values
    expect_error(quantile_normalisation(matrix(c(1, NaN, 3, 4), nrow = 2)), "Input 'x' contains missing values", fixed = TRUE)
    # Test empty numeric matrix
    expect_error(quantile_normalisation(matrix(numeric(), nrow = 0, ncol = 3)), "Input 'x' must have at least one row and one column.")
    expect_error(quantile_normalisation(matrix(numeric(), nrow = 3, ncol = 0)), "Input 'x' must have at least one row and one column.")
})


test_that("quantile_normalization alias works correctly", {
    mat <- matrix(c(1, 3, 2, 4, 6, 5, 7, 9, 8), nrow = 3, byrow = TRUE)
    expect_equal(quantile_normalization(mat), quantile_normalisation(mat))
})
