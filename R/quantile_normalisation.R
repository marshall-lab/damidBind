#' Quantile Normalisation
#'
#' Performs quantile normalisation of a numeric matrix in native R, matching
#' the algorithm used by `preprocessCore` (including its tie-handling rule).
#'
#' @param x A numeric matrix; rows are features (e.g., genes), columns are samples/arrays.
#'
#' @return A numeric matrix of the same dimensions as \code{x}, quantile normalised.
#'
#' @details
#' This function is a native R implementation of the standard quantile
#' normalisation algorithm. It is designed to be a drop-in replacement for,
#' and produce identical results to, the function of the same name in the
#' `preprocessCore` package.
#'
#' This native R version is provided within `damidBind` to avoid known
#' issues where the `preprocessCore` package can lead to errors or cause R to crash on some
#' Linux systems due to conflicts with OpenMP and/or BLAS/LAPACK library
#' configurations. By providing this native R implementation, `damidBind` ensures
#' it works reliably for all users without requiring them to recompile
#' dependencies or manage system environment variables.
#'
#' This implementation exactly mirrors the behaviour of the `preprocessCore` library’s
#' classic quantile normalisation, including its specific handling of ties:
#' average ranks are computed for ties, and if the fractional part of a rank is greater than 0.4,
#' the output value is the average of the two adjacent quantile means; otherwise, only the lower
#' (floored) quantile mean is used.
#'
#' The function stops if any NA, Inf, or NaN values are present in \code{x}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(9), nrow = 3)
#' quantile_normalisation(x)
#'
#' @export
quantile_normalisation <- function(x) {
  # Argument checks
  if (!is.matrix(x))
    stop("Input 'x' must be a matrix.")

  if (!is.numeric(x))
    stop("Input 'x' must be a numeric matrix.")

  # Check for missing or non-finite values
  if (anyNA(x))
    stop("Input 'x' contains missing values (NA). Please remove or impute them before normalisation.")

  if (any(!is.finite(x)))
    stop("Input 'x' contains non-finite values (Inf or NaN). Please remove or replace them before normalisation.")

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (n_row == 0L || n_col == 0L)
    stop("Input 'x' must have at least one row and one column.")

  if (n_row == 1L) return(x)

  # Compute quantile means
  sorted_cols <- apply(x, 2, sort, method = "quick")
  quantile_means <- rowMeans(sorted_cols)

  # Normalise each column using PreprocessCore's tie rule
  x_norm <- matrix(NA_real_, n_row, n_col)
  for (j in seq_len(n_col)) {
    col <- x[, j]
    ranks <- rank(col, ties.method = "average")
    r_floor <- floor(ranks)
    frac <- ranks - r_floor
    idx1 <- pmax(1L, r_floor)
    idx2 <- pmin(n_row, r_floor + 1L)

    values <- quantile_means[idx1]
    w <- frac > 0.4
    if (any(w))
      values[w] <- 0.5 * (quantile_means[idx1[w]] + quantile_means[idx2[w]])

    x_norm[, j] <- values
  }

  colnames(x_norm) <- colnames(x)
  rownames(x_norm) <- rownames(x)

  x_norm
}

#' Quantile Normalization (US spelling alias)
#'
#' This is a direct alias for \code{quantile_normalisation()}; it is provided for users who prefer US English spelling.
#' See \code{?quantile_normalisation} for documentation.
#'
#' @param x A numeric matrix; rows are features (e.g., genes), columns are samples/arrays.
#'
#' @aliases quantile_normalization
#' @export
quantile_normalization <- quantile_normalisation
