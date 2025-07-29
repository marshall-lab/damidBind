#' Generic plot method for DamIDResults objects
#'
#' A generic plot method that creates a default visualisation (a volcano plot)
#' for a `DamIDResults` object. For more advanced plotting options or different
#' plot types, see the specific functions `plot_volcano()`, `plot_venn()`, and
#' `analyse_go_terms()`.
#'
#' @param x A `DamIDResults` object.
#' @param y (Missing) Not used.
#' @param ... Additional arguments passed to `plot_volcano()`.
#'
#' @seealso [plot_volcano()], [plot_venn()], [analyse_go_terms()] for more
#' powerful and specific plotting functions.
#'
#' @rdname DamIDResults-class
#' @export
setMethod(
    "plot", signature(x = "DamIDResults", y = "missing"),
    function(x, y, ...) {
        # This method simply calls the main volcano plot function by default.
        plot_volcano(diff_results = x, ...)
    }
)
