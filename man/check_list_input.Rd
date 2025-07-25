% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_volcano.R
\name{check_list_input}
\alias{check_list_input}
\title{Check and Coerce List Input for Configuration}
\usage{
check_list_input(default_list, config_input)
}
\arguments{
\item{default_list}{A `list` representing the default configuration keys and
their values. This will be the base for merging or the fallback if
`config_input` is invalid.}

\item{config_input}{The user-provided configuration input. This is expected
to be a `list`, but can also be `NULL`, `FALSE`, or `0`.
\itemize{
  \item If `config_input` is a `list`, it is merged with `default_list` using
    `modifyList`, and the merged list is returned.
  \item If `config_input` is `NULL`, `FALSE`, or `0`, the function returns `NULL`,
    indicating that no configuration should be applied (e.g., to turn off a feature).
  \item If `config_input` is any other non-list single value, a message is printed
    informing the user that the input was not a list, and the `default_list`
    is returned, allowing the plot to proceed with defaults.
}}
}
\value{
A `list` (the merged configuration or `default_list` if `config_input`
  was a non-list invalid type) or `NULL` (if `config_input` was `NULL`, `FALSE`, or `0`).
}
\description{
This helper function validates and processes user-provided configuration inputs.
It ensures the input is a list and, if so, merges it with a set of default
values. It also handles non-list inputs gracefully, treating specific values
(e.g., `FALSE`, `0`) as `NULL` (meaning "no configuration"), and warning the
user if an unexpected non-list value is provided, while still merging valid
list inputs with defaults.
}
\keyword{internal}
