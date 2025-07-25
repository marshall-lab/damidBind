% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/granges_functions.R
\name{all_overlaps_to_original}
\alias{all_overlaps_to_original}
\title{Find best gene overlap(s) for each query interval}
\usage{
all_overlaps_to_original(query, subject, maxgap = 0)
}
\arguments{
\item{query}{GRanges object for which to annotate overlaps.}

\item{subject}{GRanges of gene annotation (must have metadata 'gene_name').}

\item{maxgap}{Integer. Maximum gap allowed for considered overlap (default 0=abutting).}
}
\value{
Character vector (one per query) of gene names, or "" if none.
}
\description{
Annotates each input region with gene(s) it overlaps (comma-separated list).
}
