% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyse_go_terms.R
\name{analyse_go_enrichment}
\alias{analyse_go_enrichment}
\title{Perform Gene Ontology (GO) Enrichment Analysis for Differentially Bound/Expressed Regions}
\usage{
analyse_go_enrichment(
  diff_results,
  direction = "cond1",
  org_db = org.Dm.eg.db::org.Dm.eg.db,
  ontology = "BP",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  plot_title = NULL,
  show_category = 12,
  label_format_width = 30,
  save = NULL,
  save_results_path = NULL,
  maxGSSize = 1000,
  minGSSize = 10
)
}
\arguments{
\item{diff_results}{A `DamIDResults` object, as returned by
`differential_binding()` or `differential_accessibility()`.}

\item{direction}{Character string. Specifies which set of genes to analyse, either using condition names,
"cond1" or "cond2", or "all" (for all significantly enriched genes from either direction).
Default is "cond1".}

\item{org_db}{An OrgDb object specifying the organism's annotation database.
For Drosophila, use `org.Dm.eg.db::org.Dm.eg.db`.}

\item{ontology}{Character string. The GO ontology to use: "BP" (Biological Process),
"MF" (Molecular Function), or "CC" (Cellular Component). Default is "BP".}

\item{pvalue_cutoff}{Numeric. Adjusted p-value cutoff for significance. Default: 0.05.}

\item{qvalue_cutoff}{Numeric. Q-value cutoff for significance. Default: 0.2.}

\item{plot_title}{Character string. Title for the generated dot plot.}

\item{show_category}{Integer. Number of top enriched GO categories to display in the plot. Default: 12.}

\item{label_format_width}{Integer. Max character length for GO term labels on the plot. Default: 30.}

\item{save}{List or `NULL`. Controls saving the plot to a file (dot plot).
If `NULL`, `FALSE`, or `0`, the plot is not saved.
If a `list`, it specifies saving parameters:
\itemize{
  \item \code{filename} (character): The path and base name for the output file. If not specified, the default name "damidBind_GSEA_dotplot" is used.
  \item \code{format} (character): File format ("pdf", "svg", or "png"). Default is "pdf".
  \item \code{width} (numeric): Width of the plot in inches. Default is 6.
  \item \code{height} (numeric): Height of the plot in inches. Default is 6.
}}

\item{save_results_path}{Character string or NULL. If a path is provided (e.g., "go_results.csv"), the
enrichment results table will be saved to this CSV file.}

\item{maxGSSize}{Integer. Maximum size of gene sets to consider. Default: 1000.}

\item{minGSSize}{Integer. Minimum size of gene sets to consider. Default: 10.}
}
\value{
A list containing:
  \item{enrich_go_object}{`enrichResult` object from `clusterProfiler`.}
  \item{results_table}{Data frame of enrichment results.}
  \item{dot_plot}{`ggplot` object of the dot plot.}
  NULL if no significant enrichment is found or if input validation fails.
}
\description{
This function performs Gene Ontology (GO) enrichment analysis using `clusterProfiler`
for either the upCond1ulated or upCond2ulated regions/genes identified by
`differential_binding()` or `differential_accessibility()`. It automatically extracts the relevant Flybase IDs (FBgnIDs)
and the background universe from the input `DamIDResults` object.
}
\details{
This function assumes that the `analysis` slot in the `diff_results`
object contains a `gene_ids` column.
If this column is not present, or cannot be processed, the function will return NULL.

The function includes an internal helper `clean_gene_symbols` which filters common ambiguous
gene symbols (snoRNA, snRNA, tRNA) that may not be useful for GO enrichment.
}
