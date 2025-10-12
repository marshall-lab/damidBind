#' damidBind: Differential Binding Analysis for DamID-seq Data
#'
#' @description
#' The damidBind package provides a streamlined workflow for determining
#' differential protein binding, RNA polymerase occupancy, or chromatin
#' accessibility from DamID-based sequencing experiments. It handles data loading,
#' processing, statistical analysis, and provides a suite of visualisation tools
#' for interpretation and exploration of the results.
#'
#' @details
#' The package is designed for three main experimental types:
#' \itemize{
#'   \item \strong{Transcription Factor Binding}: Analysis of conventional DamID or TaDa data to find differential binding sites for a protein of interest. Uses log-ratio data and the limma backend via \code{\link{differential_binding}}.
#'   \item \strong{Gene Transcription}: Analysis of RNA Polymerase II TaDa data to infer differential gene expression, analysed over gene bodies with \code{\link{load_data_genes}}.
#'   \item \strong{Chromatin Accessibility}: Analysis of CATaDa data to find differential accessibility. Uses count-based data and the NOIseq backend via \code{\link{differential_accessibility}}.
#' }
#'
#' \strong{Core Workflow:}
#' \enumerate{
#'   \item \strong{Load Data:} Use \code{\link{load_data_peaks}} for TF binding/accessibility or \code{\link{load_data_genes}} for RNA Pol II occupancy. These functions read bedGraph and peak files, calculate occupancy scores, and annotate regions with nearby genes.
#'   \item \strong{Perform Differential Analysis:} Use \code{\link{differential_binding}} for conventional DamID log-ratio data (limma) or \code{\link{differential_accessibility}} for CATaDa count data (NOIseq). These return a \code{\linkS4class{DamIDResults}} object.
#'   \item \strong{Visualise and Explore:} Use the plotting functions on the `DamIDResults` object: \code{\link{plot_volcano}}, \code{\link{plot_venn}} and  \code{\link{analyse_go_enrichment}}; and the interactive \code{\link{browse_igv_regions}}  to browse differentially-bound regions in an interactive IGV browser window..
#' }
#'
#' For a complete walkthrough, please see the package vignette by running:
#' `browseVignettes("damidBind")`
#'
#' @seealso
#' \strong{Primary Functions:}
#' \itemize{
#'   \item \code{\link{load_data_peaks}}: Load binding data and associated peak regions.
#'   \item \code{\link{load_data_genes}}: Load binding data summarised over gene bodies.
#'   \item \code{\link{differential_binding}}: Perform differential analysis for log-ratio data.
#'   \item \code{\link{differential_accessibility}}: Perform differential analysis for count-based data.
#'   \item \code{\linkS4class{DamIDResults}}: The main results object returned by analysis functions.
#' }
#'
#' \strong{Useful links:}
#' \itemize{
#'   \item The damidBind vignette: `browseVignettes("damidBind")`
#'   \item Report bugs at \href{https://github.com/marshall-lab/damidBind/issues}{https://github.com/marshall-lab/damidBind/issues}
#' }
#'
#' @keywords internal
"_PACKAGE"

## Function imports

## From dplyr: for data manipulation with pipes
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr rename
#'
## From ggplot2: for plotting
#' @import ggplot2
#'
## From forcats: for factor manipulation in plots
#' @importFrom forcats fct_reorder
#'
## From stringr: for text manipulation and matching
#' @importFrom stringr str_wrap
#' @importFrom stringr str_match
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#'
## From rtracklayer: for importing genomic file formats
#' @importFrom rtracklayer import
#'
## From tools: for file path manipulation
#' @importFrom tools file_path_sans_ext
#'
## From AnnotationHub: for fetching genomic annotations
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom AnnotationHub query
#'
## From ensembldb: for working with Ensembl databases
#' @importFrom ensembldb genes
#' @importFrom ensembldb metadata
#' @importFrom ensembldb dbconn
#'
## From limma: for differential analysis of log-ratio data
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#'
## From NOISeq: for differential analysis of count data
#' @importFrom NOISeq readData
#' @importFrom NOISeq noiseq
#' @importFrom NOISeq degenes
#'
## From GenomicRanges: for handling genomic intervals
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges trim
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#'
## From IRanges: for handling genomic ranges
#' @importFrom IRanges IRanges
#'
## From S4Vectors: for advanced vector and DataFrame types
#' @importFrom S4Vectors Rle
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#'
## From BiocParallel: for parallel computation
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#'
## From stats: for statistical functions used widely
#' @importFrom stats weighted.mean
#' @importFrom stats na.omit
#' @importFrom stats setNames
#'
## Imports for Shiny/IGV functionality
#' @importFrom igvShiny igvShiny
#' @importFrom igvShiny renderIgvShiny
#' @importFrom igvShiny igvShinyOutput
#' @importFrom igvShiny loadBedTrack
#' @importFrom igvShiny loadBedGraphTrack
#' @importFrom igvShiny parseAndValidateGenomeSpec
#' @importFrom igvShiny showGenomicRegion
#' @importFrom shiny shinyApp
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny sidebarLayout
#' @importFrom shiny sidebarPanel
#' @importFrom shiny mainPanel
#' @importFrom shiny observeEvent
#' @importFrom shiny p
#' @importFrom shiny h4
#' @importFrom shiny hr
#' @importFrom shiny runApp
#' @importFrom DT DTOutput
#' @importFrom DT renderDT
#' @importFrom DT datatable
#' @importFrom DT formatRound
#'
## From BioVenn: for Venn diagrams
#' @importFrom BioVenn draw.venn
#'
## From DBI: for database connections (ensdb)
#' @importFrom DBI dbDisconnect
#' @importFrom utils modifyList
#' @importFrom methods is
#' @importFrom methods new
#'
NULL
