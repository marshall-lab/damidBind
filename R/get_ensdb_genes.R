#' Extract gene annotation from Ensembl via AnnotationHub EnsDb
#'
#' Retrieves gene information for a given organism from the most appropriate
#' Ensembl database hosted via Bioconductor's AnnotationHub and ensembldb.
#'
#' @param organism_keyword Character. Unique non-case-senstive string to search
#'        for the organism (e.g., "drosophila melanogaster").
#' @param genome_build Optional character. Genome build identifier to further
#'        restrict the EnsDb selection (e.g., "BDGP6").
#' @param ensembl_version Optional integer. Specific Ensembl version to fetch.
#'        If NULL, the latest available version is used.
#' @param exclude_biotypes Character vector. Gene biotypes to exclude from the
#'        result (default: c("transposable_element", "pseudogene")).
#' @param include_gene_metadata Character vector. Metadata columns to keep for
#'        each gene (default: c("gene_id", "gene_name")).
#'
#' @return List with:
#'   \item{genes}{A GRanges object of genes (metadata columns per argument).}
#'   \item{ensembl_version}{Character. The Ensembl version string.}
#'   \item{genome_build}{Character. Genome build identifier.}
#'   \item{species}{Character. Latin binomial species name.}
#'   \item{common_name}{Character. Species common name.}
#'
#' @details
#' This function queries AnnotationHub for EnsDb objects matching a supplied
#' organism keyword, with optional filtering by genome build and Ensembl version.
#' Genes matching excluded biotypes are filtered out. Only user-selected
#' metadata fields are retained in the genes output.
#'
#' @examples
#' \dontrun{
#' result <- get_ensdb_genes("drosophila melanogaster", ensembl_version = 110)
#' }
#'
#' @import AnnotationHub
#' @import ensembldb
#' @importFrom GenomeInfoDb genome
#' @export
get_ensdb_genes <- function(
    organism_keyword = "drosophila melanogaster",
    genome_build = NULL,
    ensembl_version = NULL,
    exclude_biotypes = c("transposable_element", "pseudogene"),
    include_gene_metadata = c("gene_id", "gene_name")
) {
  message("Finding genome versions ...")

  ah <- AnnotationHub()

  # Query EnsDb databases matching the organism keyword (case-insensitive)
  query_res <- query(ah, c("EnsDb", organism_keyword))

  if (length(query_res) == 0) {
    stop("No EnsDb found for organism: ", organism_keyword)
  }

  # Extract metadata DataFrame for all query results
  meta_df <- mcols(query_res)

  # If genome_build is specified, filter query_res by genome_build (case-insensitive match)
  if (!is.null(genome_build)) {
    matched_genomes <- tolower(meta_df$genome) == tolower(genome_build)
    if (!any(matched_genomes)) {
      stop(sprintf("No EnsDb found for organism '%s' with genome build '%s'. Available builds are: %s",
                   organism_keyword, genome_build, paste(unique(meta_df$genome), collapse = ", ")))
    }
    filtered_query_res <- query_res[matched_genomes]
  } else {
    filtered_query_res <- query_res
  }

  # Extract Ensembl versions from titles, e.g. "Ensembl 113 EnsDb for ..."
  ens_versions <- as.integer(sub("Ensembl ([0-9]+) EnsDb.*", "\\1", filtered_query_res$title))

  if (all(is.na(ens_versions))) {
    # fallback: pick last resource if version parsing fails
    index <- length(filtered_query_res)
  } else if (!(is.null(ensembl_version))) {
    if (ensembl_version %in% ens_versions) {
      index <- which(ens_versions == ensembl_version)
    } else {
      stop(sprintf("Version '%s' is not available. Available builds are: %s",
                   ensembl_version, paste(ens_versions, collapse = ", ")))
    }
  } else {
    index <- which.max(ens_versions)
  }
  message(sprintf("Loading Ensembl genome version '%s'",filtered_query_res$title[[index]]))
  ensdb <- filtered_query_res[[index]]

  # Extract genes object and filter
  genes_gr <- genes(ensdb)
  genes_gr <- genes_gr[!(mcols(genes_gr)$gene_biotype %in% exclude_biotypes),]
  genes_gr <- genes_gr[,colnames(mcols(genes_gr)) %in% include_gene_metadata]

  # Extract genome build
  genome_build <- unique(genome(genes_gr))[1]

  # Extract metadata
  ensdb_meta_df <- metadata(ensdb)
  ensdb_metav <- setNames(as.character(ensdb_meta_df[, 2]), ensdb_meta_df[, 1])
  ensembl_version <- if ("ensembl_version" %in% names(ensdb_metav)) ensdb_metav[["ensembl_version"]] else NA
  species <- if ("species" %in% names(ensdb_metav)) ensdb_metav[["species"]] else NA
  common_name <- if ("common_name" %in% names(ensdb_metav)) ensdb_metav[["common_name"]] else NA

  # Close the DB connection to prevent warning
  dbcon <- dbconn(ensdb)
  dbDisconnect(dbcon)

  return(list(
    genes = genes_gr,
    ensembl_version = ensembl_version,
    genome_build = genome_build,
    species = species,
    common_name = common_name
  ))
}
