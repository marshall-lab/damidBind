# This file is sourced by all test files in tests/testthat/

# Dummy EnsDb object for mocking get_ensdb_genes
# This significantly speeds up tests that use gene annotation

make_dummy_ensdb_genes <- function(...) {
  # Mimic essential GRanges structure.
  genes_gr <- GRanges(
    seqnames = Rle("2L", 7),
    ranges = IRanges(
      start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
      end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
    ),
    strand = Rle(strand(c("+", "-", "+", "+", "-", "-", "+"))),
    gene_id = c("FBgn0000001", "FBgn0000002", "FBgn0000003", "FBgn0000004", "FBgn0000005", "FBgn0000006", "FBgn0000007"),
    gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene"),
    gene_biotype = c("protein_coding", "snoRNA", "protein_coding", "protein_coding", "tRNA", "protein_coding", "protein_coding")
  )
  genome(genes_gr) <- "dm6"

  list(
    genes = genes_gr,
    ensembl_version = "113",
    genome_build = "BDGP6",
    species = "Drosophila melanogaster",
    common_name = "Fruit fly"
  )
}


# This helper creates a list that pretends to be an EnsDb object.
# It only needs to respond to the functions that get_ensdb_genes() calls on it.
create_mock_ensdb_object <- function() {
  mock_metadata_df <- data.frame(
    name = c("ensembl_version", "species", "common_name"),
    value = c("113", "Drosophila melanogaster", "Fruit fly"),
    stringsAsFactors = FALSE
  )

  # A very simple mock a GRanges object for genes
  mock_genes_gr <- GRanges(
    seqnames = Rle(c("2L", "3R")),
    ranges = IRanges(start = c(1000, 2000), end = c(1500, 2500)),
    strand = Rle(strand(c("+", "-"))),
    gene_id = c("FBgn001", "FBgn002"),
    gene_name = c("geneA", "geneB"),
    gene_biotype = c("protein_coding", "snoRNA")
  )
  genome(mock_genes_gr) <- "BDGP6.49" # Match the genome build for consistency

  ensdb_obj <- list(
    # Functions that the mock object must respond to
    genes = function(...) mock_genes_gr,
    metadata = function(...) mock_metadata_df,
    dbconn = function(...) list() # Mock connection object
  )
  # Give it the class "EnsDb" so 'inherits' checks will pass
  class(ensdb_obj) <- "EnsDb"
  ensdb_obj
}


setClass("MockHubResult",
  slots = c(data = "DataFrame")
)
setMethod("length", "MockHubResult", function(x) {
  nrow(x@data)
})
setMethod("mcols", "MockHubResult", function(x, use.names = FALSE, ...) {
  x@data
})
setMethod("[", "MockHubResult", function(x, i, j, ..., drop = TRUE) {
  x@data <- x@data[i, , drop = FALSE]
  x
})
setMethod("$", "MockHubResult", function(x, name) {
  x@data[[name]]
})
setMethod("[[", "MockHubResult", function(x, i, j, ...) {
  create_mock_ensdb_object()
})
