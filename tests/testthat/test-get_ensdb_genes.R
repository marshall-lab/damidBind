# FILE: tests/testthat/test-get_ensdb_genes.R

library(testthat)
library(damidBind)
library(AnnotationHub)
library(ensembldb)
library(GenomicRanges)
library(S4Vectors)

context("get_ensdb_genes")

source(test_path("_test_helper.R"))

# --- Mock AnnotationHub and related functions for isolated testing ---

# This helper function creates a mock DataFrame that mimics the result of an AnnotationHub query.
create_mock_query_result <- function() {
    df <- data.frame(
        ah_id = c("AH_Dm_113", "AH_Dm_112", "AH_Hs_113"),
        title = c(
            "Ensembl 113 EnsDb for Drosophila melanogaster",
            "Ensembl 112 EnsDb for Drosophila melanogaster",
            "Ensembl 113 EnsDb for Homo sapiens"
        ),
        species = c("Drosophila melanogaster", "Drosophila melanogaster", "Homo sapiens"),
        genome = c("BDGP6.49", "BDGP6.48", "GRCh38.p14"),
        stringsAsFactors = FALSE
    )
    # It's crucial that this is an S4 DataFrame, not a base data.frame
    S4Vectors::DataFrame(df)
}

test_that("get_ensdb_genes can retrieve Drosophila EnsDb and filter biotypes", {
    # Mock the necessary functions using testthat's local_mocked_bindings
    local_mocked_bindings(
        # Mock the main AnnotationHub constructor to return a dummy object
        AnnotationHub = function(...) new("AnnotationHub"),
        # Mock the query function
        query = function(ah_obj, search_terms) {
            # This part creates the raw data.frame
            df <- S4Vectors::DataFrame(
                ah_id = c("AH_Dm_113", "AH_Dm_112", "AH_Hs_113"),
                title = c(
                    "Ensembl 113 EnsDb for Drosophila melanogaster",
                    "Ensembl 112 EnsDb for Drosophila melanogaster",
                    "Ensembl 113 EnsDb for Homo sapiens"
                ),
                species = c("Drosophila melanogaster", "Drosophila melanogaster", "Homo sapiens"),
                genome = c("BDGP6.49", "BDGP6.48", "GRCh38.p14")
            )
            if (any(grepl("drosophila melanogaster", search_terms, ignore.case = TRUE))) {
                # Filter the DataFrame
                res_df <- df[grepl("Drosophila melanogaster", df$species), ]
                # Return the DataFrame wrapped in our mock S4 object
                return(new("MockHubResult", data = res_df))
            }
            # For no matches, return the S4 object with a zero-row DataFrame
            return(new("MockHubResult", data = df[FALSE, ]))
        },
        # Mock the ensembldb functions to return mock data directly
        # This avoids the need to mock the `[[` operator
        genes = function(x, ...) {
            # This mock now simply returns the genes from our helper, ignoring the input 'x'
            create_mock_ensdb_object()$genes()
        },
        metadata = function(x, ...) {
            create_mock_ensdb_object()$metadata()
        },
        dbconn = function(x, ...) {
            create_mock_ensdb_object()$dbconn()
        },
        # Mock the DBI function
        dbDisconnect = function(conn) invisible(NULL)
    )

    result <- suppressMessages({
        get_ensdb_genes(
            organism_keyword = "drosophila melanogaster",
            exclude_biotypes = "snoRNA"
        )
    })

    expect_type(result, "list")
    expect_named(result, c("genes", "ensembl_version", "genome_build", "species", "common_name"))
    expect_s4_class(result$genes, "GRanges")

    # Check that the biotype was correctly excluded
    expect_equal(length(result$genes), 1)
    expect_equal(result$genes$gene_name, "geneA")
    expect_false("snoRNA" %in% result$genes$gene_biotype)

    # Check that metadata was correctly extracted from the mock
    expect_equal(result$ensembl_version, "113")
    expect_equal(result$genome_build, "BDGP6.49")
})

test_that("get_ensdb_genes handles errors for non-existent organisms and invalid parameters", {
    # The same mocks as above will be in scope here.
    # The key change is REMOVING the mock for `[[` and simplifying the others.
    local_mocked_bindings(
        # Mock the main AnnotationHub constructor to return a dummy object
        AnnotationHub = function(...) new("AnnotationHub"),

        # Mock the query function
        query = function(ah_obj, search_terms) {
            # This part creates the raw data.frame
            df <- S4Vectors::DataFrame(
                ah_id = c("AH_Dm_113", "AH_Dm_112", "AH_Hs_113"),
                title = c(
                    "Ensembl 113 EnsDb for Drosophila melanogaster",
                    "Ensembl 112 EnsDb for Drosophila melanogaster",
                    "Ensembl 113 EnsDb for Homo sapiens"
                ),
                species = c("Drosophila melanogaster", "Drosophila melanogaster", "Homo sapiens"),
                genome = c("BDGP6.49", "BDGP6.48", "GRCh38.p14")
            )
            if (any(grepl("drosophila melanogaster", search_terms, ignore.case = TRUE))) {
                # Filter the DataFrame
                res_df <- df[grepl("Drosophila melanogaster", df$species), ]
                # Return the DataFrame wrapped in our mock S4 object
                return(new("MockHubResult", data = res_df))
            }
            # For no matches, return the S4 object with a zero-row DataFrame
            return(new("MockHubResult", data = df[FALSE, ]))
        },

        # Mock ensembldb functions to return mock data directly, avoiding `[[`
        genes = function(x, ...) create_mock_ensdb_object()$genes(),
        metadata = function(x, ...) create_mock_ensdb_object()$metadata(),
        dbconn = function(x, ...) create_mock_ensdb_object()$dbconn(),

        # Mock the DBI function
        dbDisconnect = function(conn) invisible(NULL)
    )

    # Test for non-existent organism
    expect_error(
        suppressMessages(get_ensdb_genes("nonexistentus organismus")),
        "No EnsDb found for organism:"
    )

    # Test for invalid genome_build
    expect_error(
        suppressMessages(get_ensdb_genes("drosophila melanogaster", genome_build = "INVALID_BUILD")),
        "No EnsDb found for organism .* with genome build"
    )

    # Test for invalid ensembl_version
    expect_error(
        suppressMessages(get_ensdb_genes("drosophila melanogaster", ensembl_version = 9999)),
        "Version '9999' is not available."
    )
})
