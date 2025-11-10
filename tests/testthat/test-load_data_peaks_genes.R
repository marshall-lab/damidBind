library(testthat)
library(damidBind)
library(GenomicRanges)
library(S4Vectors)
library(dplyr)

# Source helper function
source(test_path("_test_helper.R"))

# Mock get_ensdb_genes for all load_* tests
local_mocked_bindings(get_ensdb_genes = make_dummy_ensdb_genes)

context("Data Loading: load_data_peaks")

test_that("load_data_peaks loads and processes extdata files correctly", {
    data_dir <- system.file("extdata", package = "damidBind")
    binding_profiles_path <- data_dir
    peaks_path <- data_dir

    # Call the function
    dl <- load_data_peaks(
        binding_profiles_path = binding_profiles_path,
        peaks_path = peaks_path,
        quantile_norm = FALSE, # Test without quantile normalisation first
        BPPARAM = BiocParallel::SerialParam()
    )

    # Check overall structure of the returned list
    expect_type(dl, "list")
    expect_named(dl, c("binding_profiles_data", "peaks", "pr", "occupancy", "test_category"))

    # Check binding_profiles_data
    expect_s4_class(dl$binding_profiles_data, "GRanges")
    expect_true(length(dl$binding_profiles_data) > 0)
    # Check for metadata columns in the GRanges object
    expect_true(all(c(
        "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm"
    ) %in%
        names(mcols(dl$binding_profiles_data))))

    # Check peaks (list of GRanges)
    expect_type(dl$peaks, "list")
    expect_length(dl$peaks, 4) # Should be 4 peak files
    expect_s4_class(dl$peaks[[1]], "GRanges")

    # Check pr (reduced peaks)
    expect_s4_class(dl$pr, "GRanges")
    expect_true(length(dl$pr) > 0)
    expect_true("name" %in% names(mcols(dl$pr))) # Should have names like chr:start-end

    # Check occupancy data frame
    expect_s3_class(dl$occupancy, "data.frame")
    expect_true(nrow(dl$occupancy) > 0)
    expect_true(all(c(
        "name", "nfrags",
        "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm",
        "gene_name", "gene_id"
    ) %in%
        colnames(dl$occupancy)))
    expect_true(all(rownames(dl$occupancy) == dl$occupancy$name))

    # Check test_category
    expect_equal(dl$test_category, "bound")

    # Check for gene annotations (from dummy ensdb_genes)
    expect_true(any(dl$occupancy$gene_name != ""))
    expect_true(any(dl$occupancy$gene_id != ""))
})

test_that("load_data_peaks works with quantile_norm = TRUE", {
    data_dir <- system.file("extdata", package = "damidBind")
    binding_profiles_path <- data_dir
    peaks_path <- data_dir

    dl_qnorm <- load_data_peaks(
        binding_profiles_path = binding_profiles_path,
        peaks_path = peaks_path,
        quantile_norm = TRUE,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Check for renamed, quantile-normalised metadata columns
    expect_true(all(c(
        "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm_qnorm",
        "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm_qnorm",
        "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm_qnorm",
        "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm_qnorm"
    ) %in%
        names(mcols(dl_qnorm$binding_profiles_data))))

    # Check that values are indeed normalised
    # Correctly extract the matrix from the GRanges mcols
    signal_cols_qnorm <- as.matrix(mcols(dl_qnorm$binding_profiles_data))

    col_means <- colMeans(signal_cols_qnorm, na.rm = TRUE)
    # Quantile normalization makes distributions similar; mean/median should converge
    expect_true(sd(col_means) < 0.5) # Using noisy real data
    # Also check that sample medians are similar
    sample_medians <- apply(signal_cols_qnorm, 2, median, na.rm = TRUE)
    expect_true(sd(sample_medians) < 0.5)
})

test_that("load_data_peaks throws error if no binding files found", {
    temp_dir <- tempdir()
    expect_error(
        load_data_peaks(
            binding_profiles_path = file.path(temp_dir, "*.nonexistent"),
            peaks_path = file.path(temp_dir, "*.nonexistent"),
            ensdb_genes = make_dummy_ensdb_genes()$genes
        ),
        "No binding profile files"
    )
})

context("Data Loading: load_data_genes")

# This test uses the provided Bsh TF binding data in extdata
test_that("load_data_genes loads and processes Bsh extdata (log2 ratio) files correctly", {
    data_dir <- system.file("extdata", package = "damidBind")
    binding_profiles_path <- data_dir

    # Use a dummy ensdb_genes, which would be the reference for transcripts
    dummy_ensdb <- make_dummy_ensdb_genes()$genes

    dl <- load_data_genes(
        binding_profiles_path = binding_profiles_path,
        ensdb_genes = dummy_ensdb,
        quantile_norm = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Check overall structure
    expect_type(dl, "list")
    expect_named(dl, c("binding_profiles_data", "occupancy", "test_category"))

    # Check binding_profiles_data
    expect_s4_class(dl$binding_profiles_data, "GRanges")
    expect_true(length(dl$binding_profiles_data) > 0)
    expect_true(all(c(
        "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm"
    ) %in%
        names(mcols(dl$binding_profiles_data))))

    # Check occupancy
    expect_s3_class(dl$occupancy, "data.frame")
    expect_true(nrow(dl$occupancy) > 0)
    # Filenames in occupancy should not have ".gatc"
    expect_true(all(c(
        "name", "nfrags",
        "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
        "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm",
        "gene_name", "gene_id"
    ) %in% colnames(dl$occupancy)))

    # Verify some expected values from the occupancy dataframe
    expect_true(any(grepl("gene", dl$occupancy$gene_name, ignore.case = TRUE)))
    expect_true(any(grepl("FBgn", dl$occupancy$gene_id, ignore.case = TRUE)))

    # Check test_category
    expect_equal(dl$test_category, "expressed")
})

test_that("load_data_genes works with quantile_norm = TRUE", {
    data_dir <- system.file("extdata", package = "damidBind")
    binding_profiles_path <- data_dir

    dl_qnorm <- load_data_genes(
        binding_profiles_path = binding_profiles_path,
        ensdb_genes = make_dummy_ensdb_genes()$genes,
        quantile_norm = TRUE,
        BPPARAM = BiocParallel::SerialParam()
    )

    # Check quantile normalisation in the binding_profiles_data object
    qnorm_mcols <- mcols(dl_qnorm$binding_profiles_data)
    expect_true(all(grepl("_qnorm$", names(qnorm_mcols))))

    # Check that occupancy values were calculated from the normalised data
    # The occupancy df colnames should have the _qnorm suffix
    expect_true(all(grepl("_qnorm", grep("Bsh_Dam_L", colnames(dl_qnorm$occupancy), value = TRUE))))

    # Check that values are indeed normalised
    signal_cols_qnorm <- dl_qnorm$occupancy %>%
        select(ends_with("_qnorm")) %>%
        as.matrix()

    col_means <- colMeans(signal_cols_qnorm, na.rm = TRUE)
    expect_true(sd(col_means) < 0.5)
    sample_medians <- apply(signal_cols_qnorm, 2, median, na.rm = TRUE)
    expect_true(sd(sample_medians) < 0.5)
})
