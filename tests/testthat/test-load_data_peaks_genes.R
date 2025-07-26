library(testthat)
library(damidBind)
library(GenomicRanges) # For GRanges object checks
library(S4Vectors) # For Rle if needed in comparisons

# Source helper functio
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
    quantile_norm = FALSE, # Test without quantile normalization first
    BPPARAM = BiocParallel::SerialParam()
  )

  # Check overall structure of the returned list
  expect_type(dl, "list")
  expect_named(dl, c("binding_profiles_data", "peaks", "pr", "occupancy", "test_category"))

  # Check binding_profiles_data
  expect_s3_class(dl$binding_profiles_data, "data.frame")
  expect_true(nrow(dl$binding_profiles_data) > 0)
  expect_true(all(c(
    "chr", "start", "end",
    "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm"
  ) %in%
    colnames(dl$binding_profiles_data)))

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
    "gene_names", "gene_ids"
  ) %in%
    colnames(dl$occupancy)))
  expect_true(all(rownames(dl$occupancy) == dl$occupancy$name))

  # Check test_category
  expect_equal(dl$test_category, "bound")

  # Check for gene annotations (from dummy ensdb_genes)
  # Actual extdata has genes on chr2L. Dummy ensdb_genes also has chr2L.
  # So these should largely overlap if your ranges are correctly defined in the dummy.
  expect_true(any(dl$occupancy$gene_names != ""))
  expect_true(any(dl$occupancy$gene_ids != ""))
})

test_that("load_data_peaks works with quantile_norm = TRUE", {
  # temp_extdata_dir <- copy_extdata_to_temp()
  data_dir <- system.file("extdata", package = "damidBind")
  binding_profiles_path <- data_dir
  peaks_path <- data_dir

  dl_qnorm <- load_data_peaks(
    binding_profiles_path = binding_profiles_path,
    peaks_path = peaks_path,
    quantile_norm = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(all(c(
    "chr", "start", "end",
    "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm_qnorm",
    "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm_qnorm",
    "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm_qnorm",
    "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm_qnorm"
  ) %in%
    colnames(dl_qnorm$binding_profiles_data)))

  # Check that values are indeed normalized (e.g., all means are similar)
  # This implies the quantile_normalisation function itself is working,
  # and that the effect is seen in the signal columns.
  signal_cols_qnorm <- dl_qnorm$binding_profiles_data %>%
    dplyr::select(dplyr::starts_with("Bsh_Dam_L")) %>%
    as.matrix()

  col_means <- colMeans(signal_cols_qnorm)
  # Quantile normalization makes distributions similar; mean/median should converge
  expect_true(sd(col_means) < 0.5) # Adjusted for potentially noisy real data
  # Also check that sample medians are similar
  sample_medians <- apply(signal_cols_qnorm, 2, median)
  expect_true(sd(sample_medians) < 0.5)

  # unlink(temp_extdata_dir, recursive = TRUE)
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

# This test will now use the Bsh TF binding data for load_data_genes
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
  expect_s3_class(dl$binding_profiles_data, "data.frame")
  expect_true(nrow(dl$binding_profiles_data) > 0)
  expect_true(all(c(
    "chr", "start", "end",
    "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm"
  ) %in%
    colnames(dl$binding_profiles_data)))

  # Check occupancy
  expect_s3_class(dl$occupancy, "data.frame")
  expect_true(nrow(dl$occupancy) > 0)
  expect_true(all(c(
    "name", "nfrags",
    "Bsh_Dam_L4_r1-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L4_r2-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L5_r1-ext300-vs-Dam.kde-norm",
    "Bsh_Dam_L5_r2-ext300-vs-Dam.kde-norm",
    "gene_names", "gene_ids"
  ) %in% colnames(dl$occupancy)))

  # Verify some expected values from the occupancy dataframe, e.g., that gene names were joined
  expect_true(any(grepl("gene", dl$occupancy$gene_names, ignore.case = TRUE)))
  expect_true(any(grepl("FBgn", dl$occupancy$gene_ids, ignore.case = TRUE)))

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

  # Check that values are indeed normalized in the occupancy dataframe
  signal_cols_qnorm <- dl_qnorm$occupancy %>%
    dplyr::select(dplyr::starts_with("Bsh_Dam_L")) %>%
    as.matrix()

  col_means <- colMeans(signal_cols_qnorm)
  expect_true(sd(col_means) < 0.5)
  sample_medians <- apply(signal_cols_qnorm, 2, median)
  expect_true(sd(sample_medians) < 0.5)
})
