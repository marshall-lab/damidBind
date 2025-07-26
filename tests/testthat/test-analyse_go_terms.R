library(testthat)
library(damidBind)
library(clusterProfiler) # Required for mocking its functions
library(ggplot2)       # Required for mocking ggsave

context("GO Enrichment Analysis: analyse_go_enrichment")

dummy_org_db <- org.Dm.eg.db::org.Dm.eg.db

# Helper to create a dummy diff_results object for GO enrichment tests
make_dummy_diff_results_for_go <- function() {
  analysis_table <- data.frame(
    logFC = c(2.5, -3.0, 0.5, 1.8, -1.2, 0.1, 0, 0, 0, 0, 0),
    adj.P.Val = c(0.001, 0.005, 0.1, 0.002, 0.08, 0.5, 1, 1, 1, 1, 1),
    gene_names = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "geneG", "geneH", "geneI", "geneJ", "geneK"),
    gene_ids = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007", "FBgn008", "FBgn009", "FBgn010", "FBgn011"),
    row.names = c("Locus1", "Locus2", "Locus3", "Locus4", "Locus5", "Locus6", "Locus7", "Locus8", "Locus9", "Locus10", "Locus11")
  )

  upCond1 <- analysis_table[c("Locus1", "Locus4"), ] # Significant for Cond1
  upCond2 <- analysis_table[c("Locus2"), ]           # Significant for Cond2

  new("DamIDResults",
    upCond1 = upCond1,
    upCond2 = upCond2,
    analysis = analysis_table,
    cond = c(Cond1 = "Condition_X", Cond2 = "Condition_Y"),
    data = list(test_category = "bound")
  )
}

test_that("analyse_go_enrichment runs without error and returns expected structure", {
  diff_res <- make_dummy_diff_results_for_go()

  # Mock clusterProfiler functions
  # bitr mock: map some FBgnIDs to SYMBOLS
  mock_bitr_map <- data.frame(
    FLYBASE = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007", "FBgn008", "FBgn009", "FBgn010", "FBgn011"),
    SYMBOL = c("SYM1", "SYM2", "SYM3", "SYM4", "SYM5", "SYM6", "SYM7", "SYM8", "SYM9", "SYM10", "SYM11"),
    stringsAsFactors = FALSE
  )

  # enrichGO mock: return a dummy enrichResult object for successful run
  mock_enrich_go_result <- new("enrichResult",
    readable = TRUE,
    result = data.frame(
      ID = c("GO:0001", "GO:0002"),
      Description = c("Biological Process 1", "Biological Process 2"),
      GeneRatio = c("2/5", "3/5"),
      BgRatio = c("10/100", "15/100"),
      pvalue = c(0.001, 0.005),
      p.adjust = c(0.002, 0.006),
      qvalue = c(0.01, 0.02),
      geneID = c("SYM1/SYM4", "SYM1/SYM2/SYM4"),
      Count = c(2, 3),
      stringsAsFactors = FALSE
    ),
    ontology = "BP",
    universe = c("SYM1", "SYM2", "SYM4", "SYM5", "SYM6", "SYM7", "SYM8", "SYM9", "SYM10", "SYM11"),
    gene = c("SYM1", "SYM4"), # For Cond1
    keytype = "SYMBOL",
    organism = "Drosophila melanogaster"
  )

  local_mocked_bindings(
    bitr = function(geneID, fromType, toType, OrgDb) {
      # Filter mock_bitr_map based on input geneID
      subset(mock_bitr_map, FLYBASE %in% geneID)
    },
    enrichGO = function(gene, OrgDb, universe, keyType, ont, pAdjustMethod, maxGSSize, minGSSize, pvalueCutoff, qvalueCutoff, readable) {
      if (length(gene) == 0 || length(universe) == 0) return(NULL) # Simulate no enrichment
      mock_enrich_go_result
    }
  )
  local_mocked_bindings(
    ggsave = function(filename, plot, width, height, units, device) {
      # Simulate file creation to check if ggsave was called
      file.create(filename)
    },
    .package = "ggplot2"
  )

  # Test Cond1 (default direction)
  tmp_go_file <- tempfile(pattern = "go_report", fileext = ".csv")
  tmp_plot_file <- tempfile(pattern = "go_plot", fileext = ".pdf")

  res <- analyse_go_enrichment(
    diff_results = diff_res,
    org_db = dummy_org_db,
    save_results_path = tmp_go_file,
    save = list(filename = tools::file_path_sans_ext(tmp_plot_file), format = "pdf")
  )

  # Check overall output structure
  expect_type(res, "list")
  expect_named(res, c("enrich_go_object", "results_table", "dot_plot"))
  expect_s4_class(res$enrich_go_object, "enrichResult")
  expect_s3_class(res$results_table, "data.frame")
  expect_s3_class(res$dot_plot, "ggplot")

  # Check results table content based on mock
  expect_equal(res$results_table$ID, c("GO:0001", "GO:0002"))
  expect_equal(res$results_table$Description, c("Biological Process 1", "Biological Process 2"))
  expect_equal(res$results_table$Count, c(2, 3))
  expect_equal(res$results_table$GeneRatio, c(2 / 5, 3 / 5))

  # Check if results table and plot were 'saved' (i.e., files created by mock)
  expect_true(file.exists(tmp_go_file))
  expect_true(file.exists(tmp_plot_file))

  # Clean up temp files
  unlink(c(tmp_go_file, tmp_plot_file))

  # --- Test 'all' direction ---
  # Reset mock_enrich_go_result to simulate enrichment for combined set
  mock_enrich_go_result@gene <- c("SYM1", "SYM2", "SYM4") # All significant (Locus1, Locus2, Locus4)
  res_all <- analyse_go_enrichment(
    diff_results = diff_res,
    direction = "all",
    org_db = dummy_org_db,
    save = FALSE # Do not save for this specific test
  )
  expect_s4_class(res_all$enrich_go_object, "enrichResult")
  expect_equal(res_all$enrich_go_object@gene, c("SYM1", "SYM2", "SYM4")) # Check the genes passed to enrichGO

})


test_that("analyse_go_enrichment handles cases with no significant genes", {
  diff_res_no_sig <- make_dummy_diff_results_for_go()
  diff_res_no_sig@upCond1 <- diff_res_no_sig@upCond1[FALSE, ] # No significant genes
  diff_res_no_sig@upCond2 <- diff_res_no_sig@upCond2[FALSE, ]

  # Mock enrichGO to return NULL (no enrichment found, or no genes passed)
  local_mocked_bindings(
    bitr = function(geneID, fromType, toType, OrgDb) data.frame(), # No mapping
    enrichGO = function(...) NULL, # No enrichment
    .package = "clusterProfiler"
  )

  expect_message(res <- analyse_go_enrichment(diff_res_no_sig, org_db = dummy_org_db),
    "No significant genes found for direction 'cond1'") # First message
  expect_null(res)
})


test_that("analyse_go_enrichment validates org_db input", {
  diff_res <- make_dummy_diff_results_for_go()
  expect_error(
    analyse_go_enrichment(diff_results = diff_res, org_db = "not_an_OrgDb_object"),
    "'org_db' must be a valid OrgDb object"
  )
})

test_that("analyse_go_enrichment handles invalid direction", {
  diff_res <- make_dummy_diff_results_for_go()
  expect_error(
    analyse_go_enrichment(diff_results = diff_res, direction = "invalid", org_db = dummy_org_db),
    "Invalid 'direction' specified"
  )
})

test_that("analyse_go_enrichment alias analyze_go_enrichment works", {
  diff_res <- make_dummy_diff_results_for_go()

  # Apply the same mocks as the main test to ensure consistency
  mock_bitr_map <- data.frame(
    FLYBASE = c("FBgn001", "FBgn002", "FBgn004"),
    SYMBOL = c("SYM1", "SYM2", "SYM4"),
    stringsAsFactors = FALSE
  )
  mock_enrich_go_result <- new("enrichResult",
    readable = TRUE,
    result = data.frame(
      ID = c("GO:0001"),
      Description = c("Test Process"),
      GeneRatio = c("1/2"),
      BgRatio = c("10/100"),
      pvalue = c(0.01),
      p.adjust = c(0.01),
      qvalue = c(0.01),
      geneID = c("SYM1"),
      Count = c(1),
      stringsAsFactors = FALSE
    ),
    ontology = "BP",
    universe = c("SYM1", "SYM2", "SYM3", "SYM4", "SYM5", "SYM6"),
    gene = c("SYM1"),
    keytype = "SYMBOL",
    organism = "Drosophila melanogaster"
  )

  local_mocked_bindings(
    bitr = function(geneID, fromType, toType, OrgDb) subset(mock_bitr_map, FLYBASE %in% geneID),
    enrichGO = function(...) mock_enrich_go_result
  )
  local_mocked_bindings(
    ggsave = function(...) NULL,
    .package = "ggplot2"
  )

  res_alias <- analyze_go_enrichment(diff_results = diff_res, org_db = dummy_org_db)
  expect_s4_class(res_alias$enrich_go_object, "enrichResult")
  expect_equal(res_alias$results_table$Description, "Test Process")
})

# Test for clean_gene_symbols (moved from previous test-analyse_go_terms.R to clarify grouping)
test_that("clean_gene_symbols correctly filters out specified gene types", {
  # Direct call to the internal function
  clean_gs <- damidBind:::clean_gene_symbols

  # Basic filtering
  symbols <- c("gene1", "snoRNA_abc", "tRNA_XYZ", "mRNA_def", "snRNA_123")
  expect_equal(clean_gs(symbols), c("gene1", "mRNA_def"))

  # Case insensitivity
  symbols_case <- c("SnoRNA", "tRnA", "protein_xyz")
  expect_equal(clean_gs(symbols_case), "protein_xyz")

  # No symbols to filter
  symbols_none <- c("geneA", "geneB")
  expect_equal(clean_gs(symbols_none), c("geneA", "geneB"))

  # All symbols filtered
  symbols_all <- c("snoRNA", "tRNA")
  expect_equal(clean_gs(symbols_all), character(0))

  # Handling NULL or empty input
  expect_equal(clean_gs(NULL), character(0))
  expect_equal(clean_gs(character(0)), character(0))

  # With custom extra regex
  symbols_extra <- c("gene1", "mirna_ab", "rRNA_123")
  expect_equal(clean_gs(symbols_extra, clean_extra = "miRNA|rRNA"), c("gene1"))
})
