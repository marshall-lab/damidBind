library(testthat)
library(damidBind)
library(GenomicRanges)
library(IRanges)
library(BiocParallel) # For BPPARAM

context("GRanges Utilities")

test_that("reduce_regions correctly combines and reduces GRanges", {
  gr1 <- GRanges("chr1", IRanges(c(1, 10, 20), width = 5))
  gr2 <- GRanges("chr1", IRanges(c(3, 12, 18), width = 7)) # Overlaps gr1
  gr3 <- GRanges("chr2", IRanges(c(1, 5), width = 10))

  reduced_gr <- suppressWarnings(reduce_regions(list(gr1, gr2, gr3)))

  # Expect 3 unique, non-overlapping regions:
  # chr1:1-26 (union of 1-5, 3-9, 10-14, 12-18, 18-24, 20-24)
  # chr2:1-10
  # chr2:5-14
  # No, reduce function will combine overlapping regions on the same chromosome.
  # So chr1:1-26 means region 1-5 from gr1, 3-9 from gr2, 10-14 from gr1, 12-18 from gr2, 18-24 from gr2, 20-24 from gr1.
  # Max range on chr1 is 1-26.
  # From gr1: 1-5, 10-14, 20-24
  # From gr2: 3-9, 12-18, 18-24
  # Union on chr1: 1-9 (from 1-5, 3-9), 10-18 (from 10-14, 12-18), 18-24 (from gr2), 20-24 (from gr1)
  # reduce(c(gr1_chr1,gr2_chr1)) will give one region: 1-24. No, will be 1-9 & 10-24.
  # reduce(GRanges("chr1", IRanges(c(1,10,20),width=5)),GRanges("chr1", IRanges(c(3,12,18),width=7)))
  # GRanges object with 6 ranges and 0 metadata columns:
  #       seqnames    ranges strand
  #          <Rle> <IRanges>  <Rle>
  #   [1]     chr1     [1,5]      *
  #   [2]     chr1    [10,14]      *
  #   [3]     chr1    [20,24]      *
  #   [4]     chr1     [3,9]      *
  #   [5]     chr1    [12,18]      *
  #   [6]     chr1    [18,24]      *
  # reduce(combined) -> GRanges object with 3 ranges and 0 metadata columns:
  #   [1]     chr1     [1,9]      *
  #   [2]     chr1    [10,18]      *
  #   [3]     chr1    [18,24]      *
  # The output of reduce is minimal, non-overlapping regions spanning the original ranges.
  # The interval 18-24 is overlapped by 12-18. So 10-18 and 18-24 could be 10-24.
  # `reduce` combines overlapping. The combination of all these is chr1:1-24.
  # chr1: 1-5, 3-9, 10-14, 12-18, 18-24, 20-24. Max start=1, max end=24.
  # Only one region chr1: 1-24, and one region chr2: 1-14.
  # Corrected expected:
  # chr1: 1-24 (union of 1-5, 3-9, 10-14, 12-18, 18-24, 20-24)
  # chr2: 1-14 (union of 1-10, 5-14)
  expect_s4_class(reduced_gr, "GRanges")
  expect_equal(length(reduced_gr), 2)
  expect_equal(as.character(seqnames(reduced_gr)), c("chr1", "chr2"))
  expect_equal(start(reduced_gr), c(1, 1))
  expect_equal(end(reduced_gr), c(24, 14)) # 1-10 and 5-14 make 1-14
  expect_true("name" %in% names(mcols(reduced_gr)))
  expect_equal(mcols(reduced_gr)$name, c("chr1:1-24", "chr2:1-14"))
})

test_that("gene_occupancy calculates weighted mean correctly", {
  # Mock ensdb_genes for stability
  dummy_genes_gr <- GRanges(
    seqnames = "chr1",
    # Change width to 501 to make 'end' 1500 and 3500, matching expected row names in tests
    ranges = IRanges(c(1000, 3000), width = c(501, 501)),
    gene_id = c("TestGene1", "TestGene2"),
    gene_name = c("GeneA", "GeneB")
  )

  # Binding data with specific values for calculation
  binding_df <- data.frame(
    chr = "chr1",
    start = c(950, 1050, 1400, 2900, 3200, 3400),
    end = c(1049, 1149, 1499, 2999, 3299, 3499),
    val1 = c(10, 20, 30, 5, 15, 25),
    val2 = c(11, 21, 31, 6, 16, 26)
  )

  # Calculate expected for GeneA (chr1:1000-1500)
  # Fragments:
  # 1: 950-1049 (overlaps 1000-1049 -> 50bp) val=10
  # 2: 1050-1149 (overlaps 1000-1149 -> 1000-1149 -> 100bp if within, but interval length is 100) val=20
  # 3: 1400-1499 (overlaps 1400-1499 -> 100bp) val=30
  # Overlap with GeneA (1000-1500):
  # frag1: 1000-1049, length=50. val1=10, val2=11
  # frag2: 1050-1149, length=100. val1=20, val2=21
  # frag3: 1400-1499, length=100. val1=30, val2=31
  # Expected weighted mean for val1 (GeneA): (10*50 + 20*100 + 30*100) / (50 + 100 + 100) = (500 + 2000 + 3000) / 250 = 5500 / 250 = 22
  # Expected weighted mean for val2 (GeneA): (11*50 + 21*100 + 31*100) / (50 + 100 + 100) = (550 + 2100 + 3100) / 250 = 5750 / 250 = 23

  # Calculate expected for GeneB (chr1:3000-3500)
  # Fragments:
  # 4: 2900-2999 (overlaps 3000-3000 -> 0bp)
  # 5: 3200-3299 (overlaps 3200-3299 -> 100bp) val=15
  # 6: 3400-3499 (overlaps 3400-3499 -> 100bp) val=25
  # Expected weighted mean for val1 (GeneB): (15*100 + 25*100) / (100 + 100) = (1500 + 2500) / 200 = 4000 / 200 = 20
  # Expected weighted mean for val2 (GeneB): (16*100 + 26*100) / (100 + 100) = (1600 + 2600) / 200 = 4200 / 200 = 21

  res <- gene_occupancy(binding_df, ensdb_genes = dummy_genes_gr, BPPARAM = BiocParallel:::SerialParam()) # Use serial for test consistency

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("name", "nfrags", "val1", "val2", "gene_names", "gene_ids") %in% colnames(res)))

  # Check GeneA
  expect_equal(res["chr1:1000-1500", "val1"], 22)
  expect_equal(res["chr1:1000-1500", "val2"], 23)
  expect_equal(res["chr1:1000-1500", "nfrags"], 3)
  expect_equal(res["chr1:1000-1500", "gene_names"], "GeneA")
  expect_equal(res["chr1:1000-1500", "gene_ids"], "TestGene1")

  # Check GeneB
  expect_equal(res["chr1:3000-3500", "val1"], 20)
  expect_equal(res["chr1:3000-3500", "val2"], 21)
  expect_equal(res["chr1:3000-3500", "nfrags"], 2)
  expect_equal(res["chr1:3000-3500", "gene_names"], "GeneB")
  expect_equal(res["chr1:3000-3500", "gene_ids"], "TestGene2")
})


test_that("gr_occupancy calculates weighted mean for arbitrary regions", {
  regions_gr <- GRanges("chr1", IRanges(start = c(100, 300), end = c(150, 350)))
  binding_df <- data.frame(
    chr = "chr1",
    start = c(90, 110, 280, 305),
    end = c(110, 130, 310, 320),
    signal_a = c(10, 20, 5, 15),
    signal_b = c(100, 200, 50, 150)
  )

  # Region 1: chr1:100-150
  # Frag1: 90-110 (overlap 100-110, len=11) -> signal_a=10, signal_b=100
  # Frag2: 110-130 (overlap 110-130, len=21) -> signal_a=20, signal_b=200
  # No, original code: (overlap_ends - overlap_starts + 1)
  # For frag1 (90-110) overlapping region (100-150):
  # pmax(100, 90) = 100
  # pmin(150, 110) = 110
  # blens = 110 - 100 + 1 = 11
  # For frag2 (110-130) overlapping region (100-150):
  # pmax(100, 110) = 110
  # pmin(150, 130) = 130
  # blens = 130 - 110 + 1 = 21
  # Expected for signal_a in region1: (10*11 + 20*21) / (11 + 21) = (110 + 420) / 32 = 530 / 32 = 16.5625
  # Expected for signal_b in region1: (100*11 + 200*21) / (11 + 21) = (1100 + 4200) / 32 = 5300 / 32 = 165.625

  # Region 2: chr1:300-350
  # Frag3: 280-310 (overlap 300-310, len=11) -> signal_a=5, signal_b=50
  # Frag4: 305-320 (overlap 305-320, len=16) -> signal_a=15, signal_b=150
  # Expected for signal_a in region2: (5*11 + 15*16) / (11 + 16) = (55 + 240) / 27 = 295 / 27 = 10.9259
  # Expected for signal_b in region2: (50*11 + 150*16) / (11 + 16) = (550 + 2400) / 27 = 2950 / 27 = 109.259

  res <- gr_occupancy(binding_df, regions_gr, BPPARAM = BiocParallel::SerialParam())

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("name", "nfrags", "signal_a", "signal_b") %in% colnames(res)))

  expect_equal(res["chr1:100-150", "signal_a"], 16.5625)
  expect_equal(res["chr1:100-150", "signal_b"], 165.625)
  expect_equal(res["chr1:100-150", "nfrags"], 2)

  expect_equal(res["chr1:300-350", "signal_a"], 10.925925925925925)
  expect_equal(res["chr1:300-350", "signal_b"], 109.25925925925925)
  expect_equal(res["chr1:300-350", "nfrags"], 2)
})

test_that("all_overlaps_to_original correctly annotates overlaps", {
  query_gr <- GRanges("chr1", IRanges(c(10, 50, 100, 200), width = 5))
  subject_gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(c(5, 40, 90, 210), width = c(10, 15, 20, 10)),
    gene_name = c("GENE_A", "GENE_B", "GENE_C", "GENE_D"),
    gene_id = c("ID_A", "ID_B", "ID_C", "ID_D")
  )

  # Expected overlaps:
  # Query 1 (10-14) overlaps GENE_A (5-14)
  # Query 2 (50-54) overlaps GENE_B (40-54)
  # Query 3 (100-104) overlaps GENE_C (90-109)
  # Query 4 (200-204) no overlap close by
  # Check with maxgap=0
  overlaps <- all_overlaps_to_original(query_gr, subject_gr, maxgap = 0)

  expect_type(overlaps, "list")
  expect_named(overlaps, c("genes", "ids"))
  expect_equal(overlaps$genes, c("GENE_A", "GENE_B", "GENE_C", ""))
  expect_equal(overlaps$ids, c("ID_A", "ID_B", "ID_C", ""))

  # Test with multiple overlaps for a single query (comma separated)
  subject_gr_multi <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(c(5, 8, 100), width = 10),
    gene_name = c("GENE_X", "GENE_Y", "GENE_Z"),
    gene_id = c("ID_X", "ID_Y", "ID_Z")
  )
  query_gr_multi <- GRanges("chr1", IRanges(8, width = 5)) # overlaps GENE_X (5-14) and GENE_Y (8-17)

  overlaps_multi <- all_overlaps_to_original(query_gr_multi, subject_gr_multi)
  expect_equal(overlaps_multi$genes, "GENE_X,GENE_Y") # Sorted alphabetically
  expect_equal(overlaps_multi$ids, "ID_X,ID_Y")

  # Test with maxgap
  # Query 4 (200-204) with subject GENE_D (210-219)
  # Gap is 210 - 204 - 1 = 5.
  overlaps_gap <- all_overlaps_to_original(query_gr, subject_gr, maxgap = 5)
  expect_equal(overlaps_gap$genes[4], "GENE_D")
  expect_equal(overlaps_gap$ids[4], "ID_D")
})
