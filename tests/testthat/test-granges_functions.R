library(testthat)
library(damidBind)
library(GenomicRanges)
library(IRanges)
library(BiocParallel) # For BPPARAM

context("GRanges Utilities")

test_that("reduce_regions correctly combines and reduces GRanges", {
    gr1 <- GRanges("chr1", IRanges(c(1, 10, 20), width = 5))
    gr2 <- GRanges("chr1", IRanges(c(3, 12, 35), width = 7))
    gr3 <- GRanges("chr2", IRanges(c(1, 20), width = 10))
    gr4 <- GRanges("chr2", IRanges(14,22))

    reduced_gr <- suppressWarnings(reduce_regions(list(gr1, gr2, gr3, gr4)))

    # Expect 5 unique, non-overlapping regions:
    # GRanges object with 5 ranges and 1 metadata column:
    #     seqnames    ranges strand |        name
    #        <Rle> <IRanges>  <Rle> | <character>
    # [1]     chr1      1-18      * |   chr1:1-18
    # [2]     chr1     20-24      * |  chr1:20-24
    # [3]     chr1     35-41      * |  chr1:35-41
    # [4]     chr2      1-10      * |   chr2:1-10
    # [5]     chr2     14-29      * |  chr2:14-29

    expect_s4_class(reduced_gr, "GRanges")
    expect_equal(length(reduced_gr), 5)
    expect_equal(as.character(seqnames(reduced_gr)), c(rep("chr1",3), rep("chr2",2)))
    expect_equal(start(reduced_gr), c(1, 20, 35, 1, 14))
    expect_equal(end(reduced_gr), c(18, 24, 41, 10, 29))
    expect_true("name" %in% names(mcols(reduced_gr)))
    expect_equal(mcols(reduced_gr)$name, c("chr1:1-18", "chr1:20-24", "chr1:35-41",
                                           "chr2:1-10", "chr2:14-29"))
})

test_that("gene_occupancy calculates weighted mean correctly", {
    # Mock ensdb_genes for stability
    dummy_genes_gr <- GRanges(
        seqnames = "chr1",
        ranges = IRanges(c(1000, 3000), width = c(501, 501)),
        gene_id = c("TestGene1", "TestGene2"),
        gene_name = c("GeneA", "GeneB")
    )

    # Binding data with specific values for calculation
    binding_gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = c(950, 1050, 1400, 2900, 3200, 3400),
            end = c(1049, 1149, 1499, 2999, 3299, 3499)
        ),
        val1 = c(10, 20, 30, 5, 15, 25),
        val2 = c(11, 21, 31, 6, 16, 26)
    )

    res <- calculate_occupancy(binding_gr, dummy_genes_gr, BPPARAM = BiocParallel:::SerialParam()) # Use serial for test consistency

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 2)
    expect_true(all(c("name", "nfrags", "val1", "val2", "gene_name", "gene_id") %in% colnames(res)))

    # Check GeneA
    expect_equal(res["chr1:1000-1500", "val1"], 22)
    expect_equal(res["chr1:1000-1500", "val2"], 23)
    expect_equal(res["chr1:1000-1500", "nfrags"], 3)
    expect_equal(res["chr1:1000-1500", "gene_name"], "GeneA")
    expect_equal(res["chr1:1000-1500", "gene_id"], "TestGene1")

    # Check GeneB
    expect_equal(res["chr1:3000-3500", "val1"], 20)
    expect_equal(res["chr1:3000-3500", "val2"], 21)
    expect_equal(res["chr1:3000-3500", "nfrags"], 2)
    expect_equal(res["chr1:3000-3500", "gene_name"], "GeneB")
    expect_equal(res["chr1:3000-3500", "gene_id"], "TestGene2")
})


test_that("calculate_occupancy calculates weighted mean for arbitrary regions", {
    regions_gr <- GRanges("chr1", IRanges(start = c(100, 300), end = c(150, 350)))
    binding_gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = c(90, 110, 280, 305),
            end = c(110, 130, 310, 320)
        ),
        signal_a = c(10, 20, 5, 15),
        signal_b = c(100, 200, 50, 150)
    )

    res <- calculate_occupancy(binding_gr, regions_gr, BPPARAM = BiocParallel::SerialParam())

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

    # Check with maxgap=0
    # Query 1 (10-14) overlaps GENE_A (5-14)
    # Query 2 (50-54) overlaps GENE_B (40-54)
    # Query 3 (100-104) overlaps GENE_C (90-109)
    # Query 4 (200-204) no overlap
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
    expect_equal(overlaps_multi$genes, "GENE_X,GENE_Y")
    expect_equal(overlaps_multi$ids, "ID_X,ID_Y")

    # Test with maxgap = 5
    # Query 4 (200-204) with subject GENE_D (210-219)
    overlaps_gap <- all_overlaps_to_original(query_gr, subject_gr, maxgap = 5)
    expect_equal(overlaps_gap$genes[4], "GENE_D")
    expect_equal(overlaps_gap$ids[4], "ID_D")
})
