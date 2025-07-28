#' Compute average occupancy for each gene, using DamID binding profiles
#'
#' For each gene in the provided annotation, this function calculates the average
#' log2 score by taking a weighted mean of all overlapping fragments from the
#' binding data. The weighting is based on the width of the overlap between
#' each fragment and the gene body.
#'
#' @param binding_data A data.frame as produced by `build_dataframes()`. It must
#'   contain columns 'chr', 'start', 'end', and then one numeric column per sample.
#' @param ensdb_genes A GRanges object of gene annotations. A custom GRanges
#'   object can be provided, but it must contain metadata columns 'gene_name'
#'   and 'gene_id'. If not provided, it is retrieved via `get_ensdb_genes`.
#' @param buffer Integer. Number of bases to extend gene boundaries on both sides
#'   before calculating overlaps. Default is 0.
#' @param BPPARAM A BiocParallel parameter object for parallel computation.
#'   Default: `BiocParallel::bpparam()`.
#'
#' @return A data.frame where rows are genes and columns include one per sample
#'   (containing the weighted mean occupancy), plus annotation columns such as
#'   gene name, number of overlapping fragments ('nfrags'), and gene ID.
#' @export
#' @examples
#' # 1. Create a GRanges object for gene annotations
#' genes_gr <- GenomicRanges::GRanges(
#'   "chrY",
#'   IRanges::IRanges(c(1000, 3000), width = 500),
#'   gene_name = c("gene1", "gene2"),
#'   gene_id = c("ID1", "ID2")
#' )
#'
#' # 2. Create mock binding data
#' binding_df <- data.frame(
#'   chr = "chrY",
#'   start = c(950, 1200, 3100),
#'   end = c(1050, 1300, 3200),
#'   sampleA = c(1.5, 2.0, 0.5),
#'   sampleB = c(1.8, 2.2, 0.4)
#' )
#'
#' # 3. Calculate average gene occupancy
#' # Use BiocParallel::SerialParam() for deterministic execution in examples
#' if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'   gene_occ <- gene_occupancy(binding_df, genes_gr,
#'                              BPPARAM = BiocParallel::SerialParam())
#'   print(gene_occ)
#' }
gene_occupancy <- function(
    binding_data,
    ensdb_genes = get_ensdb_genes()$genes,
    buffer = 0,
    BPPARAM = BiocParallel::bpparam()) {
  if (!is(binding_data, "data.frame")) stop("'binding_data' must be a data.frame as returned by build_dataframes().")
  if (!is(ensdb_genes, "GRanges")) stop("'ensdb_genes' must be a GRanges object.")
  if (is.null(mcols(ensdb_genes)$gene_name)) stop("ensdb_genes must have a metadata column 'gene_name' with gene names.")

  message("Calculating average occupancy per gene ...")

  # Convert fragments to GRanges
  frag_gr <- GRanges(
    seqnames = binding_data$chr,
    IRanges(start = binding_data$start, end = binding_data$end)
  )

  # Extend genes by buffer
  gene_gr <- suppressWarnings(
    trim(
      resize(
        ensdb_genes,
        width = width(ensdb_genes) + buffer * 2,
        fix = "center"
      )
    )
  )

  resolve_dups <- function(gr) {
    key <- sprintf("%s:%d-%d", seqnames(gr), start(gr), end(gr))
    dup_pos <- 1
    while (any(duplicated(key))) {
      dups <- duplicated(key)
      key[dups] <- sprintf("%s:%d-%d", seqnames(gr[dups]), start(gr[dups]), end(gr[dups]) + dup_pos)
      dup_pos <- dup_pos + 1
    }
    return(key)
  }

  locs <- resolve_dups(gene_gr)
  mcols(gene_gr)$gene_loc <- locs

  # Find overlaps between genes and fragments
  overlaps <- findOverlaps(gene_gr, frag_gr)
  # For each gene, aggregate overlapping fragment scores (weighted)
  gene_names <- as.character(mcols(gene_gr)$gene_name)
  gene_id <- as.character(mcols(gene_gr)$gene_id)
  gene_loc <- as.character(mcols(gene_gr)$gene_loc)
  sample_cols <- setdiff(colnames(binding_data), c("chr", "start", "end"))

  # For parallel aggregation
  calc_occupancy <- function(i) {
    gene_hits <- subjectHits(overlaps)[queryHits(overlaps) == i]
    if (length(gene_hits) == 0) {
      return(rep(NA, length(sample_cols) + 4))
    }

    # Get covered fragment intervals and intersect size with gene
    gene_range <- ranges(gene_gr)[i]
    # chr <- seqnames(gene_gr)[i]
    # gene_loci <- mcols(gene_gr)$gene_loc[i]

    frag_starts <- binding_data$start[gene_hits]
    frag_ends <- binding_data$end[gene_hits]

    # Trim frags to gene body
    overlap_starts <- pmax(start(gene_range), frag_starts)
    overlap_ends <- pmin(end(gene_range), frag_ends)

    # Fragment lengths for weighted mean
    blens <- overlap_ends - overlap_starts + 1
    if (any(blens < 1)) blens[blens < 1] <- 0

    res <- vapply(sample_cols, function(col) {
      vals <- binding_data[[col]][gene_hits]
      weighted.mean(vals, w = blens, na.rm = TRUE)
    }, FUN.VALUE = numeric(1))

    # loc/name, nfrags, Averages,  name, gene_id
    c(gene_loc[i], length(gene_hits), res, gene_names[i], gene_id[i])
  }

  results <- BiocParallel::bplapply(seq_along(gene_gr), calc_occupancy, BPPARAM = BPPARAM)
  results <- do.call(rbind, results)
  colnames(results)[1:2] <- c("name", "nfrags")
  colnames(results)[(ncol(results) - 1):ncol(results)] <- c("gene_names", "gene_ids")
  results_df <- as.data.frame(results)
  results_df[, 2:(ncol(results_df) - 2)] <- apply(results_df[, 2:(ncol(results_df) - 2)], 2, as.numeric)
  results_df <- na.omit(results_df)

  # Rownames
  rownames(results_df) <- results_df$name

  results_df
}


#' Compute occupancy (average) per region over a set of GRanges
#'
#' For each interval in the `regions` GRanges object, this function finds all
#' overlapping fragments in `binding_data` and computes a weighted mean of their
#' signal values. The mean is weighted by the length of the overlap between each
#' fragment and the region.
#'
#' @param binding_data A data.frame as returned by `build_dataframes()`. It must
#'   contain columns 'chr', 'start', 'end', followed by numeric sample columns.
#' @param regions A GRanges object of genomic intervals (e.g., reduced peaks)
#'   over which to calculate occupancy.
#' @param buffer Optional integer. Number of base pairs to expand each interval
#'   in `regions` on both sides before calculating occupancy. Default is 0.
#' @param BPPARAM A BiocParallel parameter object for parallel computation.
#'   Default is `BiocParallel::bpparam()`.
#'
#' @return A data.frame with one row per region from the input `regions` object.
#'   Columns include 'name', 'nfrags' (number of overlapping fragments), and one
#'   column for the calculated weighted mean occupancy for each sample.
#' @export
#' @examples
#' # 1. Create a set of regions (e.g., reduced peaks)
#' regions_gr <- GenomicRanges::GRanges("chrX", IRanges::IRanges(start = c(100, 500), width = 100))
#' S4Vectors::mcols(regions_gr)$name <- paste0(
#'   GenomicRanges::seqnames(regions_gr), ":",
#'   GenomicRanges::start(regions_gr), "-", GenomicRanges::end(regions_gr)
#' )
#'
#' # 2. Create a mock binding data data.frame
#' binding_df <- data.frame(
#'   chr = "chrX",
#'   start = c(90, 150, 480, 550),
#'   end = c(110, 170, 520, 580),
#'   sampleA = c(1.2, 0.8, 2.5, 3.0),
#'   sampleB = c(1.0, 0.9, 2.8, 2.9)
#' )
#'
#' # 3. Calculate occupancy over the regions
#' # Use BiocParallel::SerialParam() for deterministic execution in examples
#' if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'   occupancy_data <- gr_occupancy(binding_df, regions_gr,
#'                                  BPPARAM = BiocParallel::SerialParam())
#'   print(occupancy_data)
#' }
gr_occupancy <- function(
    binding_data,
    regions,
    buffer = 0,
    BPPARAM = BiocParallel::bpparam()) {
  if (!is(binding_data, "data.frame")) stop("'binding_data' must be a data.frame as from build_dataframes().")
  if (!is(regions, "GRanges")) stop("'regions' must be a GRanges object.")

  message("Calculating average occupancy per region ...")
  if (buffer != 0) {
    regions <- suppressWarnings(
      trim(
        resize(
          regions,
          width = width(regions) + buffer * 2,
          fix = "center"
        )
      )
    )
  }
  frag_gr <- GRanges(
    seqnames = binding_data$chr,
    IRanges(start = binding_data$start, end = binding_data$end)
  )
  sample_cols <- setdiff(colnames(binding_data), c("chr", "start", "end"))

  # Name regions if they don't have names
  if (!"name" %in% names(mcols(regions))) {
    mcols(regions)$name <- paste0(seqnames(regions), ":", start(regions), "-", end(regions))
  }

  # Perform a single findOverlaps call - this is the key efficiency gain
  overlaps <- findOverlaps(regions, frag_gr)

  calc_occupancy <- function(i) {
    region_hits <- subjectHits(overlaps)[queryHits(overlaps) == i]
    if (length(region_hits) == 0) {
      return(rep(NA, length(sample_cols) + 2))
    }

    region_range <- ranges(regions)[i]
    frag_starts <- binding_data$start[region_hits]
    frag_ends <- binding_data$end[region_hits]

    overlap_starts <- pmax(start(region_range), frag_starts)
    overlap_ends <- pmin(end(region_range), frag_ends)

    blens <- overlap_ends - overlap_starts + 1
    if (any(blens < 1)) blens[blens < 1] <- 0

    avg_scores <- vapply(sample_cols, function(col) {
      vals <- binding_data[[col]][region_hits]
      weighted.mean(vals, w = blens, na.rm = TRUE)
    }, FUN.VALUE = numeric(1))

    return(c(mcols(regions)$name[i], length(region_hits), avg_scores))
  }

  results <- BiocParallel::bplapply(seq_along(regions), calc_occupancy, BPPARAM = BPPARAM)
  results <- do.call(rbind, results)

  # Handle cases where no regions had overlaps
  if (is.null(results) || nrow(results) == 0) {
    # Return an empty dataframe with correct structure
    outdf <- data.frame(
      name = character(0),
      nfrags = integer(0),
      stringsAsFactors = FALSE
    )
    for (col in sample_cols) outdf[[col]] <- numeric(0)
    return(outdf)
  }

  results_df <- as.data.frame(results, stringsAsFactors = FALSE)
  colnames(results_df) <- c("name", "nfrags", sample_cols)

  # Convert appropriate columns to numeric
  numeric_cols <- c("nfrags", sample_cols)
  results_df[, numeric_cols] <- lapply(results_df[, numeric_cols], as.numeric)

  results_df <- na.omit(results_df)
  rownames(results_df) <- results_df$name
  return(results_df)
}

#' Reduce a list of GRanges to unique, non-overlapping regions
#'
#' Takes a list of GRanges objects (e.g., peak sets from multiple samples),
#' combines them, and merges any overlapping or adjacent regions into a single,
#' minimal set of genomic intervals.
#'
#' @param peaks A list of GRanges objects.
#' @return A GRanges object containing the reduced (union) regions, with a
#'   `name` metadata column in the format "chr:start-end".
#' @export
#' @examples
#' # Create a list of GRanges objects with overlapping regions
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 200), width = 50))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(120, 300), width = 50))
#' gr_list <- list(gr1,gr2)
#'
#' # Reduce the list to a single set of non-overlapping regions
#' reduced <- reduce_regions(gr_list)
#' print(reduced)
#' # The result combines overlapping regions [100-149] and [120-169] into [100-169].
reduce_regions <- function(peaks) {
  if (!is.list(peaks) || !all(vapply(peaks, function(x) is(x, "GRanges"), FUN.VALUE = logical(1)))) {
    stop("'peaks' must be a list of GRanges objects.")
  }
  combined <- do.call(c, unname(peaks))
  pr <- reduce(combined)
  # Add names for intervals
  mcols(pr)$name <- paste0(seqnames(pr), ":", start(pr), "-", end(pr))
  pr
}

#' Find best gene overlap(s) for each query interval
#'
#' Annotates each input region with the gene(s) it overlaps. A gene is considered
#' overlapping if its body is within the specified `maxgap` of the query region.
#' If a query region overlaps multiple genes, their names and IDs are returned as
#' a comma-separated string.
#'
#' @param query A GRanges object containing the regions to be annotated.
#' @param subject A GRanges object of gene annotations. It must have metadata
#'   columns named `gene_name` and, optionally, `gene_id`.
#' @param maxgap Integer. The maximum number of base pairs between the query and
#'   subject for them to be considered overlapping. Default is 0 (must be touching).
#'
#' @return A list containing two character vectors of the same length as `query`:
#'   \item{genes}{A character vector where each element contains a comma-separated
#'   list of `gene_name` values from subject regions overlapping the corresponding
#'   query region. An empty string `""` indicates no overlap.}
#'   \item{ids}{A character vector with the corresponding `gene_id` values, if
#'   the `gene_id` column exists in the subject.}
#'
#' @export
#' @examples
#' # Create a query GRanges object with regions of interest
#' query_regions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 500), width = 50))
#'
#' # Create a subject GRanges object with gene annotations
#' gene_annotations <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(90, 200, 525), width = c(30, 50, 50)),
#'   gene_name = c("geneA", "geneB", "geneC"),
#'   gene_id = c("FBgn01", "FBgn02", "FBgn03")
#' )
#'
#' # Find overlaps (query 1 overlaps geneA; query 2 overlaps geneC)
#' overlaps <- all_overlaps_to_original(query_regions, gene_annotations, maxgap = 0)
#' print(overlaps)
#'
#' # With a larger gap, query 1 now also overlaps geneB
#' overlaps_gapped <- all_overlaps_to_original(query_regions, gene_annotations, maxgap = 50)
#' print(overlaps_gapped)
all_overlaps_to_original <- function(query, subject, maxgap = 0) {
  if (!is(query, "GRanges")) stop("'query' must be a GRanges object.")
  if (!is(subject, "GRanges")) stop("'subject' must be a GRanges object.")
  if (is.null(mcols(subject)$gene_name)) stop("Subject GRanges must have a 'gene_name' metadata column.")

  ol <- findOverlaps(query, subject, maxgap = maxgap)
  hits <- as.data.frame(ol)
  by_query <- split(subject$gene_name[hits$subjectHits], hits$queryHits)
  # For queries with no overlap, return ""
  query_idx <- seq_along(query)
  out <- character(length(query))
  for (i in query_idx) {
    if (as.character(i) %in% names(by_query)) {
      out[i] <- paste(sort(unique(by_query[[as.character(i)]])), collapse = ",")
    } else {
      out[i] <- ""
    }
  }


  gene_id_out <- character(length(query))
  if (!is.null(mcols(subject)$gene_id)) {
    by_query <- split(subject$gene_id[hits$subjectHits], hits$queryHits)
    query_idx <- seq_along(query)
    names_out <- character(length(query))
    for (i in query_idx) {
      if (as.character(i) %in% names(by_query)) {
        gene_id_out[i] <- paste(sort(unique(by_query[[as.character(i)]])), collapse = ",")
      } else {
        gene_id_out[i] <- ""
      }
    }
  }

  return(list(genes = out, ids = gene_id_out))
}
