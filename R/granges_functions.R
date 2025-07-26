#' Compute average occupancy for each gene, using DamID binding profiles
#'
#' For each gene, averages log2 scores of all overlapping fragments (weighted by coverage).
#'
#' @param binding_data A data.frame as returned by build_dataframes(): Must have columns 'chr', 'start', 'end', and then sample columns.
#' @param ensdb_genes A GRanges object of gene annotations.  Retrieved by `get_ensdb_genes` by default, but a custom GRanges object of genes can be provided.  Must contain metadata columns 'gene_name' and 'gene_id'.
#' @param buffer Integer. Optional number of bases to extend gene boundaries both upstream and downstream.
#' @param BPPARAM A BiocParallel parameter object. Default: BiocParallel::bpparam(). On Windows, will use serial for safety.
#'
#' @return data.frame: rows = genes, columns = one per sample, plus annotation columns (gene name, nfrags, etc.)
#' @export
gene_occupancy <- function(
    binding_data,
    ensdb_genes = get_ensdb_genes()$genes,
    buffer = 0,
    BPPARAM = BiocParallel::bpparam()
    ) {
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
    if (length(gene_hits) == 0) return(rep(NA, length(sample_cols) + 4))

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

    res <- sapply(sample_cols, function(col) {
      vals <- binding_data[[col]][gene_hits]
      weighted.mean(vals, w = blens, na.rm = TRUE)
    })

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
#' For each region in 'regions', compute weighted mean of sample values from fragments overlapping that region.
#'
#' @param binding_data Data.frame as from build_dataframes().
#' @param regions GRanges of intervals (e.g., peaks).
#' @param buffer Optional integer to expand interval bounds (default 0).
#' @param BPPARAM BiocParallel parameter object.
#' @return data.frame with one row per region and columns as for binding_data samples.
#' @export
gr_occupancy <- function(
    binding_data,
    regions,
    buffer = 0,
    BPPARAM = BiocParallel::bpparam()
    ) {
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
    if (length(region_hits) == 0) return(rep(NA, length(sample_cols) + 2))

    region_range <- ranges(regions)[i]
    frag_starts <- binding_data$start[region_hits]
    frag_ends <- binding_data$end[region_hits]

    overlap_starts <- pmax(start(region_range), frag_starts)
    overlap_ends <- pmin(end(region_range), frag_ends)

    blens <- overlap_ends - overlap_starts + 1
    if (any(blens < 1)) blens[blens < 1] <- 0

    avg_scores <- sapply(sample_cols, function(col) {
      vals <- binding_data[[col]][region_hits]
      weighted.mean(vals, w = blens, na.rm = TRUE)
    })

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
#' @param peaks List of GRanges objects (e.g., from multiple samples)
#' @return GRanges: Reduced (union) regions.
#' @export
reduce_regions <- function(peaks) {
  if (!is.list(peaks) || !all(sapply(peaks, function(x) is(x, "GRanges")))) {
    stop("'peaks' must be a list of GRanges objects.")
  }
  combined <- do.call(c, peaks)
  pr <- reduce(combined)
  # Add names for intervals
  mcols(pr)$name <- paste0(seqnames(pr), ":", start(pr), "-", end(pr))
  pr
}

#' Find best gene overlap(s) for each query interval
#'
#' Annotates each input region with gene(s) it overlaps (comma-separated list).
#'
#' @param query GRanges object for which to annotate overlaps.
#' @param subject GRanges of gene annotation (must have metadata 'gene_name').
#' @param maxgap Integer. Maximum gap allowed for considered overlap (default 0=abutting).
#' @return Character vector (one per query) of gene names, or "" if none.
#' @export
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
