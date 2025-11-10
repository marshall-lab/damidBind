#' Compute occupancy for genomic regions
#'
#' For each interval in the `regions` GRanges object, this function finds all
#' overlapping fragments in `binding_data` and computes a weighted mean of their
#' signal values. Any metadata columns present in the input `regions`
#' object are preserved in the output data.frame.
#'
#' @param binding_data A data.frame as produced by `build_dataframes()`. It must
#'   contain columns 'chr', 'start', 'end', followed by numeric sample columns.
#' @param regions A GRanges object of genomic intervals (e.g., genes or reduced peaks)
#'   over which to calculate occupancy.
#' @param buffer Optional integer. Number of base pairs to expand each interval
#'   in `regions` on both sides before calculating occupancy. Default is 0.
#' @param BPPARAM A BiocParallel parameter object for parallel computation.
#'   Default is `BiocParallel::bpparam()`.
#'
#' @return A data.frame with one row per region from the input `regions` object.
#'   The output includes the weighted mean occupancy for each sample, `nfrags`
#'   (number of overlapping fragments), and all original metadata columns from `regions`.
#'   Rownames are generated from region coordinates to ensure uniqueness.
#' @export
#' @examples
#' # Create a set of regions with metadata
#' regions_gr <- GenomicRanges::GRanges(
#'     "chrX", IRanges::IRanges(start = c(100, 500), width = 100),
#'     gene_name = c("MyGene1", "MyGene2"), score = c(10, 20)
#' )
#'
#' # Create a mock binding data GRanges object
#' binding_gr <- GenomicRanges::GRanges(
#'     seqnames = "chrX",
#'     ranges = IRanges::IRanges(
#'         start = c(90, 150, 480, 550),
#'         end = c(110, 170, 520, 580)
#'     ),
#'     sampleA = c(1.2, 0.8, 2.5, 3.0),
#'     sampleB = c(1.0, 0.9, 2.8, 2.9)
#' )
#'
#' # Calculate occupancy over the regions
#' # Use BiocParallel::SerialParam() for deterministic execution in examples
#' if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'     occupancy_data <- calculate_occupancy(binding_gr, regions_gr,
#'         BPPARAM = BiocParallel::SerialParam()
#'     )
#'     print(occupancy_data)
#' }
calculate_occupancy <- function(
        binding_data,
        regions,
        buffer = 0,
        BPPARAM = BiocParallel::bpparam()) {

    # Input Validation
    if (!is(binding_data, "GRanges")) {
        stop("'binding_data' must be a GRanges object.")
    }
    if (!is(regions, "GRanges")) {
        stop("'regions' must be a GRanges object.")
    }
    message(sprintf("Calculating average occupancy for %d regions...", length(regions)))

    # Apply buffer if specified
    if (buffer != 0) {
        regions <- suppressWarnings(
            trim(resize(regions, width = width(regions) + buffer * 2, fix = "center"))
        )
    }

    # Resolve potential duplicate coordinates to create unique keys
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
    region_keys <- resolve_dups(regions)
    mcols(regions)$`.internal_name` <- region_keys

    # Get sample column names
    sample_cols <- colnames(mcols(binding_data))

    overlaps <- findOverlaps(regions, binding_data)

    # Worker function to calculate occupancy over one region
    calc_one_region <- function(i) {
        # Find which fragments overlap the current region 'i'
        region_hits <- subjectHits(overlaps)[queryHits(overlaps) == i]

        # If no fragments overlap, return NA
        if (length(region_hits) == 0) {
            return(rep(NA, length(sample_cols) + 1))
        }

        # Get the genomic range for the current region
        region_range <- ranges(regions)[i]

        # Get the coordinates of all overlapping fragments using robust accessors
        overlapping_fragments <- binding_data[region_hits]
        frag_starts <- start(overlapping_fragments)
        frag_ends <- end(overlapping_fragments)

        # Calculate the width of the precise overlap for each fragment
        overlap_starts <- pmax(start(region_range), frag_starts)
        overlap_ends <- pmin(end(region_range), frag_ends)
        overlap_widths <- overlap_ends - overlap_starts + 1
        overlap_widths[overlap_widths < 1] <- 0

        # Calculate the weighted mean for each sample column
        avg_scores <- vapply(sample_cols, function(col) {
            # Use the correct 'binding_data' object
            vals <- mcols(overlapping_fragments)[[col]]
            if (sum(overlap_widths, na.rm = TRUE) == 0) return(NA_real_)
            weighted.mean(vals, w = overlap_widths, na.rm = TRUE)
        }, FUN.VALUE = numeric(1))

        # Return the number of fragments and the calculated scores
        return(c(nfrags = length(region_hits), avg_scores))
    }

    # Run the calculation in parallel for all regions
    results_matrix <- do.call(rbind, bplapply(seq_along(regions), calc_one_region, BPPARAM = BPPARAM))

    # Convert the matrix of results to a data.frame
    results_df <- as.data.frame(results_matrix)
    colnames(results_df) <- c("nfrags", sample_cols)
    cols_to_convert <- colnames(results_df) # All are numeric
    results_df[cols_to_convert] <- lapply(results_df[cols_to_convert], as.numeric)

    # Get all original metadata from the input regions
    original_metadata_df <- as.data.frame(mcols(regions))

    # Combine the calculated results with the original metadata
    final_df <- cbind(original_metadata_df, results_df)

    final_df <- na.omit(final_df)
    rownames(final_df) <- final_df$`.internal_name`

    # Rename the internal name column to 'name' for consistency
    if(".internal_name" %in% colnames(final_df)){
        colnames(final_df)[colnames(final_df) == ".internal_name"] <- "name"
    }

    return(final_df)
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
#' gr_list <- list(gr1, gr2)
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
#' @noRd
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
