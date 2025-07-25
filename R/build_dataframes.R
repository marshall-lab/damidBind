#' Internal Helper: Process Binding Profile Inputs
#'
#' This internal function encapsulates the logic for loading binding profiles,
#' whether from file paths or from a list of GRanges objects. It performs
#' validation and returns a unified data.frame.
#'
#' @param binding_profiles_path Character vector. Path(s) to directories or file globs.
#' @param binding_profiles List of GRanges objects.
#' @return A dataframe containing the merged binding profiles.
#' @keywords internal
process_binding_profiles <- function(binding_profiles_path = NULL, binding_profiles = NULL) {
  # Validate that one and only one input type is provided
  if (is.null(binding_profiles_path) && is.null(binding_profiles)) {
    stop("Must supply one of binding_profiles_path or binding_profiles")
  }
  if (!is.null(binding_profiles_path) && !is.null(binding_profiles)) {
    stop("Provide only one of binding_profiles_path or binding_profiles, not both")
  }

  if (!is.null(binding_profiles_path)) {
    # Load from file paths
    message("Locating binding profile files")
    binding_files <- locate_files(binding_profiles_path, pattern = "\\.(bedgraph|bg)(\\.gz)?$")
    if (length(binding_files) == 0) {
      stop("No binding profile files (.bedgraph) found in the specified path(s).")
    }
    binding_profiles_data <- build_dataframes(binding_files)
  } else {
    # Load from a list of GRanges objects
    if (!is.list(binding_profiles) || is.null(names(binding_profiles))) {
      stop("binding_profiles must be a named list of GRanges objects")
    }
    if (!all(sapply(binding_profiles, function(x) is(x, "GRanges")))) {
      stop("All binding_profiles list elements must be GRanges objects")
    }
    message("Building binding profile dataframe from supplied GRanges objects ...")
    binding_profiles_data <- build_dataframes_from_granges(binding_profiles)
  }
  return(binding_profiles_data)
}

#' Internal Helper: Apply Quantile Normalisation
#'
#' This internal function applies quantile normalisation to a binding profile
#' data.frame if the user requests it. It is not exported.
#'
#' @param binding_profiles_data The data.frame of binding profiles.
#' @param quantile_norm Logical. If TRUE, normalisation is applied.
#' @return A data.frame, which is quantile-normalised if requested.
#' @keywords internal
apply_quantile_normalisation <- function(binding_profiles_data, quantile_norm) {
  if (isTRUE(quantile_norm)) {
    message("Applying quantile normalisation")
    coord_cols <- binding_profiles_data[, 1:3, drop = FALSE]
    signal_mat <- as.matrix(binding_profiles_data[, -(1:3), drop = FALSE])

    if (ncol(signal_mat) == 0) {
      warning("No signal columns found for quantile normalisation. Returning original data.")
      return(binding_profiles_data)
    }

    qnorm_mat <- quantile_normalisation(signal_mat)
    colnames(qnorm_mat) <- paste0(colnames(signal_mat), "_qnorm")

    binding_profiles_data <- cbind(coord_cols, as.data.frame(qnorm_mat))
  }
  return(binding_profiles_data)
}


#' Load genome-wide binding data and associated peak files or GRanges objects
#'
#' Reads DamID-seq log2 ratio binding data either from bedGraph files or
#' directly from a list of GRanges objects, and associated peak regions either
#' from GFF/bed files or from a list of GRanges objects.
#' This function is suitable for transcription factor binding analyses.
#' For peak discovery, use an external peak caller (e.g. 'find_peaks').
#'
#' One of `binding_profiles_path` or `binding_profiles` must be provided.
#' Similarly, one of `peaks_path` or `peaks` must be provided.
#'
#' When supplying GRanges lists, each GRanges should contain exactly one numeric
#' metadata column representing the binding signal, and all GRanges should be supplied
#' as a named list, with element names used as sample names.
#'
#' @param binding_profiles_path Character vector. Path(s) to directories or file
#'   globs containing log2 ratio binding tracks in bedGraph format. Wildcards ('*') supported.
#' @param peaks_path Character vector. Path(s) to directories or file globs containing
#'   the peak calls in GFF or BED format.
#' @param binding_profiles List of GRanges objects with binding profiles, one per sample.
#' @param peaks List of GRanges objects representing peak regions.
#' @param quantile_norm Logical (default: FALSE). If TRUE, quantile-normalise the signal columns across all datasets.
#' @param organism Organism string (lower case) to obtain genome annotation from (if not providing a custom `ensdb_genes` object)  Defautls to "drosophila melanogaster".
#' @param ensdb_genes GRanges object: gene annotation. Automatically obtained from `organism` if NULL.
#' @param BPPARAM BiocParallel function (defaults to BiocParallel::bpparam())
#'
#' @return A list with components:
#'   \item{binding_profiles_data}{data.frame: Signal matrix for all regions, with columns chr, start, end, sample columns.}
#'   \item{peaks}{list(GRanges): All loaded peak regions from input files or directly supplied.}
#'   \item{pr}{GRanges: Reduced (union) peak regions across samples.}
#'   \item{occupancy}{data.frame: Binding values summarised over reduced peaks, with overlap annotations.}
#'   \item{test_category}{Character scalar; will be "bound".}
#'
#' @examples
#' \dontrun{
#' # Using file paths:
#' res <- load_data_peaks(
#'   binding_profiles_path = "data/binding/",
#'   peaks_path = "data/peaks/",
#'   organism = "drosophila melanogaster")
#'
#' # Using GRanges lists:
#' res <- load_data_peaks(
#'   binding_profiles = list(
#'     cond1_n1 = gr1,
#'     cond1_n2 = gr2,
#'     cond2_n1 = gr3,
#'     cond2_n2 = gr4),
#'   peaks = list(
#'     cond1_n1 = peaks_gr1,
#'     cond1_n2 = peaks_gr2,
#'     cond2_n1 = peaks_gr3,
#'     cond2_n2 = peaks_gr4),
#'   organism = "drosophila melanogaster")
#' }
#'
#' @export
load_data_peaks <- function(
    binding_profiles_path = NULL,
    peaks_path = NULL,
    binding_profiles = NULL,
    peaks = NULL,
    quantile_norm = FALSE,
    organism = "drosophila melanogaster",
    ensdb_genes = NULL,
    BPPARAM = BiocParallel::bpparam()
) {

  if (is.null(ensdb_genes)) {
    ensdb_genes = get_ensdb_genes(organism_keyword = organism)$genes
  }
  if (!is(ensdb_genes, "GRanges")) {
    stop("ensdb_genes must be supplied as a GRanges object.")
  }

  binding_profiles_data <- process_binding_profiles(binding_profiles_path, binding_profiles)
  binding_profiles_data <- apply_quantile_normalisation(binding_profiles_data, quantile_norm)

  # Validate and load peaks
  if (is.null(peaks_path) && is.null(peaks)) {
    stop("Must supply one of peaks_path or peaks")
  }
  if (!is.null(peaks_path) && !is.null(peaks)) {
    stop("Provide only one of peaks_path or peaks, not both")
  }

  if (!is.null(peaks_path)) {
    message("Locating peak files")
    peaks_files <- locate_files(peaks_path, pattern = "\\.(gff|bed)(\\.gz)?$")
    if (length(peaks_files) == 0) {
      stop("No peak files (.gff or .bed) found in the specified path(s).")
    }
    peaks <- lapply(peaks_files, import_peaks)
  } else {
    if (!is.list(peaks) || is.null(names(peaks))) {
      stop("peaks must be a named list of GRanges objects")
    }
    if (!all(sapply(peaks, function(x) is(x, "GRanges")))) {
      stop("All peaks list elements must be GRanges objects")
    }
    message("Using supplied peaks GRanges list.")
  }

  # Process peaks and calculate occupancy
  pr <- reduce_regions(peaks)
  message("Calculating occupancy over peaks")
  occupancy <- gr_occupancy(binding_profiles_data, pr, BPPARAM = BPPARAM)

  gene_overlaps <- all_overlaps_to_original(pr, ensdb_genes, maxgap = 1000)
  occupancy$gene_names <- gene_overlaps$genes
  if (!is.null(gene_overlaps$ids)) {
    occupancy$gene_ids <- gene_overlaps$ids
  }

  list(
    binding_profiles_data = binding_profiles_data,
    peaks = peaks,
    pr = pr,
    occupancy = occupancy,
    test_category = "bound"
  )
}



#' Load genome-wide binding data for gene expression (RNA polymerase occupancy)
#'
#' Reads RNA Polymerase DamID binding profiles either from bedGraph files or
#' directly from a named list of GRanges objects. Calculates binding occupancy
#' summarised over genes.
#'
#' One of `binding_profiles_path` or `binding_profiles` must be provided.
#'
#' When supplying GRanges lists, each GRanges should contain exactly one numeric
#' metadata column representing the signal, and `binding_profiles` must be a
#' named list, with element names used as sample names.
#'
#' @param binding_profiles_path Character vector of directories or file globs
#'   containing log2 ratio binding tracks in bedGraph format. Wildcards ('*') supported.
#' @param binding_profiles Named list of GRanges objects representing binding profiles.
#' @param quantile_norm Logical (default: FALSE) quantile-normalise across all signal columns if TRUE.
#' @param organism Organism string (lower case) to obtain genome annotation from (if not providing a custom `ensdb_genes` object)
#'   Defautls to "drosophila melanogaster".
#' @param ensdb_genes GRanges object: gene annotation. Automatically obtained from `organism` if NULL.
#' @param BPPARAM BiocParallel function (defaults to BiocParallel::bpparam())
#'
#' @return List with elements:
#'   \item{binding_profiles_data}{data.frame of merged binding profiles, with chr, start, end, sample columns.}
#'   \item{occupancy}{data.frame of occupancy values summarised over genes.}
#'   \item{test_category}{Character scalar; will be "expressed".}
#'
#' @examples
#' \dontrun{
#' # Using file paths:
#' res <- load_data_genes(
#'   binding_profiles_path = "data/rnapol/",
#'   organism = "drosophila melanogaster")
#'
#' # Using GRanges list:
#' res <- load_data_genes(
#'   binding_profiles = list(
#'     cond1_n1 = gr1,
#'     cond1_n2 = gr2,
#'     cond2_n1 = gr3,
#'     cond2_n2 = gr4),
#'   organism = "drosophila melanogaster")
#' }
#'
#' @export
load_data_genes <- function(
    binding_profiles_path = NULL,
    binding_profiles = NULL,
    quantile_norm = FALSE,
    organism = "drosophila melanogaster",
    ensdb_genes = NULL,
    BPPARAM = BiocParallel::bpparam()
) {
  if (is.null(ensdb_genes)) {
    ensdb_genes = get_ensdb_genes(organism_keyword = organism)$genes
  }
  if (!is(ensdb_genes, "GRanges")) {
    stop("ensdb_genes must be supplied as a GRanges object.")
  }

  binding_profiles_data <- process_binding_profiles(binding_profiles_path, binding_profiles)
  binding_profiles_data <- apply_quantile_normalisation(binding_profiles_data, quantile_norm)

  # Calculate occupancy over genes
  occupancy <- gene_occupancy(binding_profiles_data, ensdb_genes, BPPARAM = BPPARAM)

  list(
    binding_profiles_data = binding_profiles_data,
    occupancy = occupancy,
    test_category = "expressed"
  )
}


#' Build data frames of binding profiles from bedGraph files
#'
#' Each file must be in bedGraph format: chr, start, end, value.
#' @param bedgraphs Character vector of file paths.
#' @return data.frame with merged intervals and all sample columns.
#' @keywords internal
build_dataframes <- function(bedgraphs) {
  message("Building binding profile dataframe from input files ...")
  if (length(bedgraphs) < 1) stop("No bedGraph files supplied.")

  # Helper function to deal with multiple periods in files
  strip_all_exts_recursive <- function(filepath) {
    lastpath <- ""
    while (filepath != lastpath) {
      lastpath <- filepath
      filepath <- file_path_sans_ext(filepath)
    }
    filepath
  }

  data.df <- NULL
  for (bf in bedgraphs) {
    # Name sample by filename (remove extension)
    sample_name <- basename(strip_all_exts_recursive(bf))
    gr <- import_bedgraph_as_df(bf, colname = sample_name)
    if (is.null(data.df)) {
      data.df <- gr
    } else {
      data.df <- merge(data.df, gr, by = c("chr", "start", "end"), all = F)
    }
    message(" - Loaded: ", sample_name)
  }
  # Order by chromosome and location
  data.df <- data.df[order(data.df$chr, data.df$start), ]
  data.df
}


#' Build a merged binding profile dataframe from a named list of GRanges objects
#'
#' Each GRanges should have exactly one numeric metadata column representing the binding signal.
#' The list must be named, with element names used as sample names.
#'
#' @param gr_list Named list of GRanges objects with binding signal in numeric metadata column.
#' @return data.frame with merged intervals and columns: chr, start, end, sample columns.
#' @keywords internal
build_dataframes_from_granges <- function(gr_list) {
  if (length(gr_list) == 0) stop("Empty GRanges list supplied to build_dataframes_from_granges.")

  df_list <- lapply(seq_along(gr_list), function(i) {
    gr <- gr_list[[i]]
    sample_name <- names(gr_list)[i]

    if (!is(gr, "GRanges")) {
      stop("Element ", i, " of gr_list is not a GRanges object.")
    }

    # Detect numeric metadata columns
    mcols_gr <- mcols(gr)
    numeric_cols <- names(mcols_gr)[sapply(mcols_gr, is.numeric)]

    if (length(numeric_cols) == 0) {
      stop("GRanges object '", sample_name, "' has no numeric metadata column for binding signal.")
    } else if (length(numeric_cols) > 1) {
      stop("GRanges object '", sample_name, "' has multiple numeric metadata columns; please provide exactly one.")
    }
    value_col <- numeric_cols[1]

    df <- data.frame(
      chr = as.character(seqnames(gr)),
      start = start(gr),
      end = end(gr),
      value = mcols_gr[[value_col]],
      stringsAsFactors = FALSE
    )
    colnames(df)[4] <- sample_name
    return(df)
  })

  # Merge profiles
  merged_df <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = FALSE), df_list)

  # Sort and return
  merged_df <- merged_df[order(merged_df$chr, merged_df$start), , drop = FALSE]
  return(merged_df)
}


#' Locate files with optional wildcard expansion
#' @param paths Character vector of directories or patterns (e.g., "/mypath/sample*").
#' @param pattern Optional regex for file extensions (e.g., "\\.gff$").  Case is ignored.
#' @return Character vector of file paths.
#' @keywords internal
locate_files <- function(paths, pattern = NULL) {
  files <- character()
  for (p in paths) {
    # Expand wildcards
    expanded <- Sys.glob(p)
    # List files in directory, if p is a bare directory
    if (dir.exists(p)) {
      dir_files <- list.files(p, pattern = pattern, full.names = TRUE)
      files <- c(files, dir_files)
    }
    # Check if expanded globs were dirs:
    for (checkexpdir in expanded) {
      if (dir.exists(checkexpdir)) {
        dir_files <- list.files(checkexpdir, pattern, full.names = TRUE)
        files <- c(files, dir_files)
      }
    }
    # Include expanded gene_names
    files <- c(files, expanded[grepl(pattern, expanded, ignore.case = TRUE)])
  }
  files <- unique(files)
  return(files[file.exists(files)])
}

#' Import a GFF file as a GRanges object
#' @param path File path (GFF/GTF)
#' @return GRanges
#' @keywords internal
import_peaks <- function(path) {
  tryCatch({
    import(path)
  }, error = function(e) {
    stop("Failed to read peaks file: ", path, "\\n", e$message)
  })
}

#' Import a bedGraph file as a data.frame (chr, start, end, value)
#' @param path File path (bedGraph)
#' @param colname Name of the value column (usually sample name).
#' @return data.frame
#' @keywords internal
import_bedgraph_as_df <- function(path, colname = "score") {
  gr <- tryCatch({
    import(path, format = "bedGraph")
  }, error = function(e) {
    stop("Failed to read bedGraph file: ", path, "\\n", e$message)
  })
  df <- as.data.frame(gr)[, c("seqnames", "start", "end", "score")]
  names(df) <- c("chr", "start", "end", colname)
  df
}
