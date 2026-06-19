#' @title Get metadata for a set of files
#' @description Internal helper that captures the absolute path, modification
#'   timestamp, and size for a vector of file paths.
#' @param files A character vector of file paths.
#' @return A data.frame with metadata for each file. If input is NULL or empty,
#'   returns a data.frame indicating user-supplied data.
#' @noRd
._get_file_metadata <- function(files) {
    if (is.null(files) || length(files) == 0) {
        return(data.frame(
            file = "User-supplied GRanges",
            path = "N/A",
            modified = NA,
            size_bytes = NA,
            stringsAsFactors = FALSE
        ))
    }

    # Resolve absolute paths and capture file system information
    abs_paths <- normalizePath(files, mustWork = FALSE, winslash = "/")
    info <- file.info(files)

    data.frame(
        file = basename(files),
        path = abs_paths,
        modified = info$mtime,
        size_bytes = info$size,
        stringsAsFactors = FALSE
    )
}

#' @title Assemble analysis metadata
#' @description Internal function that aggregates package, annotation, and
#'   file metadata into a single nested list for the DamIDResults object.
#' @param test_category Character. The type of analysis (bound, expressed, etc).
#' @param annotation_info List or NULL. The output from get_ensdb_genes().
#' @param binding_files Character vector. Paths to binding profile files.
#' @param peak_files Character vector or NULL. Paths to peak files.
#' @return A nested list of metadata.
#' @noRd
._get_analysis_metadata <- function(test_category, annotation_info,
                                    binding_files, peak_files = NULL,
                                    samples_kept = NULL,
                                    drop_samples_input = NULL) {
    # Package version and analysis timestamp
    meta <- list(
        package_version = as.character(utils::packageVersion("damidBind")),
        analysis_date = Sys.time(),
        test_category = test_category
    )

    # Record annotation metadata
    if (is.null(annotation_info)) {
        # User provided a custom GRanges object
        meta$annotation <- list(source = "Manual (GRanges)")
    } else {
        meta$annotation <- list(
            source = "AnnotationHub",
            ah_id = annotation_info$ah_id,
            ah_snapshot = annotation_info$ah_snapshot,
            ensembl_version = annotation_info$ensembl_version,
            genome_build = annotation_info$genome_build,
            species = annotation_info$species
        )
    }

    # Record datafiles
    meta$files <- list(
        binding_profiles = ._get_file_metadata(binding_files),
        peaks = if (!is.null(peak_files)) ._get_file_metadata(peak_files) else NULL
    )

    # Record filter provenance
    meta$provenance <- list(
        samples_kept = samples_kept,
        drop_samples_applied = !is.null(drop_samples_input),
        drop_samples_filter = drop_samples_input
    )

    return(meta)
}


#' @title Process binding profile inputs
#'
#' @description
#' This internal function encapsulates the logic for loading binding profiles,
#' whether from file paths, a list of GRanges objects, or a single merged
#' GRanges object. It performs validation and returns a unified list
#' containing the data and file manifest.
#'
#' @param binding_profiles_path Character vector. Path(s) to directories or file globs.
#' @param binding_profiles List of GRanges objects or a single GRanges object.
#' @return A list containing the merged binding profiles (GRanges) and the file paths (character).
#' @noRd
process_binding_profiles <- function(binding_profiles_path = NULL, binding_profiles = NULL) {
    # Validate that one and only one input type is provided
    if (is.null(binding_profiles_path) && is.null(binding_profiles)) {
        stop("Must supply one of binding_profiles_path or binding_profiles")
    }
    if (!is.null(binding_profiles_path) && !is.null(binding_profiles)) {
        stop("Provide only one of binding_profiles_path or binding_profiles, not both")
    }

    binding_files <- NULL

    if (!is.null(binding_profiles_path)) {
        # Load from file paths
        message("Locating binding profile files")
        binding_files <- locate_files(binding_profiles_path, pattern = "\\.(bedgraph|bg)(\\.gz)?$")
        if (length(binding_files) == 0) {
            stop("No binding profile files (.bedgraph) found in the specified path(s).")
        }
        binding_profiles_data <- build_dataframes(binding_files)

        # Convert back to GRanges post-merge
        binding_profiles_gr <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(binding_profiles_data$chr),
            ranges = IRanges::IRanges(start = binding_profiles_data$start, end = binding_profiles_data$end)
        )
        # Add the sample data to the metadata columns
        mcols(binding_profiles_gr) <- binding_profiles_data[, !(names(binding_profiles_data) %in% c("chr", "start", "end")), drop = FALSE]

    } else if (is(binding_profiles, "GRanges")) {
        # Input is a single GRanges object -- assume metadata columns are samples
        message("Using supplied binding profile GRanges object ...")

        # Basic validation: ensure it has numeric metadata columns
        mcols_gr <- mcols(binding_profiles)
        numeric_cols <- vapply(mcols_gr, is.numeric, FUN.VALUE = logical(1))

        if (!any(numeric_cols)) {
            stop("The supplied binding_profiles GRanges must have at least one numeric metadata column representing signal.")
        }

        binding_profiles_gr <- binding_profiles

    } else if (is.list(binding_profiles)) {
        # Input is a list of GRanges objects (one per sample)
        if (is.null(names(binding_profiles))) {
            stop("binding_profiles list must be named (names are used as sample IDs)")
        }
        if (!all(vapply(binding_profiles, function(x) is(x, "GRanges"), FUN.VALUE = logical(1)))) {
            stop("All binding_profiles list elements must be GRanges objects")
        }
        message("Building binding profile dataframe from supplied GRanges list ...")
        binding_profiles_data <- build_dataframes_from_granges(binding_profiles)

        # Convert to GRanges
        binding_profiles_gr <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(binding_profiles_data$chr),
            ranges = IRanges::IRanges(start = binding_profiles_data$start, end = binding_profiles_data$end)
        )
        mcols(binding_profiles_gr) <- binding_profiles_data[, !(names(binding_profiles_data) %in% c("chr", "start", "end")), drop = FALSE]

    } else {
        stop("binding_profiles must be a list of GRanges objects or a single GRanges object")
    }

    return(list(
        data = binding_profiles_gr,
        files = binding_files
    ))
}

#' Apply genome-wide data normalisation
#'
#' This internal function applies one of several potential normalisation methods
#' to a genome-wide binding profile dataset prior to peak subsetting/occupancy calculation.
#'
#' @param binding_profiles_data The GRanges object of binding profiles.
#' @param norm_method Character. The method to apply ("none", "scale", "quantile", "rpm").
#' @param pre_scale Logical.  Whether to apply scaling prior to normalisation.
#' @return A GRanges object, normalised if requested.
#' @noRd
apply_normalisation <- function(binding_profiles_data, norm_method = c("none", "quantile", "loess", "rpm"), pre_scale = TRUE) {
    norm_method <- match.arg(norm_method)

    if (norm_method == "none") return(binding_profiles_data)

    message(sprintf("Applying '%s' normalisation to binding profiles", norm_method))

    # Extract matrix
    signal_mat <- as.matrix(mcols(binding_profiles_data))

    # Apply basic scaling
    if (isTRUE(pre_scale)) {
        signal_mat <- scale(signal_mat, center=FALSE, scale=TRUE)
    }

    if (ncol(signal_mat) == 0) {
        warning("No signal columns found for normalisation. Returning original data.")
        return(binding_profiles_data)
    }

    norm_mat <- switch(norm_method,
                       "quantile" = {
                           quantile_normalisation(signal_mat)
                       },
                       "loess" = {
                           # Cyclic Loess
                           limma::normalizeBetweenArrays(signal_mat, method = "cyclicloess")
                       },
                       "rpm" = {
                           apply(signal_mat, 2, function(x) {
                               total <- sum(abs(x), na.rm = TRUE)
                               if (!is.na(total) && total > 0) (x / total) * 1e6 else x
                           })
                       }
    )

    # Final safety check for NaNs produced by non-finite inputs
    if (any(is.nan(norm_mat))) {
        warning("Normalisation produced NaNs. Check for non-finite values in input data. Returning un-normalised data.")
        return(binding_profiles_data)
    }

    # Append suffix to indicate normalisation method used
    suffix <- switch(norm_method,
                     "quantile" = "_qnorm",
                     "loess" = "_loess",
                     "rpm" = "_rpm"
    )

    colnames(norm_mat) <- paste0(colnames(signal_mat), suffix)
    mcols(binding_profiles_data) <- as.data.frame(norm_mat)

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
#' When supplying GRanges lists for `binding_profiles`, each GRanges should
#' contain exactly one numeric metadata column representing the binding signal,
#' and all GRanges should be supplied as a named list, with element names used
#' as sample names. Alternatively, a single GRanges object already containing
#' multiple sample columns may be provided.
#'
#' @param binding_profiles_path Character vector. Path(s) to directories or file
#'   globs containing log2 ratio binding tracks in bedGraph format. Wildcards ('*') supported.
#' @param peaks_path Character vector. Path(s) to directories or file globs containing
#'   the peak calls in GFF or BED format.
#' @param binding_profiles A list of GRanges objects (one per sample) or a
#'   single GRanges object containing multiple sample metadata columns.
#' @param peaks A list of GRanges objects (one per sample) or a single GRanges
#'   object representing the union of peaks.
#' @param drop_samples A character vector of sample names or patterns to remove.
#'   Matching samples are removed from the analysis before normalisation and
#'   occupancy calculation. This can be useful for excluding samples that fail
#'   initial quality checks. Default: `NULL` (no samples are dropped).
#' @param maxgap_loci Integer, the maximum bp distance between a peak boundary
#'   and a gene to associate that peak with the gene. Default: 1000.
#' @param norm_method Character. Determines how raw genome-wide signal is normalised
#'   prior to establishing peak specific occupancy. Options are:
#'   \itemize{
#'     \item \code{"none"} (default): Use original data. Recommended if loading CATaDa
#'       count data pre-processed via \code{damidseq_pipeline --catada} (which is natively
#'       RPM-normalised).
#'     \item \code{"loess"}: Uses the `cyclicloess` method from `limma`.
#'     \item \code{"quantile"}: Quantile normalisation. Forces identical target distributions.
#'     \item \code{"rpm"}: Reads Per Million scaling. A simple global scaling for raw
#'       un-normalised count data.  Use with caution: do not use for log2 ratio data, and avoid
#'       if CATaDa data has already been RPM normalised (the default behaviour
#'       with `damidseq_pipeline`)
#'   }
#' @param pre_scale Logical.  If TRUE, samples will be scaled without centering
#'   via base R `scale(center=FALSE, scale=TRUE)`, prior to normalisation being applied.
#'   default is FALSE; generally useful for DamID log-ratio data, but not recommended
#'   for CATaDa data.
#' @param quantile_norm \bold{Deprecated}. Please use `norm_method = "quantile"` instead.
#' @param organism Organism string (lower case) to obtain genome annotation from
#'   (if not providing a custom `ensdb_genes` object)
#'   Default: "drosophila melanogaster".
#' @param ensdb_genes GRanges object: gene annotation. Automatically obtained
#'   from `organism` if NULL.
#' @param BPPARAM BiocParallel function (defaults to BiocParallel::bpparam())
#' @param plot_diagnostics Logical. If `TRUE` (the default in interactive sessions),
#'   diagnostic plots (PCA and correlation heatmap) will be generated and
#'   displayed for both the raw binding data and the summarised occupancy data.
#'
#' @return A list with components:
#'   \item{binding_profiles_data}{GRanges: Signal matrix for all regions, with genomic coordinates and sample metadata columns.}
#'   \item{peaks}{list(GRanges): All loaded peak regions.}
#'   \item{pr}{GRanges: Reduced (union) peak regions across samples.}
#'   \item{occupancy}{data.frame: Binding values summarised over reduced peaks, with overlap annotations.}
#'   \item{test_category}{Character scalar; will be "bound".}
#'   \item{metadata}{List: Data provenance and analysis metadata.}
#'
#'
#' @examples
#' # Create a mock GRanges object for gene annotation
#' # This object, based on the package's unit tests, avoids network access
#' # and includes a very long gene to ensure overlaps with sample data.
#' mock_genes_gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle("2L", 7),
#'     ranges = IRanges::IRanges(
#'         start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'         end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'     ),
#'     strand = S4Vectors::Rle(GenomicRanges::strand(c("+", "-", "+", "+", "-", "-", "+"))),
#'     gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'     gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene")
#' )
#'
#' # Get path to sample data files included with the package
#' data_dir <- system.file("extdata", package = "damidBind")
#'
#' # Run loading function using sample files and mock gene annotations
#' # Normalising with either 'quantile' or 'loess' is generally recommended
#' # for DamID log2 ratios
#' loaded_data <- load_data_peaks(
#'     binding_profiles_path = data_dir,
#'     peaks_path = data_dir,
#'     ensdb_genes = mock_genes_gr,
#'     norm_method = "quantile"
#' )
#'
#' # View the structure of the output
#' str(loaded_data, max.level = 1)
#'
#' @export
load_data_peaks <- function(
        binding_profiles_path = NULL,
        peaks_path = NULL,
        binding_profiles = NULL,
        peaks = NULL,
        drop_samples = NULL,
        maxgap_loci = 1000,
        norm_method = c("none", "loess", "quantile", "rpm"),
        pre_scale = FALSE,
        quantile_norm = NULL,
        organism = "drosophila melanogaster",
        ensdb_genes = NULL,
        BPPARAM = BiocParallel::bpparam(),
        plot_diagnostics = interactive()) {

    # Handle deprecated quantile_norm argument carefully
    if (!is.null(quantile_norm)) {
        warning("The 'quantile_norm' argument is deprecated. Please use 'norm_method' instead.", call. = FALSE)
        if (isTRUE(quantile_norm)) {
            norm_method <- "quantile"
        }
    }

    norm_method <- match.arg(norm_method)

    # Get annotation and metadata
    anno_info <- NULL
    if (is.null(ensdb_genes)) {
        anno_info <- get_ensdb_genes(organism_keyword = organism)
        ensdb_genes <- anno_info$genes
    }
    if (!is(ensdb_genes, "GRanges")) {
        stop("ensdb_genes must be supplied as a GRanges object.")
    }

    # Load binding profiles
    binding_info <- process_binding_profiles(binding_profiles_path, binding_profiles)
    binding_profiles_data <- binding_info$data
    binding_files <- binding_info$files

    # Validate and load peaks
    if (is.null(peaks_path) && is.null(peaks)) {
        stop("Must supply one of peaks_path or peaks")
    }
    if (!is.null(peaks_path) && !is.null(peaks)) {
        stop("Provide only one of peaks_path or peaks, not both")
    }

    peaks_files <- NULL
    if (!is.null(peaks_path)) {
        message("Locating peak files")
        peaks_files <- locate_files(peaks_path, pattern = "\\.(gff|bed)(\\.gz)?$")
        if (length(peaks_files) == 0) {
            stop("No peak files (.gff or .bed) found in the specified path(s).")
        }
        # Helper to strip all extensions for consistent naming with binding profiles
        strip_all_exts_recursive <- function(filepath) {
            lastpath <- ""
            while (filepath != lastpath) {
                lastpath <- filepath
                filepath <- file_path_sans_ext(filepath)
            }
            filepath
        }
        # Name the files so the resulting list is named.
        names(peaks_files) <- vapply(peaks_files, function(p) basename(strip_all_exts_recursive(p)), character(1))
        peaks <- lapply(peaks_files, import_peaks)
    } else if (is(peaks, "GRanges")) {
        # Wrap single GRanges in a list for reduce_regions compatibility
        message("Using supplied peak GRanges object.")
        peaks <- list("User_supplied" = peaks)
    } else {
        if (!is.list(peaks) || is.null(names(peaks))) {
            stop("peaks must be a named list of GRanges objects or a single GRanges object")
        }
        if (!all(vapply(peaks, function(x) is(x, "GRanges"), FUN.VALUE = logical(1)))) {
            stop("All peaks list elements must be GRanges objects")
        }
        message("Using supplied peaks GRanges list.")
    }

    # Drop samples if requested by the user
    if (!is.null(drop_samples)) {
        temp_data_list <- list(
            binding_profiles_data = binding_profiles_data,
            peaks = peaks
        )
        filtered_data <- ._drop_input_samples(temp_data_list, drop_samples)
        binding_profiles_data <- filtered_data$binding_profiles_data
        peaks <- filtered_data$peaks
    }

    # Apply normalisation
    binding_profiles_data <- apply_normalisation(binding_profiles_data, norm_method, pre_scale)

    # Process peaks and calculate occupancy
    pr <- reduce_regions(peaks)
    message("Calculating occupancy over peaks")
    occupancy <- calculate_occupancy(binding_profiles_data, pr, BPPARAM = BPPARAM)

    gene_overlaps <- all_overlaps_to_original(pr, ensdb_genes, maxgap = maxgap_loci)
    occupancy$gene_name <- gene_overlaps$genes
    if (!is.null(gene_overlaps$ids)) {
        occupancy$gene_id <- gene_overlaps$ids
    }

    # Record metadata
    metadata <- ._get_analysis_metadata(
        test_category = "bound",
        annotation_info = anno_info,
        binding_files = binding_files,
        peak_files = peaks_files,
        samples_kept = colnames(mcols(binding_profiles_data)),
        drop_samples_input = drop_samples
    )

    result_list <- list(
        binding_profiles_data = binding_profiles_data,
        peaks = peaks,
        pr = pr,
        occupancy = occupancy,
        test_category = "bound",
        metadata = metadata
    )

    if (isTRUE(plot_diagnostics)) {
        plot_input_diagnostics(result_list)
    }
    return(result_list)
}


#' Load genome-wide binding data for gene expression (RNA polymerase occupancy)
#'
#' Reads RNA Polymerase DamID binding profiles either from bedGraph files or
#' directly from a list of GRanges objects. Calculates binding occupancy
#' summarised over genes.
#'
#' One of `binding_profiles_path` or `binding_profiles` must be provided.
#'
#' When supplying GRanges lists for `binding_profiles`, each GRanges should
#' contain exactly one numeric metadata column representing the binding signal,
#' and `binding_profiles` must be a named list, with element names used as
#' sample names. Alternatively, a single GRanges object already containing
#' multiple sample columns may be provided.
#'
#' @param binding_profiles_path Character vector of directories or file globs
#'   containing log2 ratio binding tracks in bedGraph format. Wildcards ('*') supported.
#' @param binding_profiles A list of GRanges objects (one per sample) or a
#'   single GRanges object containing multiple sample metadata columns.
#' @param drop_samples A character vector of sample names or patterns to remove.
#'   Matching samples are removed from the analysis before normalisation and
#'   occupancy calculation. This can be useful for excluding samples that fail
#'   initial quality checks. Default: `NULL` (no samples are dropped).
#' @param norm_method Character. Determines how raw genome-wide signal is normalised
#'   prior to establishing peak specific occupancy. Options are:
#'   \itemize{
#'     \item \code{"none"} (default): Use original data. Recommended if loading CATaDa
#'       count data pre-processed via \code{damidseq_pipeline --catada} (which is natively
#'       RPM-normalised).
#'     \item \code{"loess"}: Uses the `cyclicloess` method from `limma`.
#'     \item \code{"quantile"}: Quantile normalisation. Forces identical target distributions.
#'     \item \code{"rpm"}: Reads Per Million scaling. A simple global scaling for raw
#'       un-normalised count data.  Use with caution: do not use for log2 ratio data, and avoid
#'       if CATaDa data has already been RPM normalised (the default behaviour
#'       with `damidseq_pipeline`)
#'   }
#' @param pre_scale Logical.  If TRUE, samples will be scaled without centering
#'   via base R `scale(center=FALSE, scale=TRUE)`, prior to normalisation being applied.
#'   Default is FALSE.
#' @param quantile_norm \bold{Deprecated}. Please use `norm_method = "quantile"` instead.
#' @param organism Organism string (lower case) to obtain genome annotation from (if not providing a custom `ensdb_genes` object)
#'   Defautls to "drosophila melanogaster".
#' @param calculate_occupancy_pvals Calculate occupancy p-values as a proxy for gene expression status (see details).  Not used for differential expression analysis, but used when present for downstream analysis and plotting. (default: TRUE)
#' @param return_per_replicate_fdr Legacy option of returning BH-adjusted RNA Polymerase occupancy FDR values
#'   per replicate.  As of v0.99.12, unadjusted p-values are returned by defualt; these are
#'   then aggregated at the condition level during `differential_binding()` and the
#'   aggregate p-values adjusted to gain statistical power.  This option exists for legacy or unsual
#'   end-user applications.  Use with caution.  (default: FALSE)
#' @param null_model_iterations Number of iterations to use to determine null model for FDR (default: 100000)
#' @param ensdb_genes GRanges object: gene annotation. Automatically obtained
#'   from `organism` if NULL.
#' @param BPPARAM BiocParallel function (defaults to BiocParallel::bpparam())
#' @param plot_diagnostics Logical. If `TRUE` (the default in interactive sessions),
#'   diagnostic plots (PCA and correlation heatmap) will be generated and
#'   displayed for both the raw binding data and the summarised occupancy data.
#' @param occupancy_plot_diagnostics Logical. If `TRUE` (default in interactive sessions),
#'   diagnostic plots for the gene expression null model will be displayed.
#'
#' @return List with elements:
#'   \item{binding_profiles_data}{GRanges: Merged binding profiles, with genomic coordinates and sample metadata columns.}
#'   \item{occupancy}{data.frame: Occupancy values summarised over genes.}
#'   \item{test_category}{Character scalar; will be "expressed".}
#'   \item{metadata}{List: Data provenance and analysis metadata.}
#'
#' @examples
#' # Create a mock GRanges object for gene annotations
#' # This object avoids network access
#' # and includes a very long gene to ensure overlaps with sample data.
#' mock_genes_gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle("2L", 7),
#'     ranges = IRanges::IRanges(
#'         start = c(1000, 2000, 3000, 5000, 6000, 7000, 8000),
#'         end = c(1500, 2500, 3500, 5500, 6500, 7500, 20000000)
#'     ),
#'     strand = S4Vectors::Rle(GenomicRanges::strand(c("+", "-", "+", "+", "-", "-", "+"))),
#'     gene_id = c("FBgn001", "FBgn002", "FBgn003", "FBgn004", "FBgn005", "FBgn006", "FBgn007"),
#'     gene_name = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "LargeTestGene")
#' )
#'
#' # Get path to sample data files included with the package
#' data_dir <- system.file("extdata", package = "damidBind")
#'
#' # Run loading function using sample files and mock gene annotations
#' # This calculates occupancy over genes instead of peaks.
#' loaded_data_genes <- load_data_genes(
#'     binding_profiles_path = data_dir,
#'     ensdb_genes = mock_genes_gr,
#'     norm_method = "none",
#'     calculate_occupancy_pvals = FALSE
#' )
#'
#' # View the head of the occupancy table
#' head(loaded_data_genes$occupancy)
#'
#' @export
load_data_genes <- function(
        binding_profiles_path = NULL,
        binding_profiles = NULL,
        drop_samples = NULL,
        norm_method = c("none", "loess", "quantile", "rpm"),
        pre_scale = FALSE,
        quantile_norm = NULL,
        organism = "drosophila melanogaster",
        calculate_occupancy_pvals = TRUE,
        return_per_replicate_fdr = FALSE,
        occupancy_plot_diagnostics = interactive(),
        null_model_iterations = 100000,
        ensdb_genes = NULL,
        BPPARAM = BiocParallel::bpparam(),
        plot_diagnostics = interactive()) {

    # Handle deprecated quantile_norm argument
    if (!is.null(quantile_norm)) {
        warning("The 'quantile_norm' argument is deprecated. Please use 'norm_method' instead.", call. = FALSE)
        if (isTRUE(quantile_norm)) {
            norm_method <- "quantile"
        }
    }

    norm_method <- match.arg(norm_method)

    # Resolve annotation
    anno_info <- NULL
    if (is.null(ensdb_genes)) {
        anno_info <- get_ensdb_genes(organism_keyword = organism)
        ensdb_genes <- anno_info$genes
    }
    if (!is(ensdb_genes, "GRanges")) {
        stop("ensdb_genes must be supplied as a GRanges object.")
    }

    # Load binding profiles
    binding_info <- process_binding_profiles(binding_profiles_path, binding_profiles)
    binding_profiles_data <- binding_info$data
    binding_files <- binding_info$files

    # Drop samples if requested by the user.
    if (!is.null(drop_samples)) {
        temp_data_list <- list(binding_profiles_data = binding_profiles_data)
        filtered_data <- ._drop_input_samples(temp_data_list, drop_samples)
        binding_profiles_data <- filtered_data$binding_profiles_data
    }

    # Apply normalisation
    binding_profiles_data <- apply_normalisation(binding_profiles_data, norm_method, pre_scale)

    # Calculate occupancy over genes
    occupancy <- calculate_occupancy(binding_profiles_data, ensdb_genes, BPPARAM = BPPARAM)

    # Optionally, calculate and add FDR columns
    if (isTRUE(calculate_occupancy_pvals)) {
        occupancy <- calculate_and_add_occupancy_pvals(
            binding_data = binding_profiles_data,
            occupancy_df = occupancy,
            null_model_iterations = null_model_iterations,
            return_per_replicate_fdr = return_per_replicate_fdr,
            plot_diagnostics = occupancy_plot_diagnostics,
            BPPARAM = BPPARAM
        )
    }

    # Record metadata
    metadata <- ._get_analysis_metadata(
        test_category = "expressed",
        annotation_info = anno_info,
        binding_files = binding_files,
        peak_files = NULL,
        samples_kept = colnames(mcols(binding_profiles_data)),
        drop_samples_input = drop_samples
    )

    result_list <- list(
        binding_profiles_data = binding_profiles_data,
        occupancy = occupancy,
        test_category = "expressed",
        metadata = metadata
    )
    if (isTRUE(plot_diagnostics)) {
        plot_input_diagnostics(result_list)
    }
    return(result_list)
}


#' Build data frames of binding profiles from bedGraph files
#'
#' Each file must be in bedGraph format: chr, start, end, value.
#' @param bedgraphs Character vector of file paths.
#' @return data.frame with merged intervals and all sample columns.
#' @noRd
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
            data.df <- merge(data.df, gr, by = c("chr", "start", "end"), all = FALSE)
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
#' @noRd
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
        numeric_cols <- names(mcols_gr)[vapply(mcols_gr, is.numeric, FUN.VALUE = logical(1))]

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
#' @noRd
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
        # Include expanded gene_name
        files <- c(files, expanded[grepl(pattern, expanded, ignore.case = TRUE)])
    }
    files <- unique(files)
    return(files[file.exists(files)])
}

#' Import a GFF file as a GRanges object
#' @param path File path (GFF/GTF/BED)
#' @return GRanges
#' @noRd
import_peaks <- function(path) {
    tryCatch(
        {
            import(path)
        },
        error = function(e) {
            stop(sprintf("Failed to read peaks file '%s':\n%s", path, conditionMessage(e)))
        }
    )
}

#' Import a bedGraph file as a data.frame (chr, start, end, value)
#' @param path File path (bedGraph)
#' @param colname Name of the value column (usually sample name).
#' @return data.frame
#' @noRd
import_bedgraph_as_df <- function(path, colname = "score") {
    gr <- tryCatch(
        {
            import(path, format = "bedGraph")
        },
        error = function(e) {
            stop(sprintf("Failed to read bedGraph file '%s':\n%s", path, conditionMessage(e)))
        }
    )
    df <- as.data.frame(gr)[, c("seqnames", "start", "end", "score")]
    names(df) <- c("chr", "start", "end", colname)

    # Test for gapped offset (caused when loading closed rather than half-open datasets) and correct if present
    gaps <- df$start[-1] - df$end[-nrow(df)]
    tab <- table(gaps)
    if (as.integer(names(which.max(tab))) == 2) {
        # start of fragment(n+1) is always 2bp away from end of fragment(n): needs correcting
        df$end <- df$end + 1
    }

    df
}
