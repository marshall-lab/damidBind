# Data provenance for `inst/extdata`

This document describes the generation of the sample data files located in the `inst/extdata` directory of the `damidBind` package.

## 1. Source data

The data are a truncated subset of the following publicly-available dataset:

-   **Publication:** Xu C, Ramos TB, Marshall OJ & Doe CQ (2024). *Notch signaling and Bsh homeodomain activity are integrated to diversify Drosophila lamina neuron types*. eLife 12: RP90136.
-   **DOI:** [10.7554/eLife.90136](https://doi.org/10.7554/eLife.90136)
-   **GEO Accession:** [GSE247239](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247239)

## 2. Processed data generation pipeline

The files in `inst/extdata` were generated from the raw sequencing data (FASTQ files) available in the GEO repository mentioned above. These were processed using external command-line software as detailed below.

### 2.1 Generation of log2 ratio bedGraph files

The genome-wide log2(Dam-Bsh/Dam-only) binding profiles were generated from the raw FASTQ files using `damidseq_pipeline`.

-   **Software:** `damidseq_pipeline`
-   **Version:** v1.6
-   **Source:** <https://github.com/owenjm/damidseq_pipeline>
-   **Command:** The pipeline was run with default options against the *Drosophila melanogaster* dm6 reference genome assembly Bowtie2 indices, and the genome reference GATC fragment file:

``` sh
damidseq_pipeline \
    --bowtie2_genome_dir=/path/to/bowtie_indices/DmBDGP6 \
    --gatc_frag_file=/path/to/Dmel_BDGP6.GATC.gff \
    --data_dir=/path/to/fastq_raw_files/
```

producing genome-wide log2(Dam-Bsh/Dam-only) bedGraph binding profiles at GATC resolution.

### 2.2 Peak calling

Significant binding peaks were identified from each replicate's bedGraph file using `find_peaks`.

-   **Software:** `find_peaks`
-   **Version:** v1.2
-   **Source:** <https://github.com/owenjm/find_peaks>
-   **Command:** `find_peaks` was run with `--format=bed` and otherwise default options:

``` sh
find_peaks --format=bed *.bedgraph
```

generating `.peaks.bed` files for each biological replicate.

### 2.3 Compression

Both `bedGraph` and `bed` files were compressed via:

``` sh
 gzip *bedgraph *bed
```

### 2.4 Truncating for inst/data

To create a small example dataset, the full genome-wide bedGraph and peak files were truncated to include only data from chromosome 2L:

``` fish
for ext in bedgraph.gz bed.gz; \
  for i in *$ext; \
    gzip -dc $i | awk '$1 ~ /2L$/ || $0 ~ /^track/' | gzip > (basename $i $ext)2L.$ext; \
  end; \
end;
```

The resulting files are included in the `inst/extdata` directory.

The full, untruncated processed dataset is available at Zenodo ([10.5281/zenodo.16649477](https://doi.org/10.5281/zenodo.16649477)).
