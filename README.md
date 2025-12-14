# Introduction

DamID is a highly-sensitive means to profile the genome-wide association of proteins with chromatin in living eukaryotic cells, without fixation or the use of antibodies. Cell-type specific techniques such as Targeted DamID to profile protein binding, and CATaDA to profile chromatin accessibility, have made the technique an extremely powerful tool to understand the binding of transcription factors, chromatin proteins, RNA polymerase and chromatin changes during development and disease.

The `damidBind` package provides a straight-forward, statistically-robust formal analysis pipeline to analyse and explore differential DamID binding, gene transcription or chromatin accessibility between two conditions. The package imports processed data from DamID-seq experiments, in the form of binding bedgraphs and GFF peak calls. After optionally normalising data, combining peaks across replicates and determining per-replicate peak occupancy, the package links bound loci to nearby genes. It then uses either `limma` (for conventional log2 ratio DamID binding data) or `NOIseq` (for counts-based CATaDa chromatin accessibility data) to identify differentially-enriched regions between two conditions.

The package provides a number of visualisation tools (volcano plots, gene ontology enrichment plots via `clusterProfiler` and proportional Venn diagrams via `BioVenn` for downstream data exploration and analysis. An powerful, interactive IGV genome browser interface (powered by `Shiny` and `igvShiny`) allows users to rapidly and intuitively assess significant differentially-bound regions in their genomic context.

Although extensive customisation options are available if required, much of the data handling by `damidBind` is taken care of automatically, with sensible defaults assumed. To move from loading raw data to visualising differentially-enriched regions on a volcano plot or browsing enriched regions in an interactive IGV window is a simple three command procedure.

# Installation

You can install damidBind either from Bioconductor (stable) or Github (latest)

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# To install from Bioconductor (requires Bioconductor release 3.23):
BiocManager::install(version='devel')
BiocManager::install("damidBind")

# ... or to install the latest build from Github 
# (requires R>4.4, should work on most recent Bioconductor releases):
BiocManager::install("marshall-lab/damidBind")
```

# Vignette

A more complete guide to using `damidBind` can be found in the [damidBind vignette](https://marshall-lab.org/damidBind/articles/damidBind_vignette.html).

# Quick start guide

Using `damidBind`, in eight easy steps:

``` r
# Load up the data 
# (use load_data_genes() for RNA Polymerase occupancy data)
input <- load_data_peaks(
    binding_profiles_path = "path/to/binding_profile_bedgraphs",
    peaks_path = "path/to/peak_gffs_or_beds"
    # add quantile_norm = TRUE if appropriate
) 

# Determine differential binding 
# (use differential_accessibility() for CATaDa chromatin accessibility data)
input.diff <- differential_binding(
    input,
    cond = c(
        "Display name 1" = "Condition 1 identifying string in filenames",
        "Display name 2" = "Condition 2 identifying string in filenames"
    )
)

# The result 'input.diff' is a formal S4 object.
# You can see a summary by typing its name:
input.diff

# View the proportion of differentially bound loci
plot_venn(input.diff)

# Plot the differential binding, labelling associated genes with outliers
plot_volcano(input.diff)

# Analyse GO enrichment in peaks associated with one condition
analyse_go_enrichment(
    input.diff,
    direction = "Condition 1 identifier set with differential_binding() above"
)

# View the differentially bound loci in an Shiny/IGV browser web app, 
# with an interactive, searchable, sortable table of bound regions
browse_igv_regions(input.diff)

# Apply additional functions on the differential binding results
my_custom_function(analysisTable(input.diff))
```
