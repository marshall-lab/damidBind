# Introduction

DamID [@Steensel2000; @Steensel2001] is a highly-sensitive means to profile the genome-wide association of proteins with chromatin in living eukaryotic cells, without fixation or the use of antibodies. Cell-type specific techniques such as Targeted DamID [@Southall2013; @Marshall2016] to profile protein binding, and CATaDA [@Aughey2018] to profile chromatin accessibility, have made the technique an extremely powerful tool to understand the binding of transcription factors, chromatin proteins, RNA polymerase and chromatin changes during development and disease.

Despite the technique's growing popularity and adoption, no formal analysis pipeline or R package exists to analyse and explore differential DamID binding, gene transcription or chromatin accessibility between two conditions. The `damidBind` package provides this functionality.

`damidbind` imports processed data from DamID-seq experiments in the form of binding bedgraphs and GFF peak calls. After optionally normalising data, combining peaks across replicates and determining per-replicate peak occupancy, the package links bound loci to nearby genes. It then uses either `limma` (for conventional log2 ratio DamID binding data) or `NOIseq` (for counts-based CATaDa chromatin accessibility data) to identify differentially-enriched regions between two conditions. The package provides a number of visualisation tools (volcano plots, GSEA plots, Venn diagrams) for downstream data exploration and analysis. An interactive IGV genome browser interface (powered by `Shiny` and `igvShiny`) allows users to rapidly and intuitively assess significant differentially-bound regions.

Although extensive customisation options are available if required, much of the data handling by `damidBind` is taken care of automatically, with sensible defaults assumed. To move from loading raw data to visualising differentially-enriched regions on a volcano plot or browsing enriched regions in an interactive IGV window is a simple three command procedure.

# Quick start guide

Using `damidBind`, in seven easy steps:

```
## Example code only, not run:

# Load up the data (use load_data_genes() for RNA Polymerase occupancy data)
input <- load_data_peaks(
  binding_profiles_path = "path/to/binding_profile_bedgraphs",
  peaks_path = "path/to/peak_gffs_or_beds"
) # add quantile_norm = T if appropriate

# Deterimine differential binding  (use differential_accesibility() for CATaDa chromatin accessibility data)
input.diff <- differential_binding(
  input,
  cond = c(
    "Condition 1 identifying string in filenames",
    "Condition 2 identifying string in filenames"
  )
)

# The result 'input.diff' is a formal S4 object.
# You can see a summary by simply typing its name:
input.diff

# View the proporition of differentially bound loci
plot_venn(input.diff)

# Plot the differential binding, labelling associated genes with outliers
plot_volcano(input.diff)

# Analyse GO enrichment in peaks associated with one condition
analyse_go_enrichment(
  input.diff, 
  direction = "Condition 1 identifier set with differential_binding() above"
)

# View the differentially bound regions in an IGV browser window, with an interactive table of bound regions
browse_igv_regions(input.diff)

# Apply additional functions on the differential binding results
# Use the '@' accessor to get data 'slots' from the S4 object
my_custom_function(input.diff@analysis)
```

