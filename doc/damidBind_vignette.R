## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  comment = '##',
  warning = FALSE
)
library(DT)

## ----quickstart, eval=FALSE, echo=TRUE----------------------------------------
# ## Example code only, not run:
# 
# # Load up the data (use load_data.genes() for RNA Polymerase occupancy data)
# input <- load_data.peaks(
#   binding_profiles_path = "path/to/binding_profile_bedgraphs",
#   peaks_path = "path/to/peak_gffs_or_beds"
# ) # add quantile_norm = T if appropriate
# 
# # Deterimine differential binding  (use differential_accesibility() for CATaDa chromatin accessibility data)
# input.diff <- differential_binding(
#   input,
#   cond = c(
#     "Condition 1 identifying string in filenames",
#     "Condition 2 identifying string in filenames"
#   )
# )
# 
# # The result 'input.diff' is a formal S4 object.
# # You can see a summary by simply typing its name:
# input.diff
# 
# # View the proporition of differentially bound loci
# plot_venn(input.diff)
# 
# # Plot the differential binding, labelling associated genes with outliers
# plot_volcano(input.diff)
# 
# # Analyse GO enrichment in peaks associated with one condition
# analyse_go_enrichment(
#   input.diff,
#   direction = "Condition 1 identifier set with differential_binding() above"
# )
# 
# # View the differentially bound regions in an IGV browser window, with an interactive table of bound regions
# browse_igv_regions(input.diff)
# 
# # Apply additional functions on the differential binding results
# # Use the '@' accessor to get data 'slots' from the S4 object
# my_custom_function(input.diff@analysis)

## ----load_library-------------------------------------------------------------
library(damidBind)

## ----data_files---------------------------------------------------------------
data_dir = system.file("extdata", package = "damidBind")

# Show the files present for clarity in this vignette example:
files <- list.files(data_dir)
print(files)

## ----load_data----------------------------------------------------------------
input.bsh = load_data.peaks(
  binding_profiles_path = data_dir,
  peaks_path = data_dir, 
  quantile_norm = T
)

## ----differential_binding-----------------------------------------------------
diff.bsh = differential_binding(
  input.bsh,
  cond = c("L4","L5"),
  cond_names = c("L4 neurons","L5 neurons")
)

## ----show_object--------------------------------------------------------------
# See a summary of the results object
diff.bsh

## ----venn, fig.cap="A venn diagram of significantly bound loci by Bsh in L4 and L5 neurons.  The exclusive parts of each set represent regions that are differentially bound between the two conditions.", fig.small=TRUE, message=FALSE----
plot_venn(diff.bsh)

## ----volcanoPlotSimple, fig.cap="Differential binding of Bsh in L4 and L5 neuronal subtypes.  Genes associated with differentially bound peaks are displayed; the limitations of label overlaps means that only outliers are labelled.  (Dataset from chromosome 2L only)"----
plot_volcano(
  diff.bsh
)

## ----volcanoPlotCleanNames, fig.cap="Differential binding of Bsh in L4 and L5 neuronal subtypes.  Genes associated with differentially bound peaks are displayed, after some common, but less useful, gene label classes are removed.  (Dataset from chromosome 2L only)"----
plot_volcano(
  diff.bsh,
  label_config = list(clean_names=T)
)

## ----volcanoPlotHighlightedGenes, fig.cap="Differential binding of Bsh in L4 and L5 neuronal subtypes.  Genes that are specifically expressed in L4 neurons are highlighted.  (Dataset from chromosome 2L only)"----
L4_only_genes = c("Mp", "tnc", "grn", "rut", "mtd", "rdgB", "Octbeta2R", "msi", "Octbeta3R", "beat-IIIb", "ap", "Fili", "LRP1", "CG7378", "CG13698", "twit", "CG9336", "tok", "CG12991", "dpr1", "CG42339", "beat-IIb", "mav", "CG34377", "alpha-Man-IIb", "Pli", "CG32428", "osp", "Pka-R2", "CG15202", "CG8916", "CG15894", "side", "CG42258", "CHES-1-like", "SP2353", "CG44838", "Atg1", "Traf4", "DIP-beta", "KCNQ", "metro", "nAChRalpha1", "path", "CG10527", "Pde8", "CG30116", "CG7985", "CG1688", "dpr12", "pigs", "Eip63F-1", "CG14795", "2mit", "CG42340", "BicD", "CG18265", "hppy", "5-HT1A", "Chd64", "CG33090", "Dyb", "Btk29A", "Apc", "Rox8", "nAChRalpha5", "CG42748", "CG3257", "CG2269", "beat-IV", "CG8086", "glec", "CG31688", "oaf", "Drl-2", "CG8188", "aos", "CG31676", "REPTOR", "RabX4", "alt", "Pura", "DIP1", "ewg", "side-VIII", "nAChRalpha7", "Alh", "kug", "Ca-Ma2d", "bru2", "CG43737", "lncRNA:CR44024", "lncRNA:CR46006", "Had1", "CG3961", "comm", "Toll-6", "CG13685", "tow", "CG10019")

plot_volcano(
  diff.bsh,
  label_config = NULL,
  highlight = list(
    "L4 specific" = L4_only_genes
  ),
  highlight_config = list(
    size=1.5,
    label=T
  )
)

## ----volcannoPlotMultipleHighlights, fig.cap="Differential binding of Bsh in L4 and L5 neuronal subtypes.  Genes that are specifically expressed in each subtype are highlighted.   (Dataset from chromosome 2L only)"----
L5_only_genes = c("Ptth", "Nep2", "kek1", "CG4168", "kek3", "CG6959", "Dtg", "ND-23", "Scp2", "Octalpha2R", "Hs6st", "CG16791", "SKIP", "LpR1", "RpL34a", "Ald1", "CG10011", "heph", "nolo", "Act42A","Fkbp12", "Pkc53E", "AstC-R1", "Muc14A", "CG33543", "ChAT", "Act5C", "Ptpmeg2", "fabp", "CG31221", "Octbeta1R", "CG14669", "sdk", "Shawl", "side-V", "NaCP60E", "sif", "OtopLc", "side-II", "kuz", "CG42540", "Dscam3", "haf", "CG42673", "pdm3", "tinc", "CG42750", "sdt", "Nuak1", "Hk", "scrib", "tsr", "dpr20", "GluRIB", "CG43902", "CG44242", "Dscam2", "CG44422", "lncRNA:CR45312", "Scsalpha1", "Rop", "Con", "Hsc70-3", "dpr8", "eag", "ND-18", "Nrt", "CG17839", "fz", "CG32137", "Rh7", "Sod1", "CG32052", "dpr6", "Hsp67Ba", "axed", "GluRIA", "robo2")


plot_volcano(
  diff.bsh,
  label_config = NULL,
  highlight = list(
    "L4 specific" = L4_only_genes,
    "L5 specific" = L5_only_genes
  ),
  highlight_config = list(
    size=1.5,
    label=F
  )
)

## ----igv_shiny----------------------------------------------------------------
## Interactive, blocking session, uncomment to run
# browse_igv_regions(diff.bsh)

## ----gsea, fig.cap="Enriched GO terms for genes associated with differential Bsh binding in L4 neuron.  Dataset from chromosome 2L only."----
go.bsh_l4 = analyse_go_enrichment(
  diff.bsh, 
  direction="L4", 
  org_db = org.Dm.eg.db::org.Dm.eg.db
)

## ----paged_table, echo=FALSE, results='asis'----------------------------------
# Display the data frame as a paged table
DT::datatable(
  head(diff.bsh@analysis,n=20),
  options = list(
    pageLength = 10, # Number of rows to show per page
    scrollX = TRUE   # Enable horizontal scrolling for wide tables
  )
) %>% DT::formatRound(columns=c(1:8,11),digits=2)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

