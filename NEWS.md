# damidBind 1.1.1
* NEW: Added additional diagnostic plots (mostly used to generate supplementary figures in the damidBind manuscript): `plot_catada_mean_variance()` and `test_weighting_vs_bias_artifact()`
* NEW: `load_data_peaks()` and `load_data_genes()` now support multiple normalisation options: quantile normalisation ("quantile"), cyclic LOESS norm ("loess"), and RPM (reads per million) norm ("rpm".  The default is "none" (no normalisation).
* NEW: a `pre_scale` option has been added to both `load_data_peaks()` and `load_data_genes()`, which will perform zero-preserving (uncentred) scaling (via base R `scale(center=FALSE, scale=TRUE)`) prior to normalisation.  The default is FALSE.
* NEW: all input metadata and genome annotation metadata are saved in the `DamIDResults` S4 object.  Basic metadata are displayed with the default `show()` method, and all are accessible via the `metadata()` accessor.
* NEW: `plot_volcano()` now takes a new `labels` option, values of `"all"`, `"highlight"`, or a character vector containing one or more group names from the `highlight` list.  Only this group will be labelled on the plot.  The original label configuration options still remain.
* NEW: `plot_volcano()` supports highlight groups of either loci or gene names, including different categories in the same plot.  The default "auto" behaviour matches a highlight group against whichever definition has the most matches, but use the new `highlight_by` option to enforce matching by either "gene_name" or "id" (i.e. locus).
* CHANGED: deprecated the `quantile_norm` option for `load_data_peaks()` and `load_data_genes()`.
* Added a citation to the current bioRxiv preprint on damidBind

# damidBind 1.0.0
* damidBind was automatically bumped to version 1.0.0 via Biocondutor release 3.23.

# damidBind 0.99.14
* NEW: Added `average_tracks` option to `browse_igv_regions()`.  When `TRUE`, this feature displays averaged binding tracks per condition rather than replicates.
* NEW: Added `export_data_archive` option to `browse_igv_regions()`.  When set to a valid file path, this parameter will export the tracks that would be displayed in the Shiny IGV window as bedGraph or BED files in a zip archive, instead of launching the Shiny app.  These files can be used to render final publication figures via external utilities such as `pyGenomeTracks` if desired.  This option can of course be combined with `average_tracks` above.

# damidBind 0.99.12
* CHANGED: major rewrite of the occupancy FDR logic for gene expression status.  Briefly:
    - Model regression fits have been improved with weighted natural spline fits for Tier 2 regressions (relationship is not log-linear)
    - Threshold scores for null model fitting are now determined from input data
    - Diagnostic plots for Tier 2 regressions have been added
    - Empirical p-values now aggregated on a condition level by default, using either Stouffer's (default) or Fisher's method
    - Occupancy p-values are now calculated by default during data loading, and merged during `differential_binding()`, before BH adjustment is applied to the aggregated p-values per condition.  See the relevant man pages for more details.
* CHANGED: internal logic for handling condition names and matches
* FIXED: Venn diagram universe logic now correctly determines the expressed gene universe when handling RNA Polymerase datasets

# damidBind 0.99.10
* NEW: plot diagnostics for sample loading (PCA, clustered correlation heatmap)
* NEW: plot diagnostics for limma functions (eBayes moderation / SA plots)
* NEW: limma functions now use `trend` and `robust` by default to better fit heteroscedastic DamID data
* NEW: new function `extract_unique_sample_ids` generates simplifed, unique sample names from complex filenames for display (in diagnostic plots, IGV/Shiny browser)
* CHANGED: dense point labelling function now uses exact `dbscan` kNN functions, not approximate HNSW functions
* CHANGED: differential threshold filtering can now be set to non-zero values
* CHANGED: default replicate filtering for differential analyses is set to the minimum number of the two condition replicates

# damidBind 0.99.9
* FIXED: igvShiny code now correctly handles the new internal GRanges binding profile data objects

# damidBind 0.99.8
* Unicode removed from documentation to allow latex generation :/

# damidBind 0.99.6
* NEW: Differential analyses now optionally screen out loci with negative signal in a specified number of replicates per sample (default is 2).  These loci are removed before differential analysis.
* NEW: Geneset universe calculations now take into account the low-signal screening; universe is determined only on the loci passing filter.
* NEW: Loci labels are optionally sampled to allow labels even in highly-dense plot regions.  Sampling algorithm determines a KNN graph from points, then uses this to select only a subset of total labelled points.  A default parameter is provided, but all algorithm parameters are fully customisable.
* NEW: Test conditions specified via the `cond` parameter can now use regexes, when `regex = TRUE` is set.
* CHANGED: Test conditions specified via `cond` now use a simple named vector, merging `cond` and `cond_names`.
* CHANGED: Improved and fixed underlying logic and tests to handle the `cond` parameter change
* CHANGED: Substantially improved legend guide handling parameters and logic for volcano plots
* FIXED: a number of minor bugfixes, code refactoring and tidying.

# damidBind 0.99.4
* NEW: Added in FDR modelling and calculation routines for RNA Polymerase occupancy
* NEW: plot_volcano and plot_venn functions can now filter for FDR when plotting
* NEW: expressed() accessor method for obtaining all genes passing an FDR threshold for a specific condition
* CHANGED: binding profile data is now stored as a GRanges object
* CHANGED: analyse_go_terms now returns the clusterProfiler results object for downstream use
* Begun incorporating gseGO functionality into analyse_go_terms
* Improved documentation of DamIDResults class and accessors
* Added more details to the package vignette
* Many small bug fixes

# damidBind 0.99.2
* NEW: Added button in IGV viewer to save SVG (useful when viewing window is horizontally compressed and button in the IGV interface is not displayed)
* NEW: New legend customisation options for volcano plots (including a legend internal to the plot frame, now default)
* NEW: More GO term plot customisation options
* FIXED: All volcano plot point labels are now rendered as a single geom layer, preventing ugly overlaps
* Substantial internal refactoring
* Many small bug fixes
* Updated unit tests

# damidBind 0.99.0
* Initial submission to Bioconductor.
