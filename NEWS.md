# damidBind 0.99.7
*   Unicode removed from documentation to allow latex generation :/

# damidBind 0.99.6
*   NEW: Differential analyses now optionally screen out loci with negative signal in a specified number of replicates per sample (default is 2).  These loci are removed before differential analysis.
*   NEW: Geneset universe calculations now take into account the low-signal screening; universe is determined only on the loci passing filter.
*   NEW: Loci labels are optionally sampled to allow labels even in highly-dense plot regions.  Sampling algorithm determines a KNN graph from points, then uses this to select only a subset of total labelled points.  A default parameter is provided, but all algorithm parameters are fully customisable.
*   NEW: Test conditions specified via the `cond` parameter can now use regexes, when `regex = TRUE` is set.
*   CHANGED: Test conditions specified via `cond` now use a simple named vector, merging `cond` and `cond_names`.
*   CHANGED: Improved and fixed underlying logic and tests to handle the `cond` parameter change
*   CHANGED: Substantially improved legend guide handling parameters and logic for volcano plots
*   FIXED: a number of minor bugfixes, code refactoring and tidying.

# damidBind 0.99.4
*   NEW: Added in FDR modelling and calculation routines for RNA Polymerase occupancy
*   NEW: plot_volcano and plot_venn functions can now filter for FDR when plotting
*   NEW: expressed() accessor method for obtaining all genes passing an FDR threshold for a specific condition
*   CHANGED: binding profile data is now stored as a GRanges object
*   CHANGED: analyse_go_terms now returns the clusterProfiler results object for downstream use
*   Begun incorporating gseGO functionality into analyse_go_terms
*   Improved documentation of DamIDResults class and accessors
*   Added more details to the package vignette
*   Many small bug fixes

# damidBind 0.99.2
*   NEW: Added button in IGV viewer to save SVG (useful when viewing window is horizontally compressed and button in the IGV interface is not displayed)
*   NEW: New legend customisation options for volcano plots (including a legend internal to the plot frame, now default)
*   NEW: More GO term plot customisation options
*   FIXED: All volcano plot point labels are now rendered as a single geom layer, preventing ugly overlaps
*   Substantial internal refactoring
*   Many small bug fixes
*   Updated unit tests

# damidBind 0.99.0
*   Initial submission to Bioconductor.
