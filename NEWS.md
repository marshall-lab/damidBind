# damidBind 0.99.4
*   NEW: Added in FDR modelling and calculation routines for RNA Polymerase occupancy
*   NEW: plot_volcano and plot_venn functions can now filter for FDR when plotting
*   NEW: expressed() accessor method for obtaining all genes passing an FDR threshold for a specific condition
*   CHANGED: binding profile data is now stored as a GRanges object
*   CHANGED: analyse_go_terms now returns the clusterProfiler results object for downsteam use
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
