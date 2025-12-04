# Code Utilities (Plain English)

- **tempdir_utils.R**: Forces R’s temp directory to a fixed path for packages that need a stable temp location.
- **seurat_utils.R**: Encapsulates the routine Seurat pipeline: create objects, QC, normalize, find variable genes, cluster, embed, and label cell types.
- **cellchat_utils.R**: Runs the standard CellChat preprocessing to prepare communication networks from a Seurat object.
- **cellchat_plot_utils.R**: Provides quick plotting for CellChat outputs: overviews, pathway aggregates, and bubble charts.
- **limma_utils.R**: Runs limma differential testing (design + contrast) and returns tidy results.
- **kegg_utils.R**: Performs human KEGG enrichment from gene symbols and produces a dot plot.
- **enrichment_utils.R**: GSEA-style KEGG on ranked genes plus helpers to retrieve top pathways.
- **microbiome_utils.R**: Core microbiome analyses: Adonis, distance-based MDS with ellipses, and composition boxplots.
- **microbiome_prediction_utils.R**: Workflow pieces for microbiome prediction: Adonis with FDR, TableOne, stepwise logistic, and OR/CI data for forest-style plots.
- **iobr_utils.R**: Thin wrappers around IOBR for signature scoring and count-to-TPM conversion.
- **lefse_utils.R**: Wraps the microeco LEfSe workflow and its visualizations (bars, cladograms).
- **metabolism_utils.R**: Cleans and clusters lipid/metabolite matrices, scales/clips for heatmaps, and plots correlations with marginals.
- **multiomics_utils.R**: General multi-omics helpers for quantile binning, merging features with metadata, and subsetting by group.
- **plot_utils.R**: A ready-made boxplot with significance labels for long-format data.
- **heatmap_utils.R**: Scales/clips matrices and draws heatmaps with optional annotations.
- **theme_utils.R**: A reusable ggplot theme tuned for long x-axis labels and a clean top legend.
