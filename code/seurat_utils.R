# Seurat workflow helpers extracted from repeated patterns in the scripts.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# Create Seurat object from raw count matrix and metadata.
build_seurat_object <- function(raw_data,
                                metadata,
                                min_cells = 0,
                                min_features = 0) {
  stopifnot(all(colnames(raw_data) %in% rownames(metadata)))
  CreateSeuratObject(counts = raw_data,
                     metadata = metadata,
                     min.cells = min_cells,
                     min.features = min_features)
}

# Basic QC: add mitochondrial and ribosomal percentages and apply cutoffs.
qc_filter_basic <- function(sce,
                            min_count = 0,
                            min_feature = 0,
                            mt_pattern = "^MT-",
                            rp_pattern = "^RP[SL][[:digit:]]",
                            percent_mt_max = Inf,
                            percent_rp_max = Inf) {
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = mt_pattern)
  sce[["percent.rp"]] <- PercentageFeatureSet(sce, pattern = rp_pattern)
  subset(sce,
         subset = nCount_RNA > min_count &
                  nFeature_RNA > min_feature &
                  percent.mt < percent_mt_max &
                  percent.rp < percent_rp_max)
}

# Normalize, find variable features, scale, and run PCA.
normalize_variable_pca <- function(sce,
                                   nfeatures = 2000,
                                   scale_factor = 1e6,
                                   npcs_plot = 20) {
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = scale_factor)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = nfeatures)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce, npcs = npcs_plot, verbose = FALSE)
  sce
}

# Neighbor graph, clustering, and TSNE/UMAP embedding.
cluster_and_embed <- function(sce,
                              dims_use = 1:20,
                              resolution = 1.2,
                              method = c("tsne", "umap")) {
  method <- match.arg(method)
  sce <- FindNeighbors(sce, dims = dims_use)
  sce <- FindClusters(sce, resolution = resolution)
  if (method == "tsne") {
    sce <- RunTSNE(sce, dims.use = dims_use)
  } else {
    sce <- RunUMAP(sce, dims = dims_use)
  }
  sce
}

# Map cluster IDs to cell type labels using a named list.
label_celltypes <- function(sce, mapping) {
  stopifnot("seurat_clusters" %in% colnames(sce@meta.data))
  cluster_ids <- as.character(sce@meta.data$seurat_clusters)
  sce@meta.data$celltype <- vapply(cluster_ids, function(cl) {
    hit <- names(mapping)[vapply(mapping, function(x) cl %in% x, logical(1))]
    ifelse(length(hit) > 0, hit[1], "Unlabeled")
  }, character(1))
  sce
}

# DotPlot for marker panels.
marker_dotplot <- function(sce, markers, assay = "RNA") {
  DotPlot(sce, features = markers, assay = assay)
}
