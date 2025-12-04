# Metabolomics / lipidomics helpers based on patterns in the scripts.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggExtra)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
})

# Prepare a lipid matrix: set rownames, drop duplicates, transpose to samples x features.
prepare_lipid_matrix <- function(df, id_col = 1, drop_cols = NULL) {
  x <- df
  if (!is.null(drop_cols)) {
    x <- x[, -drop_cols, drop = FALSE]
  }
  rn <- x[[id_col]]
  x <- x[, -id_col, drop = FALSE]
  x <- x[!duplicated(rn), , drop = FALSE]
  rownames(x) <- rn[!duplicated(rn)]
  t(as.data.frame(x, check.names = FALSE))
}

# K-means clustering on scaled data, returning matrix with cluster labels.
lipid_kmeans_clusters <- function(mat, centers = 2, transpose = TRUE) {
  x <- mat
  if (transpose) x <- t(x)
  km <- kmeans(x, centers = centers)
  out <- cbind(x, cluster = km$cluster)
  list(labeled_matrix = out, clustering = km)
}

# Heatmap-ready scaling and clipping (reuse if heatmap_utils is not sourced).
scale_clip_matrix <- function(mat, clip_min = -2.5, clip_max = 2.5, transpose = FALSE) {
  x <- mat
  if (transpose) x <- t(x)
  x <- scale(x, center = TRUE, scale = TRUE)
  x[x > clip_max] <- clip_max
  x[x < clip_min] <- clip_min
  x
}

# Draw a pheatmap with optional annotation for clusters.
lipid_heatmap <- function(mat_scaled, annotation = NULL, annotation_colors = NULL,
                          cluster_cols = FALSE, cluster_rows = TRUE,
                          fontsize_row = 5, show_colnames = TRUE,
                          bk_breaks = c(seq(-2.5, -0.1, by = 0.01), seq(0, 2.5, by = 0.01))) {
  pheatmap(mat_scaled,
           annotation = annotation,
           annotation_colors = annotation_colors,
           cluster_cols = cluster_cols,
           cluster_rows = cluster_rows,
           fontsize_row = fontsize_row,
           show_colnames = show_colnames,
           color = c(colorRampPalette(c("blue", "white"))(length(bk_breaks) / 2),
                     colorRampPalette(c("white", "red"))(length(bk_breaks) / 2)),
           breaks = bk_breaks)
}

# Correlation scatter with linear fit and marginal boxplots.
cor_scatter_marginal <- function(df, x, y, group = NULL, method = "spearman") {
  aes_map <- aes_string(x = x, y = y, color = group)
  p <- ggplot(df, aes_map) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "gray") +
    theme_bw()
  if (!is.null(group)) {
    p <- p + guides(color = guide_legend(title = NULL))
  }
  cor_val <- cor.test(df[[x]], df[[y]], method = method)
  p <- ggMarginal(p, type = "boxplot", groupColour = !is.null(group), groupFill = !is.null(group))
  list(plot = p, cor = cor_val)
}
