# Heatmap helpers (scale/clip then pheatmap).
suppressPackageStartupMessages({
  library(pheatmap)
})

# Scale (center/scale), transpose optional, and clip to range for visualization.
scale_and_clip <- function(mat, transpose = FALSE, clip_min = -2.5, clip_max = 2.5) {
  x <- mat
  if (transpose) x <- t(x)
  x <- scale(x, center = TRUE, scale = TRUE)
  x[x > clip_max] <- clip_max
  x[x < clip_min] <- clip_min
  x
}

# Basic pheatmap wrapper with optional annotation colors.
draw_heatmap <- function(mat,
                         annotation = NULL,
                         annotation_colors = NULL,
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         fontsize_row = 5,
                         show_colnames = TRUE) {
  pheatmap(mat,
           annotation = annotation,
           annotation_colors = annotation_colors,
           cluster_cols = cluster_cols,
           cluster_rows = cluster_rows,
           fontsize_row = fontsize_row,
           show_colnames = show_colnames)
}
