# Multi-omics helpers reflecting repeated steps (binning, merging by group).
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Quantile binning into Q1-Q4 (or custom labels).
quantile_bin <- function(x, probs = c(0.25, 0.5, 0.75, 1), labels = NULL, na.rm = TRUE) {
  qs <- quantile(x, probs = probs, na.rm = na.rm)
  labs <- labels
  if (is.null(labels)) {
    labs <- paste0("Q", seq_along(probs))
  }
  cut(x,
      breaks = c(-Inf, qs),
      labels = labs,
      include.lowest = TRUE,
      right = TRUE)
}

# Merge feature matrix with metadata by row/column names (samples).
merge_features_metadata <- function(feature_mat, metadata, feature_id_col = NULL, meta_id_col = NULL) {
  if (!is.null(feature_id_col)) {
    rownames(feature_mat) <- feature_mat[[feature_id_col]]
    feature_mat <- feature_mat[, setdiff(colnames(feature_mat), feature_id_col), drop = FALSE]
  }
  if (!is.null(meta_id_col)) {
    rownames(metadata) <- metadata[[meta_id_col]]
  }
  stopifnot(all(colnames(feature_mat) %in% rownames(metadata)))
  data.frame(sample = colnames(feature_mat)) %>%
    left_join(metadata %>% tibble::rownames_to_column("sample"), by = "sample")
}

# Filter matrix columns by a metadata group.
subset_by_group <- function(feature_mat, metadata, group_col, group_value) {
  stopifnot(group_col %in% colnames(metadata))
  keep <- rownames(metadata)[metadata[[group_col]] %in% group_value]
  feature_mat[, colnames(feature_mat) %in% keep, drop = FALSE]
}
