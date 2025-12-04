# Microbiome analysis helpers reflecting existing patterns.
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(reshape2)
})

# Run Adonis on microbiome data (samples x features) with full metadata.
run_adonis_full <- function(feature_mat, metadata) {
  stopifnot(nrow(feature_mat) == nrow(metadata))
  adonis(feature_mat ~ ., data = metadata)
}

# Classical MDS (PCA on distance) with ellipse plot.
mds_ellipse_plot <- function(feature_mat, groups, k = 2, colors = c("red", "blue")) {
  mds <- cmdscale(dist(feature_mat), k = k, eig = TRUE)
  pts <- data.frame(mds$points)
  colnames(pts) <- paste0("X", seq_len(ncol(pts)))
  pts$group <- groups
  ggplot(pts, aes_string(x = "X1", y = "X2", color = "group", fill = "group")) +
    geom_point(size = 4, alpha = 0.6) +
    stat_ellipse(geom = "polygon", alpha = 0.5, level = 0.9) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_bw()
}

# Long-format boxplot for compositional data (genus/species) already melted.
composition_boxplot <- function(df_long, palette = NULL, angle = 60) {
  p <- ggplot(df_long, aes(x = variable, y = value, fill = groups)) +
    geom_boxplot(notch = FALSE, alpha = 0.95,
                 outlier.colour = "black",
                 outlier.shape = 16,
                 outlier.size = 0.65) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = angle, hjust = 1),
          legend.position = "top")
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }
  p
}
