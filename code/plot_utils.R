# Plot helpers for common long-format comparisons.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

plot_box_compare <- function(df_long,
                             x = "variable",
                             y = "value",
                             group = "groups",
                             palette = NULL,
                             notch = FALSE,
                             angle = 60) {
  p <- ggplot(df_long, aes_string(x = x, y = y, fill = group)) +
    geom_boxplot(notch = notch, alpha = 0.95,
                 outlier.colour = "black",
                 outlier.shape = 16,
                 outlier.size = 0.65) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = angle, hjust = 1)) +
    stat_compare_means(label = "p.signif")
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }
  p
}
