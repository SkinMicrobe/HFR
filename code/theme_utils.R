# Reusable ggplot themes extracted from repeated usage.
suppressPackageStartupMessages({
  library(ggplot2)
})

# Theme similar to the boxplot styling used across scripts.
theme_long_xlabels <- function(angle = 60) {
  theme_bw() +
    theme(plot.title = element_text(size = rel(2), hjust = 0.5),
          axis.title = element_text(size = rel(2.5)),
          axis.text = element_text(size = rel(1.9)),
          axis.text.x = element_text(face = "plain", size = 15, angle = angle, color = "black"),
          axis.text.y = element_text(face = "plain", size = 15, angle = 0, color = "black"),
          panel.grid.major = element_line(color = "white"),
          panel.grid.minor = element_line(color = "white"),
          panel.border = element_rect(color = "white"),
          axis.line = element_line(color = "black", size = 1.0),
          legend.key.size = unit(.5, "inches"),
          legend.title = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.box.just = "top",
          legend.text = element_text(colour = "black", size = 10, face = "plain"))
}
