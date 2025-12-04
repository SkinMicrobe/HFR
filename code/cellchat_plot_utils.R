# CellChat plotting helpers reflecting repeated usage patterns.
suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

# Circle overview of interaction counts and weights.
plot_cellchat_overview <- function(cellchat_obj) {
  group_size <- as.numeric(table(cellchat_obj@idents))
  p1 <- netVisual_circle(cellchat_obj@net$count,
                         vertex.weight = group_size,
                         weight.scale = TRUE,
                         label.edge = FALSE,
                         title.name = "Number of interactions")
  p2 <- netVisual_circle(cellchat_obj@net$weight,
                         vertex.weight = group_size,
                         weight.scale = TRUE,
                         label.edge = FALSE,
                         title.name = "Interaction weights/strength")
  p1 + p2
}

# Aggregate pathway visualizations (hierarchy, circle, heatmap).
plot_cellchat_pathways <- function(cellchat_obj,
                                   pathways,
                                   vertex_receiver = NULL,
                                   layout = c("circle", "hierarchy", "heatmap")) {
  layout <- match.arg(layout)
  if (layout == "circle") {
    netVisual_aggregate(cellchat_obj, signaling = pathways, layout = "circle")
  } else if (layout == "hierarchy") {
    netVisual_aggregate(cellchat_obj, signaling = pathways, vertex.receiver = vertex_receiver)
  } else {
    netVisual_heatmap(cellchat_obj, signaling = pathways, color.heatmap = "Reds")
  }
}

# Bubble plot for selected sources/targets or pathways.
plot_cellchat_bubble <- function(cellchat_obj,
                                 sources = NULL,
                                 targets = NULL,
                                 pathways = NULL,
                                 pairLR = NULL,
                                 remove_isolate = FALSE) {
  netVisual_bubble(cellchat_obj,
                   sources.use = sources,
                   targets.use = targets,
                   signaling = pathways,
                   pairLR.use = pairLR,
                   remove.isolate = remove_isolate)
}
