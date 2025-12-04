# CellChat pipeline helper.
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
})

# Run a streamlined CellChat workflow on a Seurat object.
process_cellchat <- function(seurat_obj,
                             group_by = "bigcelltype",
                             db = CellChatDB.human,
                             min_cells = 10,
                             pathways_use = NULL,
                             assay_slot = "data") {
  data.input <- GetAssayData(seurat_obj, slot = assay_slot)
  meta <- seurat_obj@meta.data
  stopifnot(group_by %in% colnames(meta))
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = group_by)
  cellchat@DB <- db
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat, pathways.use = pathways_use)
  cellchat <- aggregateNet(cellchat)
  cellchat
}
