library(CellChat)

data.input <- GetAssayData(integrated, assay = "RNA", slot = "data")
meta <- data.frame(group = Idents(integrated), row.names = colnames(integrated))
cellchat <- createCellChat(data.input, meta = meta, group.by = "group")
cellchat@DB <- CellChatDB.human

# Analysis workflow
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualize communication circle plot
netVisual_circle(cellchat@net$count)