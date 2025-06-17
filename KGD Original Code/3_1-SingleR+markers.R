library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(tidyverse)

#Cell type
#SingleR for cell annotation
#install SingleR
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
##load the merged Seurat object
TN.combined <- readRDS(file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/TN.combined_dim30-DF-35.rds")


#SingleR
counts <- GetAssayData(TN.combined_DF)

library(ExperimentHub)
eh <- ExperimentHub()
removeCache(eh)




#hpca database
hpca.se  <- HumanPrimaryCellAtlasData()
hpca.se
pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
pred.hpca
table(pred.hpca$pruned.labels)
clustering.table_hpca <- table(pred.hpca@listData[["pruned.labels"]], TN.combined@active.ident)
clustering.table_hpca
write.csv(clustering.table_hpca , file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_SingleR_hpca.csv", col.names = TRUE)

#bpe database
bpe.se <- BlueprintEncodeData()
bpe.se
pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test=1,
                    labels = bpe.se$label.main)
pred.bpe
table(pred.bpe$pruned.labels)
clustering.table_bpe <- table(pred.bpe@listData[["pruned.labels"]], TN.combined@active.ident)
clustering.table_bpe
write.csv(clustering.table_bpe , file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_SingleR_bpe.csv", col.names = TRUE)

#load hpca analysis
clustering.table_hpca <- read.csv("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_SingleR_hpca.csv")
rownames(clustering.table_hpca) <- clustering.table_hpca[,1]
clustering.table_hpca <- clustering.table_hpca[,-1]
clustering.table_hpca["annotation",] <- rownames(clustering.table_hpca)[apply(clustering.table_hpca,2,which.max)]
write.csv(clustering.table_hpca , file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_SingleR_hpca_summary.csv", col.names = TRUE)


#load bpe analysis
clustering.table_bpe <- read.csv("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_SingleR_bpe.csv")
rownames(clustering.table_bpe) <- clustering.table_bpe[,1]
clustering.table_bpe <- clustering.table_bpe[,-1]
clustering.table_bpe["annotation",] <- rownames(clustering.table_bpe)[apply(clustering.table_bpe,2,which.max)]
write.csv(clustering.table_bpe , file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_SingleR_bpe_summary.csv", col.names = TRUE)

# #mouse database
# mouse.se <- MouseRNAseqData()
# mouse.se
# pred.mouse <- SingleR(test = counts, ref = mouse.se, assay.type.test=1,
#                       labels = mouse.se$label.main)
# pred.mouse
# table(pred.mouse$pruned.labels)
# clustering.table_mouse <- table(pred.mouse@listData[["pruned.labels"]], TN.combined@active.ident)
# clustering.table_mouse
# write.csv(clustering.table_mouse , file = "C:/Users/user/Desktop/SingleR_mouse.csv", col.names = TRUE)


##Vlnplot for checking cell type
DefaultAssay(TN.combined) <- "RNA"

##Classical markers
#Epithelial cells
markers.to.plot <- c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP")
Classical_markers_Epi <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Epithelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Epi.pdf", width = 15, height = 15)
Classical_markers_Epi
dev.off()

#Sweat gland cell
markers.to.plot <- c("MUCL1","PIP","AQP5")
Classical_markers_Sweat_gland <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Sweat gland cell") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Sweat_gland.pdf", width = 15, height = 15)
Classical_markers_Sweat_gland
dev.off()

#SMC and pili muscle
markers.to.plot <- c("MCAM","ACTA2","MYL9","TAGLN","MYH11")
Classical_markers_SMC <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="SMC and pili muscle") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_SMC.pdf", width = 15, height = 15)
Classical_markers_SMC
dev.off()


#pericyte
markers.to.plot <- c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11")
Classical_markers_peicyte <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Pericyte") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_pericyte.pdf", width = 15, height = 15)
Classical_markers_peicyte
dev.off()


#fibroblasts
markers.to.plot <- c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248")
Classical_markers_fibroblasts <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_fibroblasts.pdf", width = 15, height = 15)
Classical_markers_fibroblasts
dev.off()


#Vascular endothelial cells
markers.to.plot <- c("PECAM1","VWF")
Classical_markers_vEC <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Vascular endothelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_vEC.pdf", width = 15, height = 15)
Classical_markers_vEC
dev.off()


#lymphatic endothelial cells
markers.to.plot <- c("PROX1","LYVE1")
Classical_markers_lEC <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Lymphatic endothelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_lEC.pdf", width = 15, height = 15)
Classical_markers_lEC
dev.off()

#T cells
markers.to.plot <- c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7")
Classical_markers_TC <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="T cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_TC.pdf", width = 15, height = 15)
Classical_markers_TC
dev.off()


#NK cells
markers.to.plot <- c("GNLY","NKG7")
Classical_markers_NK <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="NK cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_NK.pdf", width = 15, height = 15)
Classical_markers_NK
dev.off()


#B cells
markers.to.plot <- c("MS4A1","CD79A","SEC11C","CD79B")
Classical_markers_BC <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="B cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_BC.pdf", width = 15, height = 15)
Classical_markers_BC
dev.off()


#Plasma cells
markers.to.plot <- c("IGJ","MZB1","XBP1","CD79A","CD79B")
Classical_markers_Plasma <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Plasma cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Plasma.pdf", width = 15, height = 15)
Classical_markers_Plasma
dev.off()


#Myeloid cells (monocyte)
markers.to.plot <- c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ")
Classical_markers_monocyte <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Monocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_monocyte.pdf", width = 15, height = 15)
Classical_markers_monocyte
dev.off()

#Myeloid cells (dendritic cell)
markers.to.plot <- c("IRF7","HLA-DRA","LYZ","S100B","CD1C")
Classical_markers_dendritic <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Dendritic cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Dendritic cells.pdf", width = 15, height = 15)
Classical_markers_monocyte
dev.off()

#Granulocytes&Neutrophils
markers.to.plot <- c("ITGAX","ITGAM","FCGR2A","ANPEP")
Classical_markers_Neutrophils <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Granulocytes and Neutrophils") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Neutrophils.pdf", width = 15, height = 15)
Classical_markers_Neutrophils
dev.off()


#Mast cells
markers.to.plot <- c("ADCYAP1","CPA3","TPSAB1","VWA5A")
Classical_markers_Mast <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Mast cells") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Mast.pdf", width = 15, height = 15)
Classical_markers_Mast
dev.off()


#Melanocytes
markers.to.plot <- c("DCT","MLANA")
Classical_markers_Melanocytes <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Melanocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))


pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Melanocytes.pdf", width = 15, height = 15)
Classical_markers_Melanocytes
dev.off()


#Neuronal cells
markers.to.plot <- c("NRXN1","SCN7A","CDH19", "S100B", "IGFBP5", "MIA", "EGFL8", "NGFR", "TYR")
Classical_markers_Neuronal <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Neuronal cells") + theme(plot.title = element_text(hjust = 0.5, size=24))


pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Neuronal.pdf", width = 15, height = 15)
Classical_markers_Neuronal
dev.off()

#Schwann cells
markers.to.plot <- c("NRXN1","CCN3","MPZ","PTN","S100B")
Classical_markers_Schwann <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Schwann cells") + theme(plot.title = element_text(hjust = 0.5, size=24))


pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Schwann.pdf", width = 15, height = 15)
Classical_markers_Schwann
dev.off()

#Adipocytes
markers.to.plot <- c("ADIPOQ","PLIN1","FABP4","LEP","GPD1","CD36")
Classical_markers_Adipocytes <- DotPlot(TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Adipocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))


pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Adipocytes.pdf", width = 15, height = 15)
Classical_markers_Adipocytes
dev.off()
