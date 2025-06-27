
#### Load Packages ####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('patchwork')) {install.packages('patchwork'); library(patchwork)}


#### Create new folder ####
try({
  ## Create folder for KGD Lab
  Name_ExportFolder_KGD_CTAnnot <- paste0(Name_ExportFolder,"/","KGD_Lab_CTAnnot")
  # Create export folder if it does not exist
  if (!dir.exists(Name_ExportFolder_KGD_CTAnnot)){dir.create(Name_ExportFolder_KGD_CTAnnot)}
})


#### Set marker list ####
# 設定預設 Assay
DefaultAssay(seuratObject_Sample) <- "RNA"
Idents(seuratObject_Sample) <- "seurat_clusters"

# 建立各細胞類型與其 marker 的清單
marker_sets_KGD <- list(
  "Epithelial cells"          = c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP"),
  "Sweat gland cell"          = c("MUCL1","PIP","AQP5"),
  "SMC and pili muscle"       = c("MCAM","ACTA2","MYL9","TAGLN","MYH11"),
  "Pericyte"                  = c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11"),
  "Fibroblasts"               = c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248"),
  "Vascular endothelial cells"= c("PECAM1","VWF"),
  "Lymphatic endothelial cells" = c("PROX1","LYVE1"),
  "T cells"                   = c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7"),
  "NK cells"                  = c("GNLY","NKG7"),
  "B cells"                   = c("MS4A1","CD79A","SEC11C","CD79B"),
  "Plasma cells"              = c("IGJ","MZB1","XBP1","CD79A","CD79B"),
  "Monocytes"                 = c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ"),
  "Dendritic cells"           = c("IRF7","HLA-DRA","LYZ","S100B","CD1C"),
  "Granulocytes / Neutrophils"= c("ITGAX","ITGAM","FCGR2A","ANPEP"),
  "Mast cells"                = c("ADCYAP1","CPA3","TPSAB1","VWA5A"),
  "Melanocytes"               = c("DCT","MLANA"),
  "Neuronal cells"            = c("NRXN1","SCN7A","CDH19","S100B","IGFBP5","MIA","EGFL8","NGFR","TYR"),
  "Schwann cells"             = c("NRXN1","CCN3","MPZ","PTN","S100B"),
  "Adipocytes"                = c("ADIPOQ","PLIN1","FABP4","LEP","GPD1","CD36")
)


#### Visualization ####
# 逐一產生每個細胞類型的 DotPlot
dotplots <- lapply(names(marker_sets_KGD), function(cell_type) {
  DotPlot(
    object    = seuratObject_Sample,
    features  = marker_sets_KGD[[cell_type]],
    cols      = c("white", "darkred"),
    scale = FALSE,
    dot.scale = 8
  ) +
    RotatedAxis() +
    labs(title = cell_type) +
    theme(plot.title = element_text(hjust = 0.5, size = 24))
})

# 以 patchwork 合併所有 DotPlot；可依資料量調整 ncol 或 nrow
Plot_combined <- wrap_plots(dotplots, ncol = 3)



#### Export pdf ####
# 輸出成一張 PDF
pdf(paste0(Name_ExportFolder_KGD_CTAnnot, "/", Name_Export, "_Classical_markers_Combined_Jojie.pdf"), 
    width = 30, height = 60)
print(Plot_combined)
dev.off()

Plot_combined <- wrap_plots(dotplots, ncol = 4)
png(
  filename = paste0(Name_ExportFolder_KGD_CTAnnot, "/", Name_Export, "_Classical_markers_Combined_Jojie.png"),
  width = 30,        # 寬度 (英吋)
  height = 30,      # 高度 (英吋)
  units = "in",                 # 設定單位為英吋
  res = 600                     # 解析度 (dpi)，可依需求調整
)
print(Plot_combined)
dev.off()


# rm(dotplots, Plot_combined)


################################################################################
## Create marker dataframe
# 建立長格式資料框
marker_df_scMRMA_KGD <- data.frame(
  gene      = unlist(marker_sets_KGD, use.names = FALSE),
  CellType  = rep(names(marker_sets_KGD), times = lengths(marker_sets_KGD)),
  stringsAsFactors = FALSE
)

# 檢視結果
head(marker_df_scMRMA_KGD)





################################################################################
## # Old Version
# ##Vlnplot for checking cell type
# DefaultAssay(seuratObject_Sample) <- "RNA"
# 
# ##Classical markers
# #Epithelial cells
# markers.to.plot <- c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP")
# Classical_markers_Epi <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Epithelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Epi.pdf", width = 15, height = 15)
# Classical_markers_Epi
# dev.off()
# 
# #Sweat gland cell
# markers.to.plot <- c("MUCL1","PIP","AQP5")
# Classical_markers_Sweat_gland <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Sweat gland cell") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Sweat_gland.pdf", width = 15, height = 15)
# Classical_markers_Sweat_gland
# dev.off()
# 
# #SMC and pili muscle
# markers.to.plot <- c("MCAM","ACTA2","MYL9","TAGLN","MYH11")
# Classical_markers_SMC <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="SMC and pili muscle") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_SMC.pdf", width = 15, height = 15)
# Classical_markers_SMC
# dev.off()
# 
# 
# #pericyte
# markers.to.plot <- c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11")
# Classical_markers_peicyte <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Pericyte") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_pericyte.pdf", width = 15, height = 15)
# Classical_markers_peicyte
# dev.off()
# 
# 
# #fibroblasts
# markers.to.plot <- c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248")
# Classical_markers_fibroblasts <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_fibroblasts.pdf", width = 15, height = 15)
# Classical_markers_fibroblasts
# dev.off()
# 
# 
# #Vascular endothelial cells
# markers.to.plot <- c("PECAM1","VWF")
# Classical_markers_vEC <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Vascular endothelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_vEC.pdf", width = 15, height = 15)
# Classical_markers_vEC
# dev.off()
# 
# 
# #lymphatic endothelial cells
# markers.to.plot <- c("PROX1","LYVE1")
# Classical_markers_lEC <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Lymphatic endothelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_lEC.pdf", width = 15, height = 15)
# Classical_markers_lEC
# dev.off()
# 
# #T cells
# markers.to.plot <- c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7")
# Classical_markers_TC <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="T cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_TC.pdf", width = 15, height = 15)
# Classical_markers_TC
# dev.off()
# 
# 
# #NK cells
# markers.to.plot <- c("GNLY","NKG7")
# Classical_markers_NK <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="NK cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_NK.pdf", width = 15, height = 15)
# Classical_markers_NK
# dev.off()
# 
# 
# #B cells
# markers.to.plot <- c("MS4A1","CD79A","SEC11C","CD79B")
# Classical_markers_BC <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="B cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_BC.pdf", width = 15, height = 15)
# Classical_markers_BC
# dev.off()
# 
# 
# #Plasma cells
# markers.to.plot <- c("IGJ","MZB1","XBP1","CD79A","CD79B")
# Classical_markers_Plasma <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Plasma cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Plasma.pdf", width = 15, height = 15)
# Classical_markers_Plasma
# dev.off()
# 
# 
# #Myeloid cells (monocyte)
# markers.to.plot <- c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ")
# Classical_markers_monocyte <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Monocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_monocyte.pdf", width = 15, height = 15)
# Classical_markers_monocyte
# dev.off()
# 
# #Myeloid cells (dendritic cell)
# markers.to.plot <- c("IRF7","HLA-DRA","LYZ","S100B","CD1C")
# Classical_markers_dendritic <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Dendritic cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Dendritic cells.pdf", width = 15, height = 15)
# Classical_markers_monocyte
# dev.off()
# 
# #Granulocytes&Neutrophils
# markers.to.plot <- c("ITGAX","ITGAM","FCGR2A","ANPEP")
# Classical_markers_Neutrophils <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Granulocytes and Neutrophils") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Neutrophils.pdf", width = 15, height = 15)
# Classical_markers_Neutrophils
# dev.off()
# 
# 
# #Mast cells
# markers.to.plot <- c("ADCYAP1","CPA3","TPSAB1","VWA5A")
# Classical_markers_Mast <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Mast cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Mast.pdf", width = 15, height = 15)
# Classical_markers_Mast
# dev.off()
# 
# 
# #Melanocytes
# markers.to.plot <- c("DCT","MLANA")
# Classical_markers_Melanocytes <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Melanocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Melanocytes.pdf", width = 15, height = 15)
# Classical_markers_Melanocytes
# dev.off()
# 
# 
# #Neuronal cells
# markers.to.plot <- c("NRXN1","SCN7A","CDH19", "S100B", "IGFBP5", "MIA", "EGFL8", "NGFR", "TYR")
# Classical_markers_Neuronal <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Neuronal cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Neuronal.pdf", width = 15, height = 15)
# Classical_markers_Neuronal
# dev.off()
# 
# #Schwann cells
# markers.to.plot <- c("NRXN1","CCN3","MPZ","PTN","S100B")
# Classical_markers_Schwann <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Schwann cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Schwann.pdf", width = 15, height = 15)
# Classical_markers_Schwann
# dev.off()
# 
# #Adipocytes
# markers.to.plot <- c("ADIPOQ","PLIN1","FABP4","LEP","GPD1","CD36")
# Classical_markers_Adipocytes <- DotPlot(seuratObject_Sample, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
#   RotatedAxis() + labs(title="Adipocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))
# 
# 
# pdf("D:/Jojie/Analysis_Outputs/scRNA-11072024/071724_Classical_markers_Adipocytes.pdf", width = 15, height = 15)
# Classical_markers_Adipocytes
# dev.off()
