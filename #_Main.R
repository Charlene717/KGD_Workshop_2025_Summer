##### Presetting #####
rm(list = ls()) # 清空變數
memory.limit(150000)

#### Load Packages ####
if(!require('Seurat'))         { install.packages('Seurat');         library(Seurat) }
if(!require('tidyverse'))      { install.packages('tidyverse');      library(tidyverse) }
if(!require('dplyr'))          { install.packages('dplyr');          library(dplyr) }
if(!require('ggplot2'))        { install.packages('ggplot2');        library(ggplot2) }

library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)

#### Set Parameter ####
## Set Export 
# Generate unique export parameters
Set_Project <- "Demo"

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate a unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

## Construct Set_note
Set_note <- paste0(Name_FileID, "_", Set_Project)

Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_",Set_note)
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}



#### Data Loading and Preprocessing ####
source("Sample_GSM6111845.R")


source("Sample_2901.R")
source("Sample_3080.R")


#### Merge multiple samples ####
seurat_all <- merge(seurat_2901, y = list(seurat_3080),
                    add.cell.ids = c("HC2901","HC3080"))

seurat_all <- JoinLayers(seurat_all)  # 🔥 在 merge 後統一執行 JoinLayers

seurat_list <- SplitObject(seurat_all, split.by = "orig.ident")


# cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Normalize and score cell cycle for each sample
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = TRUE)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
})


# Data integration
features <- SelectIntegrationFeatures(seurat_list)
anchors <- FindIntegrationAnchors(seurat_list, dims = 1:30)
integrated <- IntegrateData(anchors, dims = 1:30)

# Dimensionality reduction and clustering
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
integrated <- RunPCA(integrated, npcs = 50) %>% RunUMAP(dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30) %>% FindClusters(resolution = 0.3)

DimPlot(integrated, reduction = "umap")


#### Cell Type Annotation ####
library(SingleR); library(celldex)
hpca <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(integrated, slot = "data")

pred <- SingleR(test = counts, ref = hpca, labels = hpca$label.main)

# View results
table(pred$pruned.labels)

# DotPlot to display marker genes
DefaultAssay(integrated) <- "RNA"
DotPlot(integrated, features = c("KRT14", "CD3D", "PECAM1")) + RotatedAxis()













# --------------------------
# 7. 物件儲存
# --------------------------
# ※請先於腳本上方或執行前設定 Name_ExportFolder / Name_Export 變數※
# Name_ExportFolder <- "C:/Charlene/Output"  # 匯出資料夾路徑
# Name_Export       <- "GSM6111845_preprocessed"  # 檔名前綴


# 7A. 匯出 Seurat 物件 (RDS) — 將物件儲存為 .rds 供後續流程載入
saveRDS(seurat_GSM6111845, file = paste0(Name_ExportFolder, "/", Name_Export, ".rds"))  # 單一 Seurat 物件；適合日後 readRDS() 載入


# 7B. 匯出完整 R 工作空間 (RData) — 若需重現所有物件、變數
save.image(paste0(Name_ExportFolder, "/", Name_Export, ".RData"))  # 儲存 .RData 以保留所有記憶體物件

# 7C. 紀錄當前 Session 資訊 (套件版本、平台等)
writeLines(capture.output(sessionInfo()), paste0(Name_ExportFolder, "/", Name_Export, "_session_info.txt"))  # 輸出 sessionInfo 至文字檔
sessionInfo()  # 亦於 Console 顯示，方便即時檢查


