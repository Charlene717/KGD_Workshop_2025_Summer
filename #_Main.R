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
source("Sample_GSM6111844.R")
source("Sample_GSM6111845.R")
source("Sample_GSM6111847.R")



#### Data Integration ####
source("RUN_Integration.R")


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


