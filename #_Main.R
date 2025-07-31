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

## KGD Lab 套件版本
source("Install_required_packages_KGD_Lab.R")
## 自動安裝最新版本套件
source("Install_required_packages.R")


## Speed up
if(!require('future')) install.packages('future'); library(future)
## https://github.com/immunogenomics/presto
if(!require('presto')) devtools::install_github("immunogenomics/presto"); library(presto) # Speeds up FindAllMarkers
plan(multicore, workers = 20)
options(future.globals.maxSize = 2048*100 * 1024^2) # Set memory limit to ~204.8 GB



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
source("RUN_CellType_Annotation_SingleR.R")

# DotPlot to display marker genes
DefaultAssay(seurat_all_integrated) <- "RNA"
DotPlot(seurat_all_integrated, features = c("KRT14", "CD3D", "PECAM1")) + RotatedAxis()

source("KGD_CTAnnot_MarkerList.R")


#### Differential Expression Gene (DEG) Analysis ####
source("RUN_DEG.R")


#### Enrichment Analysis ####
source("RUN_Enrichment_Analysis_ORA.R")
source("RUN_Enrichment_Analysis_GSEA.R")


#### CellChat Analysis ####
source("RUN_CellChat.R")


#### Trajectory Analysis ####
source("RUN_Trajectory_Analysis_Monocle3.R")



#### Export ####
# 匯出 Seurat 物件 (RDS) — 將物件儲存為 .rds 供後續流程載入
saveRDS(seurat_all_integrated, file = paste0(Name_ExportFolder, "/", Name_Export, ".rds"))  # 單一 Seurat 物件；適合日後 readRDS() 載入


# 匯出完整 R 工作空間 (RData) — 若需重現所有物件、變數
save.image(paste0(Name_ExportFolder, "/", Name_Export, ".RData"))  # 儲存 .RData 以保留所有記憶體物件

# 紀錄當前 Session 資訊 (套件版本、平台等)
writeLines(capture.output(sessionInfo()), paste0(Name_ExportFolder, "/", Name_Export, "_session_info.txt"))  # 輸出 sessionInfo 至文字檔
sessionInfo()  # 亦於 Console 顯示，方便即時檢查


#### Export loupeR ####
Set_Run_loupeR <- TRUE
if(Set_Run_loupeR){
  if (!require('loupeR')) {remotes::install_github("10xGenomics/loupeR"); library(loupeR)}
  loupeR::setup()
  
  try({seurat_all_integrated$CELL <- NULL})
  create_loupe_from_seurat(seurat_all_integrated, 
                           output_name = paste0(Name_ExportFolder,"/", Name_Export,""))
  
}
