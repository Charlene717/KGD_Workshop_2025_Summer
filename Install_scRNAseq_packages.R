### 📦 安裝與載入 scRNA-seq Workflow 所需套件 ###
### Version: 2025-06 ###

# CRAN 套件
if(!require('Seurat'))         { install.packages('Seurat');         library(Seurat) }
if(!require('dplyr'))          { install.packages('dplyr');          library(dplyr) }
if(!require('ggplot2'))        { install.packages('ggplot2');        library(ggplot2) }

# 為 tidyverse 使用者提供一體化選項（可選）
if(!require('tidyverse'))      { install.packages('tidyverse');      library(tidyverse) }

# Bioconductor 安裝器
if(!require('BiocManager'))    { install.packages('BiocManager');    library(BiocManager) }

# Bioconductor 套件
if(!require('SingleR'))        { BiocManager::install('SingleR');    library(SingleR) }
if(!require('celldex'))        { BiocManager::install('celldex');    library(celldex) }
if(!require('monocle'))        { BiocManager::install('monocle');    library(monocle) }

# 其他CRAN/外部套件
if(!require('DoubletFinder'))  { install.packages('DoubletFinder');  library(DoubletFinder) }
if(!require('CellChat'))       { install.packages('CellChat');       library(CellChat) }

# ✅ CellChat 資料庫 (human)
if(!"CellChatDB.human" %in% ls("package:CellChat")) {
  CellChat::CellChatDB <- CellChat::CellChatDB.human
  message("✅ CellChatDB.human 已載入")
}

### 📍 所有套件已載入完成 ###
