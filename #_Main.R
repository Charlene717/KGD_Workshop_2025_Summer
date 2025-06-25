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


#### Data Loading and Preprocessing ####
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


