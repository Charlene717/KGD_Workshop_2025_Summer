# KGD_Workshop_2025_Summer


# 🧬 scRNA-seq Workshop 實作範例流程

本實作流程涵蓋單細胞轉錄體分析的整體流程，包含資料前處理、雙細胞去除、資料整合、細胞註解、細胞通訊與發展軌跡分析等，皆以簡潔 R 語法示範。

---

## 📦 1. 資料讀取與前處理（Seurat）

```r
library(Seurat)
library(dplyr)
library(ggplot2)

# 讀取資料並建立 Seurat 物件
data_2901 <- Read10X("2901/")
seurat_2901 <- CreateSeuratObject(counts = data_2901, min.cells = 3, min.features = 200)

# 計算粒線體基因比例與品質控制
seurat_2901[["percent.mt"]] <- PercentageFeatureSet(seurat_2901, pattern = "^MT-")
seurat_2901 <- subset(seurat_2901, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
```

---

## ❌ 2. 雙細胞過濾（DoubletFinder）

```r
library(DoubletFinder)

# Seurat 常規前處理流程
seurat_2901 <- NormalizeData(seurat_2901) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
seurat_2901 <- FindNeighbors(seurat_2901, dims = 1:20)
seurat_2901 <- FindClusters(seurat_2901)
seurat_2901 <- RunUMAP(seurat_2901, dims = 1:20)

# 找出最佳 pK
sweep_res <- paramSweep(seurat_2901, PCs = 1:20)
sweep_stats <- summarizeSweep(sweep_res)
best_pk <- find.pK(sweep_stats)
pK <- as.numeric(as.character(best_pk[which.max(best_pk$BCmetric), "pK"]))

# 計算期望雙細胞數
annotations <- seurat_2901$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp <- round(0.08 * nrow(seurat_2901@meta.data) * (1 - homotypic.prop))

# 雙細胞預測並過濾
seurat_2901 <- doubletFinder(seurat_2901, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp)
seurat_2901 <- subset(seurat_2901, subset = DF.classifications == "Singlet")
```

---

## 🔗 3. 多樣本資料整合

```r
# 合併多個樣本
seurat_all <- merge(seurat_2901, y = list(seurat_3080, seurat_3116, seurat_3138),
                    add.cell.ids = c("HC2901","HC3080","DSAP3116","DSAP3138"))

seurat_list <- SplitObject(seurat_all, split.by = "orig.ident")

# 各樣本正規化與細胞週期
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData()
  CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
})

# 資料整合
features <- SelectIntegrationFeatures(seurat_list)
anchors <- FindIntegrationAnchors(seurat_list, dims = 1:30)
integrated <- IntegrateData(anchors, dims = 1:30)

# 降維與分群
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
integrated <- RunPCA(integrated, npcs = 50) %>% RunUMAP(dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30) %>% FindClusters(resolution = 0.3)
```

---

## 🧾 4. 細胞註解（SingleR + 標記基因）

```r
library(SingleR); library(celldex)
hpca <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(integrated, slot = "data")

pred <- SingleR(test = counts, ref = hpca, labels = hpca$label.main)

# 查看結果
table(pred$pruned.labels)

# DotPlot 顯示標記基因
DefaultAssay(integrated) <- "RNA"
DotPlot(integrated, features = c("KRT14", "CD3D", "PECAM1")) + RotatedAxis()
```

---

## 🔁 5. 細胞通訊分析（CellChat）

```r
library(CellChat)

data.input <- GetAssayData(integrated, assay = "RNA", slot = "data")
meta <- data.frame(group = Idents(integrated), row.names = colnames(integrated))
cellchat <- createCellChat(data.input, meta = meta, group.by = "group")
cellchat@DB <- CellChatDB.human

# 分析流程
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 繪製通訊圈圖
netVisual_circle(cellchat@net$count)
```

---

## ⏳ 6. 發展軌跡分析（Monocle2）

```r
library(monocle)

# Seurat 轉 Monocle
data <- as(as.matrix(integrated@assays$RNA@data), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", data = integrated@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(data)))
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds) %>% estimateDispersions()
cds <- reduceDimension(cds, method = "DDRTree") %>% orderCells()
plot_cell_trajectory(cds, color_by = "seurat_clusters")
```
