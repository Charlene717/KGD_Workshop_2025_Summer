
# 🧬 scRNA-seq Workshop Practical Workflow

This practical workflow demonstrates key steps in single-cell RNA sequencing (scRNA-seq) analysis, including data preprocessing, doublet removal, sample integration, cell annotation, cell-cell communication, and pseudotime analysis. The R code examples are simplified for hands-on training.

---

## 📦 1. Data Loading and Preprocessing (Seurat)

```r
library(Seurat)
library(dplyr)
library(ggplot2)

# Read data and create Seurat object
data_2901 <- Read10X("2901/")
seurat_2901 <- CreateSeuratObject(counts = data_2901, min.cells = 3, min.features = 200)

# Calculate mitochondrial gene percentage and perform quality control
seurat_2901[["percent.mt"]] <- PercentageFeatureSet(seurat_2901, pattern = "^MT-")
seurat_2901 <- subset(seurat_2901, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
```

---

## ❌ 2. Doublet Removal (DoubletFinder)

```r
library(DoubletFinder)

# Standard Seurat preprocessing steps
seurat_2901 <- NormalizeData(seurat_2901) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
seurat_2901 <- FindNeighbors(seurat_2901, dims = 1:30)
seurat_2901 <- FindClusters(seurat_2901)
seurat_2901 <- RunUMAP(seurat_2901, dims = 1:30)

# Determine optimal pK
sweep_res <- paramSweep(seurat_2901, PCs = 1:30)
sweep_stats <- summarizeSweep(sweep_res)
best_pk <- find.pK(sweep_stats)
pK <- as.numeric(as.character(best_pk[which.max(best_pk$BCmetric), "pK"]))

# Estimate expected doublet count
annotations <- seurat_2901$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp <- round(0.08 * nrow(seurat_2901@meta.data) * (1 - homotypic.prop))

# Predict and remove doublets
seurat_2901 <- doubletFinder(seurat_2901, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp)
seurat_2901 <- subset(seurat_2901, subset = DF.classifications == "Singlet")
```

---

## 🔗 3. Integration of Multiple Samples

```r
# Merge multiple samples
seurat_all <- merge(seurat_2901, y = list(seurat_3080, seurat_3116, seurat_3138),
                    add.cell.ids = c("HC2901","HC3080","DSAP3116","DSAP3138"))

seurat_list <- SplitObject(seurat_all, split.by = "orig.ident")

# Normalize and score cell cycle for each sample
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData()
  CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
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
```

---

## 🧾 4. Cell Type Annotation (SingleR + Marker Genes)

```r
library(SingleR); library(celldex)
hpca <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(integrated, slot = "data")

pred <- SingleR(test = counts, ref = hpca, labels = hpca$label.main)

# View results
table(pred$pruned.labels)

# DotPlot to display marker genes
DefaultAssay(integrated) <- "RNA"
DotPlot(integrated, features = c("KRT14", "CD3D", "PECAM1")) + RotatedAxis()
```

---

## 🔁 5. Cell-Cell Communication Analysis (CellChat)

```r
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
```

---

## ⏳ 6. Pseudotime Analysis (Monocle2)

```r
library(monocle)

# Convert Seurat object to Monocle format
data <- as(as.matrix(integrated@assays$RNA@data), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", data = integrated@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(data)))
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds) %>% estimateDispersions()
cds <- reduceDimension(cds, method = "DDRTree") %>% orderCells()
plot_cell_trajectory(cds, color_by = "seurat_clusters")
```
