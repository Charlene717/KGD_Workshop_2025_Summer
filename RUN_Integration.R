#### Data Integration ####

## Merge multiple samples
seurat_all_merge <- merge(seurat_GSM6111844, y = list(seurat_GSM6111845,seurat_GSM6111847),
                    add.cell.ids = c("GSM6111844","GSM6111845","GSM6111847"))

seurat_all_merge <- JoinLayers(seurat_all_merge)  # 🔥 在 merge 後統一執行 JoinLayers

seurat_list <- SplitObject(seurat_all_merge, split.by = "orig.ident")


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
seurat_all_integrated <- IntegrateData(anchors, dims = 1:30)

# Dimensionality reduction and clustering
DefaultAssay(seurat_all_integrated) <- "integrated"
seurat_all_integrated <- ScaleData(seurat_all_integrated, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
seurat_all_integrated <- RunPCA(seurat_all_integrated, npcs = 50) %>% RunUMAP(dims = 1:30)
seurat_all_integrated <- FindNeighbors(seurat_all_integrated, dims = 1:30) %>% FindClusters(resolution = 0.3)

DimPlot(seurat_all_integrated, reduction = "umap")
