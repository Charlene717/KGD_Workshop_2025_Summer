

#### Merge multiple samples ####
seurat_all_merge <- merge(seurat_GSM6111844, y = list(seurat_GSM6111845,seurat_GSM6111847),
                    add.cell.ids = c("GSM6111844","GSM6111845","GSM6111847"))

# 合併層（僅限 Seurat ≥ 5）
if (packageVersion("Seurat") >= "5.0.0") {
  seurat_all_merge <- JoinLayers(seurat_all_merge) # 在 merge 後統一執行 JoinLayers
}



Run_Demo_Merge <- TRUE
if(Run_Demo_Merge){
  seurat_all_merge <- NormalizeData(seurat_all_merge,normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_all_merge <- FindVariableFeatures(seurat_all_merge, selection.method = "vst", nfeatures = 2000)
  seurat_all_merge <- ScaleData(seurat_all_merge, features = rownames(seurat_all_merge), verbose = TRUE)
  seurat_all_merge <- CellCycleScoring(seurat_all_merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  
  seurat_all_merge <- ScaleData(seurat_all_merge, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
  seurat_all_merge <- RunPCA(seurat_all_merge, npcs = 50) %>% RunUMAP(dims = 1:30)
  seurat_all_merge <- FindNeighbors(seurat_all_merge, dims = 1:30) %>% FindClusters(resolution = 0.3)
  DimPlot(seurat_all_merge, reduction = "umap") +
    DimPlot(seurat_all_merge, reduction = "umap", group.by = "orig.ident")
  
}



#### Data Integration ####
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



features <- SelectIntegrationFeatures(seurat_list)
anchors <- FindIntegrationAnchors(seurat_list, dims = 1:30)
seurat_all_integrated <- IntegrateData(anchors, dims = 1:30)

# Dimensionality reduction and clustering
DefaultAssay(seurat_all_integrated) <- "integrated"
seurat_all_integrated <- ScaleData(seurat_all_integrated, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
seurat_all_integrated <- RunPCA(seurat_all_integrated, npcs = 50) %>% RunUMAP(dims = 1:30)
seurat_all_integrated <- FindNeighbors(seurat_all_integrated, dims = 1:30) %>% FindClusters(resolution = 0.3)

DimPlot(seurat_all_integrated, reduction = "umap") +
  DimPlot(seurat_all_integrated, reduction = "umap", group.by = "orig.ident")
