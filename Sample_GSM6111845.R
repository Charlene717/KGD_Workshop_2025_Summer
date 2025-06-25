library(Seurat)
library(dplyr)
library(ggplot2)

data_GSM6111845 <- Read10X("C:/Charlene/Dataset_KGD_Lab/scRNA-seq/10x/sample_filtered_feature_bc_matrix/GSM6111845_Normal_Sole/")

seurat_GSM6111845 <- CreateSeuratObject(counts = data_GSM6111845, project = "GSM6111845",
                                  min.cells = 3, min.features = 200)

# Calculate mitochondrial gene percentage and perform quality control
seurat_GSM6111845[["percent.mt"]] <- PercentageFeatureSet(seurat_GSM6111845, pattern = "^MT-")
seurat_GSM6111845 <- subset(seurat_GSM6111845, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)



library(DoubletFinder)

# Standard Seurat preprocessing steps
seurat_GSM6111845 <- NormalizeData(seurat_GSM6111845, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_GSM6111845 <- FindVariableFeatures(seurat_GSM6111845, selection.method = "vst", nfeatures = 2000)
seurat_GSM6111845 <- ScaleData(seurat_GSM6111845, features = rownames(seurat_GSM6111845), verbose = TRUE)
seurat_GSM6111845 <- RunPCA(seurat_GSM6111845)
ElbowPlot(seurat_GSM6111845)
seurat_GSM6111845 <- FindNeighbors(seurat_GSM6111845, dims = 1:30)
seurat_GSM6111845 <- FindClusters(seurat_GSM6111845)
seurat_GSM6111845 <- RunUMAP(seurat_GSM6111845, dims = 1:30)
DimPlot(seurat_GSM6111845, reduction = "umap", label = T)

# Determine optimal pK
sweep_res <- paramSweep(seurat_GSM6111845, PCs = 1:30)
sweep_stats <- summarizeSweep(sweep_res)
best_pk <- find.pK(sweep_stats)
pK <- as.numeric(as.character(best_pk[which.max(best_pk$BCmetric), "pK"]))

# Estimate expected doublet count
annotations <- seurat_GSM6111845$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp <- round(0.08 * nrow(seurat_GSM6111845@meta.data) * (1 - homotypic.prop))

# Predict and remove doublets
seurat_GSM6111845 <- doubletFinder(seurat_GSM6111845, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp) # Run DoubletFinder with specified parameters
df_col <- grep("^DF\\.classifications", colnames(seurat_GSM6111845@meta.data), value = TRUE)[1] # Identify the metadata column that starts with "DF.classifications"
DimPlot(seurat_GSM6111845, reduction = 'umap', group.by = df_col) # Visualize UMAP colored by doublet classification
table(seurat_GSM6111845@meta.data[[df_col]]) # Display count of Singlet vs Doublet classifications
seurat_GSM6111845 <- seurat_GSM6111845[, seurat_GSM6111845@meta.data[[df_col]] == "Singlet"] # Subset the object to retain only singlet cells
DimPlot(seurat_GSM6111845, reduction = "umap", group.by = df_col) # Re-visualize UMAP after singlet filtering
seurat_GSM6111845 # View the resulting Seurat object
