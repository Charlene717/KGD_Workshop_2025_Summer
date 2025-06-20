library(Seurat)
library(dplyr)
library(ggplot2)

data_2901 <- Read10X("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/KGD_Workshop_2025_Summer/Input_dataset/TN125_Keloid/")

seurat_2901 <- CreateSeuratObject(counts = data_2901, project = "HC2901",
                                  min.cells = 3, min.features = 200)

# Calculate mitochondrial gene percentage and perform quality control
seurat_2901[["percent.mt"]] <- PercentageFeatureSet(seurat_2901, pattern = "^MT-")
seurat_2901 <- subset(seurat_2901, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)



library(DoubletFinder)

# Standard Seurat preprocessing steps
seurat_2901 <- NormalizeData(seurat_2901, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_2901 <- FindVariableFeatures(seurat_2901, selection.method = "vst", nfeatures = 2000)
seurat_2901 <- ScaleData(seurat_2901, features = rownames(seurat_2901), verbose = TRUE)
seurat_2901 <- RunPCA(seurat_2901)
ElbowPlot(seurat_2901)
seurat_2901 <- FindNeighbors(seurat_2901, dims = 1:30)
seurat_2901 <- FindClusters(seurat_2901)
seurat_2901 <- RunUMAP(seurat_2901, dims = 1:30)
DimPlot(seurat_2901, reduction = "umap", label = T)

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
seurat_2901 <- doubletFinder(seurat_2901, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp) # Run DoubletFinder with specified parameters
df_col <- grep("^DF\\.classifications", colnames(seurat_2901@meta.data), value = TRUE)[1] # Identify the metadata column that starts with "DF.classifications"
DimPlot(seurat_2901, reduction = 'umap', group.by = df_col) # Visualize UMAP colored by doublet classification
table(seurat_2901@meta.data[[df_col]]) # Display count of Singlet vs Doublet classifications
seurat_2901 <- seurat_2901[, seurat_2901@meta.data[[df_col]] == "Singlet"] # Subset the object to retain only singlet cells
DimPlot(seurat_2901, reduction = "umap", group.by = df_col) # Re-visualize UMAP after singlet filtering
seurat_2901 # View the resulting Seurat object
