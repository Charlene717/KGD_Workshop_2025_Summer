library(Seurat)
library(dplyr)
library(ggplot2)

data_3080 <- Read10X("C:/Charlene/Dataset_KGD_Lab/scRNA-seq/10x/sample_filtered_feature_bc_matrix/GSM6111845_Normal_Sole/")

seurat_3080 <- CreateSeuratObject(counts = data_3080, project = "HC3080",
                                  min.cells = 3, min.features = 200)

# Calculate mitochondrial gene percentage and perform quality control
seurat_3080[["percent.mt"]] <- PercentageFeatureSet(seurat_3080, pattern = "^MT-")
seurat_3080 <- subset(seurat_3080, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)



library(DoubletFinder)

# Standard Seurat preprocessing steps
seurat_3080 <- NormalizeData(seurat_3080, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_3080 <- FindVariableFeatures(seurat_3080, selection.method = "vst", nfeatures = 2000)
seurat_3080 <- ScaleData(seurat_3080, features = rownames(seurat_3080), verbose = TRUE)
seurat_3080 <- RunPCA(seurat_3080)
ElbowPlot(seurat_3080)
seurat_3080 <- FindNeighbors(seurat_3080, dims = 1:30)
seurat_3080 <- FindClusters(seurat_3080)
seurat_3080 <- RunUMAP(seurat_3080, dims = 1:30)
DimPlot(seurat_3080, reduction = "umap", label = T)

# Determine optimal pK
sweep_res <- paramSweep(seurat_3080, PCs = 1:30)
sweep_stats <- summarizeSweep(sweep_res)
best_pk <- find.pK(sweep_stats)
pK <- as.numeric(as.character(best_pk[which.max(best_pk$BCmetric), "pK"]))

# Estimate expected doublet count
annotations <- seurat_3080$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp <- round(0.08 * nrow(seurat_3080@meta.data) * (1 - homotypic.prop))

# Predict and remove doublets
seurat_3080 <- doubletFinder(seurat_3080, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp) # Run DoubletFinder with specified parameters
df_col <- grep("^DF\\.classifications", colnames(seurat_3080@meta.data), value = TRUE)[1] # Identify the metadata column that starts with "DF.classifications"
DimPlot(seurat_3080, reduction = 'umap', group.by = df_col) # Visualize UMAP colored by doublet classification
table(seurat_3080@meta.data[[df_col]]) # Display count of Singlet vs Doublet classifications
seurat_3080 <- seurat_3080[, seurat_3080@meta.data[[df_col]] == "Singlet"] # Subset the object to retain only singlet cells
DimPlot(seurat_3080, reduction = "umap", group.by = df_col) # Re-visualize UMAP after singlet filtering
seurat_3080 # View the resulting Seurat object
