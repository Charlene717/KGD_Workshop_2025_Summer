## Ref: https://github.com/chris-mcginnis-ucsf/DoubletFinder
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)

##Remove doublet cells by DoubletFinder package
#HC2901
#Run general flow of scRNA-seq by Seurat package
HC2901 <- NormalizeData(object = HC2901)
HC2901 <- FindVariableFeatures(object = HC2901)
HC2901 <- ScaleData(object = HC2901)
HC2901 <- RunPCA(object = HC2901)
ElbowPlot(HC2901)
HC2901 <- FindNeighbors(object = HC2901, dims = 1:50)
HC2901 <- FindClusters(object = HC2901)
HC2901 <- RunUMAP(object = HC2901, dims = 1:30)
DimPlot(HC2901, reduction = "umap", label = T)

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(HC2901, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- HC2901@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(HC2901@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
HC2901 <- doubletFinder_v3(HC2901, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(HC2901)
DimPlot(HC2901, reduction = 'umap', group.by = "DF.classifications_0.25_0.25_224")
table(HC2901@meta.data$DF.classifications_0.25_0.25_224)

#Remove doublet cells
HC2901 <- subset(HC2901, subset = DF.classifications_0.25_0.25_224 == "Singlet")
DimPlot(HC2901, reduction = "umap", group.by="DF.classifications_0.25_0.25_224")
HC2901



#HC3080
#Run general flow of scRNA-seq by Seurat package
HC3080 <- NormalizeData(object = HC3080)
HC3080 <- FindVariableFeatures(object = HC3080)
HC3080 <- ScaleData(object = HC3080)
HC3080 <- RunPCA(object = HC3080)
ElbowPlot(HC3080)
HC3080 <- FindNeighbors(object = HC3080, dims = 1:50)
HC3080 <- FindClusters(object = HC3080)
HC3080 <- RunUMAP(object = HC3080, dims = 1:30)
DimPlot(HC3080, reduction = "umap", label = T)

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(HC3080, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- HC3080@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(HC3080@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
HC3080 <- doubletFinder_v3(HC3080, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(HC3080)
DimPlot(HC3080, reduction = 'umap', group.by = "DF.classifications_0.25_0.26_101")
table(HC3080@meta.data$DF.classifications_0.25_0.26_101)

#Remove doublet cells
HC3080 <- subset(HC3080, subset = DF.classifications_0.25_0.26_101 == "Singlet")
DimPlot(HC3080, reduction = "umap", group.by="DF.classifications_0.25_0.26_101")
HC3080


#DSAP3116
#Run general flow of scRNA-seq by Seurat package
DSAP3116 <- NormalizeData(object = DSAP3116)
DSAP3116 <- FindVariableFeatures(object = DSAP3116)
DSAP3116 <- ScaleData(object = DSAP3116)
DSAP3116 <- RunPCA(object = DSAP3116)
ElbowPlot(DSAP3116)
DSAP3116 <- FindNeighbors(object = DSAP3116, dims = 1:50)
DSAP3116 <- FindClusters(object = DSAP3116)
DSAP3116 <- RunUMAP(object = DSAP3116, dims = 1:30)
DimPlot(DSAP3116, reduction = "umap", label = T)

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(DSAP3116, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- DSAP3116@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(DSAP3116@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
DSAP3116 <- doubletFinder_v3(DSAP3116, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(DSAP3116)
DimPlot(DSAP3116, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_48")
table(DSAP3116@meta.data$DF.classifications_0.25_0.3_48)

#Remove doublet cells
DSAP3116 <- subset(DSAP3116, subset = DF.classifications_0.25_0.3_48 == "Singlet")
DimPlot(DSAP3116, reduction = "umap", group.by="DF.classifications_0.25_0.3_48")
DSAP3116


#DSAP3138
#Run general flow of scRNA-seq by Seurat package
DSAP3138 <- NormalizeData(object = DSAP3138)
DSAP3138 <- FindVariableFeatures(object = DSAP3138)
DSAP3138 <- ScaleData(object = DSAP3138)
DSAP3138 <- RunPCA(object = DSAP3138)
ElbowPlot(DSAP3138)
DSAP3138 <- FindNeighbors(object = DSAP3138, dims = 1:50)
DSAP3138 <- FindClusters(object = DSAP3138)
DSAP3138 <- RunUMAP(object = DSAP3138, dims = 1:30)
DimPlot(DSAP3138, reduction = "umap", label = T)

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(DSAP3138, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- DSAP3138@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(DSAP3138@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
DSAP3138 <- doubletFinder_v3(DSAP3138, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(DSAP3138)
DimPlot(DSAP3138, reduction = 'umap', group.by = "DF.classifications_0.25_0.03_83")
table(DSAP3138@meta.data$DF.classifications_0.25_0.03_83)

#Remove doublet cells
DSAP3138 <- subset(DSAP3138, subset = DF.classifications_0.25_0.03_83 == "Singlet")
DimPlot(DSAP3138, reduction = "umap", group.by="DF.classifications_0.25_0.03_83")
DSAP3138
