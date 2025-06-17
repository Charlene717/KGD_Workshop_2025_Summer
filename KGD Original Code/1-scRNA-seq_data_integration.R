install.packages("Seurat")
#library(remotes)
remotes::install_version("Seurat", version="5.0.3")
remotes::install_version("SeuratObject", version="5.0.1")

#remove.packages("Seurat")
#remove.packages("SeuratObject")

Sys.setenv(LANGUAGE='en')

#Library
library(Seurat)
library(dplyr)
library(ggsci)
library(Matrix)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(clusterProfiler)
library(gplots)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(tidyr) 

##Load the dataset
#HC2901#Healthy skin
HC2901.data <- Read10X(data.dir = "2901/")
HC2901 <- CreateSeuratObject(counts = HC2901.data, project = "2901", min.cells = 3, min.features = 10)
HC2901 ##16067 features across 3149 samples within 1 assay 

#HC3080#Healthy skin
HC3080.data <- Read10X(data.dir = "3080/")
HC3080 <- CreateSeuratObject(counts = HC3080.data, project = "3080", min.cells = 3, min.features = 10)
HC3080 ##14409 features across 1363 samples within 1 assay 

#DSAP3116#DSAP
DSAP3116.data <- Read10X(data.dir = "3116/")
DSAP3116 <- CreateSeuratObject(counts = DSAP3116.data, project = "3116", min.cells = 3, min.features = 10)
DSAP3116 ##13093 features across 677 samples within 1 assay 

#DSAP3138#DSAP
DSAP3138.data <- Read10X(data.dir = "3138/")
DSAP3138 <- CreateSeuratObject(counts = DSAP3138.data, project = "3138", min.cells = 3, min.features = 10)
DSAP3138 ##14023 features across 1188 samples within 1 assay 

##Re-check
#HC2901 Healthy 
HC2901 #16067 genes x 3149 samples

#HC3080 Healthy 
HC3080 #14409 genes x 1363 samples

#DSAP3116 DSAP
DSAP3116 #13093 genes x 677 samples 

#DSAP3138 DSAP
DSAP3138 #14023 genes x 1188 samples

##Quality control
#HC2901
#The $ operator can add columns to object metadata.
HC2901$orig.ident1 <- "1_Healthy_skin"
HC2901$orig.ident2 <- "HC2901"
sample <- "Healthy_skin_2901"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
HC2901[["percent.mt"]] <- PercentageFeatureSet(HC2901, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(HC2901@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(HC2901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(paste0("../202404/", sample,"_.pdf"), width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
HC2901 <- subset(HC2901, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#HC3080
#The $ operator can add columns to object metadata.
HC3080$orig.ident1 <- "1_Healthy_skin"
HC3080$orig.ident2 <- "HC3080"
sample <- "Healthy_skin_3080"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
HC3080[["percent.mt"]] <- PercentageFeatureSet(HC3080, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(HC3080@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(HC3080, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(paste0("../202404/", sample,"_.pdf"), width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
HC3080 <- subset(HC3080, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#DSAP3116
#The $ operator can add columns to object metadata.
DSAP3116$orig.ident1 <- "2_DSAP_skin"
DSAP3116$orig.ident2 <- "DSAP3116"
sample <- "DSAP_skin_3116"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
DSAP3116[["percent.mt"]] <- PercentageFeatureSet(DSAP3116, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(DSAP3116@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(DSAP3116, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(paste0("../202404/", sample,"_.pdf"), width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
DSAP3116 <- subset(DSAP3116, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#DSAP3138
#The $ operator can add columns to object metadata.
DSAP3138$orig.ident1 <- "3_DSAP_skin_abrocitinib"
DSAP3138$orig.ident2 <- "DSAP3138"
sample <- "DSAP_skin_abrocitinib_3138"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
DSAP3138[["percent.mt"]] <- PercentageFeatureSet(DSAP3138, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(DSAP3138@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(DSAP3138, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(paste0("../202404/", sample,"_.pdf"), width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
DSAP3138 <- subset(DSAP3138, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

##Re-check (after Doublet remove)
HC2901  ##16067 genes x 2914 samples
HC3080  ##14409 genes x 1259 samples
DSAP3116  ##13093 genes x 624 samples
DSAP3138  ##14023 genes x 1098 samples

##scRNA-seq sample integration
#Set merge file
TN.C <- merge (HC2901, y = c(HC3080, DSAP3116, DSAP3138), 
               add.cell.ids = c("HC2901","HC3080","DSAP3116","DSAP3138"), project = "groups")
head(TN.C[[]])
table(TN.C$orig.ident)
#Set splitobject file
TN.list<- SplitObject(TN.C, split.by = "orig.ident")
TN.list

# cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# NormalizeData,FindVariableFeatures, ScaleData, CellCycleScoring
TN.list <- lapply(X = TN.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = TRUE)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
})

#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = TN.list)
#Perform integration 
#reference.list <- TN.list[c("HC2901","HC3080","HC3106","DF2795","DF3019","DF3048")]
TN.anchors <- FindIntegrationAnchors(object.list = TN.list, dims = 1:30)
saveRDS(TN.anchors, file = "C:/Users/user/Desktop/202404/TN.anchorsdim30.rds")

#this command creates an 'integrated' data assay
TN.combined <- IntegrateData(anchorset = TN.anchors, dims = 1:30)
#Perform an integrated analysis
#specify that we will perform downstream analysis on the corrected data note that the
#original unmodified data still resides in the 'RNA' assay
DefaultAssay(TN.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
#TN.combined <- ScaleData(TN.combined, verbose = TRUE)
TN.combined <- ScaleData(TN.combined, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(TN.combined), verbose = TRUE)
TN.combined <- RunPCA(TN.combined, npcs = 50, verbose = TRUE)
VizDimLoadings(TN.combined, dims = 1:2)
ElbowPlot(TN.combined, ndims= 50)
# Cluster the cells (1:30)
TN.combined <- FindNeighbors(TN.combined, reduction = "pca", dims = 1:30)
TN.combined <- FindClusters(TN.combined, resolution = 0.01)
table(TN.combined@active.ident)
table(Idents(TN.combined), TN.combined$orig.ident1)

#CellNumber
CellNumber <- table(Idents(TN.combined), TN.combined$orig.ident1)
write.csv(CellNumber , file = "C:/Users/user/Desktop/202404/CellNumber.csv") 

#TN.combined <- RunUMAP(TN.combined, reduction="pca", dims = 1:30)
#DimPlot(TN.combined,reduction = "umap",group.by = "orig.ident")
##UMAP (1:30)
TN.combined <- RunUMAP(TN.combined, reduction = "pca", dims = 1:30)
DimPlot(TN.combined, reduction = "umap",label = FALSE, pt.size = 0.8)
DimPlot(TN.combined, reduction = "umap",label = TRUE, pt.size = 0.8)
DimPlot(TN.combined, group.by = "orig.ident", pt.size = 0.8)
DimPlot(TN.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2)
DimPlot(TN.combined, reduction = "umap",label = TRUE, split.by = "TN.merge", pt.size = 0.8, ncol = 3)

##All cells cluster heatmap
Heatmapall <- subset(TN.combined, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
Heatmapall.markers <- FindAllMarkers(Heatmapall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(Heatmapall.markers)
write.csv(Heatmapall.markers , file = "C:/Users/user/Desktop/202404/Findallmarkers.csv", col.names = TRUE)
top10 <- Heatmapall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Heatmapall, features = top10$gene)

##Save the merged Seurat object
saveRDS(TN.combined, file = "C:/Users/user/Desktop/202404/TN.combined_dim30.rds")


#Cellproportion
Cellproportion <- table(Idents(TN.combined), TN.combined$orig.ident1)
Cellproportion <- round(sweep(Cellproportion,MARGIN=2, STATS=colSums(Cellproportion), FUN = "/")*100,2)
write.csv(Cellproportion, file = "C:/Users/user/Desktop/202404/Cellproportion.csv", row.names = T)
Cellproportion <- as.data.frame(Cellproportion)

nb.cols <- nrow(CellNumber)
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Paired"))(nb.cols)

proportion_plot <- ggplot(Cellproportion, aes(x = Var2, y = Freq, fill = Var1)) + theme_bw(base_size = 15) + 
  geom_col(position = "fill", width = 0.6) + xlab("Sample") + ylab("Proportion") + 
  theme(legend.title = element_blank())+ scale_fill_manual(values = mycolors)

pdf("C:/Users/user/Desktop/202404/proportion_plot.pdf", width = 15, height = 15)
proportion_plot
dev.off()

TN.combined <- readRDS(file = "C:/Users/user/Desktop/202404/TN.combined_dim30.rds")
##Rename cluster (for example_Fibroblast, Macrophage......)
new.cluster.ids <- c("C0_KC","C1_KC","C2_FB","C3_KC","C4_Tcell","C5_Mo", 
                     "C6_KC","C7_vEC","C8_MLA","C9_lEC","C10_SG","C11_SMC",
                     "C12_AT","C13_KC","C14_MC")
names(new.cluster.ids) <- levels(TN.combined)
TNname.combined <- RenameIdents(TN.combined, new.cluster.ids)
#DimPlot(TNname.combined, reduction = "umap",label = FALSE, pt.size = 0.8, cols = mycolors)
#DimPlot(TNname.combined, reduction = "umap",label = TRUE, pt.size = 0.8, repel = T, cols = mycolors)
#DimPlot(TNname.combined, group.by = "orig.ident1", pt.size = 0.8, cols = c("#00BA38","#00B9E3","#F8766D","#DB72FB"))
#DimPlot(TNname.combined, reduction = "umap",label = FALSE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, cols = mycolors)
UMAP_cluster_plot_all<- DimPlot(TNname.combined, reduction = "umap",label = TRUE, pt.size = 0.8, label.size=7, repel = T, cols = mycolors)
UMAP_plot_subset<- DimPlot(TNname.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, repel = T, cols = mycolors)

pdf("C:/Users/user/Desktop/202404/UMAP_cluster_plot_all.pdf", width = 15, height = 15)
UMAP_cluster_plot_all
dev.off()

pdf("C:/Users/user/Desktop/202404/UMAP_cluster_plot_subset.pdf", width = 15, height = 15)
UMAP_plot_subset
dev.off()

#Cell type propotion 

# rename identify from cluster to cell type
celltype.ids <- c("KC","KC","FB","KC","Tcell","Mo", 
                  "KC","vEC","MLA","lEC","SG","SMC",
                  "AT","KC","MC")

names(celltype.ids) <- levels(TN.combined)
TNtype.combined <- RenameIdents(TN.combined, celltype.ids)

nb.cols <- length(unique(celltype.ids))
typecolors <- colorRampPalette(brewer.pal(nb.cols, "Paired"))(nb.cols)

# plot UMAP 
UMAP_type_plot_all<- DimPlot(TNtype.combined, reduction = "umap",label = TRUE, pt.size = 0.8, label.size=7, repel = T, cols = typecolors)
UMAP_type_plot_subset<- DimPlot(TNtype.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, repel = T, cols = typecolors)

pdf("C:/Users/user/Desktop/202404/UMAP_celltype_plot_all.pdf", width = 15, height = 15)
UMAP_type_plot_all
dev.off()

pdf("C:/Users/user/Desktop/202404/UMAP_celltype_plot_subset.pdf", width = 15, height = 15)
UMAP_type_plot_subset
dev.off()

#Cell type proportion
celltypeproportion <- table(Idents(TNtype.combined), TN.combined$orig.ident1)
celltypeproportion <- round(sweep(celltypeproportion,MARGIN=2, STATS=colSums(celltypeproportion), FUN = "/")*100,2)
write.csv(celltypeproportion, file = "C:/Users/user/Desktop/202404/Cell_type_proportion.csv", row.names = T)
celltypeproportion <- as.data.frame(celltypeproportion)
type_proportion_plot <- ggplot(celltypeproportion, aes(x = Var2, y = Freq, fill = Var1)) + theme_bw(base_size = 15) + 
  geom_col(position = "fill", width = 0.6) + xlab("Sample") + ylab("Proportion") + labs(title="Cell type propotion") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size=24))+ scale_fill_manual(values = typecolors)

pdf("C:/Users/user/Desktop/202404/proportion_plot.pdf", width = 15, height = 15)
type_proportion_plot
dev.off()



##All cells cluster heatmap
Heatmapall <- subset(TN.combined, idents = seq(0,21)) 
Heatmapall.markers <- FindAllMarkers(Heatmapall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(Heatmapall.markers)
write.csv(Heatmapall.markers , file = "C:/Users/user/Desktop/202404/Findallmarkers_2.csv", col.names = TRUE) 
top10 <- Heatmapall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
marker_heatmap <- DoHeatmap(Heatmapall, features = top10$gene)

pdf("C:/Users/user/Desktop/202404/marker_heatmap.pdf", width = 15, height = 15)
marker_heatmap
dev.off()

##Save the merged Seurat object
saveRDS(TN.combined, file = "C:/Users/user/Desktop/202404/TN.combined_20240416.rds")


