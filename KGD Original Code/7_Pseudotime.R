#From Alice
#04012025

#install Monocle2 for trajectory
BiocManager::install("monocle")

Sys.setenv(LANGUAGE='en')


library(monocle)
library(Seurat)



#040225

#Use PD combined
PD.combined <- readRDS(file = "C:/KGD_Lab/Member/Nat/NTDA_NT_Updated 032425/PD.combined.rds")


#Fibroblast
PD.combined <- readRDS(file = "C:/KGD_Lab/Member/Nat/NTDA_NT_Updated 032425/PD.combined.rds")
PD1.subset <- subset(PD.combined, idents = c(0,9,20))
DimPlot(PD1.subset, reduction = "umap",label = TRUE, pt.size = 0.8)
FB.subgroup <- PD1.subset

#Tcells
PD.combined <- readRDS(file = "C:/KGD_Lab/Member/Nat/NTDA_NT_Updated 032425/PD.combined.rds")
PD2.subset <- subset(PD.combined, idents = c(1,3,10))
DimPlot(PD2.subset, reduction = "umap",label = TRUE, pt.size = 0.8)
CD8+Tcells.subgroup <- PD2.subset


#Epithelial Cells 

PD.subset <- subset(PD.combined, idents = c(2,8,13,19,22))
DimPlot(PD.subset, reduction = "umap",label = TRUE, pt.size = 0.8)
Epi.subgroup <- PD.subset

##Rename cluster (for example?GFibroblast?BMacrophage......)
new.cluster.ids <- c("C2","C8","C13","C19","C22")
names(new.cluster.ids) <- levels(Epi.subgroup)
Epithelial.subgroup <- RenameIdents(Epi.subgroup, new.cluster.ids)
DimPlot(Epithelial.subgroup, reduction = "umap",label = TRUE, pt.size = 0.8)

#Inserted codes 
Epithelial.subgroup@meta.data[["seurat_clusters"]]<-droplevels(Epithelial.subgroup@meta.data[["seurat_clusters"]])
levels(Epithelial.subgroup@meta.data[["seurat_clusters"]])
Epi.subgroup <- Epithelial.subgroup

#Use the original numbers
#DO not RUN
##change seurat_clusters level names
levels(Epithelial.subgroup@meta.data[["seurat_clusters"]]) <- c("2","8","13","19","22")
head(Epithelial.subgroup@active.ident)
head(Epithelial.subgroup@meta.data[["seurat_clusters"]])
#Fname.subgroup1 <- subset(Fname.subgroup, idents = c("Macro1","Macro2"))
#DO not RUN just walk char

#Monocle2 
data <- as(as.matrix(Epi.subgroup@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Epi.subgroup@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
Epi.subgroup <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.25,expressionFamily = negbinomial.size())
Epi.subgroup


#View data
pData(Epi.subgroup)
fData(Epi.subgroup)

#Estimate size factors and dispersions
Epi.subgroup <- estimateSizeFactors(Epi.subgroup)
Epi.subgroup <- estimateDispersions(Epi.subgroup)

#Filtering low-quality cells 
Epi.subgroup <- detectGenes(Epi.subgroup, min_expr = 0.1)
print(head(pData(Epi.subgroup)))
expressed_genes <- row.names(subset(fData(Epi.subgroup),num_cells_expressed >= 5))
L <- log(exprs(Epi.subgroup[expressed_genes,]))
mL <- apply(L,1,function(x){mean(x[is.finite(x)])})
sdL <- apply(L,1,function(x){sd(x[is.finite(x)])})
Lstd <- (L-mL)/sdL

library(reshape)
#melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
#                                 ^^^^^

#To check how much time it would take
#system.time(your_function(detectGenes()))

melted_dens_df <- melt(as.matrix(Lstd))
qplot(value, geom = "density", data = melted_dens_df) +stat_function(fun = dnorm, size = 0.5, color = 'red') +xlab("Standardized log(FPKM)") + ylab("Density")


#Get more RAM
gc()
memory.limit(9999999999)
gc() 


#Alternative choices for cell trajectory by ordering genes
Epi.subgroup <- detectGenes(Epi.subgroup, min_expr = 0.1)
fData(Epi.subgroup)$use_for_ordering <-
  fData(Epi.subgroup)$num_cells_expressed > 0.05 * ncol(Epi.subgroup)
plot_pc_variance_explained(Epi.subgroup, return_all = F)
Epi.subgroup <- reduceDimension(Epi.subgroup,
                                  max_components = 2,
                                  norm_method = 'log',
                                  num_dim = 3,
                                  reduction_method = 'PCA',
                                  verbose = T)
Epi.subgroup <- clusterCells(Epi.subgroup, verbose = F)
plot_cell_clusters(Epi.subgroup, color_by = 'as.factor(seurat_clusters)')

Epi.subgroup <- estimateSizeFactors(Epi.subgroup)  # Fix size factors
Epi.subgroup <- estimateDispersions(Epi.subgroup)  # Ensure proper normalization
plot_pc_variance_explained(Epi.subgroup, return_all = F)  # Retry plotting


#Epi.subgroup <- clusterCells_Density_Peak(Epi.subgroup)
#find("clusterCells_Density_Peak")

#Epi.subgroup <- clusterCells(Epi.subgroup)
#plot_cells(Epi.subgroup, color_cells_by = "partition")


#This part was skipped
plot_rho_delta(Epi.subgroup, rho_threshold = NULL, delta_threshold = NULL )
Epi.subgroup <- clusterCells(Epi.subgroup,
                               rho_threshold = 2,
                               delta_threshold = 4,
                               skip_rho_sigma = T,
                               verbose = F)
#plot_cell_clusters(Epi.subgroup, color_by = 'as.factor(seurat_clusters)')
#clustering_DEG_genes <-
#  differentialGeneTest(Epi.subgroup[expressed_genes,],
#                       fullModelFormulaStr = '~Cluster',
#                      cores = 1)

#To speed up, just choose top 100 DEGs
clustering_DEG_genes <-
  differentialGeneTest(Epi.subgroup[1:100,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)

#Changed from 1:1000 to 1:100
ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:100]

Epi.subgroup <-
  setOrderingFilter(Epi.subgroup,
                    ordering_genes = ordering_genes)

#New code to speed up
Epi.subgroup <- reduceDimension(Epi.subgroup, method = "PCA", num_dim = 50,auto_param_selection = F)
Epi.subgroup <- reduceDimension(Epi.subgroup, method = 'DDRTree')

Epi.subgroup <- reduceDimension(Epi.subgroup, method = "PCA")  # Very fast
Epi.subgroup <- reduceDimension(Epi.subgroup, method = "tSNE") # Faster than DDRTree



#Changed from numberclusters=3 to auto_param_selection=F 
Epi.subgroup <-
  reduceDimension(Epi.subgroup, method = 'tSNE',auto_param_selection = F)

#Should I use this instead?
Epi.subgroup <-
  reduceDimension(Epi.subgroup, method = 'PCA',auto_param_selection = F)

Epi.subgroup <-
  orderCells(Epi.subgroup)


plot_cell_trajectory(Epi.subgroup, color_by = "seurat_clusters")




################################################################################
#Original Code
#1
#use all merge data!!!
#Macro.subgroup <- readRDS(file = "C:/Alice/Analysis_Outputs/TN.combined.0311.rds")
TN.combined <- readRDS(file = "C:/Alice/Analysis_Outputs/TN.combined.0311.rds")
TN.subset <- subset(TN.combined, idents = c(28,29))
DimPlot(TN.subset, reduction = "umap",label = TRUE, pt.size = 0.8)
Macro.subgroup <- TN.subset

#2
##Rename cluster (for example?GFibroblast?BMacrophage......)
new.cluster.ids <- c("C28","C29")
names(new.cluster.ids) <- levels(Macro.subgroup)
Fname.subgroup <- RenameIdents(Macro.subgroup, new.cluster.ids)
DimPlot(Fname.subgroup, reduction = "umap",label = TRUE, pt.size = 0.8)

#3
##change seurat_clusters level names
levels(Fname.subgroup@meta.data[["seurat_clusters"]]) <- c("C28","C29")
head(Fname.subgroup@active.ident)
head(Fname.subgroup@meta.data[["seurat_clusters"]])
#Fname.subgroup1 <- subset(Fname.subgroup, idents = c("Macro1","Macro2"))
Macro.subgroup <- Fname.subgroup

#3
#Monocle2 
data <- as(as.matrix(Macro.subgroup@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Macro.subgroup@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
Macro.subgroup <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.25,expressionFamily = negbinomial.size())
Macro.subgroup

#4
#View data
pData(Macro.subgroup)
fData(Macro.subgroup)

#Estimate size factors and dispersions
Macro.subgroup <- estimateSizeFactors(Macro.subgroup)
Macro.subgroup <- estimateDispersions(Macro.subgroup)

#Filtering low-quality cells 
Macro.subgroup <- detectGenes(Macro.subgroup, min_expr = 0.1)
print(head(pData(Macro.subgroup)))
expressed_genes <- row.names(subset(fData(Macro.subgroup),num_cells_expressed >= 5))
L <- log(exprs(Macro.subgroup[expressed_genes,]))
mL <- apply(L,1,function(x){mean(x[is.finite(x)])})
sdL <- apply(L,1,function(x){sd(x[is.finite(x)])})
Lstd <- (L-mL)/sdL

library(reshape)
#melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
#                                 ^^^^^

melted_dens_df <- melt(as.matrix(Lstd))
qplot(value, geom = "density", data = melted_dens_df) +stat_function(fun = dnorm, size = 0.5, color = 'red') +xlab("Standardized log(FPKM)") + ylab("Density")

#Get more RAM
gc()
memory.limit(9999999999)
gc() 


#Alternative choices for cell trajectory by ordering genes
Epi.subgroup <- detectGenes(Epi.subgroup, min_expr = 0.1)
fData(Epi.subgroup)$use_for_ordering <-
  fData(Epi.subgroup)$num_cells_expressed > 0.05 * ncol(Epi.subgroup)
plot_pc_variance_explained(Epi.subgroup, return_all = F)
Epi.subgroup <- reduceDimension(Epi.subgroup,
                                  max_components = 2,
                                  norm_method = 'log',
                                  num_dim = 3,
                                  reduction_method = 'tSNE',
                                  verbose = T)
Epi.subgroup <- clusterCells(Epi.subgroup, verbose = F)
plot_cell_clusters(Epi.subgroup, color_by = 'as.factor(seurat_clusters)')

plot_rho_delta(Epi.subgroup, rho_threshold = 2, delta_threshold = 4 )
Epi.subgroup <- clusterCells(Epi.subgroup,
                               rho_threshold = 2,
                               delta_threshold = 4,
                               skip_rho_sigma = T,
                               verbose = F)
plot_cell_clusters(Epi.subgroup, color_by = 'as.factor(seurat_clusters)')
clustering_DEG_genes <-
  differentialGeneTest(Epi.subgroup[expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)

ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

Epi.subgroup <-
  setOrderingFilter(Epi.subgroup,
                    ordering_genes = ordering_genes)

Epi.subgroup <-
  reduceDimension(Epi.subgroup, method = 'DDRTree')

Epi.subgroup <-
  orderCells(Epi.subgroup)


plot_cell_trajectory(Epi.subgroup, color_by = "seurat_clusters")

################################################################################


# 確保 Monocle 已載入
library(monocle)

# 確保 Macrophage.subgroup 和 Fibroblast.subgroup 存在
if (!exists("Macrophage.subgroup")) {
  stop("Error: 'Macrophage.subgroup' does not exist. Please check your data loading step.")
}
if (!exists("Fibroblast.subgroup")) {
  stop("Error: 'Fibroblast.subgroup' does not exist. Please check your data loading step.")
}

# 檢查基因表達資料
print(dim(exprs(Macrophage.subgroup)))

# 遺失數據填補
Macrophage.subgroup <- estimateSizeFactors(Macrophage.subgroup)
Macrophage.subgroup <- estimateDispersions(Macrophage.subgroup)

# 檢測基因表達
Macrophage.subgroup <- detectGenes(Macrophage.subgroup, min_expr = 0.1)

# 確保 fData 可用
if (!"num_cells_expressed" %in% colnames(fData(Macrophage.subgroup))) {
  stop("Error: 'num_cells_expressed' not found in fData(Macrophage.subgroup). Please check preprocessing.")
}

# 設定用於排序的基因
fData(Macrophage.subgroup)$use_for_ordering <- 
  fData(Macrophage.subgroup)$num_cells_expressed > (0.05 * ncol(Macrophage.subgroup))

# 降維分析
Macrophage.subgroup <- reduceDimension(Macrophage.subgroup,
                                       max_components = 2,
                                       norm_method = 'log',
                                       num_dim = 3,
                                       reduction_method = 'tSNE',
                                       verbose = TRUE)

# 聚類分析
if (ncol(Fibroblast.subgroup) > 5) {  # 確保有足夠的細胞數
  Fibroblast.subgroup <- clusterCells(Fibroblast.subgroup, verbose = FALSE)
} else {
  stop("Error: Not enough cells in Fibroblast.subgroup for clustering.")
}

# 視覺化
plot_cell_clusters(Fibroblast.subgroup, color_by = 'as.factor(seurat_clusters)')

# 執行 Differential Expression Test
expressed_genes <- row.names(subset(fData(Fibroblast.subgroup), num_cells_expressed >= 5))
if (length(expressed_genes) == 0) {
  stop("Error: No expressed genes found for differential expression analysis.")
}

clustering_DEG_genes <- differentialGeneTest(Fibroblast.subgroup[expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 1)


library(RColorBrewer)
nb.cols <- 3
mycolors <- colorRampPalette(brewer.pal(3, "Paired"))(nb.cols)

plot_cell_trajectory(Fibroblast.subgroup, color_by = "seurat_clusters")+ scale_color_manual(values = mycolors, name = "seurat_clusters")

plot_cell_trajectory(Fibroblast.subgroup, color_by = "State")

Fibroblast.subgroup <- orderCells(Fibroblast.subgroup, root_state = 1)

plot_cell_trajectory(Fibroblast.subgroup, color_by = "Pseudotime")

plot_cell_trajectory(Fibroblast.subgroup, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 3)

plot_cell_trajectory(Fibroblast.subgroup, color_by = "seurat_clusters") +
  facet_wrap(~orig.ident1, nrow = 3)

plot_cell_trajectory(Fibroblast.subgroup, color_by = "Pseudotime") +
  facet_wrap(~orig.ident1, nrow = 3)

saveRDS(Fibroblast.subgroup, file = "C:/Alice/pseudotime/TN.combined_monocle2.rds")
Fibroblast.subgroup <- readRDS(file = "C:/Alice/pseudotime/TN.combined_monocle2.rds")

#Marker expressed genes
plot_cell_trajectory(Fibroblast.subgroup, markers = c("IL4R","IL13RA1","POSTN","COL1A1","FN1"), show_tree = FALSE, use_color_gradient = TRUE,show_branch_points = FALSE)+
  scale_colour_gradient(low="lightgrey",high="brown2") + theme(legend.position = "right")+facet_wrap(~orig.ident1, nrow = 2)


#Finding Genes that Change as a Function of Pseudotime
to_be_tested <- row.names(subset(fData(Fibroblast.subgroup),
                                 gene_short_name %in% c("POSTN","COL1A1","FN1","CTHRC1","ACTA2","APOE","APCDD1")))
cds_subset <- Fibroblast.subgroup[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters", panel_order = c("POSTN","COL1A1","FN1","CTHRC1","ACTA2","APOE","APCDD1"))
plot_genes_in_pseudotime(cds_subset, color_by = "orig.ident1", panel_order = c("POSTN","COL1A1","FN1","CTHRC1","ACTA2","APOE","APCDD1"))+ scale_color_manual(values=c("#00BA38","#F8766D"))

#Finding Genes that Change as a Function of Pseudotime splited by Two sample branches
to_be_tested <- row.names(subset(fData(Fibroblast.subgroup),
                                 gene_short_name %in% c("IL4R","IL13RA1","CCL19","COL1A1","FN1","POSTN","MFAP5","APCDD1")))
plot_genes_branched_pseudotime(Fibroblast.subgroup[to_be_tested,],
                               branch_point = 2, branch_labels = c("Normal","AK"),
                               color_by = "TN.merge",
                               ncol = 1, panel_order = c("IL4R","IL13RA1","CCL19","COL1A1","FN1","POSTN","MFAP5","APCDD1"))

#Analyzing Branches in Single-Cell Trajectories by kinetic heatmap
BEAM_res <- BEAM(Fibroblast.subgroup, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]


BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(Fibroblast.subgroup[row.names(subset(BEAM_res,
                                                                 qval < 0.5e-2)),],
                            branch_point = 1,
                            num_clusters = 7,
                            cores = 1,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)

branch_genes <- plot_genes_branched_heatmap(Fibroblast.subgroup[row.names(subset(BEAM_res,
                                                                                 qval < 0.5e-2)),],
                                            branch_point = 1,
                                            num_clusters = 7,
                                            cores = 1,
                                            use_gene_short_name = TRUE,
                                            show_rownames = TRUE,
                                            return_heatmap = T)
print(branch_genes)
write.csv(branch_genes, file = "C:/Alice/pseudotime/macrophage.csv", row.names = FALSE)


saveRDS(BEAM_res, file = "C:/Alice/pseudotime/macrophage.csv")
BEAM_res <- readRDS(file = "C:/Alice/pseudotime/macrophage.csv")
