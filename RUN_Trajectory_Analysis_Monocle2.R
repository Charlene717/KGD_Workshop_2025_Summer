library(monocle)

# Convert Seurat object to Monocle format
data <- as(as.matrix(integrated@assays$RNA@data), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", data = integrated@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(data)))
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds) %>% estimateDispersions()
cds <- reduceDimension(cds, method = "DDRTree") %>% orderCells()
plot_cell_trajectory(cds, color_by = "seurat_clusters")