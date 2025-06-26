#### Cell Type Annotation ####
library(SingleR); library(celldex)
hpca <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(seurat_all_integrated, slot = "data")

pred <- SingleR(test = counts, ref = hpca, labels = hpca$label.main)

# View results
table(pred$pruned.labels)

SingleR_hpca_clustering.table <- table(pred@listData[["pruned.labels"]], seurat_all_integrated@active.ident)
SingleR_hpca_clustering.table

seurat_all_integrated$Label_SingleR_HPCA <- pred$pruned.labels

# UMAP 使用 HPCA label
DimPlot(
  seurat_all_integrated,
  reduction = "umap",
  group.by = "Label_SingleR_HPCA",   # meta.data 中對應到 "Label_SingleR_HPCA"
  label = TRUE, repel = TRUE   # 顯示 cluster label 並避免重疊
) + ggtitle("UMAP colored by SingleR (HPCA)") -> Plot_SingleR_HPCA
print(Plot_SingleR_HPCA)

