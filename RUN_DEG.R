library(Seurat)
library(dplyr)


Idents(seurat_all_integrated) <- "seurat_clusters"

# Idents(seurat_all_integrated) <- "celltype"


# 找出所有 cluster 的 marker genes（預設為 ident 設定的 clustering）
cluster_markers <- FindAllMarkers(seurat_all_integrated, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)

# 篩選顯著 marker
cluster_markers_sig <- cluster_markers %>% filter(p_val_adj < 0.05)

# 查看 top10 marker per cluster
top10_markers <- cluster_markers_sig %>% group_by(cluster) %>% top_n(10, avg_log2FC)


# cluster_markers_sig %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10
DoHeatmap(seurat_all_integrated, features = top10_markers$gene) + NoLegend()



#################################################################################


# 設定群組標籤為 orig.ident
Idents(seurat_all_integrated) <- "orig.ident"

# 僅挑出兩組樣本
seurat_subset <- subset(seurat_all_integrated, idents = c("GSM6111844", "GSM6111847"))

# DEG 分析：GSM6111844 vs GSM6111847
deg_result <- FindMarkers(seurat_subset, 
                          ident.1 = "GSM6111844", 
                          ident.2 = "GSM6111847", 
                          logfc.threshold = 0.25, 
                          min.pct = 0.1)

# 查看 top DEG
head(deg_result[order(deg_result$p_val_adj), ])

