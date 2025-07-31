if(!require('Seurat'))         { install.packages('Seurat');         library(Seurat) }
if(!require('tidyverse'))      { install.packages('tidyverse');      library(tidyverse) }
if(!require('dplyr'))          { install.packages('dplyr');          library(dplyr) }



## Speed up
if(!require('future')) install.packages('future'); library(future)
## https://github.com/immunogenomics/presto
if(!require('presto')) devtools::install_github("immunogenomics/presto"); library(presto) # Speeds up FindAllMarkers
plan(multicore, workers = 20)
options(future.globals.maxSize = 2048*100 * 1024^2) # Set memory limit to ~204.8 GB
################################################################################

DefaultAssay(seurat_all_integrated) <- "integrated"

# 設定群組標籤為 seurat_clusters
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
if(!require('Seurat'))         { install.packages('Seurat');         library(Seurat) }
if(!require('tidyverse'))      { install.packages('tidyverse');      library(tidyverse) }
if(!require('patchwork'))      { install.packages('patchwork');      library(patchwork) }
if(!require('ggplot2'))        { install.packages('ggplot2');        library(ggplot2) }
if(!require('dplyr'))          { install.packages('dplyr');          library(dplyr) }



DefaultAssay(seurat_all_integrated) <- "RNA"
# 確保 cluster 編號存在
Idents(seurat_all_integrated) <- "seurat_clusters"
clusters <- levels(seurat_all_integrated)

# 只保留 GSM6111844 和 GSM6111847
seurat_filtered <- subset(seurat_all_integrated, subset = orig.ident %in% c("GSM6111844", "GSM6111847"))

# 依 cluster 逐一做 DEG
deg_list <- list()

for (cl in clusters) {
  try({
    message("Analyzing cluster ", cl)
    
    # 篩出該 cluster 中的細胞
    cells_in_cluster <- WhichCells(seurat_filtered, idents = cl)
    seurat_sub <- subset(seurat_filtered, cells = cells_in_cluster)
    
    # 設定為 GSM ID 作為比較群
    Idents(seurat_sub) <- "orig.ident"
    
    # 做 DEG 分析
    deg <- FindMarkers(seurat_sub,
                       ident.1 = "GSM6111844",
                       ident.2 = "GSM6111847",
                       logfc.threshold = 0.25,
                       min.pct = 0.1)
    deg$gene <- rownames(deg)
    deg$cluster <- cl
    deg_list[[cl]] <- deg
    
  })
}




library(ggrepel)
library(pheatmap)

for (cl in names(deg_list)) {
  deg <- deg_list[[cl]]
  
  # --- 火山圖 ---
  deg$significance <- "Not Sig"
  deg$significance[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.25] <- "Up"
  deg$significance[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.25] <- "Down"
  
  p <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("blue", "gray", "red")) +
    ggtitle(paste("Volcano Plot - Cluster", cl)) +
    theme_minimal() +
    geom_text_repel(data = head(deg[order(deg$p_val_adj), ], 10),
                    aes(label = gene),
                    size = 3, max.overlaps = 20)
  
  print(p)

}

