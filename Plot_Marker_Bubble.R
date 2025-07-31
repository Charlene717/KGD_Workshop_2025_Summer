###############################################################################
## 0. 準備
###############################################################################
library(Seurat)      # >= 4.3
library(tidyverse)   # >= 2.0

DefaultAssay(seurat_all_integrated) <- "RNA"
Idents(seurat_all_integrated)       <- "seurat_clusters"

###############################################################################
## 1. 檢查哪些 marker 存在
###############################################################################
gene_vec_all <- unlist(marker_sets_KGD, use.names = FALSE)          # 原始列表
gene_present <- intersect(gene_vec_all, rownames(seurat_all_integrated))   # 找得到的
gene_missing <- setdiff(gene_vec_all, gene_present)                          # 找不到的

if (length(gene_missing) > 0){
  message("⚠️  下列 ", length(gene_missing),
          " 個 marker 在 RNA assay 裡找不到，將會略過：\n",
          paste(gene_missing, collapse = ", "))
}

# 和 gene_present 對應的 CellType 向量（保持順序）
celltype_present <- rep(names(marker_sets_KGD), lengths(marker_sets_KGD))[gene_vec_all %in% gene_present]

###############################################################################
## 2. 抓表達矩陣＋cluster
###############################################################################
expr_mat <- FetchData(
  object = seurat_all_integrated,
  vars   = gene_present,       # 只抓存在的基因
  slot   = "data"
)

expr_mat$cluster <- Idents(seurat_all_integrated)

###############################################################################
## 3. 轉長格式並計算 Avg / %Exp
###############################################################################
expr_long <- expr_mat |>
  rownames_to_column("cell") |>
  pivot_longer(
    cols      = all_of(gene_present),
    names_to  = "gene",
    values_to = "value"
  )

plot_df <- expr_long |>
  group_by(cluster, gene) |>
  summarise(
    avg_exp = mean(value),
    pct_exp = mean(value > 0) * 100,
    .groups = "drop"
  ) |>
  mutate(
    CellType = celltype_present[match(gene, gene_present)],
    gene_ct  = paste(CellType, gene, sep = " : ")
  )

# y 軸層級：依 CellType 區段由上到下
plot_df$gene_ct <- factor(plot_df$gene_ct,
                          levels = rev(unique(plot_df$gene_ct)))

###############################################################################
## 4. 畫泡泡圖
###############################################################################
library(ggplot2)
p <- ggplot(plot_df, aes(x = cluster, y = gene_ct)) +
  geom_point(aes(size = pct_exp, colour = avg_exp)) +
  scale_colour_gradient(low = "white", high = "darkred") +
  scale_size(range = c(0, 8)) +
  facet_grid(CellType ~ ., scales = "free_y", space = "free_y") +
  labs(
    x      = "Seurat cluster",
    y      = NULL,
    size   = "% cells",
    colour = "Avg expr",
    title  = "KGD 皮膚細胞 Marker Bubble Plot"
  ) +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 24)
  )

p
###############################################################################
## 5. 匯出
###############################################################################
ggsave(
  filename = file.path(Name_ExportFolder_KGD_CTAnnot,
                       paste0(Name_Export, "_Markers_BubblePlot_AllInOne.pdf")),
  plot     = p,
  width    = 22, height = 40, units = "in"
)
