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
## 4. 畫泡泡圖（x = gene, facet = CellType；y = Seurat cluster）
###############################################################################
library(ggplot2)

# 為了在 facet 中保持 CellType 順序，可先設定 factor levels
plot_df$CellType <- factor(plot_df$CellType, levels = names(marker_sets_KGD))

# gene 順序直接沿用在各 CellType 內的原始順序即可
plot_df$gene <- factor(plot_df$gene,
                       levels = unique(gene_present))   # 固定整體順序

# cluster 改成 factor，確保 y 軸順序由上而下 0,1,2…
plot_df$cluster <- factor(plot_df$cluster,
                          levels = sort(unique(plot_df$cluster)))

p <- ggplot(plot_df, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_exp, colour = avg_exp)) +
  scale_colour_gradient(low = "white", high = "darkred") +
  scale_size(range = c(0, 8)) +
  facet_grid(. ~ CellType, scales = "free_x", space = "free_x", switch = "x") +
  labs(
    x      = NULL,                 # x 軸留白；CellType 會顯示在 facet strip
    y      = "Seurat cluster",
    size   = "% cells",
    colour = "Avg expr",
    title  = "KGD 皮膚細胞 Marker Bubble Plot"
  ) +
  theme_bw() +
  theme(
    strip.placement   = "outside",             # strip 移到外側（上方）
    strip.text.x.top  = element_text(size = 9, face = "bold"),
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y       = element_text(size = 7),
    panel.spacing.x   = unit(0.2, "lines"),    # 壓縮 CellType 之間間距
    plot.title        = element_text(hjust = 0.5, size = 24)
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
