###############################################################################
# 0. 安裝並載入必要 R 套件  -----------------------------------------------------
#    ◆ 原則：若套件不存在才安裝 (if(!require()))，安裝完立即 library() 載入
#    ◆ Presto 不在 CRAN，須透過 devtools 從 GitHub 安裝
###############################################################################
if (!require("Seurat"))      { install.packages("Seurat");      library(Seurat)      }  # 單細胞分析核心
if (!require("tidyverse"))   { install.packages("tidyverse");   library(tidyverse)   }  # dplyr / ggplot2…整合套件
if (!require("patchwork"))   { install.packages("patchwork");   library(patchwork)   }  # ggplot2 版面拼圖
if (!require("ggrepel"))     { install.packages("ggrepel");     library(ggrepel)     }  # 文字防重疊標籤
if (!require("pheatmap"))    { install.packages("pheatmap");    library(pheatmap)    }  # 熱圖
if (!require("future"))      { install.packages("future");      library(future)      }  # 平行運算
if (!require("glue"))        { install.packages("glue");        library(glue)        }  # 字串格式化
if (!require("devtools"))    { install.packages("devtools");    library(devtools)    }  # 安裝 GitHub 套件
if (!require("presto")) {                 # Presto：加速 FindMarkers/FindAllMarkers
  devtools::install_github("immunogenomics/presto")
  library(presto)
}

###############################################################################
# 1. 設定平行化策略與記憶體上限  -----------------------------------------------
#    • Windows 不支援 multicore fork，改用 multisession
#    • workers = 20 代表同時開 20 核心/Session；視硬體調整
###############################################################################
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = 20)   # Windows platform
} else {
  plan(multicore,    workers = 20)   # macOS / Linux
}
options(future.globals.maxSize = 200 * 1024^3)  # ≈ 200 GB，防止「記憶體不足」錯誤

###############################################################################
# 2. 整體叢集 (seurat_clusters) Marker 基因偵測  -------------------------------
#    • 以整合 (integrated) 資料層做差異分析，能降低批次效應
###############################################################################
DefaultAssay(seurat_all_integrated) <- "integrated"  # 指向整合矩陣
Idents(seurat_all_integrated)       <- "seurat_clusters"  # 設定目前分群

# 2-1. FindAllMarkers：找出各群特異基因 (只保留上調基因)
cluster_markers <- FindAllMarkers(
  object          = seurat_all_integrated,
  only.pos        = TRUE,   # 只看表現量較高的基因
  min.pct         = 0.25,   # 至少 25% 細胞有表現
  logfc.threshold = 0.25    # |Log2FC| ≥ 0.25
)

# 2-2. 篩選顯著基因 (FDR < 0.05)，並各取 Top10 畫熱圖
cluster_markers_sig <- cluster_markers %>% filter(p_val_adj < 0.05)

top10_markers <- cluster_markers_sig %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  pull(gene)                       # 取出基因名稱向量

DoHeatmap(
  object   = seurat_all_integrated,
  features = top10_markers
) + NoLegend()

###############################################################################
# 3. 逐群差異表達 (DEG)：比較 GSM6111844 vs GSM6111847 ------------------------
#    • 只在指定兩個樣本內比較，同群不同樣本，減少其他批次影響
###############################################################################
DefaultAssay(seurat_all_integrated) <- "RNA"          # 切回原始計數
Idents(seurat_all_integrated)       <- "seurat_clusters"
clusters <- levels(seurat_all_integrated)             # 取得所有 cluster 標籤

# 3-1. 僅保留目標樣本
seurat_filtered <- subset(
  x      = seurat_all_integrated,
  subset = orig.ident %in% c("GSM6111844", "GSM6111847")
)

deg_list <- vector("list", length(clusters))          # 預先建立 list 容器
names(deg_list) <- clusters

# 3-2. 逐群執行 DEG
for (cl in clusters) {
  message(glue(">>> 處理 Cluster {cl}"))
  
  ## (a) 抓該群細胞；若數量過少可跳過，避免統計不穩
  cells_in_cluster <- WhichCells(seurat_filtered, idents = cl)
  if (length(cells_in_cluster) < 50) {                # ←<< 修改門檻
    warning(glue("Cluster {cl} 細胞數 < 50，略過分析"))
    next
  }
  
  ## (b) 建立 cluster 子集，並用 orig.ident (樣本 ID) 當分群
  seurat_sub <- subset(seurat_filtered, cells = cells_in_cluster)
  Idents(seurat_sub) <- "orig.ident"
  
  ## (c) FindMarkers：GSM6111844 vs GSM6111847
  deg <- FindMarkers(
    object          = seurat_sub,
    ident.1         = "GSM6111844",
    ident.2         = "GSM6111847",
    logfc.threshold = 0.25,  # fold-change 門檻
    min.pct         = 0.10   # 基因表現百分比門檻
  ) %>%
    rownames_to_column("gene") %>%        # 轉成顯式 gene 欄位
    mutate(cluster = cl)                  # 記錄來源 cluster
  
  deg_list[[cl]] <- deg                   # 收進 list
}

###############################################################################
# 4. 火山圖 (Volcano plot)：直觀呈現 DEG 結果  -------------------------------
#    • 每個 cluster 各輸出一張火山圖；標註最顯著前 10 genes
###############################################################################
for (cl in names(deg_list)) {
  deg <- deg_list[[cl]]
  if (is.null(deg) || nrow(deg) == 0) next          # 若無分析結果跳過
  
  ## 4-1. 標記顯著性 (Up / Down / Not Sig)
  deg <- deg %>%
    mutate(significance = case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
      TRUE                                   ~ "Not Sig"
    ))
  
  ## 4-2. 繪圖
  ggplot(
    deg,
    aes(x = avg_log2FC, y = -log10(p_val_adj), colour = significance)
  ) +
    geom_point(alpha = 0.6) +
    scale_colour_manual(
      values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey70")
    ) +
    geom_text_repel(                                           # 標註前 10 顯著基因
      data         = slice_head(arrange(deg, p_val_adj), n = 10),
      aes(label    = gene),
      size         = 3,
      max.overlaps = 20
    ) +
    labs(
      title = glue("Volcano Plot • Cluster {cl} (GSM6111844 vs GSM6111847)"),
      x     = "log2 Fold-Change",
      y     = expression(-log[10]("adjusted p-value"))
    ) +
    theme_minimal(base_size = 12) -> p
  
  print(p)           # RStudio / R GUI 會即時顯示圖
}

###############################################################################
# End of Script  --------------------------------------------------------------
# • 如需輸出結果到 CSV / Excel，可在 deg_list 迴圈內用 write.csv() 等函數
# • 若要將圖存檔：ggsave(filename = ..., plot = p, width = 6, height = 5)
###############################################################################
