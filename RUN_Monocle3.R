################################################################################
## 🔧 載入必要套件
################################################################################
if (!require("Seurat"))            { install.packages("Seurat");                                  library(Seurat) }     # ↪ 單細胞分析主套件
if (!require("tidyverse"))         { install.packages("tidyverse");                               library(tidyverse) }  # ↪ tidyverse: dplyr/ggplot2 等整合工具，方便資料處理

if (!require("monocle3"))          { devtools::install_github('cole-trapnell-lab/monocle3');      library(monocle3) }  # ↪ Pseudotime 分析與軌跡重建
if (!require("SeuratWrappers"))    { remotes::install_github('satijalab/seurat-wrappers');        library(SeuratWrappers) }  # ↪ Seurat 和 Monocle 等工具間的轉換橋接


################################################################################
## 📥 讀取 Seurat 物件，並設定 clustering 身分
################################################################################
# seurat_path <- "path/to/your/seurat_all_integrated.rds"
# seurat_all_integrated <- readRDS(seurat_path)             # 載入儲存好的 Seurat RDS 物件
Idents(seurat_all_integrated) <- "seurat_clusters"        # 指定 Seurat 用來分群的欄位

################################################################################
## 🔁 將 Seurat 轉換為 Monocle3 的 CellDataSet (CDS)
################################################################################
cds <- as.cell_data_set(seurat_all_integrated)            # Seurat ➜ Monocle 格式

# ➕ 把 metadata 合併到 colData（避免重複欄位）
meta_to_add <- seurat_all_integrated@meta.data
meta_to_add <- meta_to_add[, !colnames(meta_to_add) %in% colnames(colData(cds))]
colData(cds) <- cbind(colData(cds), meta_to_add)

# ➕ 把原本的 Seurat UMAP 降維結果寫入 Monocle3 的 reducedDims slot
reducedDims(cds)$UMAP <- Embeddings(seurat_all_integrated, reduction = "umap")

# ➕ 將 Seurat 分群結果複製進 Monocle3（供 plot_cells 使用）
cds@clusters$UMAP$clusters <- factor(Idents(seurat_all_integrated))

# ⚠️ 手動指定所有細胞都屬於同一個 partition，否則後續會報錯
cds@clusters$UMAP$partitions <- factor(rep(1, length(Cells(cds))))
names(cds@clusters$UMAP$partitions) <- Cells(cds)  # 確保有正確 names

# ❗不需執行 preprocess_cds() 或 reduce_dimension()，因為已從 Seurat 導入 UMAP


####################################################################################################
cds_Ori <- cds

plot_cells(cds, color_cells_by = "seurat_clusters")

############################################################
##  產生僅含指定 Seurat cluster 的 CellDataSet 子集        ##
##  • 來源物件：cds（Monocle 3 的 CellDataSet）           ##
##  • 欄位名稱：假設 Seurat cluster 已存為                 ##
##                colData(cds)$seurat_clusters             ##
############################################################

## 1. 指定要保留的 Seurat cluster 編號 --------------------
keep_clusters <- c(0, 1, 2, 3, 4, 6, 9, 10)   # 想保留的群；用向量列出

## 2. 取得符合條件的細胞 barcodes ---------------------------
cells_use <- colnames(cds)[                    # 取出所有細胞名稱（欄名）
  colData(cds)$seurat_clusters %in%            # 檢查該細胞的 cluster
    keep_clusters                              # 若在 keep_clusters 之中 → TRUE
]                                              # 產生布林向量後回傳符合者

## 3. 建立子集 CellDataSet -------------------------------
cds <- cds[, cells_use]                 # 只保留篩選出來的細胞

## 4. （選擇性）檢查結果 -------------------------------
table(colData(cds)$seurat_clusters)      # 應只出現 0,1,2,3,4,6,9,10
plot_cells(cds, color_cells_by = "seurat_clusters")


############################################################
##  重新計算 cluster 與 partition，再學習 principal graph ##
############################################################

cds <- cluster_cells(              # 1️⃣ 重新計算 k-NN → Leiden → partition
  cds,
  reduction_method = "UMAP",       # 與當前 UMAP embedding 一致
  resolution = 1e-3                # 視需要調整；只是為了產生 partition
)

cds <- learn_graph(                # 2️⃣ 現在 partitions 長度吻合 → OK
  cds,
  use_partition = TRUE,            # 預設；確保不同 partition 不互連
  close_loop    = FALSE            # 避免額外閉環
)



################################################################################
## 📈 進行 Monocle3 的 graph 重建與 Pseudotime 計算
################################################################################
# cds <- learn_graph(cds)     # 建構細胞之間的拓樸結構
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE) # 建構細胞之間的拓樸結構
cds <- order_cells(cds)     # 排定 pseudotime（可互動式選擇 root cell）


################################################################################
## 🎨 Pseudotime + Cluster/Group 表現圖
################################################################################

# ➤ pseudotime 著色
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE)

# ➤ Seurat Cluster 著色
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "seurat_clusters",
           label_cell_groups = TRUE,
           group_label_size = 5,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           show_trajectory_graph = TRUE)

# ➤ 樣本來源 orig.ident 著色（例如不同患者）
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "orig.ident",
           label_cell_groups = TRUE,
           group_label_size = 5,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           show_trajectory_graph = TRUE)




####################################################################################################


################################################################################
## 🧬 為基因繪圖做準備（補上 gene_short_name）
################################################################################
rowData(cds)$gene_short_name <- rownames(cds)   # 若 rowData 尚未有 gene_short_name，需建立

################################################################################
## 📦 整理關鍵 KC 基因，分模組記錄於 ciliated_genes
################################################################################
ciliated_genes <- list(
  ## ▶ Basal KC
    "COL17A1","KRT15","KRT14","KRT5","DST","CDH3","ITGB1","ITGB4",
    "TP63","POSTN","CXCL14","S100A2","SYT8","CYR61",
    
  ## ▶ Spinous KC
    "KRT1","KRT10","KRTDAP","KRT6A","KRT6B","KRT6C","KRT16",
    "DSG1","CDH1","SBSN","LY6D",
  
  ## ▶ Granular KC
    "LOR","FLG","SPINK5","CDSN","DSC1","SLURP1","KLK7","IVL",
    "KRT1","KRT10","TGM3","FLG2","C10orf68","H0PX","CNFN",
    "CALML5","KRT2"
) %>% unique()

# ciliated_genes <- c("KRT15", "KRT14", "POSTN", "CXCL14", "S100A2", "KRT1", "KRT10")
# ➤ 畫出這些基因在 UMAP 上的分佈情形（可加 show_trajectory_graph = TRUE）
plot_cells(cds,
           genes = ciliated_genes,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)



####################################################################################################

################################################################################
## 🧪 測試基因集：任意子集的 pseudotime 表現
################################################################################
genes_Test <-  c("KRT15", "KRT14", "POSTN", "CXCL14", "S100A2", "KRT1", "KRT10")  # ← "TGGFB1" 更正為 "TGFB1"
Test_lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes_Test, ]
Test_lineage_cds <- order_cells(Test_lineage_cds)

plot_genes_in_pseudotime(Test_lineage_cds,
                         cell_size = 3,
                         color_cells_by = "seurat_clusters",
                         min_expr = NULL)

plot_genes_in_pseudotime(
  Test_lineage_cds,
  color_cells_by = "seurat_clusters",
  min_expr       = 0,           # 或小門檻
  vertical_jitter= 0.05,
  cell_size      = 0.8
)


####################################################################################################


#### 改用ggplot作圖 ####
###############################################################################
##  0. 套件 --------------------------------------------------------------------
###############################################################################
library(monocle3)      # pseudotime()、exprs()
library(tidyverse)     # tibble / dplyr / tidyr / ggplot2

###############################################################################
##  1. 準備 meta 資訊（細胞層級）----------------------------------------------
###############################################################################
meta_df <- tibble(
  cell           = colnames(Test_lineage_cds),                # 細胞名稱
  pseudotime     = pseudotime(Test_lineage_cds),              # 直接抓向量
  seurat_cluster = colData(Test_lineage_cds)$seurat_clusters  # Seurat 分群
) %>%
  filter(!is.na(pseudotime))                                  # 可選：拿掉 NA root

###############################################################################
##  2. 整理表達矩陣 → long format ---------------------------------------------
###############################################################################
expr_long <- exprs(Test_lineage_cds)[genes_Test, ] |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  pivot_longer(
    -gene,
    names_to  = "cell",
    values_to = "expr"
  )

###############################################################################
##  3. 合併 meta ＆ 表達 -------------------------------------------------------
###############################################################################
plot_df <- expr_long |>
  inner_join(meta_df, by = "cell") |>
  mutate(
    gene = factor(gene, levels = genes_Test)   # 控制 facet 順序
  )

###############################################################################
##  4. ggplot 畫圖 -------------------------------------------------------------
###############################################################################
ggplot(plot_df,
       aes(x      = pseudotime,
           y      = expr + 1e-6,               # 避免 log10(0)
           colour = factor(seurat_cluster))) +
  geom_point(size = 0.6, alpha = 0.65) +       # 散點
  geom_smooth(aes(group = 1),                  # 每 facet 一條趨勢線
              method    = "loess",
              span      = 0.8,
              colour    = "black",
              linewidth = 0.5,
              se        = FALSE) +
  scale_y_log10() +                            # 與 monocle3 相同的 log10
  scale_colour_brewer(palette = "Set1",
                      name    = "seurat_clusters") +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(strip.text.y     = element_text(angle = 0, hjust = 0),
        legend.position  = "right") +
  labs(x = "pseudotime", y = "Expression (log10)")

###############################################################################
##  5. 可選參數 ----------------------------------------------------------------
# - 想改線性軸：把 scale_y_log10() 換成 scale_y_continuous()
# - 想要 vertical_jitter: 在 geom_point 加 position = position_jitter(width = 0, height = 0.05)
###############################################################################

