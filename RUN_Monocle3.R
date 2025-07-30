################################################################################
## 🔧 載入必要套件
################################################################################
library(Seurat)          # 單細胞分析主套件
library(monocle3)        # Pseudotime 分析與軌跡重建
library(SeuratWrappers)  # Seurat 和 Monocle 等工具間的轉換橋接

################################################################################
## 📥 讀取 Seurat 物件，並設定 clustering 身分
################################################################################
seurat_path <- "path/to/your/seurat_obj.rds"
seurat_obj <- readRDS(seurat_path)             # 載入儲存好的 Seurat RDS 物件
Idents(seurat_obj) <- "seurat_clusters"        # 指定 Seurat 用來分群的欄位

################################################################################
## 🔁 將 Seurat 轉換為 Monocle3 的 CellDataSet (CDS)
################################################################################
cds <- as.cell_data_set(seurat_obj)            # Seurat ➜ Monocle 格式

# ➕ 把 metadata 合併到 colData（避免重複欄位）
meta_to_add <- seurat_obj@meta.data
meta_to_add <- meta_to_add[, !colnames(meta_to_add) %in% colnames(colData(cds))]
colData(cds) <- cbind(colData(cds), meta_to_add)

# ➕ 把原本的 Seurat UMAP 降維結果寫入 Monocle3 的 reducedDims slot
reducedDims(cds)$UMAP <- Embeddings(seurat_obj, reduction = "umap")

# ➕ 將 Seurat 分群結果複製進 Monocle3（供 plot_cells 使用）
cds@clusters$UMAP$clusters <- factor(Idents(seurat_obj))

# ⚠️ 手動指定所有細胞都屬於同一個 partition，否則後續會報錯
cds@clusters$UMAP$partitions <- factor(rep(1, length(Cells(cds))))
names(cds@clusters$UMAP$partitions) <- Cells(cds)  # 確保有正確 names

# ❗不需執行 preprocess_cds() 或 reduce_dimension()，因為已從 Seurat 導入 UMAP


####################################################################################################


################################################################################
## 📈 進行 Monocle3 的 graph 重建與 Pseudotime 計算
################################################################################
cds <- learn_graph(cds)     # 建構細胞之間的拓樸結構
# cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE) # 建構細胞之間的拓樸結構
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
           color_cells_by = "orig.ident1",
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
## 📦 整理關鍵 fibroblast 基因，分模組記錄於 ciliated_genes
################################################################################
ciliated_genes <- c(
  ## ▶ ECM 組成與膠原蛋白沉積
  "COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "VCAN", "SPARC",
  
  ## ▶ 肌成纖維細胞標誌（Myofibroblast markers）
  "ACTA2", "TAGLN", "POSTN", "LOX", "PLOD2",
  
  ## ▶ 細胞增生與訊號傳導（Growth factor signaling）
  "PDGFRB", "PDGFRA", "IGF1", "FGF2", "EGF", "TGFBI",
  
  ## ▶ TGF-β pathway
  "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
  
  ## ▶ Wnt/β-catenin pathway
  "WNT5A", "CTNNB1", "AXIN2",
  
  ## ▶ 細胞黏附與遷移相關
  "ITGB1",
  
  ## ▶ 基質重塑與降解（ECM remodeling enzymes）
  "MMP2", "MMP9", "TIMP1",
  
  ## ▶ 免疫與發炎因子
  "IL6", "IL11", "CCL2", "CXCL12", "CXCL14",
  
  ## ▶ 血管新生與環境重塑
  "VEGFA", "ANGPTL4",
  
  ## ▶ 幹細胞樣與纖維母細胞 progenitor 標誌
  "THY1", "CD34"
)


# ciliated_genes <- c("PDGFRA", "LUM", "DCN", "COL1A1", "COL3A1", "COL5A1", "COL6A3")
# ➤ 畫出這些基因在 UMAP 上的分佈情形（可加 show_trajectory_graph = TRUE）
plot_cells(cds,
           genes = ciliated_genes,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)





####################################################################################################



################################################################################
## ⏳ 特定模組的 pseudotime 表現動態：以 ECM 基因為例
################################################################################

# ➤ 提取 ECM gene 對應的 CDS
genes_ECM <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "VCAN", "SPARC")
ECM_lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes_ECM, ]
ECM_lineage_cds <- order_cells(ECM_lineage_cds)

# ➤ 畫出 ECM 模組在 pseudotime 上的表現變化
plot_genes_in_pseudotime(ECM_lineage_cds,
                         color_cells_by = "seurat_clusters",
                         min_expr = 0.5)



################################################################################
## 🧪 測試基因集：任意子集的 pseudotime 表現（debug 用）
################################################################################
genes_Test <- c("IL6", "PLOD2", "TGFB1", "COL5A1")  # ← "TGGFB1" 更正為 "TGFB1"
Test_lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes_Test, ]
Test_lineage_cds <- order_cells(Test_lineage_cds)

plot_genes_in_pseudotime(Test_lineage_cds,
                         color_cells_by = "seurat_clusters",
                         min_expr = 0.5)

