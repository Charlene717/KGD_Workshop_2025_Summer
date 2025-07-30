###############################################################
# 📘 單細胞軌跡分析教學：使用 Seurat + Monocle2 建立 Pseudotime
# 作者：ChatGPT 改寫，適合初學者參考
###############################################################

## -----------------------------------------------
## 1. 安裝與載入必要套件
## -----------------------------------------------

# 安裝 BiocManager（如尚未安裝）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 安裝 monocle（Monocle2）套件
if (!requireNamespace("monocle", quietly = TRUE)) {
  BiocManager::install("monocle")
}

# 安裝其餘所需套件（如尚未安裝）
packages <- c("Seurat", "reshape", "ggplot2", "RColorBrewer")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# 設定語系為英文，避免中文錯誤訊息
Sys.setenv(LANGUAGE = 'en')

# 載入套件
library(Seurat)
library(monocle)
library(reshape)
library(ggplot2)
library(RColorBrewer)


## -----------------------------------------------
## 2. 讀取 Seurat 物件（請修改為自己的路徑）
## -----------------------------------------------

# 請將以下路徑改成你本地的 rds 檔案路徑
seurat_obj_path <- "your_path/PD.combined.rds"
PD.combined <- readRDS(file = seurat_obj_path)


## -----------------------------------------------
## 3. 選取 Macrophage 或 Fibroblast 子群
## -----------------------------------------------

# 以 Fibroblast 為例（假設 cluster 0, 9, 20 為 Fibroblast）
FB.subgroup <- subset(PD.combined, idents = c(0, 9, 20))

# 使用 UMAP 可視化 Fibroblast 分群
DimPlot(FB.subgroup, reduction = "umap", label = TRUE, pt.size = 0.8)


## -----------------------------------------------
## 4. 轉換 Seurat object 為 Monocle2 的 CellDataSet
## -----------------------------------------------

# 將 Seurat 中的 RNA 表達矩陣轉為 sparseMatrix 格式
expr_matrix <- as(as.matrix(FB.subgroup@assays$RNA@data), "sparseMatrix")

# 建立細胞與基因註解資料
pd <- new("AnnotatedDataFrame", data = FB.subgroup@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix)))

# 建立 CellDataSet 物件
FB.cds <- newCellDataSet(expr_matrix,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.25,
                         expressionFamily = negbinomial.size())


## -----------------------------------------------
## 5. 前處理：標準化與過濾
## -----------------------------------------------

# 估計大小因子與離散度（normalization）
FB.cds <- estimateSizeFactors(FB.cds)
FB.cds <- estimateDispersions(FB.cds)

# 偵測表達的基因（需大於指定閾值）
FB.cds <- detectGenes(FB.cds, min_expr = 0.1)

# 篩選在 ≥5 個細胞中有表達的基因
expressed_genes <- rownames(subset(fData(FB.cds), num_cells_expressed >= 5))


## -----------------------------------------------
## 6. 表達值標準化與分布可視化
## -----------------------------------------------

# 計算 log 表達值與 Z-score 標準化
L <- log(exprs(FB.cds[expressed_genes,]))
mL <- apply(L, 1, function(x) mean(x[is.finite(x)]))
sdL <- apply(L, 1, function(x) sd(x[is.finite(x)]))
Lstd <- (L - mL) / sdL

# 轉換成長格式並繪製密度圖
melted_dens_df <- melt(as.matrix(Lstd))
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") + ylab("Density")


## -----------------------------------------------
## 7. 設定用於排序的基因並進行降維
## -----------------------------------------------

# 設定哪些基因將用於 Cell ordering（重要）
fData(FB.cds)$use_for_ordering <- fData(FB.cds)$num_cells_expressed > 0.05 * ncol(FB.cds)

# 初始降維分析（PCA → DDRTree）
FB.cds <- reduceDimension(FB.cds, method = "DDRTree")

# 建立 pseudotime 軌跡
FB.cds <- orderCells(FB.cds)


## -----------------------------------------------
## 8. Pseudotime 與群集視覺化
## -----------------------------------------------

# 視覺化 pseudotime 軌跡（用 seurat_clusters 著色）
plot_cell_trajectory(FB.cds, color_by = "seurat_clusters")

# 如果 meta.data 中有 orig.ident 或樣本名稱，也可以使用：
# plot_cell_trajectory(FB.cds, color_by = "orig.ident")

# 顯示 pseudotime
plot_cell_trajectory(FB.cds, color_by = "Pseudotime")


## -----------------------------------------------
## 9. 繪製特定 marker genes 在 pseudotime 上的表現
## -----------------------------------------------

# 選擇要檢視的 marker genes（可替換成自己的）
marker_genes <- c("POSTN", "COL1A1", "FN1", "CTHRC1", "ACTA2", "APOE", "APCDD1")

# 篩選 gene_short_name 符合 marker 的 subset
to_be_tested <- rownames(subset(fData(FB.cds), gene_short_name %in% marker_genes))
cds_subset <- FB.cds[to_be_tested, ]

# 使用 spline 模型檢定這些基因是否隨 pseudotime 改變
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
print(diff_test_res[, c("gene_short_name", "pval", "qval")])

# 畫出基因在 pseudotime 中的表現趨勢
plot_genes_in_pseudotime(cds_subset, 
                         color_by = "seurat_clusters", 
                         panel_order = marker_genes)

# 若資料有 orig.ident，可依樣本來源著色
plot_genes_in_pseudotime(cds_subset, 
                         color_by = "orig.ident", 
                         panel_order = marker_genes) +
  scale_color_manual(values = c("#00BA38", "#F8766D"))  # 可自訂樣本顏色

## -----------------------------------------------
## 10. 分支特異性分析：分支 pseudotime 上的基因表現
## -----------------------------------------------

# 使用 plot_genes_branched_pseudotime 觀察 marker genes 的分支表現
branched_genes <- c("IL4R", "IL13RA1", "CCL19", "COL1A1", "FN1", "POSTN", "MFAP5", "APCDD1")

# 篩選符合的基因
branched_subset <- FB.cds[rownames(subset(fData(FB.cds),
                                          gene_short_name %in% branched_genes)), ]

# 分支點必須為已存在的節點，例如 branch_point = 1
plot_genes_branched_pseudotime(branched_subset,
                               branch_point = 1,
                               branch_labels = c("Branch1", "Branch2"),
                               color_by = "orig.ident",
                               ncol = 1,
                               panel_order = branched_genes)

## -----------------------------------------------
## 11. BEAM 分析：分支差異基因的熱圖與儲存
## -----------------------------------------------

# 執行 BEAM（Branch Expression Analysis Modeling）分析
BEAM_res <- BEAM(FB.cds, branch_point = 1, cores = 1)

# 排序並顯示前幾個顯著基因
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
head(BEAM_res)

# 繪製分支基因的 heatmap，視覺化不同分支中基因表現差異
branch_genes <- plot_genes_branched_heatmap(
  FB.cds[rownames(subset(BEAM_res, qval < 0.01)), ],
  branch_point = 1,
  num_clusters = 6,
  cores = 1,
  use_gene_short_name = TRUE,
  show_rownames = TRUE,
  return_heatmap = TRUE  # 可用來儲存為變數
)

# 將結果輸出為 CSV 檔案（請自行更改路徑）
write.csv(branch_genes, file = "output/branch_genes_heatmap.csv", row.names = FALSE)

# 儲存 Monocle2 物件
saveRDS(FB.cds, file = "output/Fibroblast_monocle2.rds")
