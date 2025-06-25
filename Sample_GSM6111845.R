## =======================================================================
## Script: GSM6111845_scRNAseq_preprocessing_annotated.R
## Description: End‑to‑end preprocessing pipeline for single‑cell RNA‑seq
##              sample GSM6111845 (10× Genomics v3 chemistry).  每一步都
##              以中文註解說明背後邏輯與參數設定，方便後續教學與維護。
## =======================================================================

# --------------------------
# 0. 環境準備 ─ 讀取套件
# --------------------------
# ▸ Seurat:  單細胞 RNA‑seq 分析的主力工具
# ▸ dplyr:    資料處理 & 管線語法 ("%>%")
# ▸ ggplot2:  繪圖系統，用於可視化 QC 與降維結果
# 套件若尚未安裝，可先執行 install.packages("<pkg>") 再 library(<pkg>)
# -----------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------
# 1. 載入 10× 格式的矩陣
# --------------------------
# ▸ Read10X() 會一次讀取 barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz
#   並回傳稀疏矩陣 (genes × cells)。
# ▸ 這裡直接給定 Cell Ranger output 之 filtered_feature_bc_matrix 資料夾
# -----------------------------------------------------------------------
# <※請依實際路徑調整※>
data_dir <- "C:/Charlene/Dataset_KGD_Lab/scRNA-seq/10x/sample_filtered_feature_bc_matrix/GSM6111845_Normal_Sole/"

data_GSM6111845 <- Read10X(data.dir = data_dir)

# --------------------------
# 2. 建立 Seurat 物件
# --------------------------
# ▸ min.cells    : 至少出現在 ≥3 個細胞才保留該基因，避免低複雜度雜訊。
# ▸ min.features : 單一細胞須表現 ≥200 個基因，濾除空液滴或低品質細胞。
# -----------------------------------------------------------------------
seurat_GSM6111845 <- CreateSeuratObject(
  counts       = data_GSM6111845,
  project      = "GSM6111845",
  min.cells    = 3,
  min.features = 200
)

# --------------------------
# 3. 品質控制 (QC)
# --------------------------
# 3A. 計算粒線體基因比例 (percent.mt)
#     ▸ ^MT- 用於人類資料；若是小鼠請改成 ^Mt-*
# 3B. 依下列條件過濾細胞：
#     ▸ nFeature_RNA: 200–5,000 之間 (基因數過少可能是破裂細胞；過多可能是雙細胞)
#     ▸ percent.mt : <30% (粒線體比例過高表示細胞壓力或破裂)
# -----------------------------------------------------------------------
seurat_GSM6111845[["percent.mt"]] <- PercentageFeatureSet(seurat_GSM6111845, pattern = "^MT-")

seurat_GSM6111845 <- subset(
  seurat_GSM6111845,
  subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30
)

# --------------------------
# 4. 引入 DoubletFinder 套件 (偵測雙細胞)
# --------------------------
# ▸ 安裝： remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
# -----------------------------------------------------------------------
library(DoubletFinder)

# --------------------------
# 5. Seurat 標準前處理流程
# --------------------------
# 5A. NormalizeData : LogNormalize → CPM × 1e4 後取 log1p
# 5B. FindVariableFeatures : VST 方法選 2,000 高變基因
# 5C. ScaleData : Z‑score normalization (可指定 features；此處全基因)
# 5D. RunPCA : PCA 降維 (自動選前 50 PCs；後續取 1:30)
# 5E. ElbowPlot : 視覺化變異解釋量，決定合理 PC 數
# 5F. FindNeighbors / FindClusters : 建圖與社群偵測 (預設 Louvain)
# 5G. RunUMAP : 非線性降維作視覺化
# -----------------------------------------------------------------------
seurat_GSM6111845 <- NormalizeData(
  seurat_GSM6111845,
  normalization.method = "LogNormalize",
  scale.factor         = 1e4
)

seurat_GSM6111845 <- FindVariableFeatures(
  seurat_GSM6111845,
  selection.method = "vst",
  nfeatures        = 2000
)

seurat_GSM6111845 <- ScaleData(
  seurat_GSM6111845,
  features = rownames(seurat_GSM6111845),
  verbose  = TRUE
)

seurat_GSM6111845 <- RunPCA(seurat_GSM6111845)

# 觀察肘部圖，以挑選可解釋較多變異的 PC 數
ElbowPlot(seurat_GSM6111845)

# 這裡先使用 30 個 PC，後續視情況微調
seurat_GSM6111845 <- FindNeighbors(seurat_GSM6111845, dims = 1:30)
seurat_GSM6111845 <- FindClusters(seurat_GSM6111845)
seurat_GSM6111845 <- RunUMAP(seurat_GSM6111845, dims = 1:30)

# 標註叢集編號於 UMAP
DimPlot(seurat_GSM6111845, reduction = "umap", label = TRUE)

# --------------------------
# 6. 使用 DoubletFinder 預測雙細胞
# --------------------------
# DoubletFinder 需要三個重要參數：
# ▸ pN:   人工雙細胞產生比例 (預設 0.25 即 25%)
# ▸ pK:   最近鄰數量正規化參數 — 需掃描找最佳值 (下方示範)
# ▸ nExp: 預期雙細胞數量 ＝ 全細胞數 × 約略雙細胞比例 × (1 − homotypic.prop)
#         homotypic.prop 估計同型 (同叢集) 雙細胞佔比，可降低過度濾除風險
# -----------------------------------------------------------------------

## 6A. 尋找最佳 pK 值 (常見 0.005–0.3)
sweep_res   <- paramSweep(seurat_GSM6111845, PCs = 1:30, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res)

best_pk <- find.pK(sweep_stats)
# 取具有最高 BCmetric 的 pK — 表示分群品質最佳
pK <- as.numeric(as.character(best_pk[which.max(best_pk$BCmetric), "pK"]))

## 6B. 計算預期雙細胞數量 (nExp)
annotations     <- seurat_GSM6111845$seurat_clusters
homotypic.prop  <- modelHomotypic(annotations)        # 模擬同型雙細胞比例
nExp            <- round(0.08 * nrow(seurat_GSM6111845@meta.data) * (1 - homotypic.prop))
# 0.08 (~8%) 為 10× 官方建議的人類樣本雙細胞比例，可依實驗加以微調

## 6C. 執行 DoubletFinder
seurat_GSM6111845 <- doubletFinder(
  seurat_GSM6111845,
  PCs = 1:30,
  pN  = 0.25,
  pK  = pK,
  nExp = nExp,
  reuse.pANN = FALSE
)

# DoubletFinder 會在 metadata 產生多個欄位，
# 名稱格式示例： DF.classifications_0.25_0.005_150
# 以下自動尋找第一個符合 DF.classifications 開頭的欄位

df_col <- grep("^DF\\.classifications", colnames(seurat_GSM6111845@meta.data), value = TRUE)[1]

# 6D. 先檢視雙細胞於 UMAP 分布與數量概況
DimPlot(seurat_GSM6111845, reduction = "umap", group.by = df_col)

table(seurat_GSM6111845@meta.data[[df_col]])

# 6E. 只保留 Singlet 細胞進入後續分析
seurat_GSM6111845 <- seurat_GSM6111845[, seurat_GSM6111845@meta.data[[df_col]] == "Singlet"]

# 再次確認 UMAP 是否仍合理 (叢集形狀、密度)
DimPlot(seurat_GSM6111845, reduction = "umap")

# --------------------------
# 7. 結果概覽
# --------------------------
# 輸出 Seurat 物件摘要，包含細胞數、基因數、metadata 欄位等
# -----------------------------------------------------------------------
seurat_GSM6111845

# (可選) 將物件儲存為 .rds 供後續流程載入
# saveRDS(seurat_GSM6111845, file = "seurat_GSM6111845_singlets.rds")

## =========================  End of Script  ==============================
