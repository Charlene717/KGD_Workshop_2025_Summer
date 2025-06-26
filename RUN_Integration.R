###############################################################################
# 0. 前置作業 ── 套件與資料 ----------------------------------------------------
###############################################################################
# (此段按需要自行載入 Seurat 與其他相依套件，故僅留空白。)

###############################################################################
# 1. 合併多個單細胞樣本 (Merge) ------------------------------------------------
###############################################################################
#### Merge multiple samples ####
# merge() 會在「細胞 barcodes」前自動加上 add.cell.ids 指定的前綴，
# 以免不同樣本的 barcode 重複；同時把 metadata 與 Assay 內容一併整合。
seurat_all_merge <- merge(
  seurat_GSM6111844,
  y            = list(seurat_GSM6111845, seurat_GSM6111847),   # ← 其餘樣本
  add.cell.ids = c("GSM6111844", "GSM6111845", "GSM6111847")   # ← 各樣本專屬前綴
)

# 合併層（僅限 Seurat ≥ 5）
# 在 Seurat v5 之後，Read10X() 或其他輸入可能存成「多層 (Layer)」結構
# （counts.X / data.X 等），JoinLayers() 會把同名層合併成單一矩陣。
if (packageVersion("Seurat") >= "5.0.0") {
  seurat_all_merge <- JoinLayers(seurat_all_merge)             # ← 建議合併後立即執行
}

###############################################################################
# 2. Demo：直接在「合併後物件」上做 QC + 降維 + 分群 ----------------------------
###############################################################################
Run_Demo_Merge <- TRUE
if (Run_Demo_Merge) {
  
  # -------- 2-1 NormalizeData -------------------------------------------------
  # LogNormalize：每細胞先除以 total counts（此處 scale.factor = 1e4），
  # 再乘以 scale.factor 並取 log1p；可校正捕獲深度差異。
  seurat_all_merge <- NormalizeData(
    seurat_all_merge,
    normalization.method = "LogNormalize",
    scale.factor         = 10000
  )
  
  # -------- 2-2 FindVariableFeatures -----------------------------------------
  # VST 法找出跨細胞變異度最高的 2,000 基因，供後續 PCA 使用。
  seurat_all_merge <- FindVariableFeatures(
    seurat_all_merge,
    selection.method = "vst",
    nfeatures        = 2000
  )
  
  # -------- 2-3 ScaleData (首次) ---------------------------------------------
  # 對所有基因進行 Z-score 標準化 (mean = 0, sd = 1)；
  # 此階段 **不** 做迴歸，只是單純 scale，方便之後展示「未回歸」的結果。
  seurat_all_merge <- ScaleData(
    seurat_all_merge,
    features = rownames(seurat_all_merge),   # ← 全基因
    verbose  = TRUE
  )
  
  # -------- 2-4 CellCycleScoring ---------------------------------------------
  # 依預載的 s.genes / g2m.genes 打分，並將結果存入 meta.data$S.Score、$G2M.Score。
  seurat_all_merge <- CellCycleScoring(
    seurat_all_merge,
    s.features   = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes,
    set.ident    = FALSE                       # ← 不改變物件的 active.ident
  )
  
  # -------- 2-5 ScaleData (第二次) with vars.to.regress -----------------------
  # 以線性模型回歸掉 S.Score、G2M.Score 與粒線體讀數百分比，並重新 Z-scale。
  seurat_all_merge <- ScaleData(
    seurat_all_merge,
    vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
  )
  
  # -------- 2-6 RunPCA / RunUMAP ---------------------------------------------
  # PCA → UMAP：npcs = 50 先計算 50 個 PC，UMAP 只取前 30 個維度。
  seurat_all_merge <- RunPCA(seurat_all_merge, npcs = 50) %>%
    RunUMAP(dims = 1:30)
  
  # -------- 2-7 FindNeighbors / FindClusters ----------------------------------
  # 以 PCA 前 30 維建 KNN (Shared Nearest Neighbor) 圖，再用 Louvain/SLM 分群。
  seurat_all_merge <- FindNeighbors(seurat_all_merge, dims = 1:30) %>%
    FindClusters(resolution = 0.3)               # ← 解析度可酌調
  
  # -------- 2-8 可視化 --------------------------------------------------------
  # (A) 全樣本 UMAP；(B) 依 orig.ident 著色，檢視批次效應。
  DimPlot(seurat_all_merge, reduction = "umap") +
    DimPlot(seurat_all_merge, reduction = "umap", group.by = "orig.ident") %>% print()
}

###############################################################################
# 3. 樣本分割後再整合 (Integration) --------------------------------------------
###############################################################################
#### Data Integration ####

#### 3-0 先依 orig.ident 拆回各自 Seurat 物件 ####
# SplitObject() 只是「把大物件切片」，不會改變 data/matrix 內容。
seurat_list <- SplitObject(
  seurat_all_merge,
  split.by = "orig.ident"          # ← 依樣本來源欄位切分，可改為其他欄位(看要對什麼做批次效應校正)
)

#### 3-1 每個樣本單獨前處理 (Normalize → HVG → Scale → CellCycle) ####
# -- 此處仍採 LogNormalize 與 VST，保持與前面 Merge Demo 一致。
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_list <- lapply(
  X = seurat_list,
  FUN = function(x) {
    x <- NormalizeData(
      x,
      normalization.method = "LogNormalize",   # ← 同上一段
      scale.factor         = 10000             # ← 每細胞歸一化目標總計數
    )
    x <- FindVariableFeatures(
      x,
      selection.method = "vst",                # ← 變異度穩定轉換法
      nfeatures        = 2000                  # ← 固定 2,000 基因，不做修改
    )
    x <- ScaleData(
      x,
      features = rownames(x),                  # ← 全基因 Z-score
      verbose  = TRUE
    )
    x <- CellCycleScoring(
      x,
      s.features   = s.genes,                  # ← 預載 S 期基因
      g2m.features = g2m.genes,                # ← 預載 G2/M 期基因
      set.ident    = FALSE                     # ← 不改 active.ident
    )
    return(x)                                  # ← 回傳處理完成的子物件
  }
)

#### 3-2 選特徵 → 建錨點 ####
# (1) SelectIntegrationFeatures()：跨樣本挑「共同且高變異」基因
features <- SelectIntegrationFeatures(
  seurat_list                      # ← 傳入物件清單
  # nfeatures 預設 2000 → 沒明示 = 使用預設值
)

# (2) FindIntegrationAnchors()：計算「錨點」，用作批次校正基礎
anchors <- FindIntegrationAnchors(
  object.list     = seurat_list,   # ← 放入剛才處理好的列表
  anchor.features = features,      # ← 用上一步篩出的特徵向量
  dims            = 1:30           # ← 取 PCA 前 30 維作相似度空間
  # reduction     = "cca"          # ← 預設 Canonical Correlation；v5 會自動選
  # k.anchor      = 5             # ← 近鄰數；保留預設以維持原本設定
)

#### 3-3 整合 (IntegrateData) ####
# 將所有子物件依 anchors 計算的轉換矩陣整合成「integrated」Assay
seurat_all_integrated <- IntegrateData(
  anchorset = anchors,             # ← 必填；用剛產出的 Anchorset
  dims      = 1:30                 # ← 與上一步保持一致
  # features / features.to.integrate 省略 → 保留 anchors 內建資訊與預設行為
)

# ★ 執行完畢後：Seurat 物件內會新增一個名為 "seurat_all_integrated" 的 Assay，
#   其 counts/data/scale.data 均已做 batch-correction，可直接用於後續降維。
###############################################################################


###############################################################################
# 4. 降維、分群 (整合後) -------------------------------------------------------
###############################################################################
DefaultAssay(seurat_all_integrated) <- "integrated"  # ← 後續 PCA / UMAP 用整合表達

# *再次* ScaleData，以 integrated 資料為基礎迴歸掉 cell cycle & percent.mt
seurat_all_integrated <- ScaleData(
  seurat_all_integrated,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
)

# PCA → UMAP → KNN → Cluster
seurat_all_integrated <- RunPCA(seurat_all_integrated, npcs = 50) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.3)

# 整合後的 UMAP：檢查不同批次是否成功重疊
DimPlot(seurat_all_integrated, reduction = "umap") +
  DimPlot(seurat_all_integrated, reduction = "umap", group.by = "orig.ident") %>% print()

# 若有開啟 Demo，可同時把「合併/整合前後」放在一張圖方便比較
if (Run_Demo_Merge) {
  DimPlot(seurat_all_merge, reduction = "umap") +
    DimPlot(seurat_all_merge, reduction = "umap", group.by = "orig.ident") +
    DimPlot(seurat_all_integrated, reduction = "umap") +
    DimPlot(seurat_all_integrated, reduction = "umap", group.by = "orig.ident") %>% print()
}

###############################################################################
# 5. 輸出整合後的 Seurat 物件 --------------------------------------------------
###############################################################################
# saveRDS()：僅壓縮單一 R 物件；後續可用 readRDS() 直接載入。
saveRDS(
  seurat_all_integrated,
  file = paste0(Name_ExportFolder, "/", Name_Export, "_Integration.rds")
)
