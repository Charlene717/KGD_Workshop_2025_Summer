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
    s.features   = s.genes,
    g2m.features = g2m.genes,
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
    DimPlot(seurat_all_merge, reduction = "umap", group.by = "orig.ident")
}

###############################################################################
# 3. 樣本分割後再整合 (Integration) --------------------------------------------
###############################################################################
#### Data Integration ####
# SplitObject() 依 metadata$orig.ident 把合併後物件拆回多個子物件
seurat_list <- SplitObject(seurat_all_merge, split.by = "orig.ident")

# -------- 3-1 每個樣本各自 Normalize / VariableFeature / Scale / CellCycle --
# 注意：此處「先 scale、後回歸」是可靠的做法；若資料量大可考慮只 scale HVGs。
seurat_list <- lapply(
  X   = seurat_list,
  FUN = function(x) {
    x <- NormalizeData(x,
                       normalization.method = "LogNormalize",
                       scale.factor         = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x, features = rownames(x), verbose = TRUE)
    x <- CellCycleScoring(x,
                          s.features   = s.genes,
                          g2m.features = g2m.genes,
                          set.ident    = FALSE)
    return(x)
  }
)

# -------- 3-2 SelectIntegrationFeatures / FindIntegrationAnchors -------------
# 挑選跨樣本的整合特徵 (HVG union) → 建立「錨點」。
features <- SelectIntegrationFeatures(seurat_list)   # ← 預設 top 2k
anchors  <- FindIntegrationAnchors(
  seurat_list,
  anchor.features = features,
  dims = 1:30                                         # ← 與後續 IntegrateData 對齊
)

# -------- 3-3 IntegrateData ---------------------------------------------------
# 基於 anchors 做 batch-correction，生成「integrated」Assay。
seurat_all_integrated <- IntegrateData(anchors, dims = 1:30)

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
  DimPlot(seurat_all_integrated, reduction = "umap", group.by = "orig.ident")

# 若有開啟 Demo，可同時把「合併/整合前後」放在一張圖方便比較
if (Run_Demo_Merge) {
  DimPlot(seurat_all_merge, reduction = "umap") +
    DimPlot(seurat_all_merge, reduction = "umap", group.by = "orig.ident") +
    DimPlot(seurat_all_integrated, reduction = "umap") +
    DimPlot(seurat_all_integrated, reduction = "umap", group.by = "orig.ident")
}

###############################################################################
# 5. 輸出整合後的 Seurat 物件 --------------------------------------------------
###############################################################################
# saveRDS()：僅壓縮單一 R 物件；後續可用 readRDS() 直接載入。
saveRDS(
  seurat_all_integrated,
  file = paste0(Name_ExportFolder, "/", Name_Export, "_Integration.rds")
)
