###############################################################################
## 0. 套件與物件準備 -----------------------------------------------------------
###############################################################################
library(Seurat)               # ≥ 5.0.0
# seurat_obj 需已經載入／建立，並且 RNA assay 為 DefaultAssay()

###############################################################################
## 1. 定義三組 marker 基因 ------------------------------------------------------
###############################################################################
marker_list <- list(
  Basal = c(
    "COL17A1","KRT15","KRT14","KRT5","DST","CDH3","ITGB1","ITGB4",
    "TP63","POSTN","CXCL14","S100A2","SYT8","CYR61"
  ),
  Spinous = c(
    "KRT1","KRT10","KRTDAP","KRT6A","KRT6B","KRT6C","KRT16",
    "DSG1","CDH1","SBSN","LY6D"
  ),
  Granular = c(
    "LOR","FLG","SPINK5","CDSN","DSC1","SLURP1","KLK7","IVL",
    "KRT1","KRT10","TGM3","FLG2","C10orf68","H0PX","CNFN",
    "CALML5","KRT2"
  )
)

## 若擔心有些基因不在資料中，可先過濾 ──
marker_list <- lapply(marker_list, \(g) g[g %in% rownames(seurat_obj)])

###############################################################################
## 2. 依序計算 Module Score -----------------------------------------------------
###############################################################################
for (nm in names(marker_list)) {
  seurat_obj <- AddModuleScore(
    object   = seurat_obj,
    features = list(marker_list[[nm]]),  # 必須包在 list() 裡
    name     = paste0(nm, "_Score"),     # 產生的 meta 資料欄前綴
    assay    = "RNA",                    # 如有多重 assay 記得指定
    search   = TRUE                      # 若有 Ensembl ID 也能自動對應 SYMBOL
  )
}

###############################################################################
## 3. 檢查結果 -----------------------------------------------------------------
###############################################################################
# 產生的欄位名稱格式：<前綴>1，例如 Basal_Score1
grep("_Score1$", colnames(seurat_obj@meta.data), value = TRUE)
#> [1] "Basal_Score1"    "Spinous_Score1"  "Granular_Score1"

head(seurat_obj@meta.data[, c("Basal_Score1",
                                       "Spinous_Score1",
                                       "Granular_Score1")])


###############################################################################
## 1. FeaturePlot ‒ 各畫一張 -----------------------------------------------------
###############################################################################
library(Seurat)      # ≥ 5.0
library(patchwork)   # 排版用，可選

p1 <- FeaturePlot(
  object     = seurat_obj,
  features   = "Basal_Score1",
  reduction  = "umap",
  slot       = "data",          # 預設即可
  pt.size    = 0.8
) + ggtitle("Basal keratinocyte score")

p2 <- FeaturePlot(
  seurat_obj, "Spinous_Score1", reduction = "umap", pt.size = 0.8
) + ggtitle("Spinous keratinocyte score")

p3 <- FeaturePlot(
  seurat_obj, "Granular_Score1", reduction = "umap", pt.size = 0.8
) + ggtitle("Granular keratinocyte score")

# -- 直接顯示；或另存 PDF / PNG ----------------------------------------
p1 | p2 | p3                         # 三張排成一列（需要 patchwork）
# ggsave("Keratinocyte_ModuleScore_UMAP.pdf", width = 12, height = 4)
