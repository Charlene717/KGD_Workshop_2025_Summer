#### 🔍 自動偵測與安裝 scRNA-seq Workflow 所需套件 ####

# 基本套件清單
packages_needed <- c(
  "Seurat", "dplyr", "ggplot2",         # Data preprocessing & plotting
  "DoubletFinder",                      # Doublet removal
  "SingleR", "celldex",                 # Cell type annotation
  "CellChat",                           # Cell-cell communication
  "monocle",                            # Pseudotime analysis
  "BiocManager"                         # 用於安裝Bioconductor套件
)

# 安裝未安裝的套件
for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing:", pkg))
    if (pkg %in% c("SingleR", "celldex", "monocle")) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  library(pkg, character.only = TRUE)
}

# ✅ 若需確認 Bioconductor 已安裝且為最新版本：
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 檢查 CellChat 資料庫是否存在（若不存在則安裝）
if (!"CellChatDB.human" %in% rownames(installed.packages())) {
  message("Installing CellChatDB.human...")
  CellChat::CellChatDB <- CellChat::CellChatDB.human
}

# 如果你用的是 monocle2，而不是 monocle3，請加入：
if (!"monocle" %in% installed.packages()[,"Package"]) {
  BiocManager::install("monocle", version = "3.16")
}

#### 📢 所有必要套件已載入完畢 ####
