##########################################################
## 讀取 Seurat 物件（示範路徑，先註解起來以免誤執行）
##########################################################
# seurat_all_integrated <- readRDS(
#   "C:/Charlene/Code_GitHub_BioInport2025/KGD_Workshop_2025_Summer/Export_2025062617YDG_Demo/2025062617YDG_Integration.rds"
# )  # ← 從硬碟載入先前整合好的 Seurat 物件，包含所有細胞表達與 metadata


##########################################################
## Step 1：安裝與載入套件
##########################################################
if (!require("Seurat"))        { install.packages("Seurat");         library(Seurat) }     # ↪ 若尚未安裝 Seurat 就自動安裝；安裝後載入
if (!require("tidyverse"))     { install.packages("tidyverse");      library(tidyverse) }  # ↪ tidyverse: dplyr/ggplot2 等整合工具，方便資料處理
if (!require("ggplot2"))       { install.packages("ggplot2");        library(ggplot2) }    # ↪ 畫圖核心：CellChat 某些函式需依賴 ggplot2
if (!require("patchwork"))     { install.packages("patchwork");      library(patchwork) }  # ↪ 組合多圖用，常見於 Seurat/CellChat 視覺化
if (!require("future"))        { install.packages("future");         library(future) }     # ↪ 提供平行運算 (multisession) ，加速大型計算

if (!require("remotes"))       { install.packages("remotes");        library(remotes) }    # ↪ 從 GitHub 安裝套件的輔助工具
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat")                           # ↪ 從 GitHub 安裝最新版 CellChat
  library(CellChat)                                                   # ↪ 載入 CellChat 套件
}


## 可選：設定平行運算參數以提速 (根據硬體調整) --------------------------
plan("multisession", workers = 4)                 # ↪ 使用 4 個本地 CPU 執行緒並行；Windows/macOS 無需額外套件
options(future.globals.maxSize = 8 * 1024^3)      # ↪ 允許 future 物件佔用最大 8 GB RAM，避免大型矩陣拋錯


##########################################################
## Step 2：準備 Seurat 資料
##########################################################


#Bug# seurat_all_integrated <- JoinLayers(object = seurat_all_integrated, layers = "data")

# ✅ 若使用 Seurat v5、RNA assay 會預設多層 (layers) 儲存；先合併到單層 “data”
# 將 RNA assay 的多層合併為一層（預設合併到 "data" slot）
seurat_all_integrated[["RNA"]] <- JoinLayers(object = seurat_all_integrated[["RNA"]])  # ↪ 將 raw / normalized / scaled 圖層合併，留方便給 CellChat 用的 data slot

# 取出「gene × cell」表現矩陣 (log-normalized counts) -------------------
data.input <- GetAssayData(                       # ↪ 從 Seurat 物件中取資料
  object = seurat_all_integrated,
  slot   = "data",                                # ↪ 取 log-normalized 表達值；CellChat 以此估算配體/受體活性
  assay  = "RNA"
)

# 取出細胞 metadata（每細胞一列，包括分群資訊） --------------------------
meta <- seurat_all_integrated@meta.data           # ↪ 用於後續指定群組 (group.by)

## 建議：若群組欄位名稱含特殊字元，先轉成簡潔 ID -----------------------
# meta$cellchat_clusters <- paste0("C", seurat_all_integrated$seurat_clusters)  # ↪ 範例：把 0/1/2… 變成 C0/C1/C2…
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_clusters")

## 在此示範以 “Cell_Type” 欄位作群組 -------------------------------
cellchat <- createCellChat(                       # ↪ 建立 CellChat 物件
  object  = data.input,                           # ↪ 表達矩陣
  meta    = meta,                               
  group.by = "Cell_Type"                          # ↪ 依照 Cell_Type 欄位 (如 Basal, Spinous, Fibroblast…)
)


##########################################################
## Step 3：指定物種資料庫 (human / mouse)
##########################################################
CellChatDB <- CellChatDB.human                    # ↪ 人類資料庫；若是 Mus musculus 改 CellChatDB.mouse
cellchat@DB <- CellChatDB                         # ↪ 將資料庫掛到 cellchat 物件
showDatabaseCategory(CellChatDB)                  # ↪ 列出資料庫中的互動分類 (Secreted, ECM, Cell–Cell Contact)


##########################################################
## Step 4：預處理與基因篩選
##########################################################

cellchat <- subsetData(cellchat)                  # ↪ 只保留資料庫內存在的配體/受體基因，降低矩陣維度
cellchat <- identifyOverExpressedGenes(cellchat)  # ↪ 依群組計算 Z-score，挑高表達基因
cellchat <- identifyOverExpressedInteractions(cellchat)  # ↪ 結合上一步基因表達結果，鎖定可能活化的 LR 對

# （選擇性）投影到 PPI 網絡以納入複合體資料 -------------------------
# cellchat <- projectData(cellchat, PPI.human)    # ↪ 用人類 PPI 網絡加強複合體關係 (mouse 可省略)


##########################################################
## Step 5：推論細胞間通訊網絡
##########################################################

cellchat <- computeCommunProb(cellchat)           # ↪ 依平均表達量 × scale factor → 算每對 cell type 間 LR 機率
cellchat <- filterCommunication(cellchat, min.cells = 10)  # ↪ 至少有 10 個細胞的群組才保留互動，避免統計不穩
df.net   <- subsetCommunication(cellchat)         # ↪ 匯出資料框查看所有顯著 LR 對 (可寫 csv)

cellchat <- computeCommunProbPathway(cellchat)    # ↪ 進一步將 LR 對匯總到訊號 pathway 層級
cellchat <- aggregateNet(cellchat)                # ↪ 將多重 LR → 單一 pathway 網絡，方便後續視覺化


##########################################################
## Step 6：網絡可視化
##########################################################

groupSize <- as.numeric(table(cellchat@idents))   # ↪ 每個 cell type 的細胞數，用於節點大小比例

netVisual_circle(                                 # ↪ 圓環網絡：粗線 = 互動次數；色塊 = cell type
  cellchat@net$count,
  vertex.weight = groupSize,                      # ↪ 節點大小依細胞數
  weight.scale  = TRUE,                           # ↪ 線粗-節點大小皆自動縮放
  label.edge    = FALSE                           # ↪ 不顯示邊標籤
)

## 指定 pathway “GAP” 多圖示範 ---------------------------------------
netVisual_aggregate(cellchat, signaling = "GAP", layout = "circle") # ↪ 圓環
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = "GAP", layout = "chord")  # ↪ Chord diagram
par(mfrow = c(1,1))
netVisual_heatmap(cellchat, signaling = "GAP", color.heatmap = "Reds") # ↪ 熱圖

netAnalysis_contribution(cellchat, signaling = "GAP")                # ↪ 各 LR 對對此 pathway 貢獻百分比
netVisual_bubble(                                                    # ↪ 氣泡圖，單向 Granular → 全細胞
  cellchat,
  sources.use   = "Granular keratinocytes",
  targets.use   = levels(cellchat@idents),
  remove.isolate = FALSE
)


##########################################################
## Step 7：中心性分析 (辨識發送者 / 接收者)
##########################################################
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # ↪ 算 betweenness、in/out degree…

netAnalysis_signalingRole_network(                                     # ↪ 在 pathway 網絡中標記主要 sender/receiver
  cellchat,
  signaling = "GAP"
)


##########################################################
## Step 8：通訊角色熱圖 (整體或分向)
##########################################################
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all")           # ↪ Sender + Receiver 混合視角
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") +    # ↪ 只看發送
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")    # ↪ 只看接收


##########################################################
## Step 9：儲存 / 讀取 CellChat 物件
##########################################################
saveRDS(cellchat, file = "cellchat_result.rds")  # ↪ 將分析全貌 (含網絡/圖) 序列化保存

# 之後可透過 readRDS() 直接載回：
# cellchat <- readRDS("cellchat_result.rds")
