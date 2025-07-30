# library(CellChat)
# 
# data.input <- GetAssayData(integrated, assay = "RNA", slot = "data")
# meta <- data.frame(group = Idents(integrated), row.names = colnames(integrated))
# cellchat <- createCellChat(data.input, meta = meta, group.by = "group")
# cellchat@DB <- CellChatDB.human
# 
# # Analysis workflow
# cellchat <- subsetData(cellchat)
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- computeCommunProb(cellchat)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# 
# # Visualize communication circle plot
# netVisual_circle(cellchat@net$count)


###############################################
## Step 1: 安裝與載入套件
###############################################

# 安裝 remotes 套件（若尚未安裝）
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 安裝 CellChat
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat")
}

# 載入必要套件
library(CellChat)
library(patchwork)
library(Seurat)
library(ggplot2)
library(future)

# 可選：加速處理（根據 CPU 核心數調整 workers）
plan("multisession", workers = 4)  
options(future.globals.maxSize = 8 * 1024^3)  # 8GB RAM 限制，可依需求調整


###############################################
## Step 2: 準備 Seurat 資料
###############################################

# 假設你已經有一個 Seurat 物件 seurat_obj，並完成分群
# 取出資料與 metadata
data.input <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
meta <- seurat_obj@meta.data

# 建立 CellChat 物件，這裡以 seurat_clusters 作為群組依據
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "seurat_clusters")


###############################################
## Step 3: 指定物種資料庫
###############################################

# 若為人類資料，使用 CellChatDB.human；小鼠請改用 CellChatDB.mouse
CellChatDB <- CellChatDB.human  # 或 CellChatDB.mouse
cellchat@DB <- CellChatDB


###############################################
## Step 4: 預處理與篩選 gene-interaction
###############################################

# 選取 database 中有用的基因
cellchat <- subsetData(cellchat)  

# 偵測高度表現的配體/受體
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 加入蛋白質交互網路（PPI network）
cellchat <- projectData(cellchat, PPI.human)  # human 為例，mouse 可忽略


###############################################
## Step 5: 推論細胞間通訊網絡
###############################################

# 推論通訊概率
cellchat <- computeCommunProb(cellchat)

# 過濾至少出現在多少細胞以上的互動（建議 10）
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 顯示互動資料框
df.net <- subsetCommunication(cellchat)

# pathway 層級通訊
cellchat <- computeCommunProbPathway(cellchat)

# 整合 network
cellchat <- aggregateNet(cellchat)


###############################################
## Step 6: 可視化 CellChat 網絡
###############################################

# 查看每個群組的細胞數量
groupSize <- as.numeric(table(cellchat@idents))

# 繪製總體通訊量的網絡圖（線粗代表互動多寡）
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE)

# 繪製某 pathway 的網絡（例：TGFb）
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "circle")


###############################################
## Step 7: pathway 特定通訊熱圖
###############################################

# 繪製所有通訊路徑的細胞通訊角色熱圖（發送與接收）
netAnalysis_signalingRole_heatmap(cellchat, pattern = "both")


###############################################
## Step 8: 推論細胞在通訊中的角色（中心性分析）
###############################################

# 計算中心性（例如傳送者、接收者）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 顯示特定 pathway 的角色（例：CXCL）
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL")


###############################################
## Step 9: 模式識別（Pattern discovery）
###############################################

# 發送訊息的 pattern 模式
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing")

# 可視化模式（Pattern 1 為例）
netAnalysis_river(cellchat, pattern = 1)

# 氣泡圖視覺化 pattern
netAnalysis_dot(cellchat, pattern = 1)


###############################################
## Step 10: 儲存結果
###############################################

# 儲存 CellChat 分析結果
saveRDS(cellchat, file = "cellchat_result.rds")

# 若要日後讀取：
# cellchat <- readRDS("cellchat_result.rds")

