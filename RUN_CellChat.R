# seurat_all_integrated <- readRDS("C:/Charlene/Code_GitHub_BioInport2025/KGD_Workshop_2025_Summer/Export_2025062617YDG_Demo/2025062617YDG_Integration.rds")             # 載入儲存好的 Seurat RDS 物件

###############################################
## Step 1: 安裝與載入套件
###############################################
if (!require("Seurat"))       { install.packages("Seurat");         library(Seurat) }
if (!require("tidyverse"))    { install.packages("tidyverse");      library(tidyverse) }
if (!require("ggplot2"))      { install.packages("ggplot2");        library(ggplot2) }
if (!require("patchwork"))    { install.packages("patchwork");      library(patchwork) }
if (!require("future"))       { install.packages("future");         library(future) }

if (!require("remotes"))      { install.packages("remotes");        library(remotes) }
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat"); library(CellChat)
}


# 可選：加速處理（根據 CPU 核心數調整 workers）
plan("multisession", workers = 4)  
options(future.globals.maxSize = 8 * 1024^3)  # 8GB RAM 限制，可依需求調整


###############################################
## Step 2: 準備 Seurat 資料
###############################################

# # ✅ 將 Seurat v5 多層 assay 合併為單一層（預設會合併為 "data"）
# seurat_all_integrated <- JoinLayers(object = seurat_all_integrated, layers = "data")

# 將 RNA assay 的多層合併為一層（預設合併到 "data" slot）
seurat_all_integrated[["RNA"]] <- JoinLayers(object = seurat_all_integrated[["RNA"]])



# 假設你已經有一個 Seurat 物件 seurat_all_integrated，並完成分群
# 取出資料與 metadata
data.input <- GetAssayData(seurat_all_integrated, slot = "data", assay = "RNA")
meta <- seurat_all_integrated@meta.data

# ## 建立 CellChat 物件，這裡以 seurat_clusters 作為群組依據
# #Bug# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "seurat_clusters")
# 
# # 把原本的 seurat_clusters（例如 "0", "1", "2"）改成 "C0", "C1", ...
# meta$cellchat_clusters <- paste0("C", as.character(seurat_all_integrated$seurat_clusters))
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_clusters")

## 建立 CellChat 物件，這裡以 Cell_Type 作為群組依據
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_Type")

###############################################
## Step 3: 指定物種資料庫
###############################################

# 若為人類資料，使用 CellChatDB.human；小鼠請改用 CellChatDB.mouse
CellChatDB <- CellChatDB.human  # 或 CellChatDB.mouse
cellchat@DB <- CellChatDB

showDatabaseCategory(CellChatDB)

###############################################
## Step 4: 預處理與篩選 gene-interaction
###############################################

# 選取 database 中有用的基因
cellchat <- subsetData(cellchat)  

# 偵測高度表現的配體/受體
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# # 加入蛋白質交互網路（PPI network）
# cellchat <- projectData(cellchat, PPI.human)  # human 為例，mouse 可忽略


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


# 繪製某 pathway 的網絡（例：GAP）
cellchat@netP$pathways # %>% head()
netVisual_aggregate(cellchat, signaling = "GAP", layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = "GAP", layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "GAP", color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = "GAP")


netVisual_bubble(cellchat, sources.use = "Granular keratinocytes", targets.use = levels(cellchat@idents), remove.isolate = FALSE)



###############################################
## Step 7: 推論細胞在通訊中的角色（中心性分析）
###############################################

# 計算中心性（例如傳送者、接收者）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 顯示特定 pathway 的角色（例：GAP）
netAnalysis_signalingRole_network(cellchat, signaling = "GAP")


###############################################
## Step 8: pathway 特定通訊熱圖
###############################################

# 繪製所有通訊路徑的細胞通訊角色熱圖（發送與接收）
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all")

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") + 
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
 

###############################################
## Step 9: 儲存結果
###############################################

# 儲存 CellChat 分析結果
saveRDS(cellchat, file = "cellchat_result.rds")

# 若要日後讀取：
# cellchat <- readRDS("cellchat_result.rds")

