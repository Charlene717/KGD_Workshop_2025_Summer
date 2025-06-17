#Library
library(CellChat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(Seurat)
options(stringsAsFactors = FALSE)

##Split object
Split_cellchat.list<- SplitObject(TNname.combined, split.by = "orig.ident1")
Split_cellchat.list
Split1 <-Split_cellchat.list$'1_PDAC'
Split2 <-Split_cellchat.list$'2_ADJ'

#Check
DefaultAssay(Split1) <- "RNA"
DefaultAssay(Split2) <- "RNA"

#Split1_1PDAC#umap
DimPlot(Split1, reduction = "umap", label = F, pt.size = 0.8)
DimPlot(Split1, reduction = "umap", label = T, pt.size = 0.8)
#Split2_1ADJ#umap
DimPlot(Split2, reduction = "umap", label = F, pt.size = 0.8)
DimPlot(Split2, reduction = "umap", label = T, pt.size = 0.8)


##Part I: Data input & processing and initialization of CellChat object
#Extract the CellChat input files from a Seurat V3 object
data.input <- GetAssayData(Split2, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Split2)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#Create a CellChat object using data matrix as input
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#Add cell information into meta slot of the object
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
##Download all interaction pathways of CellChatDB
CellChatDB$interaction
write.csv(CellChatDB$interaction, file="C:/Users/Administrator/Downloads/Nath/CellChatDB_interaction_S2.csv", row.names=FALSE)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF???
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T???
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB???
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF???
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T???
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago???
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""???
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep???
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435???
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre???
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2???

##Download all interaction pathways of CellChatDB
#CellChatDB$interaction
#write.csv(CellChatDB$interaction, file="C:/Users/Administrator/Downloads/Nath/CellChatDB_interaction.csv", row.names=FALSE)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#Isa-isa sa pag-run#
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 3)
#note: 4_Regression
#The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  C29_Mac


#Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction 
#strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#Isa-isa
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")

#Due to the complicated cell-cell communication network, 
#we can examine the signaling sent from each cell group. 
#Here we also control the parameter edge.weight.max so that
#we can compare edge weights between different networks.
mat <- cellchat@net$count
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


mat <- cellchat@net$count
mat <- mat[rowSums(mat)!=0,]
mat <- mat[,colSums(mat)!=0]
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

if (nrow(mat) != ncol(mat)) {
  stop("The input matrix 'mat' must be square.")
}

for (i in 1:nrow(mat)) {
  # Create a zero matrix with same dimensions and dimnames as mat
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  # Copy the i-th row of mat to mat2
  mat2[i, ] <- mat[i, ]
  
  # Visualize the network
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

mat <- cellchat@net$weight
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = max(mat), title.name = rownames(mat)[i] )
}

mat <- cellchat@net$weight
mat <- mat[rowSums(mat)!=0,]
mat <- mat[,colSums(mat)!=0]
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= T, edge.weight.max = max(mat), title.name = rownames(mat)[i] )
}



##Part III: Visualization of cell-cell communication network
##Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
#Here we take input of one signaling pathway as an example. 
#All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways.
write.csv(cellchat@netP$pathways, file="C:/Users/Administrator/Downloads/Nath/Split2_cellchat_netP_pathways.csv", row.names=FALSE)


## move to merge script 
pathways.show <- c("CALCR")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

##Compute the contribution of each ligand-receptor pair to the overall signaling pathway 
#and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

#We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.SN <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.SN[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

##Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#Bubble plot
# show all the significant interactions (L-R pairs) from some cell 
#groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object

#Chord diagram
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(0:9), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(0:9), signaling = c("CCL","CXCL"),legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:7), slot.name = "netP", legend.pos.x = 10)

##Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "SN")
plotGeneExpression(cellchat, signaling = "SN", enriched.only = FALSE)

##Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
viz <- netAnalysis_signalingRole_network(cellchat, width = 10, height = 5, font.size = 9)
#cellchat, signaling = pathways.show, width

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
#gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CDH"))
#> Signaling role analysis on the cell-cell communication network from user's input
#gg1 + gg2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN", "BSP"))
ht

##Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
library(NMF)

##Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat, pattern = "outgoing")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

##Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
#Cophenetic values begin to drop when the number of incoming patterns is 4.
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#skip
library("reticulate")
use_virtualenv("C:/Users/Administrator/Documents/.virtualenvs/r-reticulate")

##Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional") #errorhere #end here
#
tmp <- cellchat@netP$similarity$functional$dr$single
tmp2 <- cellchat@netP$prob  
tmpin_names <- rownames(tmp[!is.nan(tmp[,1]),])
tmp_m <- tmp[tmpin_names,]
tmp2_m <- tmp2[,,tmpin_names]
cellchat@netP$similarity$functional$dr$single <- tmp_m
cellchat@netP$prob <- tmp2_m
#
cellchat <- netClustering(cellchat, type = "functional")

#Re-run
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#
tmp <- cellchat@netP$similarity$functional$dr$single
tmp2 <- cellchat@netP$prob  
tmpin_names <- rownames(tmp[!is.nan(tmp[,1]),])
tmp_m <- tmp[tmpin_names,]
tmp2_m <- tmp2[,,tmpin_names]
cellchat@netP$similarity$functional$dr$single <- tmp_m
cellchat@netP$prob <- tmp2_m
#
cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

#resume
##Part V: Save the CellChat object
saveRDS(cellchat, file = "C:/Users/Administrator/Downloads/Nath/CellChat_2ADJ.rds")