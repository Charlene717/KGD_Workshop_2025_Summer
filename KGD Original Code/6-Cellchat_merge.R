###Comparison analysis of multiple datasets using CellChat

#Load the required libraries
library(CellChat)
library(patchwork)

#Load CellChat object of each dataset and then merge togethe
cellchat.TN1 <- readRDS(file = "C:/Users/user/Desktop/202404/Cellchat_1_healthy_skin.rds")
cellchat.TN2 <- readRDS(file = "C:/Users/user/Desktop/202404/Cellchat_2_DSAP_skin.rds")
cellchat.TN3 <- readRDS(file = "C:/Users/user/Desktop/202404/Cellchat_3_DSAP_skin_abrocitinib.rds ")

cellchat.TN1 <- updateCellChat(cellchat.TN1)
cellchat.TN2 <- updateCellChat(cellchat.TN2)
cellchat.TN3 <- updateCellChat(cellchat.TN3)


#adjust cell cluster numbers, if one sample has empty cell clusters
cellchat.TN1 <- liftCellChat(cellchat.TN1,group.new=levels(cellchat.TN2@idents)) 
cellchat.TN1 <- liftCellChat(cellchat.TN1,group.new=levels(cellchat.TN3@idents))
# Compute the network centrality scores
cellchat.TN1 <- netAnalysis_computeCentrality(cellchat.TN1, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat.TN2 <- netAnalysis_computeCentrality(cellchat.TN2, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat.TN3 <- netAnalysis_computeCentrality(cellchat.TN3, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#cellchat.TN4 <- netAnalysis_computeCentrality(cellchat.TN4, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#Because the numbers of cell population in TN71 (17 populations) is different from TN72
#(18 populations), we can use liftcellchat to skip this problem.
object.list <- list(Healthy_skin = cellchat.TN1, DSAP_skin = cellchat.TN2, DSAP_skin_abrocitinib = cellchat.TN3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##Part I: Predict general principles of cell-cell communication
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), color.use = c("#00B9E3","#F8766D","darkGreen"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight",color.use = c("#00B9E3","#F8766D","darkGreen"))
gg1 + gg2
ggsave("C:/Users/user/Desktop/202404/interactions_strength.pdf")

#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, vertex.label.cex = 0.9, weight.scale = T)
netVisual_diffInteraction(cellchat, vertex.label.cex = 0.9, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

pdf("C:/Users/user/Desktop/202404/heatmap_comparision.pdf", width = 8, height = 8)
gg1 + gg2
dev.off()


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, vertex.label.cex = 0.9, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, vertex.label.cex = 0.9, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]))
}


#INSERTED#####
pathways.show <- c("SEMA4") 

#Contribution
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,5) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
c1 <-  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 8)
p <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = pathways.show)

##END



#Differential number of interactions or interaction strength among different cell types
group.cellType <- c(rep("C1_Epithelial", 1), rep("C4_Mac", 1),rep("C9_Mac", 1),rep("C12_Mac", 1), rep("C14_Epithelial", 1), rep("C15_Epithelial", 1)) #1st compare
group.cellType <- factor(group.cellType, levels = c("C1_Epithelial","C4_Mac","C9_Mac","C12_Mac", "C14_Epithelial", "C15_Epithelial"))
group.cellType <- NULL
group.cellType <- c(rep("C1_Epithelial", 4),rep("C12_Mac", 4), rep("C14_Epithelial", 4), rep("C15_Epithelial", 4)) #1st compare
group.cellType <- factor(group.cellType, levels = c("C1_Epithelial","C12_Mac", "C14_Epithelial", "C15_Epithelial"))

group.cellType <- c(rep("C13_T", 4), rep("C14_T", 4), rep("C20_Ductal", 4), rep("C22_T", 4)) #2nd compare
group.cellType <- factor(group.cellType, levels = c("C13_T","C14_T","C20_Ductal","C22_T"))

group.cellType <- c(rep("C3_Neu", 4), rep("C8_Mac", 4), rep("C11_Mac", 4), rep("C12_Mac", 4), rep("C17_Mac", 4), rep("C20_Ductal", 4)) #3rd compare
group.cellType <- factor(group.cellType, levels = c("C3_Neu","C8_Mac","C11_Mac","C12_Mac","C17_Mac","C20_Ductal"))


object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, vertex.label.cex = 0.9, weight.scale = T, label.edge= F, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","weight", "weight.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, vertex.label.cex = 0.9, weight.scale = T, label.edge= F, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
num.link <- unlist(num.link)
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()


for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax, label.size = 4)+ scale_y_continuous(limits = c(0,15)) + scale_x_continuous(limits = c(0,20))
}

pdf("C:/Users/user/Desktop/202404/2D_interaction_comparision.pdf", width = 10, height = 10)
patchwork::wrap_plots(plots = gg)
dev.off()


##Part II: Identify the conserved and context-specific signaling pathways
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional", comparison = c(1,2,3))
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional", comparison = c(1,2,3))
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional", comparison = c(1,2,3))
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5, comparison = c(1,2,3))
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2, comparison = c(1,2,3))

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural", comparison = c(1,2,3))
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural", comparison = c(1,2,3))
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural", comparison = c(1,2,3))
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5, comparison = c(1,2,3))
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2, comparison = c(1,2,3))

#Compute and visualize the pathway distance in the learned joint manifold
#Larger distance implies larger difference of the communication networks between 
#two datasets in terms of either functional or structure similarity. 
rankSimilarity(cellchat, type = "functional", comparison1 = c(1,2,3))
rankSimilarity(cellchat, type = "structural", comparison1 = c(1,2,3))
#> Compute the distance of signaling networks between datasets 1 2

#Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#00BFC4","#F8766D"), comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#00BFC4","#F8766D"),comparison = c(1,2))
gg1 + gg2
ggsave("C:/Users/user/Desktop/202404/information_flow_1_2.pdf", width = 30, height = 30, units = "cm")

# 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#00BFC4","#F8766D"), comparison = c(2,3))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#00BFC4","#F8766D"),comparison = c(2,3))
gg1 + gg2
ggsave("C:/Users/user/Desktop/202404/information_flow_2_3.pdf", width = 30, height = 30, units = "cm")
# 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#F8766D","#00BFC4"), comparison = c(1,3))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#F8766D","#00BFC4"),comparison = c(1,3))
gg1 + gg2
ggsave("C:/Users/user/Desktop/202404/information_flow_1_3.pdf", width = 30, height = 30, units = "cm")

#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)

i=1
# combining all the identified signaling pathways from different datasets 
# outgoing pathway patterns
pathway.union <- union(object.list[[i+1]]@netP$pathways, object.list[[i]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12, font.size = 5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 12, font.size = 5)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 8, height = 12, font.size = 5)
# ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 8, height = 12, font.size = 5)
draw(ht1 + ht2  , ht_gap = unit(0.5, "cm"))
draw(ht1 + ht3  , ht_gap = unit(0.5, "cm"))
draw(ht2 + ht3  , ht_gap = unit(0.5, "cm"))
# draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

pdf("output/outgoing_patterns.pdf", width = 9, height = 6)
draw(ht1 + ht2  , ht_gap = unit(0.5, "cm"))
dev.off()


# incoming pathway patterns
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12, color.heatmap = "GnBu", font.size = 5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 12, color.heatmap = "GnBu", font.size = 5)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 8, height = 12, color.heatmap = "GnBu", font.size = 5)
# ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 8, height = 12, color.heatmap = "GnBu", font.size = 5)
draw(ht1 + ht2 , ht_gap = unit(0.5, "cm"))
draw(ht1 + ht3  , ht_gap = unit(0.5, "cm"))
draw(ht2 + ht3  , ht_gap = unit(0.5, "cm"))
# draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

pdf("output/incoming_patterns.pdf", width = 9, height = 6)
draw(ht1 + ht2  , ht_gap = unit(0.5, "cm"))
dev.off()

# all pathway patterns
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12, color.heatmap = "OrRd", font.size = 5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 12, color.heatmap = "OrRd", font.size = 5)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 8, height = 12, color.heatmap = "OrRd", font.size = 5)
# ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+3], width = 8, height = 12, color.heatmap = "OrRd", font.size = 5)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
draw(ht1 + ht3  , ht_gap = unit(0.5, "cm"))
draw(ht2 + ht3  , ht_gap = unit(0.5, "cm"))
# draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

pdf("output/all_patterns.pdf", width = 9, height = 6)
draw(ht1 + ht2  , ht_gap = unit(0.5, "cm"))
dev.off()

##Part III: Identify the up-regulated and down-regulated signaling ligand-receptor pairs
##We can compare the communication probabilities mediated by ligand-receptor pairs from some cell groups to other cell groups
netVisual_bubble(cellchat, sources.use = c(17), targets.use = c(10),  comparison = c(1,2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(1,3,6,9,11,12,14), targets.use = c(2,15,16),  comparison = c(1,2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(2,15,16), targets.use = c(5,10,13,17),  comparison = c(1,2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(5,10,13,17), targets.use = c(2,15,16),  comparison = c(1,2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(3,21,22), targets.use = c(10,11),  comparison = c(1,2,3), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(10,11), targets.use = c(3,21,22),  comparison = c(1,2,3), angle.x = 45)

interactions <- data.frame(interaction_name = "ICAM1_SPN")
netVisual_bubble(cellchat, sources.use = c(5,10,13), targets.use = c(16),  comparison = c(1,2), angle.x = 45, pairLR.use = interactions)
netVisual_bubble_ICAM1_SPN


#> Comparing communications on a merged object
gg1 <- netVisual_bubble(cellchat, sources.use = c(1,7,8,9), targets.use = c(5),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,7,8,9), targets.use = c(5),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, sources.use = c(1,7,8,9), targets.use = c(5),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,7,8,9), targets.use = c(5),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, sources.use = c(7), targets.use = c(1,9),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(7), targets.use = c(1,9),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg1 + gg2

#The ligand-receptor pairs shown in the bubble plot can be accessed via
signaling.LSIncreased = gg1$data

##We can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis. 
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Active"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#>Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Active",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Adjacentskin",ligand.logFC = -0.1, receptor.logFC = -0.1)
#Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

#We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1,9), targets.use = c(7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1,9), targets.use = c(7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1,9), targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = c(1,9), targets.use = c(7), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


##Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
#Contribution
#in this dataset (scSeq_2862vs2901), we investigated 5 pathways NECTIN, DESMOSOME, NT, EPGN, ncWNT
#NECTIN, DESMOSOME were presented in both samples 
pathways.show <- c("CALCR")
i1 <-  netAnalysis_contribution(object.list[[1]], signaling = pathways.show, title = paste(pathways.show, names(object.list)[1]))
i2 <-  netAnalysis_contribution(object.list[[2]], signaling = pathways.show, title = paste(pathways.show, names(object.list)[2]))
i3 <-  netAnalysis_contribution(object.list[[3]], signaling = pathways.show, title = paste(pathways.show, names(object.list)[3]))
# i4 <-  netAnalysis_contribution(object.list[[4]], signaling = pathways.show, title = paste(pathways.show, names(object.list)[4]))
#i1+i2+i3+i4
i1+i2+i3

pdf(paste("C:/Users/user/Desktop/202404/",pathways.show, "_PathwayContribution.pdf"), width = 10, height = 6)
i1+i2+i3
dev.off()


# Hierarchy plot
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10)
# Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pdf(paste("output/cellchat/",pathways.show, "_PathwayHierarchyPlot.pdf"), width = 9, height = 6)
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


#Circle plot 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pdf(paste("output/cellchat/",pathways.show, "_PathwayCirclePlot.pdf"), width = 9, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

#Heatmap
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
# ComplexHeatmap::draw(ht[[1]] + ht[[2]]+ht[[3]] + ht[[4]], ht_gap = unit(0.5, "cm"))
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


pdf(paste("output/cellchat/",pathways.show, "_PathwayHeatmap.pdf"), width = 9, height = 6)
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# Compute the network centrality scores
cellchat.TN1 <- netAnalysis_computeCentrality(cellchat.TN1, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat.TN2 <- netAnalysis_computeCentrality(cellchat.TN2, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat.TN3 <- netAnalysis_computeCentrality(cellchat.TN3, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# cellchat.TN4 <- netAnalysis_computeCentrality(cellchat.TN4, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.TN1, signaling = pathways.show, width = 8, height = 2.5, font.size = 8)
netAnalysis_signalingRole_network(cellchat.TN2, signaling = pathways.show, width = 8, height = 2.5, font.size = 8)
netAnalysis_signalingRole_network(cellchat.TN3, signaling = pathways.show, width = 8, height = 2.5, font.size = 8)
# netAnalysis_signalingRole_network(cellchat.TN4, signaling = pathways.show, width = 8, height = 2.5, font.size = 8)



## below are codes for pathways that only existed in one sample
## in this dataset (scSeq_2862vsS1), we investigated 2 pathways EPGN, NEGR
##Part III : Visualization of cell-cell communication network for single pathway existed in certain sample
##Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
#Here we take input of one signaling pathway as an example. 
#All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways.
pathways.show <- c("BSP","CSPG4") 

#Contribution
i1 <- netAnalysis_contribution(cellchat.TN3, signaling = pathways.show)

pdf(paste("C:/Users/user/Desktop/202404/",pathways.show, "_PathwayContribution.pdf"), width = 9, height = 6)
i1
dev.off()


# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,10) # a numeric vector. 
p <- netVisual_aggregate(cellchat.TN2, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")

pdf(paste("C:/Users/user/Desktop/202404/",pathways.show, "_PathwayHierarchyPlot.pdf"), width = 9, height = 6)
p
dev.off()


# Circle plot
par(mfrow=c(1,1))
p <- netVisual_aggregate(cellchat.TN2, signaling = pathways.show, layout = "circle")

pdf(paste("output/cellchat/",pathways.show, "_PathwayCirclePlot.pdf"), width = 9, height = 6)
p
dev.off()


# Heatmap
par(mfrow=c(1,1))
p <- netVisual_heatmap(cellchat.TN2, signaling = pathways.show, color.heatmap = "Reds")

pdf(paste("output/cellchat/",pathways.show, "_PathwayHeatmap.pdf"), width = 9, height = 6)
p
dev.off()

# Compute the network centrality scores
cellchat.TN2 <- netAnalysis_computeCentrality(cellchat.TN2, slot.name = c("netP")) # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
c1 <-  netAnalysis_signalingRole_network(cellchat.TN2, signaling = pathways.show, width = 8, height = 2.5, font.size = 6)

#Plot the signaling gene expression distribution using violin/dot plot
p <- plotGeneExpression(cellchat.TN2, signaling = pathways.show)

pdf(paste("C:/Users/user/Desktop/202404/",pathways.show, "_Pathwaygeneexpression.pdf"), width = 9, height = 6)
p
dev.off()

#We can also visualize the cell-cell communication mediated by a single ligand-receptor pair.
pairLR.CXCL <- extractEnrichedLR(cellchat.TN2, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[45,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,14) # a numeric vector
netVisual_individual(cellchat.TN2, signaling = pathways.show,  pairLR.use = LR.show, layout = "hierarchy", vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat.TN2, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


##Part V: Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Healthy_skin", "DSAP_skin","DSAP_skin_abrocitinib")) # set factor level
plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, color.use = c("#00BA38","#F8766D","#00B9E3"))
#plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, color.use = c("#00BA38","#00B9E3","#F8766D","#DB72FB"))

#Extract the raw data of cellchat
CommunicationRawdata <- subsetCommunication(cellchat)
#CommunicationRawdata_netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(CommunicationRawdata$Healthy_skin, file = "C:/Users/user/Desktop/202404/CommunicationRawdata_Cellchat_Healthy_skin.csv", sep = "\t", row.names = F)
write.csv(CommunicationRawdata$DSAP_skin, file = "C:/Users/user/Desktop/202404/CommunicationRawdata_Cellchat_Adjacent_skin.csv", sep = "\t", row.names = F)
write.csv(CommunicationRawdata$DSAP_skin_abrocitinib, file = "C:/Users/user/Desktop/202404/CommunicationRawdata_Cellchat_Active.csv", sep = "\t", row.names = F)
#write.csv(CommunicationRawdata$Inactive, file = "C:/Users/user/Desktop/CommunicationRawdata_netP_Cellchat_Inactive.csv", sep = "\t", row.names = F)

##Save the merged CellChat object
saveRDS(cellchat, file = "C:/Users/user/Desktop/202404/cellchat_Merge.rds")
cellchat <- readRDS(file = "C:/Users/user/Desktop/202404/cellchat_Merge.rds")
