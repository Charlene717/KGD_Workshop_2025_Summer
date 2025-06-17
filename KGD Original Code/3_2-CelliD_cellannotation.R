##CelliD for cell annotation
# #install
# install.packages("devtools")
#setRepositories(ind = c(1,2,3))
#if(!require("tidyverse")) install.packages("tidyverse")
#if(!require("ggpubr")) install.packages("ggpubr")
#BiocManager::install("CelliD")

#load library 
library(CelliD)
library(tidyverse) # general purpose library for data handling
library(ggpubr) #library for plotting

#load TN.combined
TN.combined <- readRDS(file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/TN.combined_dim30-DF-35.rds")
TN.combined_DF_copy <- TN.combined_DF
#DefaultAssay(TN.combined) <- "RNA"

#CelliD dimensionality reduction through MCA
Baron <- RunMCA(TN.combined_DF_copy)
#download all cell-type gene signatures from panglaoDB
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
#filter to get human specific genes
panglao_all <- panglao %>%  filter(str_detect(species,"Hs"))
# convert dataframes to a list of named vectors which is the format for CelliD input
panglao_all <- panglao_all %>%  
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))
all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)
#remove very short signatures
all_gs <- all_gs[sapply(all_gs, length) >= 10]
#RunCellHGT
HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]
Baron$all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")
#Output data
DimPlot(Baron, group.by = "all_gs_prediction_signif", reduction = "umap",label = TRUE, label.size = 2,repel = TRUE)+
  theme(legend.text = element_text(size = 2), aspect.ratio = 1)
clustering.table_CelliD <- table(Baron@meta.data[["all_gs_prediction_signif"]], Baron@active.ident)
clustering.table_CelliD
write.csv(clustering.table_CelliD , file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/082024_CelliD_PanglaoDB.csv")

#load CelliD analysis 
clustering.table_CelliD <- read.csv("D:/Jojie/Analysis_Outputs/scRNA-11072024/082024_CelliD_PanglaoDB.csv")
rownames(clustering.table_CelliD) <- clustering.table_CelliD[,1]
clustering.table_CelliD <- clustering.table_CelliD[-nrow(clustering.table_CelliD),-1]
clustering.table_CelliD["annotation",] <- rownames(clustering.table_CelliD)[apply(clustering.table_CelliD,2,which.max)]
write.csv(clustering.table_CelliD , file = "D:/Jojie/Analysis_Outputs/scRNA-11072024/082024_CelliD_PanglaoDB_summary_default.csv", col.names = TRUE) 
