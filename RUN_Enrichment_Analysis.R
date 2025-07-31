Name_ExportFolder_EA <- paste0(Name_ExportFolder,"/Enrichment Analysis") 
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder_EA)){dir.create(Name_ExportFolder_EA)}

###################################################################################
#### Enrichment Analysis ####
# # 基因名稱轉換 (SYMBOL -> ENTREZID)
# deg_upregulated_entrez <- bitr(deg_upregulated_genes_df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
# deg_downregulated_entrez <- bitr(deg_downregulated_genes_df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID


# 基因名稱轉換 (ENSEMBL -> ENTREZID)
deg_upregulated_entrez <- bitr(deg_upregulated_genes_df$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
deg_downregulated_entrez <- bitr(deg_downregulated_genes_df$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID



#### GO 富集分析 ####
try({
  go_up <- enrichGO(
    gene = deg_upregulated_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",   # 分析所有的 GO 分類 (BP, MF, CC)
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  ) 
  
  PlotDot_GO_UP <- dotplot(go_up, showCategory = 20) + ggtitle("GO Enrichment Analysis of Upregulated Genes")
  PlotDot_GO_UP
  
})


try({
  go_down <- enrichGO(
    gene = deg_downregulated_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  PlotDot_GO_DOWN <- dotplot(go_down, showCategory = 20) + ggtitle("GO Enrichment Analysis of Downregulated Genes")
  PlotDot_GO_DOWN
  
})


# 保存結果並繪製圖表
try({
  write.csv(go_up, paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_GO_Upregulated.csv"))
  
  jpeg(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_GO_Upregulated.jpg"), width = 500, height = 700)
  print(PlotDot_GO_UP)
  dev.off()
  
})


try({
  write.csv(go_down,paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_GO_Downregulated.csv"))
  
  jpeg(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_GO_Downregulated.jpg"), width = 500, height = 700)
  print(PlotDot_GO_DOWN)
  dev.off()
  
})




#### KEGG 富集分析 ####
try({
  kegg_up <- enrichKEGG(
    gene = deg_upregulated_entrez,
    organism = 'hsa',
    pvalueCutoff = 0.05
  )
  
  PlotDot_KEGG_UP <- dotplot(kegg_up, showCategory = 20) + ggtitle("KEGG Pathways of Upregulated Genes")
  PlotDot_KEGG_UP
  
})


try({
  kegg_down <- enrichKEGG(
    gene = deg_downregulated_entrez,
    organism = 'hsa',
    pvalueCutoff = 0.05
  )
  
  PlotDot_KEGG_DOWN <- dotplot(kegg_down, showCategory = 20) + ggtitle("KEGG Pathways of Downregulated Genes")
  PlotDot_KEGG_DOWN
  
  
})


# 保存結果並繪製圖表
try({
  write.csv(kegg_up, paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_KEGG_Upregulated.csv"))
  
  jpeg(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_KEGG_Upregulated.jpg"), width = 500, height = 700)
  print(PlotDot_KEGG_UP)
  dev.off()
  
})


try({
  write.csv(kegg_down, paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_KEGG_Downregulated.csv"))
  
  jpeg(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_KEGG_Downregulated.jpg"), width = 500, height = 700)
  print(PlotDot_KEGG_DOWN)
  dev.off()
  
})


#### Reactome 富集分析 ####
try({
  reactome_up <- enrichPathway(
    gene = deg_upregulated_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  PlotDot_Reactome_UP <- dotplot(reactome_up, showCategory = 20) + ggtitle("Reactome Pathways of Upregulated Genes")
  PlotDot_Reactome_UP
  
})



try({
  reactome_down <- enrichPathway(
    gene = deg_downregulated_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  PlotDot_Reactome_DOWN <- dotplot(reactome_down, showCategory = 20) + ggtitle("Reactome Pathways of Downregulated Genes")
  PlotDot_Reactome_DOWN
  
})




# 保存結果並繪製圖表
try({
  write.csv(reactome_up, paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_Reactome_Upregulated.csv"))
  
  jpeg(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_Reactome_Upregulated.jpg"), width = 500, height = 700)
  print(PlotDot_Reactome_UP)
  dev.off()
  
})


try({
  write.csv(reactome_down, paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_Reactome_Downregulated.csv"))
  
  jpeg(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_Reactome_Downregulated.jpg"), width = 500, height = 700)
  print(PlotDot_Reactome_DOWN)
  dev.off()
  
})




pdf(paste0(Name_ExportFolder_EA,"/", Name_Export,"_EA_DotPlot.pdf"), width = 8, height = 10)

try(print(PlotDot_GO_UP))
try(print(PlotDot_GO_DOWN))
try(print(PlotDot_KEGG_UP))
try(print(PlotDot_KEGG_DOWN))
try(print(PlotDot_Reactome_UP))
try(print(PlotDot_Reactome_DOWN))

dev.off()

