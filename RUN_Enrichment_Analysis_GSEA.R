#############################################################
##  Quick GSEA Pipeline – Teaching Example                 ##
#############################################################

## ---------------- (0) 套件載入 ---------------- ##
pkgs <- c("clusterProfiler", "enrichplot", "msigdbr",
          "org.Hs.eg.db", "tidyverse", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("org.Hs.eg.db", "msigdbr", "clusterProfiler", "enrichplot"))
      BiocManager::install(p)
    else
      install.packages(p)
  }
  library(p, character.only = TRUE)
}

## ---------------- (1) 基本設定 ---------------- ##
cell_type  <- "Basal keratinocytes"   # ← 換成目標群別
score_col  <- "avg_log2FC"            # ← 用哪個欄位做排序
output_dir <- "."                     # ← 輸出資料夾（"." = 工作目錄）
output_pre <- gsub("\\s+", "_", cell_type)

## ---------------- (2) 準備 ranking vector ---------------- ##
# 2-1 取該群的 DEG data.frame
deg_df <- deg_list[[cell_type]]

# 2-2 (可選) 若只想保留顯著差異基因，可在此加 p-value 篩選
# deg_df <- deg_df %>% filter(p_val_adj < 0.05)

# 2-3 基因 Symbol → ENTREZ ID
sym2ent <- bitr(deg_df$gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db,
                drop     = FALSE)

deg_df <- left_join(deg_df, sym2ent,
                    by = c("gene" = "SYMBOL")) %>% drop_na(ENTREZID)

# 2-4 建 gene_rank：named numeric 向量
gene_rank <- deframe(deg_df %>%                   # 取兩欄組成向量
                      dplyr::select(ENTREZID, !!score_col)) 
gene_rank <- sort(gene_rank, decreasing = TRUE)   # 由大到小排序

## ---------------- (3) 取得 MSigDB 基因集 ---------------- ##
species <- "Homo sapiens"

gmt_df <- bind_rows(
  msigdbr(species, category = "H"),   # Hallmark
  # msigdbr(species, category = "C2"),  # Canonical pathways
  # msigdbr(species, category = "C3"),  # TF targets
  # msigdbr(species, category = "C7")   # Immunologic
) %>% dplyr::select(gs_name, ENTREZID = entrez_gene) %>% distinct()

## ---------------- (4) 執行 GSEA ---------------- ##
gsea_res <- GSEA(geneList     = gene_rank,
                 TERM2GENE    = gmt_df,
                 pvalueCutoff = 1,          # 不過濾，後續再篩
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 verbose      = FALSE)

## ---------------- (5) 基本可視化 ---------------- ##
# 5-1 NES dotplot（前 30 條）
plot_dot <- dotplot(gsea_res, showCategory = 30) +
  ggtitle(paste("GSEA –", cell_type))

# 5-2 Ridgeplot 總覽
plot_ridge <- ridgeplot(gsea_res, showCategory = 30) +
  ggtitle(paste("GSEA ridgeplot –", cell_type))

## ---------------- (6) 輸出結果 ---------------- ##
# 6-1 CSV
write.csv(gsea_res@result,
          file = file.path(output_dir,
                           paste0(output_pre, "_GSEA_MSigDB.csv")),
          row.names = FALSE)

# 6-2 JPG 圖檔
ggsave(file.path(output_dir, paste0(output_pre, "_GSEA_dotplot.jpg")),
       plot_dot, width = 8, height = 14, dpi = 300)
ggsave(file.path(output_dir, paste0(output_pre, "_GSEA_ridgeplot.jpg")),
       plot_ridge, width = 8, height = 14, dpi = 300)

# 6-3 PDF 整合
pdf(file.path(output_dir, paste0(output_pre, "_GSEA_plots.pdf")),
    width = 8, height = 14)
print(plot_dot)
print(plot_ridge)
dev.off()

cat("🎉 GSEA finished!  結果已存於：", normalizePath(output_dir), "\n")
