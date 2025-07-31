#############################################################
##  One-click ORA for any cell type – Teaching Example     ##
#############################################################

## ---------------- (0) 套件載入 ---------------- ##
if (!require("tidyverse"))        { install.packages("tidyverse");        library(tidyverse) }
if (!require("clusterProfiler"))  { install.packages("clusterProfiler");  library(clusterProfiler) }
if (!require("org.Hs.eg.db"))     { BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db) }
if (!require("ReactomePA"))       { BiocManager::install("ReactomePA");   library(ReactomePA) }

## ---------------- (1) 基本設定 ---------------- ##
cell_type      <- "Basal keratinocytes"   # <-- 改成目標細胞類型
output_dir     <- "."                     # 儲存目錄；改成路徑亦可
output_prefix  <- gsub("\\s+", "_", cell_type)  # 檔名前綴，空白改底線

## ---------------- (2) 擷取該類型基因 ---------------- ##
# 假設 top25_lines 形如：names(top25_lines) = cell_type 字串，值為 "Type: G1, G2, ..."
line_raw <- top25_lines[grepl(paste0("^", cell_type), names(top25_lines))]

if (length(line_raw) == 0) {
  stop("❗ 在 top25_lines 找不到指定 cell_type，請確認拼寫。")
}

gene_vec <- line_raw %>%              # 取出基因字串
  sub("^[^:]+:\\s*", "", .) %>%       # 去掉「Basal keratinocytes:」
  str_split(",\\s*") %>%              # 逗號分割
  unlist() %>% unique()               # 轉向量並去掉重複

cat("✔ 讀到基因數量：", length(gene_vec), "\n")

## ---------------- (3) SYMBOL → ENTREZ 轉換 ---------------- ##
entrez_vec <- bitr(gene_vec,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db) %>%
  pull(ENTREZID) %>% unique()

## ---------------- (4) ORA：GO / KEGG / Reactome ---------- ##
# ---- 4-1 GO ----
go_res <- enrichGO(gene         = entrez_vec,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "ALL",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.05,
                   readable     = TRUE)

# ---- 4-2 KEGG ----
kegg_res <- enrichKEGG(gene         = entrez_vec,
                       organism     = "hsa",
                       pvalueCutoff = 0.05)

# ---- 4-3 Reactome ----
react_res <- enrichPathway(gene         = entrez_vec,
                           organism     = "human",
                           pvalueCutoff = 0.05,
                           readable     = TRUE)

## ---------------- (5) 產生點圖 ---------------- ##
plot_go    <- dotplot(go_res,    showCategory = 20) + ggtitle(paste(cell_type, "– GO"))
plot_kegg  <- dotplot(kegg_res,  showCategory = 20) + ggtitle(paste(cell_type, "– KEGG"))
plot_react <- dotplot(react_res, showCategory = 20) + ggtitle(paste(cell_type, "– Reactome"))

## ---------------- (6) 輸出結果 ---------------- ##
# 6-1 CSV
write.csv(go_res,    file = file.path(output_dir, paste0(output_prefix, "_GO.csv")))
write.csv(kegg_res,  file = file.path(output_dir, paste0(output_prefix, "_KEGG.csv")))
write.csv(react_res, file = file.path(output_dir, paste0(output_prefix, "_Reactome.csv")))

# 6-2 JPG 點圖
for (plt in list(GO = plot_go, KEGG = plot_kegg, Reactome = plot_react)) {
  jpeg(file.path(output_dir, paste0(output_prefix, "_", names(plt), ".jpg")),
       width = 600, height = 800)
  print(plt[[1]])
  dev.off()
}

# 6-3 PDF 整合
pdf(file.path(output_dir, paste0(output_prefix, "_DotPlots.pdf")),
    width = 7, height = 9)
print(plot_go); print(plot_kegg); print(plot_react)
dev.off()

cat("🎉 ORA finished!  相關檔案已存於：", normalizePath(output_dir), "\n")
