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
# ---------------------- 4-1  Gene Ontology (GO)  ----------------------
go_res <- enrichGO(
  gene          = entrez_vec,      # 向量：欲檢定的基因（Entrez IDs）
  OrgDb         = org.Hs.eg.db,    # 物種註解資料庫，這裡使用人類 (org.Hs.eg.db)
  keyType       = "ENTREZID",      # 指定輸入基因 ID 類型
  ont           = "ALL",           # GO 三大面向皆檢定：BP / MF / CC
  pAdjustMethod = "BH",            # 多重假設校正方法—Benjamini-Hochberg FDR
  pvalueCutoff  = 0.05,            # 校正後 p 值門檻 (q ≤ 0.05) 才列入結果
  readable      = TRUE             # 將 Entrez 轉換為基因符號，易於閱讀
)

# ---------------------- 4-2  KEGG Pathway  ---------------------------
kegg_res <- enrichKEGG(
  gene          = entrez_vec,      # 同樣輸入 Entrez ID 向量
  organism      = "hsa",           # KEGG 物種代碼：hsa = Homo sapiens
  pvalueCutoff  = 0.05             # （單尾）p 值門檻；enrichKEGG 內建 BH 校正
)

# ---------------------- 4-3  Reactome Pathway  -----------------------
react_res <- enrichPathway(
  gene          = entrez_vec,      # Entrez 基因清單
  organism      = "human",         # 物種指定—Reactome 使用英文名稱
  pvalueCutoff  = 0.05,            # 校正後 p 值門檻 (BH by default)
  readable      = TRUE             # 轉換為易讀基因符號
)


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
