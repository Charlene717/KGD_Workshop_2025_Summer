###########################################
##  Basal keratinocytes ORA – R Script   ##
###########################################

## ============================= ##
## 0. 套件安裝與載入（無 pacman） ##
## ============================= ##
if (!require("tidyverse"))     { install.packages("tidyverse");     library(tidyverse) }
if (!require("clusterProfiler")){ install.packages("clusterProfiler"); library(clusterProfiler) }
if (!require("org.Hs.eg.db"))   { BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db) }
if (!require("ReactomePA"))     { BiocManager::install("ReactomePA");   library(ReactomePA) }
if (!require("ggplot2"))        { install.packages("ggplot2");          library(ggplot2) }

## ========================== ##
## 1. 擷取 Basal keratinocyte ##
## ========================== ##
# 假設 top25_lines 已存在（字串向量）
# 取出包含 "Basal keratinocytes" 的元素
basal_line  <- top25_lines[grepl("^Basal keratinocytes", names(top25_lines))]
# 去除前綴 "Basal keratinocytes: " 並拆成基因向量
Basal_genes <- basal_line %>%
  sub("^[^:]+:\\s*", "", .) %>%         # 去掉「Basal keratinocytes:」
  str_split(",\\s*") %>%                # 逗號切分
  unlist() %>%                          # 轉成向量
  unique()                              # 移除重複

# 檢查結果
print(Basal_genes)

## ======================================== ##
## 2. SYMBOL → ENTREZID 轉換（clusterProfiler::bitr） ##
## ======================================== ##
Basal_entrez <- bitr(Basal_genes,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Hs.eg.db) %>%
  pull(ENTREZID) %>% unique()

## ======================= ##
## 3. ORA – GO / KEGG / RE ##
## ======================= ##
# 統一輸出設定
output_dir    <- "."                     # 若想存到別處，改成 "path/to/dir"
output_prefix <- "BasalKeratinocytes_ORA"

# ---------- 3.1 GO ----------
go_res  <- enrichGO(gene           = Basal_entrez,
                    OrgDb          = org.Hs.eg.db,
                    keyType        = "ENTREZID",
                    ont            = "ALL",
                    pAdjustMethod  = "BH",
                    pvalueCutoff   = 0.05,
                    qvalueCutoff   = 0.10,
                    readable       = TRUE)

plot_GO <- dotplot(go_res, showCategory = 20) +
  ggtitle("GO enrichment – Basal keratinocytes")
plot_GO

# ---------- 3.2 KEGG ----------
kegg_res <- enrichKEGG(gene         = Basal_entrez,
                       organism     = "hsa",
                       pvalueCutoff = 0.05)

plot_KEGG <- dotplot(kegg_res, showCategory = 20) +
  ggtitle("KEGG pathways – Basal keratinocytes")
plot_KEGG

# ---------- 3.3 Reactome ----------
reac_res <- enrichPathway(gene         = Basal_entrez,
                          organism     = "human",
                          pvalueCutoff = 0.05,
                          readable     = TRUE)

plot_REAC <- dotplot(reac_res, showCategory = 20) +
  ggtitle("Reactome pathways – Basal keratinocytes")
plot_REAC

## ====================== ##
## 4.  儲存結果與圖檔       ##
## ====================== ##
# -- 4.1 表格 --
write.csv(go_res,   file = file.path(output_dir, paste0(output_prefix, "_GO.csv")))
write.csv(kegg_res, file = file.path(output_dir, paste0(output_prefix, "_KEGG.csv")))
write.csv(reac_res, file = file.path(output_dir, paste0(output_prefix, "_Reactome.csv")))

# -- 4.2 點圖（JPG） --
jpeg(file.path(output_dir, paste0(output_prefix, "_GO.jpg")),    width = 600, height = 800); print(plot_GO);   dev.off()
jpeg(file.path(output_dir, paste0(output_prefix, "_KEGG.jpg")),  width = 600, height = 800); print(plot_KEGG); dev.off()
jpeg(file.path(output_dir, paste0(output_prefix, "_Reactome.jpg")), width = 600, height = 800); print(plot_REAC); dev.off()

# -- 4.3 整合 PDF --
pdf(file.path(output_dir, paste0(output_prefix, "_DotPlots.pdf")), width = 7, height = 9)
try(print(plot_GO))
try(print(plot_KEGG))
try(print(plot_REAC))
dev.off()

message("✔ ORA completed and files saved to: ", normalizePath(output_dir))
