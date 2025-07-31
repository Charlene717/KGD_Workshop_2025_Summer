############################################################
##  教學版 GSEA（含關鍵字經典曲線） – v2025-07-31        ##
############################################################
## 必裝套件版本：clusterProfiler ≥ 4.0，enrichplot ≥ 1.14 ##
############################################################

## ===== 0. 套件安裝與載入 =================================
pkgs <- c("clusterProfiler","enrichplot","msigdbr",  # 富集分析主套件
          "org.Hs.eg.db","tidyverse",                # 物種註解 & 數據整理
          "patchwork","ggplot2")                     # 圖形排版 & 圖形繪製
for (p in pkgs){
  if (!requireNamespace(p, quietly = TRUE)){         # 若尚未安裝…
    if (p %in% c("clusterProfiler","enrichplot",
                 "msigdbr","org.Hs.eg.db")){
      BiocManager::install(p, ask = FALSE)           # Bioconductor 套件
    } else install.packages(p, quiet = TRUE)         # CRAN 套件
  }
  library(p, character.only = TRUE)                  # 載入套件
}

## ===== 0.1 基本參數設定 ==================================
cell_type      <- "Basal keratinocytes"        # 目標細胞群名稱
score_col      <- "avg_log2FC"                 # 排名向量用的欄位
output_dir     <- "."                          # 輸出資料夾（"." = 工作目錄）
output_prefix  <- gsub("\\s+","_",cell_type)   # 檔名前綴，空白改底線
species        <- "Homo sapiens"               # MSigDB 物種名稱
OrgDb          <- org.Hs.eg.db                 # 對應的 Annotation 資料庫

## ===== 1. 取出目標細胞群的 DEG 表 ========================
deg_df <- deg_list[[cell_type]] |>              # 從預先算好的 deg_list 取表
  dplyr::select(gene, !!score_col) |>           # 只保留基因與分數欄
  rename(Gene = gene,                           # 改欄位名稱統一
         log2FoldChange = !!score_col)

## ===== 2. SYMBOL → ENTREZ，並建立 ranking vector =========
sym2ent <- bitr(deg_df$Gene,                    # 批次 ID 轉換
                fromType = "SYMBOL",            # 來源是基因符號
                toType   = "ENTREZID",          # 目標為 ENTREZ ID
                OrgDb    = OrgDb,               # 物種資料庫
                drop     = FALSE)               # 保留無法對應者（NA）

deg_df <- left_join(deg_df, sym2ent,            # 合併轉換結果
                    by = c("Gene" = "SYMBOL")) |>
  drop_na(ENTREZID) |>                          # 移除沒有 ENTREZID 的列
  distinct(ENTREZID, .keep_all = TRUE)          # 確保 ENTREZID 唯一

gene_rank <- deg_df$log2FoldChange              # 取分數向量
names(gene_rank) <- deg_df$ENTREZID             # 向量名稱設為 ENTREZID
gene_rank <- sort(gene_rank, decreasing = TRUE) # 由大到小排序

## ===== 3. 準備 MSigDB TERM2GENE 表（保持數值型 ENTREZ）===
gmt_df <- bind_rows(
  msigdbr(species, category = "H"),             # Hallmark gene sets
  # msigdbr(species, category = "C2"),          # Canonical pathways
  # msigdbr(species, category = "C3"),          # TF targets
  # msigdbr(species, category = "C7")           # Immunologic signatures
) |>
  dplyr::select(gs_name, ENTREZID = entrez_gene)|># 取 gene set 名與 ENTREZ
  distinct()                                    # 去重；ENTREZ 仍為 integer

## ===== 4. 執行 GSEA 分析（不設 p-value 篩選）==============
gsea_res <- GSEA(geneList     = gene_rank,      # 排名向量
                 TERM2GENE    = gmt_df,         # 基因集對照表
                 pvalueCutoff = 1,              # 不過濾，完整結果
                 minGSSize    = 10,             # 最小基因集大小
                 maxGSSize    = 500,            # 最大基因集大小
                 verbose      = FALSE)          # 不顯示進度

## ===== 5. dotplot / ridgeplot 總覽 ========================
dot_top30 <- dotplot(gsea_res, showCategory = 30) +         # NES dotplot
  ggtitle(paste("GSEA –", cell_type))
ridge_top30 <- ridgeplot(gsea_res, showCategory = 30) +      # 圓滑密度圖
  ggtitle(paste("GSEA ridgeplot –", cell_type))

ggsave(file.path(output_dir,                           # 儲存 dotplot
                 paste0(output_prefix,
                        "_GSEA_dotplot.jpg")),
       dot_top30, width = 7, height = 12, dpi = 300)
ggsave(file.path(output_dir,                           # 儲存 ridgeplot
                 paste0(output_prefix,
                        "_GSEA_ridgeplot.jpg")),
       ridge_top30, width = 7, height = 12, dpi = 300)

## ===== 6. 依關鍵字繪製經典 GSEA 曲線 ======================
kw_vec        <- c("NFKB","TNFA","WNT","TGFB")   # 關鍵字集合
max_term_plot <- 10                              # 最多曲線數

## 6-1. 篩選符合關鍵字的 pathways
hit_tbl <- gsea_res@result |>
  filter(str_detect(toupper(ID),                   # 對 ID 做大寫比對
                    paste(kw_vec, collapse = "|"))) |>
  arrange(p.adjust) |>                             # 依調整後 p 排序
  slice_head(n = max_term_plot)                    # 取前 N 條

if (nrow(hit_tbl) == 0){
  warning("⚠ 找不到符合關鍵字的 pathway；請檢查 kw_vec")
} else {
  
  ## 6-2. 定義安全版 gseaplot2（相容新版無 combine 參數）
  library(patchwork)                               # 用於排版
  safe_gseaplot <- function(gsea_obj, term_id,
                            line_col = "steelblue"){
    tryCatch({
      p <- enrichplot::gseaplot2(
        gsea_obj,
        geneSetID    = term_id,
        ES_geom      = "line",
        pvalue_table = TRUE,                       # NES 與 padj 顯示
        color        = line_col)
      if (!inherits(p, "patchwork"))               # 舊版回 ggplot
        p <- patchwork::wrap_plots(p, ncol = 1)    # 統一轉 patchwork
      p + ggtitle(term_id) +
        theme(plot.title = element_text(
          hjust = .5, face = "bold"))
    }, error = function(e){
      message("❗ 無法繪製 ", term_id, "：", e$message)
      NULL
    })
  }
  
  ## 6-3. 批次產出曲線圖
  gsea_plots <- lapply(hit_tbl$ID, safe_gseaplot,
                       gsea_obj = gsea_res)
  names(gsea_plots) <- hit_tbl$ID
  gsea_plots <- gsea_plots[!sapply(gsea_plots, is.null)]  # 去掉失敗者
  
  ## 6-4. PDF（每頁 1 圖）
  pdf(file.path(output_dir,
                paste0(output_prefix,
                       "_GSEA_keywordCurves.pdf")),
      width = 10, height = 8)
  for (p in gsea_plots) print(p)                  # 逐頁輸出
  dev.off()
  
  ## 6-5. 個別 PNG（300 dpi）
  for (nm in names(gsea_plots)){
    ggsave(file.path(output_dir,
                     paste0(output_prefix,
                            "_GSEA_", nm, ".png")),
           gsea_plots[[nm]],
           width = 10, height = 8, dpi = 300)
  }
  
  cat("✔ 已輸出", length(gsea_plots),
      "張關鍵字 GSEA 曲線圖至：",
      normalizePath(output_dir), "\n")
}
