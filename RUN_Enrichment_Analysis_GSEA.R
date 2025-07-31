############################################################
##  教學版 GSEA（含關鍵字經典曲線） – v2025-07-31        ##
############################################################
## 必裝：clusterProfiler >= 4.0，enrichplot >= 1.14
pkgs <- c("clusterProfiler","enrichplot","msigdbr",
          "org.Hs.eg.db","tidyverse","patchwork","ggplot2")
for (p in pkgs){
  if (!requireNamespace(p, quietly = TRUE)){
    if (p %in% c("clusterProfiler","enrichplot","msigdbr","org.Hs.eg.db")){
      BiocManager::install(p, ask = FALSE)
    } else install.packages(p, quiet = TRUE)
  }
  library(p, character.only = TRUE)
}

## ---------- 0. 基本參數 ---------- ##
cell_type      <- "Basal keratinocytes"        # 想分析的群
score_col      <- "avg_log2FC"                 # 排序分數欄
output_dir     <- "."                          # 輸出目錄
output_prefix  <- gsub("\\s+","_",cell_type)   # 檔名前綴
species        <- "Homo sapiens"               # 物種
OrgDb          <- org.Hs.eg.db


## ---------- 1. 取該群 DE 資料框 ---------- ##
deg_df <- deg_list[[cell_type]] |>
  dplyr::select(gene, !!score_col) |>
  rename(Gene = gene,
         log2FoldChange = !!score_col)


## ---------- 2. SYMBOL → ENTREZ，建排名向量 ---------- ##
sym2ent <- bitr(deg_df$Gene, fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = OrgDb,
                drop     = FALSE)

deg_df <- left_join(deg_df, sym2ent,
                    by = c("Gene" = "SYMBOL")) |>
  drop_na(ENTREZID) |>
  distinct(ENTREZID, .keep_all = TRUE)  # ENTREZID 唯一

gene_rank <- deg_df$log2FoldChange
names(gene_rank) <- deg_df$ENTREZID           # **保持 numeric 名稱**
gene_rank       <- sort(gene_rank, decreasing = TRUE)


## ---------- 3. MSigDB TERM2GENE（數值型 ENTREZ） ---------- ##
gmt_df <- bind_rows(
  msigdbr(species, category="H"),    # Hallmark
  # msigdbr(species, category="C2"),   # Canonical pathways
  # msigdbr(species, category="C3"),   # TF targets
  # msigdbr(species, category="C7")    # Immunologic
) |>
  dplyr::select(gs_name, ENTREZID = entrez_gene) |>
  distinct()                            # ENTREZID 仍為 integer


## ---------- 4. 執行 GSEA（不過濾） ---------- ##
gsea_res <- GSEA(geneList     = gene_rank,
                 TERM2GENE    = gmt_df,
                 pvalueCutoff = 1,
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 verbose      = FALSE)


## ---------- 5. dotplot / ridgeplot 總覽 ---------- ##
dot_top30   <- dotplot(gsea_res, showCategory = 30) +
  ggtitle(paste("GSEA –", cell_type))
ridge_top30 <- ridgeplot(gsea_res, showCategory = 30) +
  ggtitle(paste("GSEA ridgeplot –", cell_type))

ggsave(file.path(output_dir,
                 paste0(output_prefix,"_GSEA_dotplot.jpg")),
       dot_top30, width = 7, height = 12, dpi = 300)
ggsave(file.path(output_dir,
                 paste0(output_prefix,"_GSEA_ridgeplot.jpg")),
       ridge_top30, width = 7, height = 12, dpi = 300)


## ---------- 6. 關鍵字經典曲線 ---------- ##
kw_vec        <- c("NFKB","TNFA","WNT","TGFB")   # 自訂關鍵字
max_term_plot <- 10                              # 最多畫幾條

## 6-1 找出符合關鍵字的 term
hit_tbl <- gsea_res@result |>
  filter(str_detect(toupper(ID),
                    paste(kw_vec, collapse = "|"))) |>
  arrange(p.adjust) |>
  slice_head(n = max_term_plot)

if (nrow(hit_tbl) == 0){
  warning("⚠ 找不到符合關鍵字的 pathway；請檢查 kw_vec")
} else {
  
  ## 6-2 產生 gseaplot2（無 combine 參數）
  library(patchwork)
  
  safe_gseaplot <- function(gsea_obj, term_id,
                            line_col = "steelblue"){
    tryCatch({
      p <- enrichplot::gseaplot2(
        gsea_obj,
        geneSetID    = term_id,
        ES_geom      = "line",
        pvalue_table = TRUE,
        color        = line_col)
      
      ## enrichplot <1.21 回傳 ggplot / list，
      ## ≥1.21 回傳 patchwork – 用 wrap_plots 統一
      if (!inherits(p, "patchwork"))
        p <- patchwork::wrap_plots(p, ncol = 1)
      
      p + ggtitle(term_id) +
        theme(plot.title = element_text(
          hjust = .5, face = "bold"))
    },
    error = function(e){
      message("❗ 無法繪製 ", term_id, "：", e$message)
      NULL
    })
  }
  
  gsea_plots <- lapply(hit_tbl$ID,
                       safe_gseaplot,
                       gsea_obj = gsea_res)
  names(gsea_plots) <- hit_tbl$ID
  gsea_plots <- gsea_plots[!sapply(gsea_plots, is.null)]
  
  ## 6-3 輸出 PDF（每頁 1 圖）
  pdf(file.path(output_dir,
                paste0(output_prefix,
                       "_GSEA_keywordCurves.pdf")),
      width = 10, height = 8)
  for (p in gsea_plots) print(p)
  dev.off()
  
  ## 6-4 輸出個別 PNG
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
