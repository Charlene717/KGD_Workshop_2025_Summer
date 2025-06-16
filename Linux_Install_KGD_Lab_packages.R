# Ref: https://chatgpt.com/c/6805e0a5-bdf8-8002-a169-612b25c3bff1
if(!require('devtools')) {install.packages('devtools'); library(devtools)}
# sudo chown -R charlene:charlene /opt/R/4.3.1/lib/R/library
# sudo chmod -R 775 /opt/R/4.3.1/lib/R/library



library(devtools)

# 套件與版本列表
target_pkgs <- data.frame(
  Package = c(
    "abind", "annotate", "askpass", "BH", "BiocIO", "biomaRt", "bitops", "caTools",
    "cli", "colorspace", "commonmark", "cpp11", "curl", "data.table", "DESeq2",
    "digest", "dotCall64", "DoubletFinder", "fastDummies", "fields", "fitdistrplus",
    "FNN", "fs", "future.apply", "genefilter", "geneplotter", "GEOquery", "ggrepel",
    "glue", "gtable", "HDF5Array", "igraph", "jsonlite", "later", "leidenbase",
    "locfit", "LoomExperiment", "maps", "matrixStats", "metap", "multcomp",
    "openssl", "openxlsx", "parallelly", "polyclip", "progressr", "promises", "purrr",
    "R.methodsS3", "R.oo", "R.utils", "R6", "RANN", "Rcpp", "RcppArmadillo",
    "RcppEigen", "remotes", "reticulate", "rhdf5", "rhdf5filters", "Rhdf5lib",
    "rPanglaoDB", "sandwich", "spam", "spatstat.explore", "spatstat.geom",
    "spatstat.random", "spatstat.utils", "stringi", "sys", "TCGAbiolinks",
    "TCGAbiolinksGUI.data", "uwot", "withr", "xfun", "yaml", "zoo"
  ),
  Version = c(
    "1.4-8", "1.72.0", "1.1", "1.87.0-1", "1.4.0", "2.50.3", "1.0-7", "1.18.2",
    "3.6.1", "2.1-0", "1.9.0", "0.5.2", "5.0.0", "1.14.8", "1.34.0",
    "0.6.31", "1.0-2", "2.0.4", "1.7.5", "14.1", "1.2-2",
    "1.1.3.2", "1.5.2", "1.11.3", "1.76.0", "1.72.0", "2.62.2", "0.9.3",
    "1.6.2", "0.3.6", "1.22.1", "1.4.2", "1.8.4", "1.3.0", "0.1.18",
    "1.5-9.12", "1.12.0", "3.4.1", "0.63.0", "1.11", "1.4-28",
    "2.0.6", "4.2.8", "1.35.0", "1.10-4", "0.15.1", "1.2.0.1", "1.0.1",
    "1.8.2", "1.27.0", "2.13.0", "2.6.1", "2.6.1", "1.0.13-1", "0.12.2.0.0",
    "0.3.3.9.3", "2.5.0", "1.28", "2.38.1", "1.6.0", "1.16.0",
    "0.2.1", "3.1-1", "2.9-1", "3.1-0", "3.1-0",
    "3.1-4", "3.1-2", "1.8.4", "3.4.1", "2.22.4",
    "1.14.1", "0.1.14", "3.0.2", "0.39", "2.3.7", "1.8-12"
  )
)

# 指定安裝路徑
install_path <- "/opt/R/4.3.1/lib/R/library"

# 設定安裝路徑
.libPaths(c(install_path, .libPaths()))

# 安裝特定版本的套件（修正 pkg 為 package）
for (i in seq_len(nrow(target_pkgs))) {
  pkg_name <- target_pkgs$Package[i]
  pkg_version <- target_pkgs$Version[i]

  try({
    if (!requireNamespace(pkg_name, quietly = TRUE) || packageVersion(pkg_name) != pkg_version) {
      devtools::install_version(package = pkg_name, version = pkg_version, lib = install_path, dependencies = TRUE, upgrade = "never")
    }

  })
}



################################################################################
as.data.frame(installed.packages()[, c("Package", "Version")])
# Package    Version
# abind                                   abind      1.4-8
# annotate                             annotate     1.72.0
# askpass                               askpass        1.1
# BH                                         BH   1.87.0-1
# BiocIO                                 BiocIO      1.4.0
# biomaRt                               biomaRt     2.50.3
# bitops                                 bitops      1.0-7
# caTools                               caTools     1.18.2
# cli                                       cli      3.6.1
# colorspace                         colorspace      2.1-0
# commonmark                         commonmark      1.9.0
# cpp11                                   cpp11      0.5.2
# curl                                     curl      5.0.0
# data.table                         data.table     1.14.8
# DESeq2                                 DESeq2     1.34.0
# digest                                 digest     0.6.31
# dotCall64                           dotCall64      1.0-2
# DoubletFinder                   DoubletFinder      2.0.4
# fastDummies                       fastDummies      1.7.5
# fields                                 fields       14.1
# fitdistrplus                     fitdistrplus      1.2-2
# FNN                                       FNN    1.1.3.2
# fs                                         fs      1.5.2
# future.apply                     future.apply     1.11.3
# genefilter                         genefilter     1.76.0
# geneplotter                       geneplotter     1.72.0
# GEOquery                             GEOquery     2.62.2
# ggrepel                               ggrepel      0.9.3
# glue                                     glue      1.6.2
# gtable                                 gtable      0.3.6
# HDF5Array                           HDF5Array     1.22.1
# igraph                                 igraph      1.4.2
# jsonlite                             jsonlite      1.8.4
# later                                   later      1.3.0
# leidenbase                         leidenbase     0.1.18
# locfit                                 locfit   1.5-9.12
# LoomExperiment                 LoomExperiment     1.12.0
# maps                                     maps      3.4.1
# matrixStats                       matrixStats     0.63.0
# metap                                   metap       1.11
# multcomp                             multcomp     1.4-28
# openssl                               openssl      2.0.6
# openxlsx                             openxlsx      4.2.8
# parallelly                         parallelly     1.35.0
# polyclip                             polyclip     1.10-4
# progressr                           progressr     0.15.1
# promises                             promises    1.2.0.1
# purrr                                   purrr      1.0.1
# R.methodsS3                       R.methodsS3      1.8.2
# R.oo                                     R.oo     1.27.0
# R.utils                               R.utils     2.13.0
# R6                                         R6      2.6.1
# RANN                                     RANN      2.6.1
# Rcpp                                     Rcpp   1.0.13-1
# RcppArmadillo                   RcppArmadillo 0.12.2.0.0
# RcppEigen                           RcppEigen  0.3.3.9.3
# remotes                               remotes      2.5.0
# reticulate                         reticulate       1.28
# rhdf5                                   rhdf5     2.38.1
# rhdf5filters                     rhdf5filters      1.6.0
# Rhdf5lib                             Rhdf5lib     1.16.0
# rPanglaoDB                         rPanglaoDB      0.2.1
# sandwich                             sandwich      3.1-1
# spam                                     spam      2.9-1
# spatstat.explore             spatstat.explore      3.1-0
# spatstat.geom                   spatstat.geom      3.1-0
# spatstat.random               spatstat.random      3.1-4
# spatstat.utils                 spatstat.utils      3.1-2
# stringi                               stringi      1.8.4
# sys                                       sys      3.4.1
# TCGAbiolinks                     TCGAbiolinks     2.22.4
# TCGAbiolinksGUI.data     TCGAbiolinksGUI.data     1.14.1
# uwot                                     uwot     0.1.14
# withr                                   withr      3.0.2
# xfun                                     xfun       0.39
# yaml                                     yaml      2.3.7
# zoo                                       zoo     1.8-12
# abind.1                                 abind      1.4-5
# AnnotationDbi                   AnnotationDbi     1.56.2
# AnnotationHub                   AnnotationHub      3.2.2
# ape                                       ape        5.8
# aplot                                   aplot      0.2.3
# askpass.1                             askpass      1.2.0
# assertthat                         assertthat      0.2.1
# backports                           backports      1.5.0
# base                                     base      4.1.3
# base64enc                           base64enc      0.1-3
# beachmat                             beachmat     2.10.0
# beeswarm                             beeswarm      0.4.0
# bench                                   bench      1.1.3
# BH.1                                       BH   1.84.0-0
# Biobase                               Biobase     2.54.0
# BiocFileCache                   BiocFileCache      2.2.1
# BiocGenerics                     BiocGenerics     0.40.0
# BiocManager                       BiocManager    1.30.23
# BiocNeighbors                   BiocNeighbors     1.12.0
# BiocParallel                     BiocParallel     1.28.3
# BiocSingular                     BiocSingular     1.10.0
# BiocVersion                       BiocVersion     3.14.0
# biocViews                           biocViews     1.62.1
# Biostrings                         Biostrings     2.62.0
# bit                                       bit      4.0.5
# bit64                                   bit64      4.0.5
# bitops.1                               bitops      1.0-7
# blob                                     blob      1.2.4
# boot                                     boot     1.3-28
# brew                                     brew      1.0-8
# brio                                     brio      1.1.5
# broom                                   broom      1.0.6
# bslib                                   bslib      0.7.0
# cachem                                 cachem      1.1.0
# Cairo                                   Cairo      1.6-0
# callr                                   callr      3.7.6
# car                                       car      3.1-2
# carData                               carData      3.0-5
# caTools.1                             caTools     1.18.2
# CellChat                             CellChat      1.6.1
# celldex                               celldex      1.4.0
# CelliD                                 CelliD     1.12.0
# cellranger                         cellranger      1.1.0
# circlize                             circlize     0.4.16
# class                                   class     7.3-20
# cli.1                                     cli      3.6.3
# clipr                                   clipr      0.8.0
# clue                                     clue     0.3-65
# cluster                               cluster      2.1.2
# clusterProfiler               clusterProfiler      4.2.2
# coda                                     coda   0.19-4.1
# codetools                           codetools     0.2-18
# colorspace.1                       colorspace      2.1-0
# combinat                             combinat      0.0-8
# commonmark.1                       commonmark      1.9.1
# compiler                             compiler      4.1.3
# ComplexHeatmap                 ComplexHeatmap     2.10.0
# conflicted                         conflicted      1.2.0
# corrplot                             corrplot       0.92
# covr                                     covr      3.6.4
# cowplot                               cowplot      1.1.3
# cpp11.1                                 cpp11      0.4.7
# crayon                                 crayon      1.5.3
# credentials                       credentials      1.3.2
# crosstalk                           crosstalk      1.2.1
# curl.1                                   curl      5.2.1
# data.table.1                       data.table     1.15.4
# datasets                             datasets      4.1.3
# DBI                                       DBI      1.2.3
# dbplyr                                 dbplyr      2.3.4
# DDRTree                               DDRTree      0.1.5
# DelayedArray                     DelayedArray     0.20.0
# DelayedMatrixStats         DelayedMatrixStats     1.16.0
# deldir                                 deldir      2.0-4
# densityClust                     densityClust      0.3.3
# DEoptimR                             DEoptimR    1.1-3-1
# desc                                     desc      1.4.3
# devtools                             devtools      2.4.5
# diffobj                               diffobj      0.3.5
# digest.1                               digest     0.6.36
# DO.db                                   DO.db        2.9
# docopt                                 docopt      0.7.1
# doParallel                         doParallel     1.0.17
# DOSE                                     DOSE     3.20.1
# dotCall64.1                         dotCall64      1.1-1
# DoubletFinder.1                 DoubletFinder      2.0.4
# downlit                               downlit      0.4.2
# downloader                         downloader        0.4
# dplyr                                   dplyr      1.1.4
# dqrng                                   dqrng      0.4.1
# DT                                         DT       0.33
# dtplyr                                 dtplyr      1.3.1
# ellipsis                             ellipsis      0.3.2
# enrichplot                         enrichplot     1.14.2
# evaluate                             evaluate     0.24.0
# ExperimentHub                   ExperimentHub      2.2.1
# expm                                     expm    0.999-9
# fansi                                   fansi      1.0.6
# farver                                 farver      2.1.2
# fastDummies.1                     fastDummies      1.7.3
# fastICA                               fastICA      1.2-4
# fastmap                               fastmap      1.2.0
# fastmatch                           fastmatch      1.1-4
# fgsea                                   fgsea     1.20.0
# fields.1                               fields       14.1
# filelock                             filelock      1.0.3
# fitdistrplus.1                   fitdistrplus      1.2-1
# FNN.1                                     FNN      1.1.4
# fontawesome                       fontawesome      0.5.2
# forcats                               forcats      1.0.0
# foreach                               foreach      1.5.2
# foreign                               foreign     0.8-82
# formatR                               formatR       1.14
# Formula                               Formula      1.2-5
# fs.1                                       fs      1.6.4
# futile.logger                   futile.logger      1.4.3
# futile.options                 futile.options      1.0.1
# future                                 future     1.33.2
# future.apply.1                   future.apply     1.11.2
# gargle                                 gargle      1.5.2
# generics                             generics      0.1.3
# GenomeInfoDb                     GenomeInfoDb     1.30.1
# GenomeInfoDbData             GenomeInfoDbData      1.2.7
# GenomicRanges                   GenomicRanges     1.46.1
# gert                                     gert      1.9.2
# GetoptLong                         GetoptLong      1.0.5
# ggalluvial                         ggalluvial     0.12.5
# ggbeeswarm                         ggbeeswarm      0.7.1
# ggforce                               ggforce      0.4.2
# ggfun                                   ggfun      0.1.5
# ggnetwork                           ggnetwork     0.5.13
# ggnewscale                         ggnewscale     0.4.10
# ggplot2                               ggplot2      3.5.1
# ggplotify                           ggplotify      0.1.2
# ggprism                               ggprism      1.0.4
# ggpubr                                 ggpubr      0.6.0
# ggraph                                 ggraph      2.1.0
# ggrastr                               ggrastr      1.0.1
# ggrepel.1                             ggrepel      0.9.5
# ggridges                             ggridges      0.5.6
# ggsci                                   ggsci      3.2.0
# ggsignif                             ggsignif      0.6.4
# ggtree                                 ggtree      3.2.1
# gh                                         gh      1.4.0
# gitcreds                             gitcreds      0.1.2
# GlobalOptions                   GlobalOptions      0.1.2
# globals                               globals     0.16.3
# glue.1                                   glue      1.7.0
# GO.db                                   GO.db     3.14.0
# goftest                               goftest      1.2-3
# googledrive                       googledrive      2.1.1
# googlesheets4                   googlesheets4      1.1.1
# GOSemSim                             GOSemSim     2.20.0
# gplots                                 gplots    3.1.3.1
# GPTCelltype                       GPTCelltype      1.0.1
# graph                                   graph     1.72.0
# graphics                             graphics      4.1.3
# graphlayouts                     graphlayouts      1.1.1
# grDevices                           grDevices      4.1.3
# grid                                     grid      4.1.3
# gridBase                             gridBase      0.4-7
# gridExtra                           gridExtra        2.3
# gridGraphics                     gridGraphics      0.5-1
# gtable.1                               gtable      0.3.5
# gtools                                 gtools      3.9.5
# haven                                   haven      2.5.4
# hdf5r                                   hdf5r     1.3.11
# here                                     here      1.0.1
# highr                                   highr       0.11
# hms                                       hms      1.1.3
# HSMMSingleCell                 HSMMSingleCell     1.14.0
# htmltools                           htmltools    0.5.8.1
# htmlwidgets                       htmlwidgets      1.6.4
# httpuv                                 httpuv     1.6.15
# httr                                     httr      1.4.7
# httr2                                   httr2      0.2.2
# ica                                       ica      1.0-3
# ids                                       ids      1.0.1
# igraph.1                               igraph      1.4.2
# ini                                       ini      0.3.1
# interactiveDisplayBase interactiveDisplayBase     1.32.0
# IRanges                               IRanges     2.28.0
# irlba                                   irlba    2.3.5.1
# isoband                               isoband      0.2.7
# iterators                           iterators     1.0.14
# janitor                               janitor      2.2.0
# jquerylib                           jquerylib      0.1.4
# jsonlite.1                           jsonlite      1.8.8
# KEGGREST                             KEGGREST     1.34.0
# KernSmooth                         KernSmooth    2.23-20
# knitr                                   knitr       1.48
# labeling                             labeling      0.4.3
# Lahman                                 Lahman     11.0-0
# lambda.r                             lambda.r      1.2.4
# later.1                                 later      1.3.2
# lattice                               lattice    0.20-45
# lazyeval                             lazyeval      0.2.2
# leiden                                 leiden      0.4.3
# lifecycle                           lifecycle      1.0.4
# limma                                   limma     3.50.3
# listenv                               listenv      0.9.1
# lme4                                     lme4   1.1-35.5
# lmtest                                 lmtest     0.9-40
# lobstr                                 lobstr      1.1.2
# loupeR                                 loupeR      1.0.2
# lubridate                           lubridate      1.9.3
# magrittr                             magrittr      2.0.3
# maps.1                                   maps      3.4.1
# MASS                                     MASS     7.3-55
# mathjaxr                             mathjaxr      1.6-0
# Matrix                                 Matrix      1.5-0
# MatrixGenerics                 MatrixGenerics      1.6.0
# MatrixModels                     MatrixModels      0.5-1
# matrixStats.1                     matrixStats      1.1.0
# memoise                               memoise      2.0.1
# methods                               methods      4.1.3
# mgcv                                     mgcv     1.8-39
# microbenchmark                 microbenchmark     1.4.10
# mime                                     mime       0.12
# miniUI                                 miniUI    0.1.1.1
# minqa                                   minqa      1.2.7
# mnormt                                 mnormt      2.1.1
# modelr                                 modelr     0.1.11
# monocle                               monocle     2.22.0
# multtest                             multtest     2.50.0
# munsell                               munsell      0.5.1
# mutoss                                 mutoss     0.1-13
# mvtnorm                               mvtnorm      1.1-3
# network                               network     1.18.2
# nlme                                     nlme    3.1-155
# nloptr                                 nloptr      2.1.1
# NMF                                       NMF       0.27
# nnet                                     nnet     7.3-17
# numDeriv                             numDeriv 2016.8-1.1
# nycflights13                     nycflights13      1.0.2
# openai                                 openai      0.4.1
# openssl.1                             openssl      2.2.0
# paletteer                           paletteer      1.5.0
# parallel                             parallel      4.1.3
# parallelly.1                       parallelly     1.37.1
# patchwork                           patchwork      1.2.0
# pbapply                               pbapply      1.7-2
# pheatmap                             pheatmap     1.0.12
# pillar                                 pillar      1.9.0
# pkgbuild                             pkgbuild      1.4.4
# pkgconfig                           pkgconfig      2.0.3
# pkgdown                               pkgdown      2.0.7
# pkgload                               pkgload      1.4.0
# plogr                                   plogr      0.2.0
# plotly                                 plotly     4.10.4
# plotrix                               plotrix      3.8-4
# plyr                                     plyr      1.8.9
# png                                       png      0.1-8
# polyclip.1                           polyclip     1.10-6
# polynom                               polynom      1.4-1
# praise                                 praise      1.0.0
# prettyunits                       prettyunits      1.2.0
# prismatic                           prismatic      1.1.1
# processx                             processx      3.8.4
# profmem                               profmem      0.6.0
# profvis                               profvis      0.3.7
# progress                             progress      1.2.3
# progressr.1                         progressr     0.14.0
# promises.1                           promises      1.3.0
# proxy                                   proxy     0.4-27
# ps                                         ps      1.7.7
# purrr.1                                 purrr      1.0.2
# qlcMatrix                           qlcMatrix      0.9.8
# qqconf                                 qqconf      1.3.2
# qs                                         qs     0.26.3
# quantreg                             quantreg       5.98
# qvalue                                 qvalue     2.26.0
# R6.1                                       R6      2.5.1
# ragg                                     ragg      1.3.2
# RANN.1                                   RANN      2.6.1
# RApiSerialize                   RApiSerialize      0.1.3
# rappdirs                             rappdirs      0.3.3
# RBGL                                     RBGL     1.70.0
# rbibutils                           rbibutils     2.2.13
# rcmdcheck                           rcmdcheck      1.4.0
# RColorBrewer                     RColorBrewer      1.1-3
# Rcpp.1                                   Rcpp     1.0.13
# RcppAnnoy                           RcppAnnoy     0.0.22
# RcppArmadillo.1                 RcppArmadillo   14.0.0-1
# RcppEigen.1                         RcppEigen  0.3.4.0.0
# RcppHNSW                             RcppHNSW      0.6.0
# RcppParallel                     RcppParallel      5.1.8
# RcppProgress                     RcppProgress      0.4.2
# RcppTOML                             RcppTOML      0.2.2
# RCurl                                   RCurl  1.98-1.12
# Rdpack                                 Rdpack      2.6.2
# readr                                   readr      2.1.5
# readxl                                 readxl      1.4.3
# reformulas                         reformulas      0.4.0
# registry                             registry      0.5-1
# rematch                               rematch      2.0.0
# rematch2                             rematch2      2.1.2
# remotes.1                             remotes      2.5.0
# reprex                                 reprex      2.1.1
# reshape                               reshape      0.8.9
# reshape2                             reshape2      1.4.4
# reticulate.1                       reticulate     1.38.0
# rex                                       rex      1.2.1
# rjson                                   rjson     0.2.21
# rlang                                   rlang      1.1.4
# rmarkdown                           rmarkdown       2.27
# RMySQL                                 RMySQL    0.10.27
# rngtools                             rngtools      1.5.2
# robustbase                         robustbase     0.95-1
# ROCR                                     ROCR     1.0-11
# roxygen2                             roxygen2      7.2.3
# rpart                                   rpart     4.1.16
# RPostgreSQL                       RPostgreSQL      0.7-6
# rprojroot                           rprojroot      2.0.4
# RSpectra                             RSpectra     0.16-2
# RSQLite                               RSQLite      2.3.7
# rstatix                               rstatix      0.7.2
# rstudioapi                         rstudioapi     0.16.0
# rsvd                                     rsvd      1.0.5
# Rtsne                                   Rtsne       0.17
# RUnit                                   RUnit     0.4.33
# rversions                           rversions      2.1.2
# rvest                                   rvest      1.0.4
# S4Vectors                           S4Vectors     0.32.4
# sandwich.1                           sandwich      3.0-2
# sass                                     sass      0.4.9
# ScaledMatrix                     ScaledMatrix      1.2.0
# scales                                 scales      1.3.0
# scater                                 scater     1.22.0
# scattermore                       scattermore        1.2
# scatterpie                         scatterpie      0.2.3
# scCATCH                               scCATCH      3.2.2
# scCustomize                       scCustomize      1.1.3
# sctransform                       sctransform      0.4.1
# scuttle                               scuttle      1.4.0
# selectr                               selectr      0.4-2
# sessioninfo                       sessioninfo      1.2.2
# Seurat                                 Seurat      4.3.0
# SeuratObject                     SeuratObject      4.1.3
# shadowtext                         shadowtext      0.1.4
# shape                                   shape    1.4.6.1
# shiny                                   shiny    1.8.1.1
# SingleCellExperiment     SingleCellExperiment     1.16.0
# SingleR                               SingleR      1.8.1
# sitmo                                   sitmo      2.0.2
# slam                                     slam     0.1-51
# sn                                         sn      2.1.1
# sna                                       sna      2.7-2
# snakecase                           snakecase     0.11.0
# snow                                     snow      0.4-4
# sourcetools                       sourcetools    0.1.7-1
# sp                                         sp      2.1-4
# spam.1                                   spam     2.10-0
# SparseM                               SparseM     1.84-2
# sparseMatrixStats           sparseMatrixStats      1.6.0
# sparsesvd                           sparsesvd      0.2-2
# spatial                               spatial     7.3-15
# spatstat.data                   spatstat.data      3.1-2
# spatstat.explore.1           spatstat.explore      3.3-1
# spatstat.geom.1                 spatstat.geom      3.3-2
# spatstat.random.1             spatstat.random      3.3-1
# spatstat.sparse               spatstat.sparse      3.1-0
# spatstat.univar               spatstat.univar      3.0-0
# spatstat.utils.1               spatstat.utils      3.0-5
# splines                               splines      4.1.3
# statnet.common                 statnet.common      4.9.0
# stats                                   stats      4.1.3
# stats4                                 stats4      4.1.3
# stringfish                         stringfish     0.16.0
# stringi.1                             stringi      1.8.4
# stringr                               stringr      1.5.1
# SummarizedExperiment     SummarizedExperiment     1.24.0
# survival                             survival     3.2-13
# svglite                               svglite      2.1.3
# sys.1                                     sys      3.4.2
# systemfonts                       systemfonts      1.1.0
# tcltk                                   tcltk      4.1.3
# tensor                                 tensor        1.5
# testthat                             testthat    3.2.1.1
# textshaping                       textshaping      0.4.0
# TFisher                               TFisher      0.2.0
# TH.data                               TH.data      1.1-2
# tibble                                 tibble      3.2.1
# tictoc                                 tictoc      1.2.1
# tidygraph                           tidygraph      1.2.3
# tidyr                                   tidyr      1.3.1
# tidyselect                         tidyselect      1.2.1
# tidytree                             tidytree      0.4.6
# tidyverse                           tidyverse      2.0.0
# timechange                         timechange      0.3.0
# tinytex                               tinytex       0.51
# tools                                   tools      4.1.3
# translations                     translations      4.1.3
# treeio                                 treeio     1.18.1
# tweenr                                 tweenr      2.0.3
# tzdb                                     tzdb      0.4.0
# umap                                     umap   0.2.10.0
# urlchecker                         urlchecker      1.0.1
# usethis                               usethis      2.1.6
# utf8                                     utf8      1.2.4
# utils                                   utils      4.1.3
# uuid                                     uuid      1.2-0
# uwot.1                                   uwot     0.1.14
# vctrs                                   vctrs      0.6.5
# VGAM                                     VGAM     1.1-11
# vipor                                   vipor      0.4.5
# viridis                               viridis      0.6.2
# viridisLite                       viridisLite      0.4.2
# vroom                                   vroom      1.6.5
# waldo                                   waldo      0.5.2
# whisker                               whisker      0.4.1
# withr.1                                 withr      3.0.0
# xfun.1                                   xfun       0.46
# XML                                       XML  3.99-0.14
# xml2                                     xml2      1.3.6
# xopen                                   xopen      1.0.0
# xtable                                 xtable      1.8-4
# XVector                               XVector     0.34.0
# yaml.1                                   yaml      2.3.9
# yulab.utils                       yulab.utils      0.1.4
# zip                                       zip      2.3.0
# zlibbioc                             zlibbioc     1.40.0


################################################################################
################################################################################



# 指定安裝路徑
install_path <- "/opt/R/4.3.1/lib/R/library"

# 設定安裝路徑
.libPaths(c(install_path, .libPaths()))

devtools::install_version(package = "Matrix", version = "1.6.4", lib = install_path, dependencies = TRUE, upgrade = "never")

devtools::install_version(package = "TFisher", version = "0.2.0", lib = install_path, dependencies = TRUE, upgrade = "never")

devtools::install_version(package = "parallelly", version = "1.43.0", lib = install_path, dependencies = TRUE, upgrade = "never")
devtools::install_version(package = "promises", version = "1.3.2", lib = install_path, dependencies = TRUE, upgrade = "never")

devtools::install_version(package = "Seurat", version = "4.3.0", lib = install_path, dependencies = TRUE, upgrade = "never")



##################################################################
devtools::install_version(package = "tidyverse", version = "2.0.0", lib = install_path, dependencies = TRUE, upgrade = "never")

devtools::install_version(package = "ggpubr", version = "0.6.0", lib = install_path, dependencies = TRUE, upgrade = "never")


##################################################################
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
#
# # 從 Bioconductor 的 GitHub 存檔安裝 fgsea 1.20.0
# remotes::install_github("ctlab/fgsea@v1.20.0")




# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# # 指定Bioconductor歷史版本安裝
# BiocManager::install("fgsea", version = "3.14")
# # 錯誤: Bioconductor version '3.14' requires R version '4.1'; use `version = '3.18'` with R
# # version 4.3; see https://bioconductor.org/install


# Sys.setenv("CXX14FLAGS"="-std=c++14")
# BiocManager::install("fgsea")  # 自動裝最新版

# BiocManager::install("fgsea", type="binary")
# BiocManager::install("fgsea", version = "3.17")

# install.packages("BH")
# Sys.setenv("CXX14FLAGS"="-std=c++14")
# BiocManager::install("fgsea")  # 自動裝最新版

# install.packages("https://bioconductor.org/packages/3.17/bioc/src/contrib/fgsea_1.26.0.tar.gz",
#                  repos=NULL, type="source")
# library(fgsea)

install.packages("BiocManager")
BiocManager::install("fgsea")
library(fgsea)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# devtools::install_version(package = "CelliD", version = "1.12.0", lib = install_path, dependencies = TRUE, upgrade = "never")
BiocManager::install("CelliD", version = "3.18", lib = install_path, dependencies = TRUE, update = FALSE)


###################################################################
as.data.frame(installed.packages()[, c("Package", "Version")])


