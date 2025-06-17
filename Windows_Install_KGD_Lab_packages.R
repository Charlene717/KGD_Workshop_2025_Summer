# Ref: https://chatgpt.com/c/6805e0a5-bdf8-8002-a169-612b25c3bff1
if(!require('devtools')) {install.packages('devtools'); library(devtools)}


#### 套件安裝 ####
# ❌ BiocIO (1.4.0)
# ❌ biomaRt (2.50.3)
# ❌ DESeq2 (1.34.0)
# ❌ DoubletFinder (2.0.4)

# ❌ genefilter (1.76.0)
# ❌ geneplotter (1.72.0)
# ❌ GEOquery (2.62.2)
# ❌ HDF5Array (1.22.1)

# 套件與版本配對
package_versions <- list(
  abind       = "1.4-8",
  askpass     = "1.1",
  BH          = "1.87.0-1",
  BiocIO      = "1.4.0",
  biomaRt     = "2.50.3",
  bitops      = "1.0-7",
  caTools     = "1.18.2",
  cli         = "3.6.1",
  colorspace  = "2.1-0",
  commonmark  = "1.9.0",
  cpp11       = "0.5.2",
  curl        = "5.0.0",
  data.table  = "1.14.8",
  DESeq2      = "1.34.0",
  digest      = "0.6.31",
  dotCall64   = "1.0-2",
  DoubletFinder = "2.0.4",
  fastDummies = "1.7.5",
  fields      = "14.1",
  fitdistrplus = "1.2-2",
  FNN         = "1.1.3.2",
  fs          = "1.5.2",
  future.apply = "1.11.3",
  genefilter  = "1.76.0",
  geneplotter = "1.72.0",
  GEOquery    = "2.62.2",
  ggrepel     = "0.9.3",
  glue        = "1.6.2",
  gtable      = "0.3.6",
  HDF5Array   = "1.22.1",
  igraph      = "1.4.2",
  jsonlite    = "1.8.4",
  later       = "1.3.0",
  leidenbase  = "0.1.18",
  locfit      = "1.5-9.12",
  LoomExperiment = "1.12.0",
  maps        = "3.4.1",
  matrixStats = "0.63.0",
  metap       = "1.11",
  multcomp    = "1.4-28",
  openssl     = "2.0.6",
  openxlsx    = "4.2.8",
  parallelly  = "1.35.0",
  polyclip    = "1.10-4",
  progressr   = "0.15.1",
  promises    = "1.2.0.1",
  purrr       = "1.0.1",
  R.methodsS3 = "1.8.2",
  R.oo        = "1.27.0",
  R.utils     = "2.13.0",
  R6          = "2.6.1",
  RANN        = "2.6.1",
  Rcpp        = "1.0.13-1",
  RcppArmadillo = "0.12.2.0.0",
  RcppEigen   = "0.3.3.9.3",
  remotes     = "2.5.0",
  reticulate  = "1.28",
  rhdf5       = "2.38.1",
  rhdf5filters = "1.6.0",
  Rhdf5lib    = "1.16.0",
  rPanglaoDB  = "0.2.1",
  sandwich    = "3.1-1",
  spam        = "2.9-1",
  spatstat.explore = "3.1-0",
  spatstat.geom    = "3.1-0",
  spatstat.random  = "3.1-4",
  spatstat.utils   = "3.1-2",
  stringi     = "1.8.4",
  sys         = "3.4.1",
  TCGAbiolinks = "2.22.4",
  TCGAbiolinksGUI.data = "1.14.1",
  uwot        = "0.1.14",
  withr       = "3.0.2",
  xfun        = "0.39",
  yaml        = "2.3.7",
  zoo         = "1.8-12"
)


# 安裝函數
install_specific_version <- function(pkg, ver) {
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE) || as.character(packageVersion(pkg)) != ver) {
      devtools::install_version(pkg, version = ver, upgrade = "never", dependencies = TRUE)
    }
  }, error = function(e) message(sprintf("❌ %s (%s)", pkg, ver)))
}

# 執行安裝
for (pkg in names(package_versions)) {
  install_specific_version(pkg, package_versions[[pkg]])
}

#### devtools 套件安裝 ####
if(!require('devtools')) {install.packages('devtools'); library(devtools)}

devtools::install_version("Matrix", version = "1.6.4")
devtools::install_version("TFisher", version = "0.2.0")
devtools::install_version("parallelly", version = "1.43.0")
devtools::install_version("promises", version = "1.3.2")
devtools::install_version("Seurat", version = "4.3.0")
devtools::install_version("tidyverse", version = "2.0.0")
devtools::install_version("ggpubr", version = "0.6.0")

#### Bioconductor 套件安裝 ####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("annotate")
BiocManager::install("annotate", version = "3.17",  dependencies = TRUE, update = FALSE)
BiocManager::install("fgsea")
BiocManager::install("CelliD", version = "3.18",  dependencies = TRUE, update = FALSE)


#### 顯示目前已安裝的版本 ####
as.data.frame(installed.packages()[, c("Package", "Version")])