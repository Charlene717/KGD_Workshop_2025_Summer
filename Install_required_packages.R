#!/usr/bin/env Rscript
# Install_required_packages.R
# ------------------------------------------------------------
# Automatic installer for all R packages needed in the 2025/07
# scRNA‑seq & Spatial Transcriptomics workshop.
# Run this script once before the class:
#   > source("Install_required_packages.R")
# ------------------------------------------------------------

# ------------ Helper --------------------------------------------------------
install_and_load <- function(pkg, install_fun) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    install_fun(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ------------ CRAN packages -------------------------------------------------
# if(!require('Seurat'))         { install.packages('Seurat');         library(Seurat) }
# if(!require('tidyverse'))      { install.packages('tidyverse');      library(tidyverse) } 


cran_pkgs <- c(
  "Seurat", "SeuratObject", "SeuratDisk",
  "tidyverse", "dplyr", "ggplot2", "patchwork",
  "viridis", "future", "future.apply", "fastDummies"
)

cran_install <- function(pkg)
  install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)


# ------------ Bioconductor packages -----------------------------------------
# if(!require('SingleR'))        { BiocManager::install('SingleR');    library(SingleR) }
# if(!require('celldex'))        { BiocManager::install('celldex');    library(celldex) }
# if(!require('monocle'))        { BiocManager::install('monocle');    library(monocle) }

bioc_pkgs <- c(
  "org.Hs.eg.db", "org.Mm.eg.db",
  "clusterProfiler", "enrichplot", "ReactomePA",
  "fgsea", "DOSE", "GSEABase"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)

bioc_install <- function(pkg)
  BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)

# ------------ GitHub packages -----------------------------------------------
# if(!require('DoubletFinder'))  { devtools::install_github('chris-mcginnis-ucsf/DoubletFinder');  library(DoubletFinder) }
# if(!require('DoubletFinder'))  { remotes::install_github('chris-mcginnis-ucsf/DoubletFinder');  library(DoubletFinder) }


github_repos <- c(
  "chris-mcginnis-ucsf/DoubletFinder",
  "sqjin/CellChat"
)

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)

github_install <- function(repo) {
  pkg <- sub(".*/", "", repo)
  remotes::install_github(repo, upgrade = "never", quiet = TRUE)
  pkg
}

# ------------ Installation pipeline -----------------------------------------
message("🛠️  Checking and installing CRAN packages ...")
for (pkg in cran_pkgs) {
  install_and_load(pkg, cran_install)
}

message("🛠️  Checking and installing Bioconductor packages ...")
for (pkg in bioc_pkgs) {
  install_and_load(pkg, bioc_install)
}

message("🛠️  Checking and installing GitHub packages ...")
for (repo in github_repos) {
  pkg <- sub(".*/", "", repo)
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    pkg_name <- github_install(repo)
    library(pkg_name, character.only = TRUE)
  }
}

# ------------ CellChat Database (Human) -------------------------------------
message("🛠️  Setting up CellChat human database ...")
if (suppressWarnings(require("CellChat", character.only = TRUE))) {
  if (!"CellChatDB.human" %in% ls("package:CellChat")) {
    message("⚠️  Current CellChat version does not include 'CellChatDB.human'; consider upgrading the package.")
  } else {
    CellChat::CellChatDB <- CellChat::CellChatDB.human
    message("✅ CellChatDB.human 已載入 (CellChat::CellChatDB)")
  }
}

message("✅ All required packages are now installed, loaded, and configured.")


## Displays the currently installed version
as.data.frame(installed.packages()[, c("Package", "Version")])