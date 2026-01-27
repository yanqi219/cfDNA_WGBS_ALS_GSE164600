# 00_setup.R
# Reproducible environment setup (publication-grade)
#
# Goal:
# - Anyone cloning the repo can reproduce results using `renv.lock`.
# - No hidden system state; dependencies are pinned and restored deterministically.
#
# Usage:
# - From cfDNA_WGBS_ALS_GSE164600/: source("scripts/00_setup.R")
# - From Take_home_task/:          source("cfDNA_WGBS_ALS_GSE164600/scripts/00_setup.R")
#
# IMPORTANT:
# - Commit `renv.lock` + `renv/` to GitHub after bootstrapping once.
#
# Notes:
# - macOS: prefer binaries to avoid toolchain issues.
# - `parallel` is a base/recommended package that ships with R.

if (Sys.info()[["sysname"]] == "Darwin") {
  options(pkgType = "binary")
}
options(repos = c(CRAN = "https://cloud.r-project.org"))

# This project is pinned to Bioconductor 3.22 (R 4.5.x).
# If you update R/Bioconductor, regenerate and recommit renv.lock.
BIOC_VERSION <- "3.22"

get_script_path <- function() {
  # 1) Rscript --file=... (non-interactive)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    p <- sub("^--file=", "", file_arg[[1]])
    p <- gsub("~+~", " ", p, fixed = TRUE) # some environments encode spaces
    return(normalizePath(p, mustWork = FALSE))
  }

  # 2) source("...") (interactive)
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NA_character_)
  if (!is.na(ofile) && nzchar(ofile)) {
    return(normalizePath(ofile, mustWork = FALSE))
  }

  NA_character_
}

get_project_root <- function() {
  script_path <- get_script_path()
  if (!is.na(script_path) && nzchar(script_path) && file.exists(script_path)) {
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE))
  }

  # Fallback: infer from common working directories
  wd <- normalizePath(getwd(), mustWork = FALSE)
  candidates <- c(wd, file.path(wd, "cfDNA_WGBS_ALS_GSE164600"))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "scripts", "00_setup.R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  wd
}

project_root <- get_project_root()
if (!dir.exists(project_root)) stop("Project root does not exist: ", project_root)

old_wd <- getwd()
setwd(project_root)
on.exit(setwd(old_wd), add = TRUE)

# --- renv bootstrap / restore ---
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (exists("consent", asNamespace("renv"), inherits = FALSE)) renv::consent(provided = TRUE)

# Ensure renv infrastructure exists (safe to re-run)
if (!dir.exists("renv") || !file.exists(file.path("renv", "activate.R"))) {
  renv::init(bare = TRUE)
}

lockfile <- "renv.lock"

if (file.exists(lockfile)) {
  # Publication mode: deterministic restore
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  options(repos = BiocManager::repositories(version = BIOC_VERSION))
  renv::settings$bioconductor.version(BIOC_VERSION)

  renv::restore(prompt = FALSE)
  renv::load(project_root)
} else {
  # Bootstrap mode: one-time environment creation
  message("No renv.lock found. Bootstrapping environment and creating renv.lock...")

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

  # Bioconductor repos
  options(repos = BiocManager::repositories(version = BIOC_VERSION))
  # Pin the Bioconductor version in renv settings (recorded in renv.lock)
  renv::settings$bioconductor.version(BIOC_VERSION)
  renv::load(project_root)

  # Packages required by scripts/01_qc_processing.R
  cran_packages <- c(
    # Core tidyverse components (avoid the tidyverse meta-package for stability)
    "ggplot2",
    "dplyr",
    "tidyr",
    "readr",
    "purrr",
    "stringr",
    "tibble",
    "caret",
    "pROC",
    "patchwork",
    "here"
  )

  bioc_packages <- c(
    "Rsamtools",
    "GenomicAlignments",
    "GenomicRanges",
    "IRanges",
    "GenomeInfoDb",
    "BiocGenerics",
    "Biostrings",
    "BSgenome.Hsapiens.UCSC.hg38"
  )

  snapshot_packages <- c(cran_packages, bioc_packages)

  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }

  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = TRUE, type = "source")
    }
  }

  # Snapshot captures exact versions into renv.lock.
  # Use an explicit package list for deterministic snapshots.
  options(repos = BiocManager::repositories(version = BIOC_VERSION))
  renv::settings$bioconductor.version(BIOC_VERSION)
  renv::snapshot(packages = snapshot_packages, prompt = FALSE, force = TRUE)
}

# --- Load packages used by analysis scripts ---
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(parallel))

message("Setup complete. renv restored and packages loaded (project root: ", project_root, ").")
