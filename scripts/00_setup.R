# 00_setup.R
# Install and load packages for cfDNA WGBS analysis
# Run this script once to set up the environment

# --- Initialize renv ---
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv if not already done
if (!file.exists("renv.lock")) {
  renv::init()
}

# --- Install packages ---

# CRAN packages
cran_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "data.table",   # Fast data processing (important for scalability)
  "caret",        # Classification and metrics
  "pROC",         # ROC curves
  "parallel"      # Parallel processing
)

# Bioconductor packages
bioc_packages <- c(
  "Rsamtools",    # BAM file handling
  "GenomicAlignments",  # Read alignments
  "Biostrings",   # Sequence manipulation
  "BSgenome.Hsapiens.UCSC.hg38"  # Reference genome for motif extraction
)

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# --- Snapshot renv ---
renv::snapshot()

# --- Load packages ---
library(tidyverse)
library(data.table)
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)

message("Setup complete. All packages installed and loaded.")
