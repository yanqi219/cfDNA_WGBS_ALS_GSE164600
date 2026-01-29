# 04_DNAm_analysis.R
#
# cfDNA WGBS DNA methylation analysis (ALS vs Control)
#
# Goals
# - Extract per-CpG methylation calls from Bismark-aligned BAMs (methylKit).
# - Aggregate into 500 bp tiling windows and test for differential methylation (ALS vs Ctrl).
# - Store reusable objects/tables + generate publication-quality figures.
#
# Key references / widely used methods
# - Bismark aligner + tags in BAM: Krueger & Andrews, Bioinformatics (2011).
# - methylKit (incl. processBismarkAln, tiling, logistic regression): Akalin et al., Genome Biol (2012).
# - DSS beta-binomial framework for BS-seq DM (context for standard practice): Feng et al., NAR (2014).
# - BSmooth/bsseq for WGBS smoothing/DMR (context for standard practice): Hansen et al., Genome Biol (2012).
# - cfDNA methylome feature engineering often uses fixed windows + MAD selection for PCA (cfMeDIP-seq): Shen et al., Nature (2018).
#
# Notes
# - This repo stores BAMs under: cfDNA_WGBS_ALS_GSE164600/data/raw/ (not tracked by git).
# - sample IDs + BAM filenames are in: cfDNA_WGBS_ALS_GSE164600/data/sample_metadata.csv
# - This dataset is downsampled to chr21; analysis will therefore focus on chr21.

# Function to install and load packages if missing
install_and_load <- function(pkg, bioc=FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cran_pkgs <- c(
  "ggplot2", "tidyverse", "purrr", "stringr", "tibble",
  "caret", "pROC", "patchwork", "here", "parallel", "ggpubr", "circlize"
)
bioc_pkgs <- c(
  "Rsamtools", "cigarillo", "GenomicRanges", "IRanges", "ChIPseeker", "biomaRt", "org.Hs.eg.db", "AnnotationDbi",
  "Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "GenomeInfoDb", "TxDb.Hsapiens.UCSC.hg38.knownGene", "GenomicFeatures",
  "ComplexHeatmap", "methylKit"
)
invisible(lapply(cran_pkgs, install_and_load, bioc = FALSE))
invisible(lapply(bioc_pkgs, install_and_load, bioc = TRUE))

# Config ----
set.seed(1106)

WIN_SIZE_BP <- 500L              # requested window size (non-overlapping tiles)
STEP_SIZE_BP <- 500L

MINCOV_PER_C <- 5L               # per-cytosine minimum coverage for methylation calling
MINQUAL <- 20L                   # phred threshold (methylKit default in docs)
HI_COV_PERC <- 99.9              # remove extreme coverage (PCR bias)

MIN_PER_GROUP <- 3L              # unite(): loci present in >= MIN_PER_GROUP samples per group
COV_BASES_PER_TILE <- 5L         # tileMethylCounts(): min cytosines per window

FDR_Q <- 0.05
DELTA_BETA <- 0.10               # 10% methylation difference threshold (common in BS-seq DMR practice)

NCORES <- max(1L, min(8L, parallel::detectCores(logical = FALSE)))

# Paths / inputs ----
PROJECT_DIR <- here::here()
DATA_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/data")
RAW_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
METH_DIR <- file.path(PROCESSED_DIR, "methylation")
CALLS_DIR <- file.path(METH_DIR, "methylKit_calls")
RDS_DIR <- file.path(METH_DIR, "rds")

RESULTS_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/results")
FIG_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")

if (!dir.exists(METH_DIR)) dir.create(METH_DIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(CALLS_DIR)) dir.create(CALLS_DIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(RDS_DIR)) dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(Group = factor(.data$Group, levels = c("Ctrl", "ALS")))

bam_paths <- file.path(RAW_DIR, metadata$Bam_file)
names(bam_paths) <- metadata$sample_id

missing_bams <- bam_paths[!file.exists(bam_paths)]
if (length(missing_bams) > 0) {
  stop(
    "Missing BAM files under `cfDNA_WGBS_ALS_GSE164600/data/raw/`.\n",
    "Expected (from sample_metadata.csv):\n",
    paste0("  - ", basename(missing_bams), collapse = "\n"),
    "\n\nAdd the BAMs and re-run this script.",
    call. = FALSE
  )
}

# Ensure BAM indices exist (needed by many downstream tools)
message("Indexing BAMs (if needed)...")
invisible(lapply(bam_paths, function(p) {
  bai1 <- paste0(p, ".bai")
  bai2 <- sub("\\.bam$", ".bai", p, ignore.case = TRUE)
  if (!file.exists(bai1) && !file.exists(bai2)) {
    Rsamtools::indexBam(p)
  }
  invisible(TRUE)
}))

# ------------------------------------------------------------------------------
# Plot style (match previous scripts)
# ------------------------------------------------------------------------------
theme_pub <- function(base_size = 14, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.25), margin = margin(b = 8)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1.07), margin = margin(b = 6)),
      axis.title.x = element_text(face = "bold", size = rel(1.08), margin = margin(t = 9)),
      axis.title.y = element_text(face = "bold", size = rel(1.08), margin = margin(r = 9)),
      axis.text = element_text(color = "black", size = rel(1.00)),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.line = element_line(color = "black", linewidth = 0.7),
      legend.title = element_text(face = "bold", size = rel(1.08)),
      legend.text = element_text(size = rel(1.00)),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      legend.justification = "center",
      legend.box.background = element_rect(color = NA, fill = NA),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.30, linetype = "dotted"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = rel(1.05), margin = margin(t = 3, b = 3)),
      plot.margin = margin(16, 16, 16, 16)
    )
}

COLORS <- c("ALS" = "#E64B35", "Ctrl" = "#4DBBD5")

# ------------------------------------------------------------------------------
# 1) Extract CpG methylation calls from Bismark BAMs (methylKit)
# ------------------------------------------------------------------------------
raw_rds <- file.path(RDS_DIR, "methylRawList_CpG.rds")
if (file.exists(raw_rds)) {
  message("Loading cached CpG methylation calls: ", raw_rds)
  meth_raw <- readRDS(raw_rds)
} else {
  message("Extracting CpG methylation calls from Bismark BAMs (this can take a while)...")

  # treatment vector: methylKit convention (0/1)
  treatment <- ifelse(metadata$Group == "ALS", 1L, 0L)

  # processBismarkAln supports multiple BAMs when 'location' and 'sample.id' are lists
  meth_raw <- methylKit::processBismarkAln(
    location = as.list(unname(bam_paths)),
    sample.id = as.list(metadata$sample_id),
    assembly = "hg38",
    save.folder = CALLS_DIR,
    save.context = "CpG",
    read.context = "CpG",
    nolap = TRUE,
    mincov = MINCOV_PER_C,
    minqual = MINQUAL,
    treatment = treatment,
    verbose = TRUE
  )

  # Coverage filtering: drop extreme high-coverage sites (PCR bias)
  meth_raw <- methylKit::filterByCoverage(
    meth_raw,
    lo.count = MINCOV_PER_C,
    hi.perc = HI_COV_PERC
  )

  saveRDS(meth_raw, raw_rds)
}

# Optional: restrict to chr21 only (works for both in-memory and DB-backed objects)
message("Restricting to: ", paste(CHROMS_TO_USE, collapse = ", "))
meth_raw <- methylKit::select(meth_raw, chr = CHROMS_TO_USE)

# ------------------------------------------------------------------------------
# 2) Unite CpGs across samples (destrand CpG) and tile into 500 bp windows
# ------------------------------------------------------------------------------
mbase_rds <- file.path(RDS_DIR, "methylBase_CpG_united.rds")
if (file.exists(mbase_rds)) {
  message("Loading cached united CpG table: ", mbase_rds)
  mbase <- readRDS(mbase_rds)
} else {
  message("Uniting CpGs across samples (destrand=TRUE, min.per.group=", MIN_PER_GROUP, ")...")
  mbase <- methylKit::unite(
    meth_raw,
    destrand = TRUE,
    min.per.group = MIN_PER_GROUP,
    mc.cores = NCORES
  )
  saveRDS(mbase, mbase_rds)
}

tiles_rds <- file.path(RDS_DIR, sprintf("methylBase_tiles_%dbp.rds", WIN_SIZE_BP))
if (file.exists(tiles_rds)) {
  message("Loading cached tiles: ", tiles_rds)
  tiles <- readRDS(tiles_rds)
} else {
  message("Tiling methylation into ", WIN_SIZE_BP, " bp windows...")
  tiles <- methylKit::tileMethylCounts(
    mbase,
    win.size = WIN_SIZE_BP,
    step.size = STEP_SIZE_BP,
    cov.bases = COV_BASES_PER_TILE,
    mc.cores = NCORES
  )
  saveRDS(tiles, tiles_rds)
}

# ------------------------------------------------------------------------------
# 3) Differential methylation (ALS vs Ctrl) on 500 bp windows
# ------------------------------------------------------------------------------
diff_rds <- file.path(RDS_DIR, sprintf("methylDiff_tiles_%dbp_ALS_vs_Ctrl.rds", WIN_SIZE_BP))
if (file.exists(diff_rds)) {
  message("Loading cached differential methylation results: ", diff_rds)
  diff_tiles <- readRDS(diff_rds)
} else {
  message("Running differential methylation testing (logistic regression; overdispersion='MN'; adjust='BH')...")
  diff_tiles <- methylKit::calculateDiffMeth(
    tiles,
    overdispersion = "MN",
    test = "F",
    adjust = "BH",
    mc.cores = NCORES
  )
  saveRDS(diff_tiles, diff_rds)
}

# Extract significant tiles as "DMRs" (fixed-width windows)
diff_tbl <- as.data.frame(diff_tiles) %>%
  as_tibble() %>%
  rename(
    chr = chr,
    start = start,
    end = end,
    strand = strand,
    pvalue = pvalue,
    qvalue = qvalue,
    meth_diff = meth.diff
  ) %>%
  mutate(
    mid = as.integer(floor((start + end) / 2)),
    direction = case_when(
      meth_diff >= 0 ~ "Hyper (ALS>Ctrl)",
      TRUE ~ "Hypo (ALS<Ctrl)"
    ),
    is_sig = is.finite(qvalue) & qvalue <= FDR_Q & abs(meth_diff) >= (DELTA_BETA * 100)
  )

readr::write_csv(diff_tbl, file.path(TABLE_DIR, sprintf("dmr_tiles_%dbp_ALS_vs_Ctrl_all.csv", WIN_SIZE_BP)))

dmr_tbl <- diff_tbl %>%
  filter(.data$is_sig) %>%
  arrange(.data$qvalue)

readr::write_csv(dmr_tbl, file.path(TABLE_DIR, sprintf("dmr_tiles_%dbp_ALS_vs_Ctrl_sig_q%.2f_delta%.0f.csv", WIN_SIZE_BP, FDR_Q, DELTA_BETA * 100)))

dmr_summary <- dmr_tbl %>%
  count(direction, name = "n") %>%
  mutate(total = sum(n))
readr::write_csv(dmr_summary, file.path(TABLE_DIR, sprintf("dmr_tiles_%dbp_ALS_vs_Ctrl_summary.csv", WIN_SIZE_BP)))

# ------------------------------------------------------------------------------
# 4) Store methylation matrices (for reuse)
# ------------------------------------------------------------------------------
tile_pct <- methylKit::percMethylation(tiles)
colnames(tile_pct) <- metadata$sample_id[match(colnames(tile_pct), metadata$sample_id)]
saveRDS(tile_pct, file.path(RDS_DIR, sprintf("tiles_%dbp_percent_methylation_matrix.rds", WIN_SIZE_BP)))

# Also store as a lightweight long table (useful for non-R consumers)
tile_pct_long <- as_tibble(tile_pct) %>%
  mutate(tile_id = row_number()) %>%
  pivot_longer(-tile_id, names_to = "sample_id", values_to = "pct_meth") %>%
  left_join(diff_tbl %>% transmute(tile_id = row_number(), chr, start, end, mid), by = "tile_id") %>%
  left_join(metadata %>% select(sample_id, Group), by = "sample_id")
readr::write_csv(tile_pct_long, file.path(TABLE_DIR, sprintf("tiles_%dbp_percent_methylation_long.csv", WIN_SIZE_BP)))

# ------------------------------------------------------------------------------
# 5) Publication-quality figures
# ------------------------------------------------------------------------------
message("Generating figures...")

# 5A) Tile methylation distributions (density)
tile_pct_df <- tile_pct_long %>%
  filter(is.finite(.data$pct_meth)) %>%
  mutate(Group = factor(Group, levels = c("Ctrl", "ALS")))

p_density <- ggplot(tile_pct_df, aes(x = pct_meth, color = Group)) +
  geom_density(linewidth = 0.9, alpha = 0.9) +
  scale_color_manual(values = COLORS) +
  labs(
    title = sprintf("cfDNA WGBS methylation distribution (%d bp tiles)", WIN_SIZE_BP),
    subtitle = "Percent methylation across tiled windows (chr21)",
    x = "Percent methylation (%)",
    y = "Density",
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

save_pdf_png(p_density, sprintf("fig6A_tile_methylation_density_%dbp", WIN_SIZE_BP), w = 7.5, h = 5.2)

# 5B) PCA on most variable tiles (MAD selection, cfDNA methylome papers commonly do this)
mat <- tile_pct
keep_rows <- rowSums(is.na(mat)) == 0
mat <- mat[keep_rows, , drop = FALSE]

mad_v <- apply(mat, 1, stats::mad, na.rm = TRUE)
top_n <- min(10000L, nrow(mat))
top_idx <- order(mad_v, decreasing = TRUE)[seq_len(top_n)]
mat_top <- mat[top_idx, , drop = FALSE]

pca <- prcomp(t(scale(t(mat_top))), center = TRUE, scale. = TRUE)
pca_df <- tibble(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>% left_join(metadata, by = "sample_id")

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.0, alpha = 0.9) +
  geom_text(aes(label = sample_id), size = 2.6, vjust = -0.9, check_overlap = TRUE) +
  scale_color_manual(values = COLORS) +
  labs(
    title = sprintf("PCA of %d most variable %d bp tiles (chr21)", top_n, WIN_SIZE_BP),
    subtitle = sprintf("PC1 %.1f%%, PC2 %.1f%% variance explained", 100 * var_expl[1], 100 * var_expl[2]),
    x = "PC1",
    y = "PC2",
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

save_pdf_png(p_pca, sprintf("fig6B_pca_tiles_%dbp", WIN_SIZE_BP), w = 7.5, h = 5.8)

# 5C) Heatmap of most variable tiles (row Z-scored)
z_mat <- t(scale(t(mat_top)))
z_mat[!is.finite(z_mat)] <- 0

sample_order <- metadata %>% arrange(Group, sample_id) %>% pull(sample_id)
sample_order <- intersect(sample_order, colnames(z_mat))
z_mat <- z_mat[, sample_order, drop = FALSE]

anno_df <- metadata %>%
  distinct(sample_id, Group) %>%
  filter(sample_id %in% colnames(z_mat)) %>%
  arrange(match(sample_id, colnames(z_mat))) %>%
  column_to_rownames("sample_id")

ha <- HeatmapAnnotation(
  df = anno_df,
  col = list(Group = COLORS),
  annotation_name_gp = grid::gpar(fontface = "bold", fontsize = 10),
  annotation_legend_param = list(title_gp = grid::gpar(fontface = "bold"))
)

ht <- Heatmap(
  z_mat,
  name = "Z",
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "Samples",
  row_title = sprintf("Top %d variable %d bp tiles (chr21)", top_n, WIN_SIZE_BP),
  col = circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
  heatmap_legend_param = list(title_gp = grid::gpar(fontface = "bold"), labels_gp = grid::gpar(fontsize = 9))
)

grDevices::cairo_pdf(file.path(FIG_DIR, sprintf("fig6C_heatmap_tiles_%dbp_top%d.pdf", WIN_SIZE_BP, top_n)), width = 8.5, height = 8.0)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

png(file.path(FIG_DIR, sprintf("fig6C_heatmap_tiles_%dbp_top%d.png", WIN_SIZE_BP, top_n)),
    width = 8.5, height = 8.0, units = "in", res = 450, type = "cairo-png")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

# 5D) Volcano plot (delta methylation vs -log10(q))
p_volcano <- ggplot(diff_tbl, aes(x = meth_diff, y = -log10(qvalue))) +
  geom_point(aes(color = is_sig), alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = c(-DELTA_BETA * 100, DELTA_BETA * 100), linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_Q), linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray70")) +
  labs(
    title = sprintf("Differential methylation (ALS vs Ctrl) — %d bp tiles", WIN_SIZE_BP),
    subtitle = sprintf("Significant: q≤%.2f and |Δmeth|≥%.0f%%", FDR_Q, DELTA_BETA * 100),
    x = "Δ methylation (%) (ALS − Ctrl)",
    y = expression(-log10(qvalue)),
    color = "Significant"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

save_pdf_png(p_volcano, sprintf("fig6D_volcano_tiles_%dbp", WIN_SIZE_BP), w = 7.2, h = 5.6)

# 5E) Manhattan-style plot along chr21
diff_chr <- diff_tbl %>%
  filter(chr %in% CHROMS_TO_USE) %>%
  mutate(
    sig_dir = case_when(
      is_sig & meth_diff >= 0 ~ "Hyper",
      is_sig & meth_diff < 0 ~ "Hypo",
      TRUE ~ "Not sig"
    )
  )

p_manhattan <- ggplot(diff_chr, aes(x = mid, y = -log10(qvalue))) +
  geom_point(aes(color = sig_dir), alpha = 0.75, size = 1.1) +
  geom_hline(yintercept = -log10(FDR_Q), linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_color_manual(values = c("Hyper" = "#E64B35", "Hypo" = "#4DBBD5", "Not sig" = "gray75")) +
  labs(
    title = sprintf("Differential methylation along %s (%d bp tiles)", paste(CHROMS_TO_USE, collapse = ","), WIN_SIZE_BP),
    subtitle = sprintf("q≤%.2f dashed", FDR_Q),
    x = sprintf("Genomic position on %s (bp)", paste(CHROMS_TO_USE, collapse = ",")),
    y = expression(-log10(qvalue)),
    color = NULL
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

save_pdf_png(p_manhattan, sprintf("fig6E_manhattan_tiles_%dbp_chr21", WIN_SIZE_BP), w = 10.0, h = 5.2)

# ------------------------------------------------------------------------------
# Final report
# ------------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("DNA methylation analysis complete")
message(strrep("=", 60))
message("Inputs:")
message("  - BAMs: ", RAW_DIR)
message("Outputs written to:")
message("  - R objects: ", RDS_DIR)
message("  - Tables:    ", TABLE_DIR)
message("  - Figures:   ", FIG_DIR)
message("Summary (significant tiles):")
print(dmr_summary)
message(strrep("=", 60))

