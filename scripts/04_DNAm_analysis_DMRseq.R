# 04_DNAm_analysis_DMRseq.R
#
# cfDNA WGBS DNA methylation analysis (ALS vs Control) using dmrseq
#
# Goals
# - Load Bismark coverage files into bsseq objects.
# - Detect differentially methylated regions (DMRs) using dmrseq.
# - Store reusable objects/tables + generate publication-quality figures.
#
# Key references
# - dmrseq: Korthauer et al., Biostatistics (2018). doi:10.1093/biostatistics/kxy007
# - BSmooth/bsseq: Hansen et al., Genome Biology (2012). doi:10.1186/gb-2012-13-10-r83
# - Bismark: Krueger & Andrews, Bioinformatics (2011).
#
# Notes
# - dmrseq is specifically designed for low-coverage WGBS data.
# - BSmooth smoothing borrows strength from neighboring CpGs (works with ~4x coverage).
# - Region-level inference with proper FDR control via permutation.

# Setup ----
install_and_load <- function(pkg, bioc = FALSE) {

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

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "tibble", "stringr",
               "patchwork", "here", "circlize")
bioc_pkgs <- c("bsseq", "dmrseq", "GenomicRanges", "ComplexHeatmap",
               "BiocParallel", "annotatr", "BSgenome.Hsapiens.UCSC.hg38")

invisible(lapply(cran_pkgs, install_and_load, bioc = FALSE))
invisible(lapply(bioc_pkgs, install_and_load, bioc = TRUE))

# Configuration ----
theme_pub <- function(base_size = 14, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.25),
                                margin = margin(b = 8)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1.07),
                                   margin = margin(b = 6)),
      axis.title.x = element_text(face = "bold", size = rel(1.08),
                                  margin = margin(t = 9)),
      axis.title.y = element_text(face = "bold", size = rel(1.08),
                                  margin = margin(r = 9)),
      axis.text = element_text(color = "black", size = rel(1.00)),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.line = element_line(color = "black", linewidth = 0.7),
      legend.title = element_text(face = "bold", size = rel(1.08)),
      legend.text = element_text(size = rel(1.00)),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.30,
                                      linetype = "dotted"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = rel(1.05)),
      plot.margin = margin(16, 16, 16, 16)
    )
}

COLORS <- c("ALS" = "#E64B35", "Ctrl" = "#4DBBD5")

genome <- BSgenome.Hsapiens.UCSC.hg38

set.seed(1106)

# dmrseq parameters
CUTOFF        <- 0.05    # Effect size cutoff for candidate region detection
MIN_NUM_REGION <- 3      # Minimum CpGs per candidate region
SIG_THRESHOLD <- 0.05    # Significance threshold for DMRs

# Parallelization
NCORES <- max(1L, min(4L, parallel::detectCores(logical = FALSE)))

# Paths ----
PROJECT_DIR <- here::here()
DATA_DIR    <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/data")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
BISMARK_COV_DIR <- file.path(PROCESSED_DIR, "methylation_extractor")

RESULTS_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/results")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
TABLE_DIR   <- file.path(RESULTS_DIR, "tables")
RDS_DIR     <- file.path(DATA_DIR, "processed/methylation/rds")

for (d in c(FIG_DIR, TABLE_DIR, RDS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Load metadata ----
metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"),
                            show_col_types = FALSE) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Ctrl", "ALS")))

cov_tbl <- tibble::tibble(
  cov_file = list.files(
    path = BISMARK_COV_DIR,
    pattern = "\\.bismark\\.cov\\.gz$",
    full.names = TRUE,
    recursive = TRUE
  )) %>%
  dplyr::mutate(sample_id = sub("\\..*", "", basename(.data$cov_file)))
metadata <- metadata %>%
  left_join(cov_tbl, by = "sample_id") %>%
  filter(!is.na(cov_file))

message(sprintf("Found %d samples with Bismark coverage files.", nrow(metadata)))

# Process data into BSseq object ----
message("Reading Bismark coverage files into BSseq object...")

# Read Bismark coverage files
# Format: chr, start, end, meth%, count_meth, count_unmeth
bs <- bsseq::read.bismark(
  files = metadata$cov_file,
  rmZeroCov = TRUE,
  strandCollapse = TRUE,
  verbose = TRUE
)

# Add sample metadata
pData(bs)$Group <- metadata$Group
pData(bs)$sample_id <- metadata$sample_id

message(sprintf("BSseq object: %d CpG loci, %d samples",
                nrow(bs), ncol(bs)))

# Save BSseq object
saveRDS(bs, file.path(RDS_DIR, "bsseq_object.rds"))

## Filter low-coverage CpGs ----
# Keep CpGs with coverage in at least 2 samples per group
cov_mat <- bsseq::getCoverage(bs, type = "Cov")

# Count samples with coverage > 0 per group
n_ctrl <- sum(metadata$Group == "Ctrl")
n_als  <- sum(metadata$Group == "ALS")

ctrl_idx <- which(metadata$Group == "Ctrl")
als_idx  <- which(metadata$Group == "ALS")

ctrl_covered <- rowSums(cov_mat[, ctrl_idx] > 0) >= 1
als_covered  <- rowSums(cov_mat[, als_idx] > 0) >= 1

keep_loci <- ctrl_covered & als_covered
bs_filt <- bs[keep_loci, ]

message(sprintf("After filtering: %d CpG loci (%.1f%% retained)",
                nrow(bs_filt), 100 * nrow(bs_filt) / nrow(bs)))

# DMR analysis ----
message("Running dmrseq DMR detection (this may take a few minutes)...")

# dmrseq performs BSmooth smoothing internally
if(file.exists(file.path(RDS_DIR, "dmrseq_results_raw.rds"))) {
  dmrs <- readRDS(file.path(RDS_DIR, "dmrseq_results_raw.rds"))
} else {
  dmrs <- dmrseq::dmrseq(
    bs = bs_filt,
    testCovariate = "Group",
    cutoff = CUTOFF,
    minNumRegion = MIN_NUM_REGION,
    bpSpan = 1000,           # Smoothing span in bp
    minInSpan = 15,          # Min CpGs in smoothing window
    maxGapSmooth = 2500,     # Max gap for smoothing
    maxGap = 1000,           # Max gap within a DMR
    verbose = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = 4)
  )
  # Save raw dmrseq output
  saveRDS(dmrs, file.path(RDS_DIR, "dmrseq_results_raw.rds"))
}

message(sprintf("dmrseq found %d candidate regions.", length(dmrs)))

## Format DMR results ----
# Convert to data frame
dmr_df <- as.data.frame(dmrs) %>%
  tibble::as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  dplyr::mutate(
    width = end - start + 1,
    mid = as.integer(floor((start + end) / 2)),
    direction = dplyr::case_when(
      beta > 0 ~ "Hyper (ALS>Ctrl)",
      beta < 0 ~ "Hypo (ALS<Ctrl)",
      TRUE ~ "None"
    ),
    is_sig = !is.na(pval) & pval <= SIG_THRESHOLD
  )

# Write all results
readr::write_csv(dmr_df, file.path(TABLE_DIR, "dmrseq_results_all.csv"))

# Significant DMRs
dmr_sig <- dmr_df %>%
  dplyr::filter(is_sig) %>%
  dplyr::arrange(pval)

# Summary statistics
dmr_summary <- dmr_sig %>%
  dplyr::count(direction, name = "n_DMRs") %>%
  dplyr::mutate(
    total = sum(n_DMRs),
    pct = round(100 * n_DMRs / total, 1)
  )

readr::write_csv(dmr_summary, file.path(TABLE_DIR, "dmrseq_DMRs_summary.csv"))

message("\nDMR Summary:")
print(dmr_summary)

## Extract methylation values for visualization ----
message("\nExtracting methylation values for visualization...")

# For visualization, we'll compute mean methylation per sample
# across genomic windows (1kb tiles)
tile_gr <- GenomicRanges::tileGenome(
  seqlens <- GenomeInfoDb::seqlengths(genome)[paste0("chr", c(21))],
  tilewidth = 1000,
  cut.last.tile.in.chrom = TRUE
)

# Compute mean methylation per tile per sample
message("Computing tiled methylation matrix...")
meth_by_tile <- bsseq::getMeth(bs_filt, regions = tile_gr, type = "raw",
                               what = "perRegion")
colnames(meth_by_tile) <- metadata$sample_id
rownames(meth_by_tile) <- sprintf("%s_%d_%d", as.character(GenomicRanges::seqnames(tile_gr)), 
                                  GenomicRanges::start(tile_gr), GenomicRanges::end(tile_gr))

# Remove tiles with all NA
keep_tiles <- rowSums(!is.na(meth_by_tile)) > 0
tile_mat <- meth_by_tile[keep_tiles, ] * 100  # Convert to percent

saveRDS(tile_mat, file.path(RDS_DIR, "dmrseq_tile_methylation_matrix.rds"))

# Figures ----
message("\nGenerating figures...")

# Methylation distribution by group ----
meth_long <- as.data.frame(tile_mat) %>%
  tibble::rownames_to_column("tile_id") %>%
  tidyr::pivot_longer(-tile_id, names_to = "sample_id", values_to = "pct_meth") %>%
  dplyr::left_join(metadata %>% dplyr::select(sample_id, Group), by = "sample_id") %>%
  dplyr::filter(!is.na(pct_meth))

p_density <- ggplot(meth_long, aes(x = pct_meth, color = Group)) +
  geom_density(linewidth = 0.9, alpha = 0.9) +
  scale_color_manual(values = COLORS) +
  labs(
    title = "cfDNA WGBS methylation distribution (1 kb tiles)",
    subtitle = "Percent methylation across tiled windows (chr21)",
    x = "Percent methylation (%)",
    y = "Density",
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(FIG_DIR, "fig6A_dmrseq_methylation_density.png"), p_density, width = 7.5, height = 5.2, dpi = 450, type = "cairo-png")

## PCA on most variable tiles ----
mat <- tile_mat
keep_rows <- rowSums(is.na(mat)) == 0
mat <- mat[keep_rows, , drop = FALSE]

if (nrow(mat) > 100) {
  mad_v <- apply(mat, 1, stats::mad, na.rm = TRUE)
  top_n <- min(5000L, nrow(mat))
  top_idx <- order(mad_v, decreasing = TRUE)[seq_len(top_n)]
  mat_top <- mat[top_idx, , drop = FALSE]
  
  pca <- prcomp(t(mat_top), center = TRUE, scale. = TRUE)
  pca_df <- tibble::tibble(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) %>% dplyr::left_join(metadata, by = "sample_id")
  
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = sample_id), size = 2.8,
                             max.overlaps = 20) +
    scale_color_manual(values = COLORS) +
    labs(
      title = sprintf("PCA of %d most variable 1 kb tiles (chr21)", top_n),
      subtitle = sprintf("PC1 %.1f%%, PC2 %.1f%% variance explained",
                         100 * var_expl[1], 100 * var_expl[2]),
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = "Group"
    ) +
    theme_pub(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave(file.path(FIG_DIR, "fig6B_dmrseq_pca.png"), p_pca, width = 7.5, height = 6.0, dpi = 450, type = "cairo-png")
} else {
  message("Insufficient tiles for PCA visualization.")
}

## Heatmap of most variable tiles ----
if (exists("mat_top") && nrow(mat_top) > 50) {
  z_mat <- t(scale(t(mat_top)))
  z_mat[!is.finite(z_mat)] <- 0
  
  # Order samples by group
  sample_order <- metadata %>%
    dplyr::arrange(Group, sample_id) %>%
    dplyr::pull(sample_id)
  sample_order <- intersect(sample_order, colnames(z_mat))
  z_mat <- z_mat[, sample_order, drop = FALSE]
  
  # Annotation
  anno_df <- metadata %>%
    dplyr::distinct(sample_id, Group) %>%
    dplyr::filter(sample_id %in% colnames(z_mat)) %>%
    dplyr::arrange(match(sample_id, colnames(z_mat))) %>%
    tibble::column_to_rownames("sample_id")
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = anno_df,
    col = list(Group = COLORS),
    annotation_name_gp = grid::gpar(fontface = "bold", fontsize = 10)
  )
  
  ht <- ComplexHeatmap::Heatmap(
    z_mat,
    name = "Z-score",
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 8),
    column_title = "Samples",
    row_title = sprintf("Top %d variable tiles", nrow(z_mat)),
    col = circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
    heatmap_legend_param = list(title_gp = grid::gpar(fontface = "bold"))
  )
  
  png(file.path(FIG_DIR, "fig6C_dmrseq_heatmap.png"),
      width = 9, height = 8, units = "in", res = 450, type = "cairo-png")
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  invisible(dev.off())
}
 
## Volcano plot (effect size vs -log10 p-value) ----
if (nrow(dmr_df) > 0) {
  p_volcano <- ggplot(dmr_df, aes(x = beta, y = -log10(pval))) +
    geom_point(aes(color = is_sig), alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-CUTOFF, CUTOFF), linetype = "dashed",
               color = "gray45", linewidth = 0.4) +
    geom_hline(yintercept = -log10(SIG_THRESHOLD), linetype = "dashed",
               color = "gray45", linewidth = 0.4) +
    scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray70"),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not sig.")) +
    labs(
      title = "Differential methylation (ALS vs Ctrl) — dmrseq",
      subtitle = sprintf("Significant: q ≤ %.2f", SIG_THRESHOLD),
      x = "Effect size (beta coefficient)",
      y = expression(-log[10](p-value)),
      color = NULL
    ) +
    theme_pub(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave(file.path(FIG_DIR, "fig6D_dmrseq_volcano.png"), p_volcano, width = 7.2, height = 5.6, dpi = 450, type = "cairo-png")
}

## Manhattan-style plot along chr21 ----
if (nrow(dmr_df) > 0) {
  dmr_chr <- dmr_df %>%
    dplyr::mutate(
      sig_dir = dplyr::case_when(
        is_sig & beta > 0 ~ "Hyper",
        is_sig & beta < 0 ~ "Hypo",
        TRUE ~ "Not sig."
      )
    )
  
  p_manhattan <- ggplot(dmr_chr, aes(x = mid / 1e6, y = -log10(pval))) +
    geom_point(aes(color = sig_dir), alpha = 0.75, size = 1.3) +
    geom_hline(yintercept = -log10(SIG_THRESHOLD), linetype = "dashed",
               color = "gray45", linewidth = 0.4) +
    scale_color_manual(values = c("Hyper" = "#E64B35", "Hypo" = "#4DBBD5",
                                  "Not sig." = "gray75")) +
    labs(
      title = "Differential methylation along chr21 — dmrseq",
      subtitle = sprintf("q ≤ %.2f threshold shown", SIG_THRESHOLD),
      x = "Genomic position (Mb)",
      y = expression(-log[10](p-value)),
      color = NULL
    ) +
    theme_pub(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave(file.path(FIG_DIR, "fig6E_dmrseq_manhattan_chr21.png"), p_manhattan, width = 10.0, height = 5.2, dpi = 450, type = "cairo-png")
}

## Plot example DMRs (top 10 significant) ----
if (nrow(dmr_sig) > 0) {
  message("\nPlotting top DMRs...")
  
  n_plot <- min(10, nrow(dmr_sig))
  pdf(file.path(FIG_DIR, "fig6F_dmrseq_top_DMRs.pdf"))
  for (i in seq_len(n_plot)) {
    region_gr <- GenomicRanges::GRanges(dmrs[which(dmr_df$is_sig)[i]])
    tryCatch({
    p_dmr <- dmrseq::plotDMRs(
      bs_filt,
      regions = region_gr,
      testCovariate = "Group",
      extend = 1000,
      qval = FALSE,
      addRegions = region_gr)
      print(p_dmr)
    }, error = function(e) {
      message(sprintf("Could not plot DMR %d: %s", i, e$message))
    })
  }
  dev.off()
}

# Muscle-specific unmethylated regions analysis ----
# These regions are expected to be unmethylated in muscle tissue.
message("Analyzing muscle-specific unmethylated regions...")

muscle_regions <- GenomicRanges::GRanges(
  seqnames = "chr21",
  ranges = IRanges::IRanges(
    start = c(18269892, 41052886, 41604878, 44462101),
    end   = c(18270174, 41053145, 41604960, 44462229)
  ),
  region_name = c("Region1", "Region2", "Region3", "Region4")
)

# Extract methylation values for these regions
# getMeth returns NA for positions without coverage
muscle_meth <- bsseq::getMeth(bs_filt, regions = muscle_regions,
                              type = "raw", what = "perRegion")

# Check coverage in these regions
muscle_cov <- bsseq::getCoverage(bs_filt, regions = muscle_regions,
                                 type = "Cov", what = "perRegionTotal")

# Report coverage status
message("Coverage in muscle-specific regions (total reads per region):")
muscle_cov_df <- as.data.frame(muscle_cov)
colnames(muscle_cov_df) <- metadata$sample_id
muscle_cov_df$region <- c("chr21:18269892-18270174", "chr21:41052886-41053145",
                          "chr21:41604878-41604960", "chr21:44462101-44462229")
print(muscle_cov_df %>% dplyr::select(region, everything()))

# Final combined plot ----
ht_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                        annotation_legend_side = "right")
)
p_ht <- patchwork::wrap_elements(full = ht_grob)

p_dmrseq_combined <- p_density / (p_pca | p_ht) / (p_volcano | p_manhattan) +
  patchwork::plot_layout(heights = c(1, 1.2, 1.2)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(FIG_DIR, "fig6_dmrseq_combined.png"),
  p_dmrseq_combined,
  width = 17, height = 20, dpi = 300
)
