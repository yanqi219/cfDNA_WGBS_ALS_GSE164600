# 03_end_motif_analysis.R
#
# End-motif analysis (5' 4-mers; 256 motifs)
#
# What this script produces (Figure-4-like panels):
# - (A) Unsupervised hierarchical clustering heatmap of all 256 4-mer motifs (row Z-scored)
# - (B) Shannon entropy-based motif diversity score (MDS) per sample
# - (D) Sequence logos of cfDNA end motifs (probability logos; pooled within group)
#
# References / best-practice inspiration:
# - Liu et al., STAR Protocols (2024): heatmap of row Z-scored 4-mer frequencies; cfDNA end motif workflows
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC11489062/
# - Jiang et al., Cancer Discovery (2020): plasma DNA end-motif profiling as fragmentomic marker
# - Helzer et al., Nat Commun (2025): MDS defined from 256 4-mer motif frequency distribution
#
# Notes:
# - This script consumes `data/processed/all_metrics.rds` produced by `01_qc_processing.R`.
# - `all_metrics[[i]]$motifs_5p` is a character vector of 4-mer sequences (genome-derived) per sample.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(patchwork))

# Paths ----
PROJECT_DIR <- here::here()
DATA_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/data")
RESULTS_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/results")
FIG_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(TABLE_DIR)) dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Load sample metadata ----
metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(Group = factor(.data$Group, levels = c("Ctrl", "ALS")))

# Theme + palette (match 02_fragmentation_analysis.R) ----
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

# Helpers ----
all_4mers <- function() {
  bases <- c("A", "C", "G", "T")
  grid <- expand.grid(b1 = bases, b2 = bases, b3 = bases, b4 = bases, stringsAsFactors = FALSE)
  apply(grid, 1, paste0, collapse = "")
}

is_valid_4mer <- function(x) {
  nchar(x) == 4 & grepl("^[ACGT]{4}$", x)
}

shannon_entropy_bits <- function(p) {
  p <- p[is.finite(p) & p > 0]
  -sum(p * log2(p))
}

logo_pwm_from_counts <- function(df_counts) {
  # df_counts: motif, count (non-negative integer / numeric)
  bases <- c("A", "C", "G", "T")
  df_counts <- df_counts %>%
    mutate(motif = toupper(.data$motif)) %>%
    filter(is_valid_4mer(.data$motif), .data$count > 0)

  total <- sum(df_counts$count)
  if (!is.finite(total) || total <= 0) {
    m <- matrix(0, nrow = 4, ncol = 4, dimnames = list(bases, 1:4))
    return(m)
  }

  m <- matrix(0, nrow = 4, ncol = 4, dimnames = list(bases, 1:4))
  for (pos in 1:4) {
    b <- substr(df_counts$motif, pos, pos)
    # ensure all bases exist
    v <- tapply(df_counts$count, b, sum, default = 0)
    m[, pos] <- unname(v[bases])
  }
  m / total
}

# Load extracted end motifs ----
all_metrics <- readRDS(file.path(DATA_DIR, "processed", "all_metrics.rds"))
motifs_all <- all_4mers()

# Build 256-motif counts per sample (include zeros) ----
motif_counts <- purrr::map_dfr(all_metrics, function(m) {
  motifs <- toupper(m$motifs_5p)
  motifs <- motifs[is_valid_4mer(motifs)]
  tab <- table(factor(motifs, levels = motifs_all))
  tibble(
    sample_id = m$sample_id,
    motif = names(tab),
    count = as.integer(tab)
  )
})

motif_long <- motif_counts %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(
    total = sum(.data$count),
    frequency = ifelse(.data$total > 0, .data$count / .data$total, 0)
  ) %>%
  ungroup()

# Wide frequency matrix (rows=motifs, cols=samples) ----
freq_mat <- motif_long %>%
  select(motif, sample_id, frequency) %>%
  pivot_wider(names_from = sample_id, values_from = frequency, values_fill = 0) %>%
  column_to_rownames("motif") %>%
  as.matrix()

sample_order <- metadata %>%
  arrange(Group, sample_id) %>%
  pull(sample_id)
sample_order <- intersect(sample_order, colnames(freq_mat))
freq_mat <- freq_mat[, sample_order, drop = FALSE]

# Save tidy tables ----
readr::write_csv(
  motif_long %>% select(sample_id, Group, motif, count, frequency),
  file.path(TABLE_DIR, "end_motif_4mer_frequencies_long.csv")
)

# =============================================================================
# Figure 4A: Unsupervised hierarchical clustering heatmap (256 motifs)
# =============================================================================

# Row Z-score (standard in end-motif heatmaps)
z_mat <- t(scale(t(freq_mat)))
z_mat[!is.finite(z_mat)] <- 0

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
  name = "Z-score",
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 6),
  column_title = "Samples",
  row_title = "5' 4-mer end motifs (n = 256)",
  col = circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9)
  )
)

grDevices::cairo_pdf(
  file.path(FIG_DIR, "fig4A_end_motif_4mer_clustering_heatmap.pdf"),
  width = 9.5,
  height = 9.5
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

png(file.path(FIG_DIR, "fig4A_end_motif_4mer_clustering_heatmap.png"),
    width = 9.5, height = 9.5, units = "in", res = 450, type = "cairo-png")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

# =============================================================================
# Figure 4B: Motif Diversity Score (MDS; Shannon entropy)
# =============================================================================

mds_tbl <- tibble(sample_id = colnames(freq_mat)) %>%
  mutate(
    entropy_bits = apply(freq_mat, 2, shannon_entropy_bits),
    mds = entropy_bits / log2(nrow(freq_mat)) # normalize to [0, 1]
  ) %>%
  left_join(metadata, by = "sample_id")

readr::write_csv(mds_tbl, file.path(TABLE_DIR, "end_motif_mds.csv"))

p_mds <- ggplot(mds_tbl, aes(x = Group, y = mds, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 2.0, alpha = 0.75) +
  scale_fill_manual(values = COLORS) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "End-motif diversity (MDS)",
    subtitle = "Shannon entropy of 256 4-mer frequencies (normalized by log2(256))",
    x = NULL,
    y = "Motif diversity score (0–1)"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig4B_end_motif_mds.pdf"),
       p_mds, width = 6.2, height = 4.8, dpi = 300, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "fig4B_end_motif_mds.png"),
       p_mds, width = 6.2, height = 4.8, dpi = 450)

# =============================================================================
# Figure 4D: Sequence logos (cfDNA end motifs)
# =============================================================================

if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
  install.packages("ggseqlogo", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(ggseqlogo))

group_counts <- motif_long %>%
  group_by(Group, motif) %>%
  summarise(count = sum(count), .groups = "drop")

pwm_by_group <- group_counts %>%
  split(.$Group) %>%
  lapply(logo_pwm_from_counts)

logo_plots <- lapply(names(pwm_by_group), function(g) {
  ggseqlogo::ggseqlogo(pwm_by_group[[g]], method = "prob", seq_type = "dna") +
    labs(title = g, x = "Position", y = "Probability") +
    theme_pub(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
})

p_logo <- wrap_plots(logo_plots, ncol = 1) +
  plot_annotation(title = "cfDNA 5' end motif sequence logo (4-mer; pooled within group)") &
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.2)))

ggsave(file.path(FIG_DIR, "fig4D_end_motif_sequence_logo.pdf"),
       p_logo, width = 7.2, height = 6.8, dpi = 300, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "fig4D_end_motif_sequence_logo.png"),
       p_logo, width = 7.2, height = 6.8, dpi = 450)

message("Done.")
message("Outputs written to:")
message(sprintf("  - %s", FIG_DIR))
message(sprintf("  - %s", TABLE_DIR))

