# 03_end_motif_analysis.R
#
# End-motif analysis (5' 4-mers; 256 motifs) from fragment BED files.
#
# Fragment BED format (from `01_qc_processing.R`):
# chrom, start0, end0, name, mapq, strand, length_bp, gc
# - start0 is 0-based
# - end0 is end-exclusive (BED). length_bp == end0 - start0.
#
# Motif extraction (standard cfDNA end-motif definition; strand-aware):
# - left genomic end:  [start0, start0+4) on strand "+"
# - right genomic end: [end0-4, end0) on strand "-" (getSeq returns reverse-complement)
# -> two 4-mers per fragment.
#
# Main outputs:
# - fig4A_end_motif_4mer_clustering_heatmap.pdf/png
# - fig4B_end_motif_mds.pdf/png
# - fig4D_end_motif_sequence_logo.pdf/png
# - fig4C_top20_end_motifs_barplot.pdf/png
# - fig4E_end_motif_pca.pdf/png
# - fig4F_end_motif_volcano.pdf/png
#
# Tables:
# - end_motif_4mer_frequencies_long.csv
# - end_motif_mds.csv
# - end_motif_pca_scores.csv
# - end_motif_differential_motifs.csv
#
# References (methods conventions):
# - Liu et al., STAR Protocols (2024): workflow + row Z-score heatmap convention
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC11489062/
# - Jiang et al., Cancer Discovery (2020): 256 4-mer end motif profiling

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
  "ComplexHeatmap"
)
invisible(lapply(cran_pkgs, install_and_load, bioc = FALSE))
invisible(lapply(bioc_pkgs, install_and_load, bioc = TRUE))

# Paths ----
PROJECT_DIR <- here::here()
DATA_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/data")
RESULTS_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/results")
FIG_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
FRAG_DIR <- file.path(DATA_DIR, "processed", "fragments")

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(TABLE_DIR)) dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Input ----
metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"), show_col_types = FALSE) %>%
  mutate(Group = factor(.data$Group, levels = c("Ctrl", "ALS")))

genome <- BSgenome.Hsapiens.UCSC.hg38
valid_seqnames <- seqnames(genome)

# Style ----
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
      panel.grid.major = element_line(color = "gray85", linewidth = 0.30, linetype = "dotted"),
      panel.grid.minor = element_blank(),
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

shannon_entropy_bits <- function(p) {
  p <- p[is.finite(p) & p > 0]
  -sum(p * log2(p))
}

FRAG_COLS <- c("chrom", "start0", "end0", "name", "mapq", "strand", "length_bp", "gc")
FRAG_COL_CLASSES <- c(
  chrom = "character",
  start0 = "integer",
  end0 = "integer",
  name = "character",
  mapq = "integer",
  strand = "character",
  length_bp = "integer",
  gc = "numeric"
)

read_frag_chunk <- function(con, chunk_size) {
  utils::read.delim(
    file = con,
    header = FALSE,
    sep = "\t",
    nrows = as.integer(chunk_size),
    col.names = FRAG_COLS,
    colClasses = FRAG_COL_CLASSES,
    stringsAsFactors = FALSE,
    quote = ""
  )
}

extract_4mer_counts_from_frag_bed <- function(frag_bed_gz, sample_id, genome, valid_seqnames, motifs_all, chunk_size = 200000L) {
  counts <- setNames(integer(length(motifs_all)), motifs_all)
  con <- gzfile(frag_bed_gz, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  total_extracted <- 0L

  repeat {
    x <- tryCatch(read_frag_chunk(con, chunk_size), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) break

    x <- x %>%
      filter(
        .data$chrom %in% valid_seqnames,
        is.finite(.data$start0),
        is.finite(.data$end0),
        .data$start0 >= 0L,
        .data$end0 > .data$start0
      )
    if (nrow(x) == 0) next

    gr_l <- GRanges(
      seqnames = x$chrom,
      ranges = IRanges(start = x$start0 + 1L, width = 4L),
      strand = "+"
    )

    r_start0 <- x$end0 - 4L
    keep_r <- is.finite(r_start0) & r_start0 >= 0L
    gr_r <- GRanges(
      seqnames = x$chrom[keep_r],
      ranges = IRanges(start = r_start0[keep_r] + 1L, width = 4L),
      strand = "-"
    )

    s_l <- getSeq(genome, gr_l)
    s_r <- if (length(gr_r) > 0) getSeq(genome, gr_r) else DNAStringSet()

    motifs <- toupper(c(as.character(s_l), as.character(s_r)))
    motifs <- motifs[nchar(motifs) == 4 & grepl("^[ACGT]{4}$", motifs)]
    if (length(motifs) > 0) {
      counts <- counts + as.integer(table(factor(motifs, levels = motifs_all)))
      total_extracted <- total_extracted + length(motifs)
    }

    rm(x)
    gc(verbose = FALSE)
  }

  if (sum(counts) == 0L) {
    stop(sprintf("No motifs extracted for %s. Check contig names and hg38 availability.", sample_id), call. = FALSE)
  }

  message(sprintf("  %s: extracted %s motifs", sample_id, format(total_extracted, big.mark = ",")))
  tibble(sample_id = sample_id, motif = names(counts), count = as.integer(counts))
}

# Build motif counts ----
motifs_all <- all_4mers()

sample_frag_paths <- metadata %>%
  transmute(
    sample_id = .data$sample_id,
    frag_bed = file.path(FRAG_DIR, paste0(.data$sample_id, ".fragments.bed.gz"))
  )

missing <- sample_frag_paths %>% filter(!file.exists(.data$frag_bed))
if (nrow(missing) > 0) stop("Missing fragment BED(s):\n", paste(missing$frag_bed, collapse = "\n"))

message("Extracting 4-mer end motifs from fragment BEDs (two 5' ends per fragment).")
motif_counts <- purrr::map2_dfr(
  sample_frag_paths$sample_id,
  sample_frag_paths$frag_bed,
  ~extract_4mer_counts_from_frag_bed(.y, .x, genome = genome, valid_seqnames = valid_seqnames, motifs_all = motifs_all)
)

motif_long <- motif_counts %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(
    total = sum(.data$count),
    frequency = ifelse(.data$total > 0, .data$count / .data$total, 0)
  ) %>%
  ungroup()

readr::write_csv(
  motif_long,
  file.path(TABLE_DIR, "end_motif_4mer_frequencies_long.csv")
)

freq_mat <- motif_long %>%
  dplyr::select(motif, sample_id, frequency) %>%
  tidyr::pivot_wider(names_from = sample_id, values_from = frequency, values_fill = 0) %>%
  tibble::column_to_rownames("motif") %>%
  as.matrix()

freq_mat <- freq_mat[, match(metadata$sample_id, colnames(freq_mat))]

# Top 20 motif frequencies by group ----
motif_summary <- motif_long %>%
  dplyr::group_by(Group, motif) %>%
  dplyr::summarise(mean_freq = mean(frequency), sd_freq = sd(frequency), .groups = "drop")

top20 <- motif_summary %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(total_mean = sum(mean_freq), .groups = "drop") %>%
  arrange(desc(total_mean)) %>%
  slice_head(n = 20) %>%
  pull(motif)

p_top20 <- motif_summary %>%
  filter(motif %in% top20) %>%
  mutate(motif = factor(motif, levels = rev(top20))) %>%
  ggplot(aes(x = motif, y = mean_freq, fill = Group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = COLORS) +
  coord_flip() +
  labs(
    title = "Top 20 cfDNA 5' end motifs (4-mer)",
    subtitle = "Mean motif frequency by group",
    x = NULL,
    y = "Mean frequency",
    fill = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(FIG_DIR, "fig5A_top20_end_motifs_barplot.png"),
       p_top20, width = 7.8, height = 5.2, dpi = 450)

# Unsupervised hierarchical clustering heatmap (256 motifs) ----
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

png(file.path(FIG_DIR, "fig5B_end_motif_4mer_clustering_heatmap.png"),
    width = 5000, height = 7000, res = 450)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Motif Diversity Score (MDS; Shannon entropy) ----
mds_tbl <- tibble(sample_id = colnames(freq_mat)) %>%
  mutate(
    entropy_bits = apply(freq_mat, 2, shannon_entropy_bits),
    mds = entropy_bits / log2(nrow(freq_mat))
  ) %>%
  left_join(metadata, by = "sample_id")

readr::write_csv(mds_tbl, file.path(TABLE_DIR, "end_motif_mds.csv"))

p_mds <- ggplot(mds_tbl, aes(x = Group, y = mds, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 2.0, alpha = 0.75) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "End-motif diversity (MDS)",
    subtitle = "Shannon entropy of 256 4-mer frequencies (normalized by log2(256))",
    x = NULL,
    y = "Motif diversity score (0–1)"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig5C_end_motif_mds.png"),
       p_mds, width = 8, height = 7, dpi = 450)

# PCA on samples using 256 motif frequencies ----
X <- t(freq_mat) # samples x motifs
X <- log10(X + 1e-6)
pca <- prcomp(X, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(metadata, by = "sample_id")

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = sample_id)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.6, alpha = 0.2) +
  geom_text(size = 2.7, hjust = 0.5, vjust = -0.7, show.legend = FALSE, alpha = 0.88, check_overlap = TRUE) +
  scale_color_manual(values = COLORS) +
  labs(
    title = "PCA of samples using 256 end-motif frequencies",
    x = pc1_lab,
    y = pc2_lab,
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(FIG_DIR, "fig5D_end_motif_pca.png"),
       p_pca, width = 8, height = 7, dpi = 450)

# Differential motifs + volcano ----
eps <- 1e-6
diff_tbl <- motif_long %>%
  dplyr::select(Group, motif, frequency) %>%
  group_by(motif) %>%
  summarise(
    mean_Ctrl = mean(frequency[Group == "Ctrl"]),
    mean_ALS = mean(frequency[Group == "ALS"]),
    log2fc = log2((mean_ALS + eps) / (mean_Ctrl + eps)),
    p = if (dplyr::n_distinct(Group) < 2) NA_real_ else wilcox.test(frequency ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p, method = "BH"),
    neglog10_padj = -log10(padj),
    neglog10_p = -log10(p)
  )

readr::write_csv(diff_tbl, file.path(TABLE_DIR, "end_motif_differential_motifs.csv"))

sig_cut <- 0.05
lfc_cut <- 0.05
diff_tbl <- diff_tbl %>%
  mutate(
    status = case_when(
      !is.na(p) & p < sig_cut & log2fc > lfc_cut ~ "Higher in ALS",
      !is.na(p) & p < sig_cut & log2fc < -lfc_cut ~ "Higher in Ctrl",
      TRUE ~ "Not significant"
    ),
    status = factor(status, levels = c("Higher in ALS", "Higher in Ctrl", "Not significant"))
  )

label_tbl <- diff_tbl %>%
  filter(status != "Not significant") %>%
  arrange(p) %>%
  slice_head(n = 10)

p_volcano <- ggplot(diff_tbl, aes(x = log2fc, y = neglog10_p, color = status)) +
  geom_hline(yintercept = -log10(sig_cut), linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_point(alpha = 0.85, size = 2.2) +
  geom_text(
    data = label_tbl,
    aes(label = motif, color = status),
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Higher in ALS" = unname(COLORS["ALS"]),
      "Higher in Ctrl" = unname(COLORS["Ctrl"]),
      "Not significant" = "gray70"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Differential end motifs (ALS vs Ctrl)",
    subtitle = "Wilcoxon test per motif",
    x = "log2 fold-change (ALS / Ctrl)",
    y = expression(-log[10]("FDR")),
    color = NULL
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(FIG_DIR, "fig5E_end_motif_volcano.png"),
       p_volcano, width = 7.2, height = 5.6, dpi = 450)

# Draw combined figure ----
ht_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
)
p_ht <- patchwork::wrap_elements(full = ht_grob)

p_end_motif_combined <- (p_top20 | p_ht) / (p_mds | p_pca | p_volcano) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(2.5, 1.6))

ggsave(
  file.path(FIG_DIR, "fig5_end_motif_analysis.png"),
  p_end_motif_combined,
  width = 18, height = 18, dpi = 300
)
