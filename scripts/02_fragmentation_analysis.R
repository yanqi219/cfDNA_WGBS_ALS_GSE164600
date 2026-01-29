# 02_fragmentation_analysis.R
# https://www.sciencedirect.com/science/article/pii/S0009898124022861?casa_token=xJj0Aav-s0IAAAAA:Epr21Bl-R3YfswQro5bbHfV_sMgJHIDpGc2ALqO5zA_4S5AKenI-dGB6gEVBbgOs_ZqzHHZuiug

# Insert size / fragmentation analysis

# What this script produces:
# - Fragment size distribution summaries (30–500 bp) + short/long fragment ratios
# - Genome-binned fragment size coverage (short vs long) and fragmentation ratio tracks
# - Nucleosome positioning proxy: TSS-centered cfDNA midpoint density profiles (Snyder et al.-style)
#
# References (conceptual inspiration; implementation is self-contained):
# - Snyder et al., Cell (2016): cfDNA nucleosome footprint / TSS-centered profiles
# - Cristiano et al., Nature (2019): genome-wide fragmentation patterns; short/long ratio tracks
# - cfDNAPro (Bioconductor): standardized fragment size metrics (prop/cdf, periodicity)
#
# Notes:
# - This script consumes gzipped fragment BEDs produced by scripts/01_qc_processing.R:
#   chrom, start0, end0, name, mapq, strand, length_bp, gc
# - start0/end0 are BED-style 0-based (end exclusive); we convert to 1-based for GRanges ops.
#
# Output:
# - results/tables/fragmentation_metrics.csv
# - results/figures/fig3_fragmentation_overview.png
# - results/figures/fig4_fragmentation_ratio_tracks.png
# - results/figures/fig5_tss_nucleosome_positioning.png

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
  "caret", "pROC", "patchwork", "here", "parallel", "ggpubr"
)
bioc_pkgs <- c(
  "Rsamtools", "cigarillo", "GenomicRanges", "IRanges", "ChIPseeker", "biomaRt", "org.Hs.eg.db", "AnnotationDbi",
  "Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "GenomeInfoDb", "TxDb.Hsapiens.UCSC.hg38.knownGene", "GenomicFeatures"
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

# Load sample metadata ----
metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"), show_col_types = FALSE)
metadata <- metadata %>%
  dplyr::mutate(Group = factor(.data$Group, levels = c("Ctrl", "ALS")))

# Genome ----
genome <- BSgenome.Hsapiens.UCSC.hg38

# Annotation ----
# Widely used Bioconductor annotation source for promoters/TSS:
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Theme + palette ----
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

# Chunked reader for gzipped fragment BEDs without readr.
# Reads TSV with no header: chrom, start0, end0, name, mapq, strand, length_bp
read_fragments_chunked <- function(path, chunk_size = 200000L, FUN) {
  chunk_size <- as.integer(chunk_size)
  con <- gzfile(path, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  repeat {
    x <- utils::read.delim(
      file = con,
      header = FALSE,
      sep = "\t",
      nrows = chunk_size,
      col.names = FRAG_COLS,
      colClasses = FRAG_COL_CLASSES,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    )
    if (nrow(x) == 0) break
    FUN(x)
  }

  invisible(TRUE)
}

stop_if_missing_fragments <- function(sample_ids, frag_dir = FRAG_DIR) {
  paths <- file.path(frag_dir, paste0(sample_ids, ".fragments.bed.gz"))
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing fragment BEDs (run scripts/01_qc_processing.R first):\n",
      paste0("  - ", missing, collapse = "\n")
    )
  }
  invisible(paths)
}

moving_average <- function(x, k = 21L) {
  k <- as.integer(k)
  if (k < 3L) return(x)
  if (k %% 2L == 0L) k <- k + 1L
  w <- rep(1 / k, k)
  as.numeric(stats::filter(x, w, sides = 2, method = "convolution"))
}

# Convert a BED chunk (0-based, end exclusive) to fragment midpoint GRanges (1bp)
chunk_to_midpoint_gr <- function(df) {
  # Convert to 1-based inclusive coords:
  start1 <- df$start0 + 1L
  end1 <- df$end0
  mid <- as.integer(floor((start1 + end1) / 2))
  GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = mid, width = 1L),
    strand = "*",
    length_bp = df$length_bp
  )
}

# 5' fragment start coordinates (strand-aware) ----
# Fragment BED definition from scripts/01_qc_processing.R:
# - start0 is 0-based
# - end0 is compatible with 1-based inclusive end (end1 = end0)
# - strand is derived from the first mate (bedpe-like convention)
# Therefore 5' coordinate (1-based) is:
# - '+' : start0 + 1
# - '-' : end0
chunk_to_5p_start_gr <- function(df) {
  start1 <- df$start0 + 1L
  end1 <- df$end0
  pos1 <- ifelse(df$strand == "-", end1, start1)

  GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = as.integer(pos1), width = 1L),
    strand = "*"
  )
}

# Inputs ----
stop_if_missing_fragments(metadata$sample_id, FRAG_DIR)
fragment_paths <- file.path(FRAG_DIR, paste0(metadata$sample_id, ".fragments.bed.gz"))
names(fragment_paths) <- metadata$sample_id
all_metrics <- readRDS(file.path(DATA_DIR, "processed", "all_metrics.rds"))

# Figure 2: Insert size (fragment length) QC ----
frag_df <- map_dfr(all_metrics, function(m) {
  data.frame(
    sample_id = m$sample_id,
    fragment_length = as.integer(m$fragment_lengths)
  )
}) %>%
  dplyr::left_join(metadata, by = "sample_id")

calc_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_integer_)
  x <- as.integer(x)
  tabulate(x, nbins = max(x)) |> which.max()
}

insert_metrics <- frag_df %>%
  dplyr::group_by(sample_id, Group) %>%
  dplyr::summarise(
    n_fragments = dplyr::n(),
    median_bp = stats::median(fragment_length),
    mode_bp = calc_mode(fragment_length),
    short_n = sum(fragment_length >= 100 & fragment_length <= 150),
    long_n = sum(fragment_length >= 151 & fragment_length <= 220),
    short_long_ratio = dplyr::if_else(long_n > 0, short_n / long_n, NA_real_),
    .groups = "drop"
  )
print(
  insert_metrics %>%
    dplyr::arrange(Group, sample_id) %>%
    dplyr::mutate(
      median_bp = round(median_bp, 1),
      short_long_ratio = round(short_long_ratio, 4)
    )
)

## plot ----
### insert size per group ----
metric_long <- insert_metrics %>%
  dplyr::select(sample_id, Group, median_bp, mode_bp, short_long_ratio) %>%
  tidyr::pivot_longer(
    cols = c(median_bp, mode_bp, short_long_ratio),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("median_bp", "mode_bp", "short_long_ratio"),
      labels = c("Median (bp)", "Mode (bp)", "Short/Long ratio")
    )
  )

is_metrics <- ggplot(metric_long, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = Group), width = 0.12, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  ggpubr::stat_compare_means(method = "t.test") +
  labs(
    title = "Insert size metrics",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "none")

### sawtooth plot (binwidth = 1) per sample ----
frag_df_saw <- frag_df %>%
  dplyr::filter(fragment_length >= 80 & fragment_length <= 450)

counts_saw <- frag_df_saw %>%
  dplyr::count(sample_id, Group, fragment_length, name = "n")

ann_df <- counts_saw %>%
  dplyr::group_by(sample_id, Group) %>%
  dplyr::summarise(ymax = max(n), .groups = "drop") %>%
  dplyr::left_join(insert_metrics, by = c("sample_id", "Group")) %>%
  dplyr::mutate(
    label = sprintf("Median: %.0f bp", median_bp),
    x = 245,
    y = 0.92 * ymax
  )

is_saw <- ggplot(counts_saw, aes(x = fragment_length, y = n)) +
  geom_line(linewidth = 0.35, color = "gray20") +
  geom_vline(xintercept = 167, linetype = "dashed", color = "#2E2E2E", linewidth = 0.4) +
  geom_vline(xintercept = 334, linetype = "dashed", color = "#2E2E2E", linewidth = 0.4) +
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.6,
    lineheight = 0.95,
    color = "gray10"
  ) +
  facet_wrap(~sample_id, scales = "free_y", ncol = 4) +
  labs(
    title = "Sawtooth insert-size pattern (binwidth = 1)",
    subtitle = "Dashed lines: mono-nucleosome (~167 bp) and di-nucleosome (~334 bp)",
    x = "Fragment length (bp)",
    y = "Count"
  ) +
  theme_pub(base_size = 10) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8)
  )

### Group-level median fragment size distributions ----
len_min <- 30L
len_max <- 500L
len_grid <- len_min:len_max

frag_df_30_500 <- frag_df %>%
  dplyr::filter(!is.na(fragment_length), fragment_length >= len_min, fragment_length <= len_max) %>%
  dplyr::select(sample_id, Group, fragment_length)

sample_groups <- frag_df_30_500 %>%
  dplyr::distinct(sample_id, Group)

counts_30_500 <- frag_df_30_500 %>%
  dplyr::count(sample_id, Group, fragment_length, name = "n")

counts_full <- tidyr::crossing(sample_groups, fragment_length = len_grid) %>%
  dplyr::left_join(counts_30_500, by = c("sample_id", "Group", "fragment_length")) %>%
  dplyr::mutate(n = dplyr::coalesce(n, 0L)) %>%
  dplyr::group_by(sample_id, Group) %>%
  dplyr::mutate(
    prop = n / sum(n),
    cum_prop = cumsum(prop)
  ) %>%
  dplyr::ungroup()

group_median_dist <- counts_full %>%
  dplyr::group_by(Group, fragment_length) %>%
  dplyr::summarise(
    median_prop = stats::median(prop, na.rm = TRUE),
    median_cum_prop = stats::median(cum_prop, na.rm = TRUE),
    .groups = "drop"
  )

is_prop <- ggplot(group_median_dist, aes(x = fragment_length, y = median_prop, color = Group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 167, linetype = "dashed", color = "#2E2E2E", linewidth = 0.4) +
  geom_vline(xintercept = 334, linetype = "dashed", color = "#2E2E2E", linewidth = 0.4) +
  scale_color_manual(values = COLORS) +
  scale_x_continuous(limits = c(len_min, len_max), breaks = seq(50, 500, by = 50)) +
  labs(
    title = "Median fragment size distribution (by group)",
    subtitle = sprintf("Proportion at each length (%dbp–%dbp)", len_min, len_max),
    x = "Fragment length (bp)",
    y = "Median proportion",
    color = "Group"
  ) +
  theme_pub(base_size = 11) +
  theme(legend.position = "top")

is_cum <- ggplot(group_median_dist, aes(x = fragment_length, y = median_cum_prop, color = Group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 167, linetype = "dashed", color = "#2E2E2E", linewidth = 0.4) +
  geom_vline(xintercept = 334, linetype = "dashed", color = "#2E2E2E", linewidth = 0.4) +
  scale_color_manual(values = COLORS) +
  scale_x_continuous(limits = c(len_min, len_max), breaks = seq(50, 500, by = 50)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Median cumulative fragment distribution (by group)",
    subtitle = sprintf("Proportion of reads with length ≤ N (%dbp–%dbp)", len_min, len_max),
    x = "Fragment length (bp)",
    y = "Median cumulative proportion",
    color = "Group"
  ) +
  theme_pub(base_size = 11) +
  theme(legend.position = "top")

is_combined <- is_metrics / is_saw / (is_prop | is_cum) +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1, 1)
  )

ggsave(file.path(FIG_DIR, "fig2_insert_size.png"),
       is_combined, width = 14, height = 18, dpi = 300)

# Figure 3: Genome-binned fragment size coverage + short/long ratio tracks ----
message("Computing genome-binned fragmentation ratio tracks (short vs long)...")

bin_size <- 2e6
chroms <- paste0("chr", c(21))
seqlens <- GenomeInfoDb::seqlengths(genome)[chroms]
seqlens <- seqlens[!is.na(seqlens)]

bins_gr <- GenomicRanges::tileGenome(
  seqlengths = seqlens,
  tilewidth = bin_size,
  cut.last.tile.in.chrom = TRUE
)
bins_df <- tibble::tibble(
  bin_id = seq_along(bins_gr),
  chrom = as.character(GenomeInfoDb::seqnames(bins_gr)),
  start = start(bins_gr),
  end = end(bins_gr),
  mid = as.integer(floor((start + end) / 2))
)

## GC content per bin ----
# Note: gc is computed on A/C/G/T only (Ns excluded); we also track n_frac for filtering.
bin_seqs <- Biostrings::getSeq(genome, bins_gr)
freq <- Biostrings::letterFrequency(bin_seqs, letters = c("A", "C", "G", "T", "N"), as.prob = FALSE)
acgt <- rowSums(freq[, c("A", "C", "G", "T"), drop = FALSE])
gc <- (freq[, "G"] + freq[, "C"]) / pmax(acgt, 1)
n_frac <- freq[, "N"] / Biostrings::width(bin_seqs)

bins_df <- bins_df %>%
  dplyr::mutate(
    gc = as.numeric(gc),
    n_frac = as.numeric(n_frac)
  )
S4Vectors::mcols(bins_gr)$gc <- bins_df$gc
S4Vectors::mcols(bins_gr)$n_frac <- bins_df$n_frac

# Filter bins for LOESS fitting: remove extreme GC and high-N bins
gc_min <- 0.10
gc_max <- 0.85
n_frac_max <- 0.45
bins_keep_for_gc <- is.finite(bins_df$gc) &
  bins_df$gc >= gc_min & bins_df$gc <= gc_max &
  is.finite(bins_df$n_frac) & bins_df$n_frac <= n_frac_max

## LOESS GC-bias correction ----
correct_gc_bias <- function(counts, gc_values, fit_mask = NULL, pseudocount = 1, span = 0.75) {
  counts <- as.numeric(counts)
  gc_values <- as.numeric(gc_values)
  pseudocount <- as.numeric(pseudocount)

  if (is.null(fit_mask)) {
    fit_mask <- is.finite(gc_values) & is.finite(counts)
  } else {
    fit_mask <- as.logical(fit_mask) & is.finite(gc_values) & is.finite(counts)
  }
  fit_mask_in <- fit_mask

  y <- log(counts + pseudocount)
  fit_mask <- fit_mask & is.finite(y)

  # Need enough bins to fit a stable loess curve
  if (sum(fit_mask) < 20L) {
    return(counts)
  }

  fit <- stats::loess(
    y[fit_mask] ~ gc_values[fit_mask],
    span = span,
    degree = 1,
    family = "symmetric",
    control = stats::loess.control(surface = "direct")
  )

  pred <- stats::predict(fit, newdata = gc_values)
  med <- stats::median(counts[fit_mask] + pseudocount, na.rm = TRUE)
  corrected <- (counts + pseudocount) / exp(pred) * med

  # For bins excluded from the GC model fit (extreme GC / high-N), keep NA to avoid artifacts.
  corrected[!fit_mask_in] <- NA_real_
  corrected[!is.finite(pred)] <- NA_real_
  corrected
}

bin_counts_one_sample <- function(path, bins_gr, short_range = c(90L, 150L), long_range = c(151L, 220L)) {
  short_counts <- integer(length(bins_gr))
  long_counts <- integer(length(bins_gr))

  cb <- function(x) {
    l <- x$length_bp
    keep <- !is.na(l) & l >= 30L & l <= 500L
    if (!any(keep)) return(invisible())
    x <- x[keep, , drop = FALSE]
    l <- x$length_bp

    mid_gr <- chunk_to_midpoint_gr(x)

    # Short
    idx_s <- which(l >= short_range[1] & l <= short_range[2])
    if (length(idx_s) > 0) {
      hits <- findOverlaps(mid_gr[idx_s], bins_gr, ignore.strand = TRUE)
      if (length(hits) > 0) {
        tab <- tabulate(S4Vectors::subjectHits(hits), nbins = length(bins_gr))
        short_counts <<- short_counts + tab
      }
    }

    # Long
    idx_l <- which(l >= long_range[1] & l <= long_range[2])
    if (length(idx_l) > 0) {
      hits <- findOverlaps(mid_gr[idx_l], bins_gr, ignore.strand = TRUE)
      if (length(hits) > 0) {
        tab <- tabulate(S4Vectors::subjectHits(hits), nbins = length(bins_gr))
        long_counts <<- long_counts + tab
      }
    }

    invisible()
  }

  read_fragments_chunked(path = path, chunk_size = 200000L, FUN = cb)

  # GC correction (fit on bins_keep_for_gc; apply per-sample)
  short_corrected <- correct_gc_bias(short_counts, bins_df$gc, fit_mask = bins_keep_for_gc, pseudocount = 1, span = 0.75)
  long_corrected <- correct_gc_bias(long_counts, bins_df$gc, fit_mask = bins_keep_for_gc, pseudocount = 1, span = 0.75)
  ratio_corrected <- short_corrected / (long_corrected + 1)

  tibble::tibble(
    bin_id = seq_along(bins_gr),
    short_raw = short_counts,
    long_raw = long_counts,
    short_gc_corrected = short_corrected,
    long_gc_corrected = long_corrected,
    short_long_ratio = ratio_corrected
  )
}

bin_tracks <- purrr::imap_dfr(fragment_paths, function(path, sample_id) {
  bc <- bin_counts_one_sample(path, bins_gr = bins_gr)
  bc %>%
    dplyr::mutate(
      sample_id = sample_id,
      Group = metadata$Group[match(sample_id, metadata$sample_id)]
    )
})

bin_tracks_summary <- bin_tracks %>%
  dplyr::group_by(Group, bin_id) %>%
  dplyr::summarise(
    median_short_long_ratio = stats::median(short_long_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(bins_df, by = "bin_id") %>%
  dplyr::arrange(match(chrom, chroms), start) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(
    chrom = factor(chrom, levels = chroms),
    genome_x = {
      # cumulative coordinate per chromosome
      chroms_present <- chroms[chroms %in% names(seqlens)]
      chr_lens <- as.numeric(seqlens[chroms_present])
      chr_offsets <- cumsum(c(0, head(chr_lens, -1L)))
      names(chr_offsets) <- chroms_present
      chr_offsets[as.character(chrom)] + mid
    }
  ) %>%
  dplyr::ungroup()

chrom_boundaries <- tibble::tibble(
  chrom = chroms[chroms %in% names(seqlens)],
  chr_len = as.numeric(seqlens[chroms[chroms %in% names(seqlens)]])
) %>%
  dplyr::mutate(
    offset = cumsum(dplyr::lag(chr_len, default = 0)),
    boundary = offset + chr_len,
    center = offset + chr_len / 2
  )

bin_label <- if (is.finite(bin_size) && bin_size >= 1e6) {
  sprintf("%.0f Mb", bin_size / 1e6)
} else if (is.finite(bin_size) && bin_size >= 1e3) {
  sprintf("%.0f kb", bin_size / 1e3)
} else {
  sprintf("%s bp", format(bin_size, scientific = FALSE, big.mark = ","))
}

p_tracks <- ggplot(bin_tracks_summary, aes(x = genome_x, y = median_short_long_ratio, color = Group)) +
  geom_hline(yintercept = 1, color = "gray70", linewidth = 0.35) +
  geom_line(linewidth = 0.65, alpha = 0.95) +
  geom_vline(
    data = chrom_boundaries,
    aes(xintercept = boundary),
    inherit.aes = FALSE,
    color = "gray88",
    linewidth = 0.3
  ) +
  scale_color_manual(values = COLORS) +
  scale_x_continuous(
    breaks = chrom_boundaries$center,
    labels = as.character(chrom_boundaries$chrom),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Genome-wide fragmentation ratio track (group median)",
    subtitle = sprintf(
      "GC-corrected (LOESS): median short/long ratio in %s bins (GC in [%.2f, %.2f], N ≤ %.2f)",
      bin_label, gc_min, gc_max, n_frac_max
    ),
    x = "Chromosome",
    y = "Median short/long ratio",
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

ggsave(
  file.path(FIG_DIR, "fig3_fragmentation_ratio_tracks.png"),
  p_tracks,
  width = 16, height = 6, dpi = 300
)

## Chromosome-style visualization with Gviz ----
if (!requireNamespace("Gviz", quietly = TRUE)) {
  message("Package `Gviz` not installed; skipping Gviz chromosome track plot.")
} else {
  suppressPackageStartupMessages(library(Gviz))

genome_build <- "hg38"
chrom <- as.character(chroms[[1]])
from_bp <- min(bins_df$start, na.rm = TRUE)
to_bp <- max(bins_df$end, na.rm = TRUE)

# One value per bin (same order as bins_df)
gr_bins <- GenomicRanges::GRanges(
  seqnames = bins_df$chrom,
  ranges = IRanges::IRanges(start = bins_df$start, end = bins_df$end)
)

# ALS/Ctrl ratio of the group-median short/long ratios (per bin)
eps <- 1e-6
ratio_wide <- bin_tracks_summary %>%
  dplyr::select(bin_id, Group, median_short_long_ratio) %>%
  tidyr::pivot_wider(names_from = Group, values_from = median_short_long_ratio) %>%
  dplyr::arrange(bin_id) %>%
  dplyr::mutate(
    als_over_ctrl = (`ALS` + eps) / (`Ctrl` + eps),
    log2_als_over_ctrl = log2(als_over_ctrl)
  )

# Prefer log2 scale for interpretability (baseline = 0), but label as ALS/Ctrl
y <- ratio_wide$log2_als_over_ctrl
y_lim <- stats::quantile(y, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
y_pad <- 0.15 * diff(y_lim)
y_lim <- c(y_lim[1] - y_pad, y_lim[2] + y_pad)

axis_track <- Gviz::GenomeAxisTrack(genome = genome_build, chromosome = chrom, name = "")

ratio_track <- Gviz::DataTrack(
  range = gr_bins,
  data = y,
  genome = genome_build,
  chromosome = chrom,
  name = "log2(ALS/Ctrl)\nshort/long",
  type = "l",
  col = "#2E2E2E",
  lwd = 2,
  ylim = y_lim,
  baseline = 0,
  col.baseline = "gray45",
  lty.baseline = 2,
  lwd.baseline = 1,
  cex.axis = 0.9,
  cex.title = 0.70
)

ideo_track <- tryCatch(
  {
    Gviz::IdeogramTrack(genome = genome_build, chromosome = chrom)
  },
  error = function(e) NULL
)

tracks <- list(axis_track, ratio_track)
if (!is.null(ideo_track)) tracks <- c(list(ideo_track), tracks)

png(
  file.path(FIG_DIR, "fig3_fragmentation_ratio_tracks_gviz.png"),
  width = 2400,
  height = if (!is.null(ideo_track)) 900 else 700,
  res = 300
)
Gviz::plotTracks(
  tracks,
  from = from_bp,
  to = to_bp,
  # main = "Fragmentation ratio track (ALS/Ctrl)",
  background.title = "white",
  col.title = "black",
  cex.title = 0.95,
  background.panel = "white",
  col.axis = "black",
  col.frame = "gray85"
)
dev.off()
}

# Figure 4: Fragment start distributions ----
extract_fragment_starts <- function(fragment_bed_gz, chroms = chroms, chunk_size = 200000L) {
  gr_list <- list()
  i <- 0L
  read_fragments_chunked(
    path = fragment_bed_gz,
    chunk_size = chunk_size,
    FUN = function(x) {
      gr <- chunk_to_5p_start_gr(x)
      gr <- gr[as.character(GenomeInfoDb::seqnames(gr)) %in% chroms]
      if (length(gr) == 0) return(invisible())
      i <<- i + 1L
      gr_list[[i]] <<- gr
      invisible()
    }
  )
  do.call(c, gr_list)
}

## Genomic annotation + stacked bar plot ----
annotate_fragment_starts_one_sample <- function(fragment_bed_gz, sample_id, group, chunk_size = 200000L) {
  counts <- c(Promoter = 0L, Exon = 0L, Intron = 0L, Distal_Intergenic = 0L, Three_UTR = 0L, Five_UTR = 0L, Downstream = 0L)
  total <- 0L

  read_fragments_chunked(
    path = fragment_bed_gz,
    chunk_size = chunk_size,
    FUN = function(x) {
      gr <- chunk_to_5p_start_gr(x)
      if (length(gr) == 0) return(invisible())
      total <<- total + length(gr)

      gr <- ChIPseeker::annotatePeak(gr, TxDb = txdb, level = "gene", addFlankGeneInfo=TRUE)
      gr_anno <- data.frame(gr@anno)
      gr_anno <- gr_anno %>%
        dplyr::mutate(annotation_simple = gsub(" \\(.*", "", annotation)) %>%
        dplyr::mutate(annotation_simple = case_when(
          annotation_simple == "5' UTR" ~ "Five_UTR",
          annotation_simple == "3' UTR" ~ "Three_UTR",
          annotation_simple == "Distal Intergenic" ~ "Distal_Intergenic",
          TRUE ~ annotation_simple
        ))

      counts_table <- table(gr_anno$annotation_simple)
      for (annot in names(counts_table)) {
        counts[[annot]] <<- counts[[annot]] + counts_table[[annot]]
      }
      invisible()
    }
  )

  tibble::tibble(
    sample_id = sample_id,
    Group = group,
    total_starts = total,
    Promoter = counts[["Promoter"]],
    Exon = counts[["Exon"]],
    Intron = counts[["Intron"]],
    Three_UTR = counts[["Three_UTR"]],
    Five_UTR = counts[["Five_UTR"]],
    Distal_Intergenic = counts[["Distal_Intergenic"]],
    Downstream = counts[["Downstream"]]
  )
}

message("Annotating fragment 5' starts (Promoter/Exon/Intron/Intergenic)...")
start_annot <- purrr::imap_dfr(fragment_paths, function(path, sample_id) {
  annotate_fragment_starts_one_sample(
    fragment_bed_gz = path,
    sample_id = sample_id,
    group = metadata$Group[match(sample_id, metadata$sample_id)]
  )
})

start_annot_long <- start_annot %>%
  tidyr::pivot_longer(cols = c(Promoter, Exon, Intron, Three_UTR, Five_UTR, Distal_Intergenic, Downstream), names_to = "annotation", values_to = "n") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(annotation = factor(annotation, levels = c("Promoter", "Exon", "Intron", "Three_UTR", "Five_UTR", "Distal_Intergenic", "Downstream")))

p_start_annot <- ggplot(start_annot_long, aes(x = sample_id, y = pct, fill = annotation)) +
  geom_col(width = 0.85) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Genomic distribution of fragment 5' start sites",
    x = NULL,
    y = "Percent of fragment",
    fill = "Annotation"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "top"
  )

# ggsave(
#   file.path(FIG_DIR, "fig4_fragment_start_annotation_chr21.png"),
#   p_start_annot,
#   width = 14, height = 5.5, dpi = 300
# )

## Aggregate TSS enrichment ----
# Best practice in Bioconductor is to derive TSS/promoters directly from a TxDb object.
# Using gene-level ranges avoids redundant TSSs from alternative transcripts.
genes_chr21 <- suppressMessages(GenomicFeatures::genes(txdb))
genes_chr21 <- genes_chr21[as.character(GenomeInfoDb::seqnames(genes_chr21)) %in% chroms]
profile_tss_enrichment_starts <- function(fragment_bed_gz, tss_windows, tss_site, tss_strand,
                                         flank = 2000L, chunk_size = 200000L, edge_width = 200L) {
  flank <- as.integer(flank)
  rel_grid <- (-flank):flank
  prof <- numeric(length(rel_grid))
  n_total <- 0L

  read_fragments_chunked(
    path = fragment_bed_gz,
    chunk_size = chunk_size,
    FUN = function(x) {
      gr <- chunk_to_5p_start_gr(x)
      if (length(gr) == 0) return(invisible())
      n_total <<- n_total + length(gr)

      hits <- GenomicRanges::findOverlaps(gr, tss_windows, ignore.strand = TRUE)
      if (length(hits) == 0) return(invisible())

      q <- S4Vectors::queryHits(hits)
      s <- S4Vectors::subjectHits(hits)
      pos <- start(gr[q])
      rel <- pos - tss_site[s]
      rel[tss_strand[s] == "-"] <- -rel[tss_strand[s] == "-"]
      rel <- rel[rel >= -flank & rel <= flank]
      if (length(rel) > 0) prof <<- prof + tabulate(rel + flank + 1L, nbins = length(rel_grid))
      invisible()
    }
  )

  cpm <- (prof / max(n_total, 1L)) * 1e6
  edge_width <- as.integer(edge_width)
  edge_idx <- c(seq_len(edge_width), (length(cpm) - edge_width + 1L):length(cpm))
  edge_mean <- mean(cpm[edge_idx], na.rm = TRUE)
  enrich <- if (is.finite(edge_mean) && edge_mean > 0) cpm / edge_mean else rep(NA_real_, length(cpm))

  tibble::tibble(rel_bp = rel_grid, enrichment = enrich)
}

message("Computing TSS enrichment of fragment 5' starts (chr21 gene-level TSS; TxDb promoters)...")
tss_points_pc <- unique(GenomicRanges::promoters(genes_chr21, upstream = 0L, downstream = 1L))
tss_windows_pc <- unique(GenomicRanges::promoters(genes_chr21, upstream = 2000L, downstream = 2000L))
tss_site_pc <- start(tss_points_pc)
tss_strand_pc <- as.character(strand(tss_points_pc))

tss_profiles <- purrr::imap_dfr(fragment_paths, function(path, sample_id) {
  profile_tss_enrichment_starts(
    fragment_bed_gz = path,
    tss_windows = tss_windows_pc,
    tss_site = tss_site_pc,
    tss_strand = tss_strand_pc,
    flank = 2000L
  ) %>%
    dplyr::mutate(
      sample_id = sample_id,
      Group = metadata$Group[match(sample_id, metadata$sample_id)]
    )
})

tss_summary <- tss_profiles %>%
  dplyr::group_by(Group, rel_bp) %>%
  dplyr::summarise(
    mean_enrich = mean(enrichment, na.rm = TRUE),
    se_enrich = stats::sd(enrichment, na.rm = TRUE) / sqrt(sum(is.finite(enrichment))),
    n = sum(is.finite(enrichment)),
    .groups = "drop"
  )

tss_smooth_k <- 40L # small running-mean to reduce jaggedness while preserving ~10 bp phasing
tss_summary_s <- tss_summary %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(rel_bp) %>%
  dplyr::mutate(
    mean_smooth = moving_average(mean_enrich, k = tss_smooth_k),
    se_smooth = moving_average(se_enrich, k = tss_smooth_k),
    mean_smooth = dplyr::if_else(is.na(mean_smooth), mean_enrich, mean_smooth),
    se_smooth = dplyr::if_else(is.na(se_smooth), se_enrich, se_smooth)
  ) %>%
  dplyr::ungroup()

p_tss_global <- ggplot(tss_summary_s, aes(x = rel_bp, color = Group, fill = Group)) +
  geom_hline(yintercept = 1, color = "gray60", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.35) +
  # Light raw mean (keeps fine-scale periodicity visible)
  geom_line(aes(y = mean_smooth), linewidth = 0.95, na.rm = TRUE) +
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "TSS enrichment of fragment 5' start sites (chr21 protein-coding genes)",
    subtitle = sprintf(
      "Normalized by mean signal at window edges (±2000 bp; baseline = 1). Display: %d-bp running mean overlaid.",
      tss_smooth_k
    ),
    x = "Distance to TSS (bp)",
    y = "Normalized fragment density",
    color = "Group",
    fill = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# ggsave(
#   file.path(FIG_DIR, "fig5_tss_enrichment_starts_chr21_protein_coding.png"),
#   p_tss_global,
#   width = 8, height = 5.5, dpi = 300
# )

## Muscle-related genes (COL6A1, COL6A2, SOD1) ----
muscle_symbols <- c("COL6A1", "COL6A2")

muscle_map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = muscle_symbols,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
) %>%
  dplyr::filter(!is.na(.data$ENTREZID)) %>%
  dplyr::mutate(ENTREZID = as.character(.data$ENTREZID))

genes_muscle <- genes_chr21[as.character(genes_chr21$gene_id) %in% unique(muscle_map$ENTREZID)]
tss_points_m <- unique(GenomicRanges::promoters(genes_muscle, upstream = 0L, downstream = 1L))
tss_windows_m <- unique(GenomicRanges::promoters(genes_muscle, upstream = 2000L, downstream = 2000L))
tss_site_m <- start(tss_points_m)
tss_strand_m <- as.character(strand(tss_points_m))

tss_profiles_m <- purrr::imap_dfr(fragment_paths, function(path, sample_id) {
  profile_tss_enrichment_starts(
    fragment_bed_gz = path,
    tss_windows = tss_windows_m,
    tss_site = tss_site_m,
    tss_strand = tss_strand_m,
    flank = 2000L
  ) %>%
    dplyr::mutate(
      sample_id = sample_id,
      Group = metadata$Group[match(sample_id, metadata$sample_id)]
    )
})

tss_summary_m <- tss_profiles_m %>%
  dplyr::group_by(Group, rel_bp) %>%
  dplyr::summarise(
    mean_enrich = mean(enrichment, na.rm = TRUE),
    se_enrich = stats::sd(enrichment, na.rm = TRUE) / sqrt(sum(is.finite(enrichment))),
    n = sum(is.finite(enrichment)),
    .groups = "drop"
  )

genes_label <- paste0(unique(muscle_map$SYMBOL), collapse = ", ")
tss_summary_m_s <- tss_summary_m %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(rel_bp) %>%
  dplyr::mutate(
    mean_smooth = moving_average(mean_enrich, k = 5),
    se_smooth = moving_average(se_enrich, k = 5),
    mean_smooth = dplyr::if_else(is.na(mean_smooth), mean_enrich, mean_smooth),
    se_smooth = dplyr::if_else(is.na(se_smooth), se_enrich, se_smooth)
  ) %>%
  dplyr::ungroup()

p_tss_muscle <- ggplot(tss_summary_m_s, aes(x = rel_bp, color = Group, fill = Group)) +
  geom_hline(yintercept = 1, color = "gray60", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.35) +
  geom_line(aes(y = mean_smooth), linewidth = 0.95, na.rm = TRUE) +
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "TSS enrichment of fragment 5' start sites (muscle genes on chr21)",
    subtitle = paste0("Genes: ", genes_label, " (baseline normalized to window edges). Display: ", 5, "-bp running mean overlaid."),
    x = "Distance to TSS (bp)",
    y = "Normalized fragment density",
    color = "Group",
    fill = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")


p_tss_combined <- p_start_annot / p_tss_global / p_tss_muscle +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.4, 1, 1))

ggsave(
  file.path(FIG_DIR, "fig4_fragmentation_analysis.png"),
  p_tss_combined,
  width = 14, height = 18, dpi = 300
)

