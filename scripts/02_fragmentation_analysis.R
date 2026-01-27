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
#   chrom, start0, end0, name, mapq, strand, length_bp
# - start0/end0 are BED-style 0-based (end exclusive); we convert to 1-based for GRanges ops.
#
# Output:
# - results/tables/fragmentation_metrics.csv
# - results/figures/fig3_fragmentation_overview.png
# - results/figures/fig4_fragmentation_ratio_tracks.png
# - results/figures/fig5_tss_nucleosome_positioning.png

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(parallel))

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
FRAG_COLS <- c("chrom", "start0", "end0", "name", "mapq", "strand", "length_bp")
FRAG_COL_CLASSES <- c(
  chrom = "character",
  start0 = "integer",
  end0 = "integer",
  name = "character",
  mapq = "integer",
  strand = "character",
  length_bp = "integer"
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

# Inputs ----
stop_if_missing_fragments(metadata$sample_id, FRAG_DIR)
fragment_paths <- file.path(FRAG_DIR, paste0(metadata$sample_id, ".fragments.bed.gz"))
names(fragment_paths) <- metadata$sample_id

# 1) Genome-binned fragment size coverage + short/long ratio tracks ----
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

  tibble::tibble(
    bin_id = seq_along(bins_gr),
    short = short_counts,
    long = long_counts,
    short_long_ratio = (short_counts) / (long_counts + 1)
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
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.35) +
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
    subtitle = sprintf("Median (short 90–150) / (long 151–220 + 1) in %s bins", bin_label),
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
  file.path(FIG_DIR, "fig4_fragmentation_ratio_tracks.png"),
  p_tracks,
  width = 16, height = 6, dpi = 300
)

# Optional: Chromosome-style visualization with Gviz ----
if (!requireNamespace("Gviz", quietly = TRUE)) {
  BiocManager::install("Gviz")
}
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
    name = "ALS/Ctrl\nshort/long",
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
    file.path(FIG_DIR, "fig4_fragmentation_ratio_tracks_gviz.png"),
    width = 2400,
    height = if (!is.null(ideo_track)) 900 else 700,
    res = 300
  )
  Gviz::plotTracks(
    tracks,
    from = from_bp,
    to = to_bp,
    main = "Fragmentation ratio track (ALS/Ctrl)",
    background.title = "white",
    col.title = "black",
    cex.title = 0.20,
    background.panel = "white",
    col.axis = "black",
    col.frame = "gray85"
  )
  dev.off()

# =============================================================================
# 3) Nucleosome positioning proxy: TSS-centered midpoint density ----
# =============================================================================

message("Computing TSS-centered nucleosome positioning profiles...")

ensure_txdb <- function() {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE) ||
      !requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    stop(
      "Missing Bioconductor annotation packages required for TSS profiles.\n",
      "Please install (once) and re-run:\n",
      "  BiocManager::install(c('GenomicFeatures','TxDb.Hsapiens.UCSC.hg38.knownGene'))\n",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  TxDb.Hsapiens.UCSC.hg38.knownGene
}

txdb <- ensure_txdb()

# TSS points (strand-aware) and windows
tss_window <- 2000L
tss_points <- suppressWarnings(GenomicFeatures::promoters(txdb, upstream = 0L, downstream = 1L)) # 1-bp at TSS
tss_windows <- suppressWarnings(GenomicFeatures::promoters(txdb, upstream = tss_window, downstream = tss_window))

# Keep standard chroms only, trim out-of-bound windows, de-duplicate exact coordinates
tss_points <- GenomeInfoDb::keepStandardChromosomes(tss_points, pruning.mode = "coarse")
tss_windows <- GenomeInfoDb::keepStandardChromosomes(tss_windows, pruning.mode = "coarse")
tss_points <- GenomeInfoDb::keepSeqlevels(tss_points, value = intersect(seqlevels(tss_points), chroms), pruning.mode = "coarse")
tss_windows <- GenomeInfoDb::keepSeqlevels(tss_windows, value = intersect(seqlevels(tss_windows), chroms), pruning.mode = "coarse")
tss_windows <- GenomicRanges::trim(tss_windows)
tss_points <- unique(tss_points)
tss_windows <- unique(tss_windows)

# Optionally cap to keep runtime deterministic (profiles stabilize quickly at ~10–20k TSS)
max_tss <- 20000L
if (length(tss_windows) > max_tss) {
  set.seed(1)
  idx <- sample(seq_along(tss_windows), max_tss)
  tss_windows <- tss_windows[idx]
  tss_points <- tss_points[idx]
}

tss_site <- start(tss_points) # width 1 -> start == end
tss_strand <- as.character(strand(tss_points))

# Profile accumulator: counts by relative position [-tss_window, +tss_window]
rel_grid <- (-tss_window):tss_window
nbins_rel <- length(rel_grid)

profile_one_sample <- function(path, tss_windows, tss_site, tss_strand, len_range = c(30L, 500L)) {
  prof_all <- numeric(nbins_rel)
  prof_short <- numeric(nbins_rel)
  prof_long <- numeric(nbins_rel)

  n_all <- 0L
  n_short <- 0L
  n_long <- 0L

  cb <- function(x) {
    l <- x$length_bp
    keep <- !is.na(l) & l >= len_range[1] & l <= len_range[2]
    if (!any(keep)) return(invisible())
    x <- x[keep, , drop = FALSE]
    l <- x$length_bp

    mid_gr <- chunk_to_midpoint_gr(x)

    # All
    hits <- findOverlaps(mid_gr, tss_windows, ignore.strand = TRUE)
    if (length(hits) > 0) {
      q <- S4Vectors::queryHits(hits)
      s <- S4Vectors::subjectHits(hits)
      mid <- start(mid_gr[q])
      rel <- mid - tss_site[s]
      rel[tss_strand[s] == "-"] <- -rel[tss_strand[s] == "-"]
      rel <- rel[rel >= -tss_window & rel <= tss_window]
      if (length(rel) > 0) prof_all <<- prof_all + tabulate(rel + tss_window + 1L, nbins = nbins_rel)
      n_all <<- n_all + length(q)
    }

    # Short / long subsets (same overlap indices computed separately for correctness)
    idx_s <- which(l >= 90L & l <= 150L)
    if (length(idx_s) > 0) {
      hits_s <- findOverlaps(mid_gr[idx_s], tss_windows, ignore.strand = TRUE)
      if (length(hits_s) > 0) {
        q <- S4Vectors::queryHits(hits_s)
        s <- S4Vectors::subjectHits(hits_s)
        mid <- start(mid_gr[idx_s][q])
        rel <- mid - tss_site[s]
        rel[tss_strand[s] == "-"] <- -rel[tss_strand[s] == "-"]
        rel <- rel[rel >= -tss_window & rel <= tss_window]
        if (length(rel) > 0) prof_short <<- prof_short + tabulate(rel + tss_window + 1L, nbins = nbins_rel)
        n_short <<- n_short + length(q)
      }
    }

    idx_l <- which(l >= 151L & l <= 220L)
    if (length(idx_l) > 0) {
      hits_l <- findOverlaps(mid_gr[idx_l], tss_windows, ignore.strand = TRUE)
      if (length(hits_l) > 0) {
        q <- S4Vectors::queryHits(hits_l)
        s <- S4Vectors::subjectHits(hits_l)
        mid <- start(mid_gr[idx_l][q])
        rel <- mid - tss_site[s]
        rel[tss_strand[s] == "-"] <- -rel[tss_strand[s] == "-"]
        rel <- rel[rel >= -tss_window & rel <= tss_window]
        if (length(rel) > 0) prof_long <<- prof_long + tabulate(rel + tss_window + 1L, nbins = nbins_rel)
        n_long <<- n_long + length(q)
      }
    }

    invisible()
  }

  read_fragments_chunked(path = path, chunk_size = 200000L, FUN = cb)

  list(
    rel = rel_grid,
    all = prof_all,
    short = prof_short,
    long = prof_long,
    n_all = n_all,
    n_short = n_short,
    n_long = n_long
  )
}

profiles <- purrr::imap(fragment_paths, function(path, sample_id) {
  profile_one_sample(
    path = path,
    tss_windows = tss_windows,
    tss_site = tss_site,
    tss_strand = tss_strand
  ) %>%
    append(list(
      sample_id = sample_id,
      Group = metadata$Group[match(sample_id, metadata$sample_id)]
    ))
})

profile_df <- purrr::map_dfr(profiles, function(p) {
  # Normalize to "per million hits" to make profiles comparable
  tibble::tibble(
    sample_id = p$sample_id,
    Group = p$Group,
    rel_bp = p$rel,
    all = if (p$n_all > 0) (p$all / p$n_all) * 1e6 else NA_real_,
    short = if (p$n_short > 0) (p$short / p$n_short) * 1e6 else NA_real_,
    long = if (p$n_long > 0) (p$long / p$n_long) * 1e6 else NA_real_
  )
})

profile_summary <- profile_df %>%
  tidyr::pivot_longer(cols = c(all, short, long), names_to = "set", values_to = "cpm") %>%
  dplyr::group_by(Group, set, rel_bp) %>%
  dplyr::summarise(median_cpm = stats::median(cpm, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    set = factor(set, levels = c("all", "short", "long"), labels = c("All (30–500)", "Short (90–150)", "Long (151–220)"))
  )

p_tss <- ggplot(profile_summary, aes(x = rel_bp, y = median_cpm, color = Group)) +
  geom_vline(xintercept = 0, color = "gray35", linewidth = 0.35) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~set, ncol = 1, scales = "free_y") +
  scale_color_manual(values = COLORS) +
  labs(
    title = "Nucleosome positioning proxy around TSS",
    subtitle = sprintf(
      "Median midpoint density in \u00b1%d bp windows around TSS (normalized CPM); %d TSS used",
      tss_window, length(tss_windows)
    ),
    x = "Position relative to TSS (bp)",
    y = "Median CPM",
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

# Short/long enrichment as a single curve (per group)
ratio_summary <- profile_summary %>%
  dplyr::select(Group, set, rel_bp, median_cpm) %>%
  tidyr::pivot_wider(names_from = set, values_from = median_cpm) %>%
  dplyr::mutate(
    log2_short_over_long = log2((`Short (90–150)` + 1e-3) / (`Long (151–220)` + 1e-3))
  )

p_ratio <- ggplot(ratio_summary, aes(x = rel_bp, y = log2_short_over_long, color = Group)) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "gray35", linewidth = 0.35) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(values = COLORS) +
  labs(
    title = "Short vs long enrichment around TSS",
    subtitle = "log2(CPM_short / CPM_long) (pseudocount 1e-3)",
    x = "Position relative to TSS (bp)",
    y = "log2(short/long)",
    color = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

fig_tss <- p_tss / p_ratio + plot_annotation(tag_levels = "A")

ggsave(
  file.path(FIG_DIR, "fig5_tss_nucleosome_positioning.png"),
  fig_tss,
  width = 12, height = 12, dpi = 300
)

message("Fragmentation analysis complete.")
