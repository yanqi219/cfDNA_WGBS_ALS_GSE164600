# 01_qc_processing.R

# Setup ----
# setup_candidates <- c(
#   "scripts/00_setup.R",
#   "cfDNA_WGBS_ALS_GSE164600/scripts/00_setup.R"
# )
# setup_path <- setup_candidates[file.exists(setup_candidates)][1]
# if (is.na(setup_path) || !nzchar(setup_path)) {
#   stop(
#     "Cannot find scripts/00_setup.R. ",
#     "Please set your working directory to either:\n",
#     "  - cfDNA_WGBS_ALS_GSE164600/\n",
#     "  - Take_home_task/\n"
#   )
# }
# source(setup_path)

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

# Paths
PROJECT_DIR <- here::here()
DATA_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/data")
RAW_DIR <- file.path(DATA_DIR, "raw")
RESULTS_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/results")
FIG_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
FRAG_DIR <- file.path(DATA_DIR, "processed", "fragments")
# Create output dir
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(FIG_DIR)) {
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(TABLE_DIR)) {
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(FRAG_DIR)) {
  dir.create(FRAG_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Load sample metadata
metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"), show_col_types = FALSE)
metadata <- metadata %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Ctrl", "ALS")),
    sample_id = .data$Sample_id
  )

# Load genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Figure theme
# Enhanced publication-quality theme for Nature-style figures
theme_pub <- function(base_size = 14, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = rel(1.25),
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = rel(1.07),
        margin = margin(b = 6)
      ),
      axis.title.x = element_text(
        face = "bold",
        size = rel(1.08),
        margin = margin(t = 9)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = rel(1.08),
        margin = margin(r = 9)
      ),
      axis.text = element_text(
        color = "black",
        size = rel(1.00)
      ),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.line = element_line(color = "black", linewidth = 0.7),
      legend.title = element_text(
        face = "bold",
        size = rel(1.08)
      ),
      legend.text = element_text(
        size = rel(1.00)
      ),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      legend.justification = "center",
      legend.box.background = element_rect(color = NA, fill = NA),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "gray85", size = 0.30, linetype = "dotted"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(
        face = "bold",
        size = rel(1.05),
        margin = margin(t = 3, b = 3)
      ),
      plot.margin = margin(16, 16, 16, 16)
    )
}

# Color palette
COLORS <- c("ALS" = "#E64B35", "Ctrl" = "#4DBBD5")

# QC Functions ----
{
  bam_path <- file.path(RAW_DIR, metadata$Bam_file[1])
  sample_id <- metadata$sample_id[1]
  genome <- BSgenome.Hsapiens.UCSC.hg38
  chunk_size = 1e6
  frag_dir = FRAG_DIR
}

## flagstat ----
get_flagstat <- function(bam_path) {
  
  param_all <- ScanBamParam(what = "flag")
  bf <- BamFile(bam_path, yieldSize = 5e6)
  open(bf)
  
  total <- 0
  mapped <- 0
  paired <- 0
  proper_pair <- 0
  duplicates <- 0
  secondary <- 0
  supplementary <- 0
  
  repeat{
    flags <- scanBam(bf, param = param_all)[[1]]$flag
    if (length(flags) == 0) break
    
    total <- total + length(flags)
    mapped <- mapped + sum(bitwAnd(flags, 4) == 0)        # bit 4 = unmapped
    paired <- paired + sum(bitwAnd(flags, 1) > 0)         # bit 1 = paired
    proper_pair <- proper_pair + sum(bitwAnd(flags, 2) > 0)  # bit 2 = proper pair
    duplicates <- duplicates + sum(bitwAnd(flags, 1024) > 0) # bit 1024 = duplicate
    secondary <- secondary + sum(bitwAnd(flags, 256) > 0)    # bit 256 = secondary
    supplementary <- supplementary + sum(bitwAnd(flags, 2048) > 0) # bit 2048 = supplementary
  }
  close(bf)
  
  data.frame(
    total_reads = total,
    mapped = mapped,
    unmapped = total - mapped,
    paired = paired,
    properly_paired = proper_pair,
    duplicates = duplicates,
    secondary = secondary,
    supplementary = supplementary,
    mapping_rate = mapped / total,
    proper_pair_rate = proper_pair / paired,
    duplicate_rate = duplicates / total
  )
}

## Extract QC metrics ----
# bam_path Path to BAM file
# sample_id Sample identifier
# genome BSgenome object for sequence extraction
# chunk_size Number of reads per chunk (memory management)
extract_bam_metrics <- function(bam_path, sample_id, genome, chunk_size = 1e6, frag_dir = FRAG_DIR) {

  message(sprintf("Processing: %s", sample_id))
  message(sprintf("BAM file: %s", bam_path))

  if (!file.exists(bam_path)) {
    warning(sprintf("BAM file not found: %s", bam_path))
    return(NULL)
  }

  # Get flag statistics
  flagstat <- get_flagstat(bam_path)

  # BAM parameters
  param <- ScanBamParam(
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE,
                       isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                       isDuplicate = FALSE),
    what = c("qname", "flag", "rname", "pos", "mapq", "isize", "seq", "qwidth"),
    tag = c("XM")  # Bismark methylation tag
  )

  # Open file
  bf <- BamFile(bam_path, yieldSize = chunk_size)
  open(bf)

  frag_lengths <- integer()
  motif_5p <- character()
  meth_calls <- list()
  total_reads <- 0
  mapq_values <- integer()
  
  # Stream fragments to BED (gzipped) for downstream fragmentation analyses
  # BED columns: chrom, start(0-based), end(0-based, exclusive), name, mapq, strand, length
  fragment_bed_path <- file.path(frag_dir, paste0(sample_id, ".fragments.bed.gz"))
  bed_con <- gzfile(fragment_bed_path, open = "wt")
  on.exit(try(close(bed_con), silent = TRUE), add = TRUE)

  # Process in chunks
  repeat {
    chunk_raw <- scanBam(bf, param = param)[[1]]

    if (length(chunk_raw$qname) == 0) break

    # mapq scores, filter if < 30
    mapq_values <- c(mapq_values, chunk_raw$mapq)
    keep_idx <- which(chunk_raw$mapq >= 30)
    chunk <- lapply(chunk_raw, function(x) {
      if (length(x) == length(chunk_raw$mapq)) {
        return(x[keep_idx])
      }
      if (is.list(x)) {
        return(lapply(x, function(tag_vector) {
          return(tag_vector[keep_idx])
        }))
      }
      return(x)
    })

    # Total reads after filtering
    total_reads <- total_reads + length(chunk$qname)

    # Fragment lengths, Flag 64 = first in pair
    is_first <- bitwAnd(chunk$flag, 64) > 0
    isize <- abs(chunk$isize[is_first])
    # isize <- isize[isize > 0 & isize < 1000]
    frag_lengths <- c(frag_lengths, isize)

    # Generate fragment BED entries
    idx_frag <- which(
      is_first &
        !is.na(chunk$isize) & chunk$isize != 0 &
        !is.na(chunk$pos) &
        !is.na(chunk$rname)
    )
    if (length(idx_frag) > 0) {
      tlen <- chunk$isize[idx_frag]
      pos <- chunk$pos[idx_frag]
      qwidth <- chunk$qwidth[idx_frag]

      frag_start <- ifelse(tlen > 0, pos, pos + qwidth + tlen)
      frag_len <- abs(tlen)
      frag_end <- frag_start + frag_len - 1

      keep_frag <- is.finite(frag_start) & is.finite(frag_end) & frag_len > 0 & frag_start > 0 & frag_end >= frag_start
      if (any(keep_frag)) {
        chrom <- as.character(chunk$rname[idx_frag][keep_frag])
        bed_start0 <- as.integer(frag_start[keep_frag] - 1)
        bed_end0 <- as.integer(frag_end[keep_frag])
        name <- as.character(chunk$qname[idx_frag][keep_frag])
        mapq <- as.integer(chunk$mapq[idx_frag][keep_frag])
        strand <- ifelse(tlen[keep_frag] > 0, "+", "-")
        length_bp <- as.integer(frag_len[keep_frag])

        bed_lines <- paste(chrom, bed_start0, bed_end0, name, mapq, strand, length_bp, sep = "\t")
        writeLines(bed_lines, con = bed_con, useBytes = TRUE)
      }
    }

    # End motifs
    if (length(chunk$rname) > 0) {
      idx <- idx_frag

      if (length(idx) > 0) {
        tlen <- chunk$isize[idx]
        pos <- chunk$pos[idx]
        qwidth <- chunk$qwidth[idx]

        frag_start <- ifelse(tlen > 0, pos, pos + qwidth + tlen)
        frag_width <- abs(tlen)
        frag_end <- frag_start + frag_width - 1

        gr_left <- GRanges(
          seqnames = chunk$rname[idx],
          ranges = IRanges(start = frag_start, width = 4),
          strand = "+"
        )

        gr_right <- GRanges(
          seqnames = chunk$rname[idx],
          ranges = IRanges(start = frag_end - 3, width = 4),
          strand = "-"
        )

        # Validate coordinates
        valid_left <- start(gr_left) > 0
        valid_right <- start(gr_right) > 0

        # get sequences
        if (any(valid_left)) {
          seqs_left <- tryCatch(
            getSeq(genome, gr_left[valid_left]),
            error = function(e) DNAStringSet()
          )
          if (length(seqs_left) > 0) {
            motif_5p <- c(motif_5p, as.character(seqs_left))
          }
        }

        if (any(valid_right)) {
          seqs_right <- tryCatch(
            getSeq(genome, gr_right[valid_right]),
            error = function(e) DNAStringSet()
          )
          if (length(seqs_right) > 0) {
            motif_5p <- c(motif_5p, as.character(seqs_right))
          }
        }
      }
    }

    # methylation calls
    if (!is.null(chunk$tag$XM)) {
      xm_tags <- chunk$tag$XM[!is.na(chunk$tag$XM)]
      if (length(xm_tags) > 0) {
        meth_calls[[length(meth_calls) + 1]] <- parse_xm_tags(xm_tags)
      }
    }

    # Memory cleanup
    rm(chunk)
    gc(verbose = FALSE)
  }

  close(bf)

  # Aggregate methylation
  meth_summary <- if (length(meth_calls) > 0) {
    do.call(rbind, meth_calls) %>%
      summarise(across(everything(), sum))
  } else {
    data.frame(Z = 0, z = 0, X = 0, x = 0, H = 0, h = 0)
  }

  list(
    sample_id = sample_id,
    fragment_bed = fragment_bed_path,
    flagstat = flagstat,
    filtered_reads = total_reads,  # Reads passing quality filters
    fragment_lengths = frag_lengths,
    motifs_5p = motif_5p,
    mapq = mapq_values,
    methylation = meth_summary
  )
}

## Parse Bismark XM methylation tags ----
parse_xm_tags <- function(xm_tags) {
  all_chars <- paste(xm_tags, collapse = "")

  data.frame(
    Z = str_count(all_chars, "Z"),  # methylated CpG
    z = str_count(all_chars, "z"),  # unmethylated CpG
    X = str_count(all_chars, "X"),  # methylated CHG
    x = str_count(all_chars, "x"),  # unmethylated CHG
    H = str_count(all_chars, "H"),  # methylated CHH
    h = str_count(all_chars, "h")   # unmethylated CHH
  )
}

## Calculate summary statistics ----
# metrics List from extract_bam_metrics
# group Sample group
calculate_summary_stats <- function(metrics) {

  if (is.null(metrics)) {
    return(data.frame(sample_id = NA))
  }

  fl <- metrics$fragment_lengths
  meth <- metrics$methylation
  fs <- metrics$flagstat

  # CpG methylation rate
  cpg_total <- meth$Z + meth$z
  cpg_meth_rate <- if (cpg_total > 0) meth$Z / cpg_total else NA

  # CHH methylation rate - bisulfite conversion efficiency
  # Should be very low (<1%) for complete conversion
  chh_total <- meth$H + meth$h
  chh_meth_rate <- if (chh_total > 0) meth$H / chh_total else NA
  bisulfite_conv <- if (!is.na(chh_meth_rate)) 1 - chh_meth_rate else NA

  data.frame(
    sample_id = metrics$sample_id,
    # Flagstat metrics (raw BAM statistics)
    total_reads = fs$total_reads,
    mapped_reads = fs$mapped,
    unmapped_reads = fs$unmapped,
    properly_paired = fs$properly_paired,
    duplicates = fs$duplicates,
    mapping_rate = fs$mapping_rate,
    proper_pair_rate = fs$proper_pair_rate,
    duplicate_rate = fs$duplicate_rate,
    # Filtered reads (passing QC)
    filtered_reads = metrics$filtered_reads,
    filter_pass_rate = metrics$filtered_reads / fs$total_reads,
    # Fragment length statistics
    n_fragments = length(fl),
    mean_frag_length = mean(fl),
    median_frag_length = median(fl),
    sd_frag_length = sd(fl),
    frag_mode = as.numeric(names(sort(table(fl), decreasing = TRUE))[1]),
    short_frag_ratio = sum(fl < 150) / length(fl),  # <150bp
    mono_nuc_ratio = sum(fl >= 150 & fl <= 220) / length(fl),  # mono-nucleosome
    di_nuc_ratio = sum(fl >= 300 & fl <= 400) / length(fl),  # di-nucleosome
    mean_mapq = mean(metrics$mapq),
    n_motifs = length(metrics$motifs_5p),
    cpg_methylation = cpg_meth_rate,
    chh_methylation = chh_meth_rate,
    bisulfite_conversion = bisulfite_conv
  )
}

# Process all samples ----
all_metrics <- lapply(seq_len(nrow(metadata)), function(i) {
  row <- metadata[i, ]
  bam_path <- file.path(RAW_DIR, row$Bam_file)

  metrics <- extract_bam_metrics(
    bam_path = bam_path,
    sample_id = row$sample_id,
    genome = genome,
    frag_dir = FRAG_DIR
  )

  return(metrics)
})

# Filter out failed samples
all_metrics <- all_metrics[!sapply(all_metrics, is.null)]

if (length(all_metrics) == 0) {
  stop("No BAM files found. Please check data/raw/ directory.")
}

# Summary statistics table ----
summary_stats <- map_dfr(all_metrics, calculate_summary_stats)
summary_stats <- summary_stats %>%
  dplyr::left_join(metadata, by = "sample_id")

# Save summary table
readr::write_csv(summary_stats, file.path(TABLE_DIR, "qc_summary_stats.csv"))

# Figure 1: basic QCs ----
key_metrics <- summary_stats %>%
  dplyr::select(
    sample_id,
    Group,
    total_reads,
    filtered_reads,
    median_frag_length,
    mean_mapq,
    cpg_methylation,
    bisulfite_conversion
  ) %>%
  dplyr::mutate(
    total_reads_m = total_reads / 1e6,
    filtered_reads_m = filtered_reads / 1e6,
    bisulfite_conversion_pct = 100 * bisulfite_conversion,
    cpg_methylation_pct = 100 * cpg_methylation
  )

key_metrics_by_group <- key_metrics %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    total_reads_m_mean = mean(total_reads_m, na.rm = TRUE),
    filtered_reads_m_mean = mean(filtered_reads_m, na.rm = TRUE),
    median_frag_length_mean = mean(median_frag_length, na.rm = TRUE),
    mean_mapq_mean = mean(mean_mapq, na.rm = TRUE),
    cpg_methylation_pct_mean = mean(cpg_methylation_pct, na.rm = TRUE),
    bisulfite_conversion_pct_mean = mean(bisulfite_conversion_pct, na.rm = TRUE),
    .groups = "drop"
  )

message("\n================ QC key metrics (per sample) ================\n")
print(
  key_metrics %>%
    dplyr::arrange(Group, dplyr::desc(filtered_reads)) %>%
    dplyr::select(sample_id, Group, total_reads, filtered_reads, median_frag_length, mean_mapq, cpg_methylation, bisulfite_conversion)
)

message("\n================ QC key metrics (by group) ================\n")
print(key_metrics_by_group)

## plot ----
key_metrics_long_sample <- key_metrics %>%
  dplyr::select(sample_id, Group, total_reads_m, filtered_reads_m, median_frag_length, mean_mapq, cpg_methylation_pct, bisulfite_conversion_pct) %>%
  tidyr::pivot_longer(
    cols = c(total_reads_m, filtered_reads_m, median_frag_length, mean_mapq, cpg_methylation_pct, bisulfite_conversion_pct),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    sample_id = factor(
      sample_id,
      levels = key_metrics %>%
        dplyr::arrange(Group, sample_id) %>%
        dplyr::pull(sample_id)
    )
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("total_reads_m", "filtered_reads_m", "median_frag_length", "mean_mapq", "cpg_methylation_pct", "bisulfite_conversion_pct"),
      labels = c("Total reads (M)", "Filtered reads (M)", "Median fragment length (bp)", "Mean MAPQ", "CpG methylation (%)", "Bisulfite conversion (%)")
    )
  )

p_qc_sample <- ggplot(
  key_metrics_long_sample,
  aes(x = value, y = sample_id, fill = Group)
) +
  geom_col(width = 0.8, alpha = 0.95) +
  facet_wrap(~metric, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "QCs by sample",
    x = NULL,
    y = "Sample",
    fill = "Group"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 7)
  )

p_qc_group <- ggplot(
  key_metrics_long_sample,
  aes(x = Group, y = value, fill = Group)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.95) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "QCs by group",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

p_qc_combined <- p_qc_group / p_qc_sample +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.4))

ggsave(file.path(FIG_DIR, "fig1_qc_key_metrics.png"),
       p_qc_combined, width = 14, height = 16, dpi = 300)

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

# Group-level median fragment size distributions ----
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


# Figure 3: End Motif Analysis ----
## Aggregate motif counts ----
motif_df <- map_dfr(all_metrics, function(m) {
  if (length(m$motifs_5p) == 0) return(NULL)
  
  # Count motifs
  motif_counts <- table(m$motifs_5p)
  
  data.frame(
    sample_id = m$sample_id,
    motif = names(motif_counts),
    count = as.integer(motif_counts)
  )
})
motif_df <- motif_df %>%
  dplyr::left_join(metadata, by = "sample_id")

# Calculate frequencies per sample
motif_freq <- motif_df %>%
  group_by(sample_id, Group) %>%
  mutate(frequency = count / sum(count)) %>%
  ungroup()

# Average frequency by group
motif_summary <- motif_freq %>%
  group_by(group, motif) %>%
  summarise(
    mean_freq = mean(frequency),
    sd_freq = sd(frequency),
    .groups = "drop"
  )

# Top 20 motifs
top_motifs <- motif_summary %>%
  group_by(motif) %>%
  summarise(total_freq = sum(mean_freq)) %>%
  arrange(desc(total_freq)) %>%
  head(20) %>%
  pull(motif)

# 2A: Bar plot of top motifs
p2a <- motif_summary %>%
  filter(motif %in% top_motifs) %>%
  ggplot(aes(x = reorder(motif, mean_freq), y = mean_freq, fill = group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = COLORS) +
  coord_flip() +
  labs(
    title = "Top 20 Fragment End Motifs (5')",
    x = "4-mer Motif",
    y = "Mean Frequency",
    fill = "Group"
  ) +
  theme_pub()

# 2B: Motif frequency comparison (ALS vs Ctrl)
motif_wide <- motif_summary %>%
  select(group, motif, mean_freq) %>%
  pivot_wider(names_from = group, values_from = mean_freq, values_fill = 0)

if ("ALS" %in% names(motif_wide) && "Ctrl" %in% names(motif_wide)) {
  motif_wide <- motif_wide %>%
    mutate(
      log2fc = log2((ALS + 1e-6) / (Ctrl + 1e-6)),
      significant = abs(log2fc) > 0.5
    )
  
  p2b <- ggplot(motif_wide, aes(x = Ctrl, y = ALS)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    geom_text(data = filter(motif_wide, significant), 
              aes(label = motif), size = 2.5, vjust = -0.5, check_overlap = TRUE) +
    scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray50")) +
    labs(
      title = "Motif Frequency: ALS vs Control",
      x = "Control Frequency",
      y = "ALS Frequency"
    ) +
    theme_pub() +
    theme(legend.position = "none")
} else {
  p2b <- ggplot() + theme_void() + 
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data")
}

# 2C: Heatmap of motif frequencies by sample
motif_matrix <- motif_freq %>%
  filter(motif %in% top_motifs) %>%
  select(sample_id, motif, frequency) %>%
  pivot_wider(names_from = motif, values_from = frequency, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# Simple heatmap with ggplot
motif_heatmap_df <- motif_freq %>%
  filter(motif %in% top_motifs) %>%
  mutate(motif = factor(motif, levels = top_motifs))

p2c <- ggplot(motif_heatmap_df, aes(x = motif, y = sample_id, fill = frequency)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "End Motif Frequencies",
    x = "Motif",
    y = "Sample",
    fill = "Frequency"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7)
  )

# Combine motif plots
p2_combined <- (p2a | p2b) / p2c +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.2))

ggsave(file.path(FIG_DIR, "fig2_end_motifs.pdf"), 
       p2_combined, width = 12, height = 10, dpi = 300)
ggsave(file.path(FIG_DIR, "fig2_end_motifs.png"), 
       p2_combined, width = 12, height = 10, dpi = 300)

# =============================================================================
# Figure 3: Methylation QC
# =============================================================================

# 3A: CpG methylation distribution
p3a <- ggplot(summary_stats, aes(x = group, y = cpg_methylation, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "CpG Methylation Level",
    x = "Group",
    y = "Methylation Rate"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 3B: Bisulfite conversion efficiency
p3b <- ggplot(summary_stats, aes(x = group, y = bisulfite_conversion * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 99, linetype = "dashed", color = "darkgreen") +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Bisulfite Conversion Efficiency",
    subtitle = "Estimated from CHH methylation (should be >99%)",
    x = "Group",
    y = "Conversion Rate (%)"
  ) +
  ylim(95, 100) +
  theme_pub() +
  theme(legend.position = "none")

# 3C: CpG vs CHH methylation (quality check)
p3c <- ggplot(summary_stats, aes(x = chh_methylation * 100, y = cpg_methylation * 100, 
                                  color = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_id), size = 2, vjust = -0.8, check_overlap = TRUE) +
  scale_color_manual(values = COLORS) +
  labs(
    title = "CpG vs CHH Methylation",
    subtitle = "Low CHH indicates good bisulfite conversion",
    x = "CHH Methylation (%)",
    y = "CpG Methylation (%)",
    color = "Group"
  ) +
  theme_pub()

# Combine methylation plots
p3_combined <- (p3a | p3b | p3c) +
  plot_annotation(
    title = "Methylation Quality Control",
    tag_levels = "A"
  )

ggsave(file.path(FIG_DIR, "fig3_methylation_qc.pdf"), 
       p3_combined, width = 14, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "fig3_methylation_qc.png"), 
       p3_combined, width = 14, height = 5, dpi = 300)

# =============================================================================
# Figure 4: Mapping & Alignment Statistics (Flagstat)
# =============================================================================

# 4A: Total reads per sample
p4a <- ggplot(summary_stats, aes(x = reorder(sample_id, total_reads), 
                                  y = total_reads / 1e6, fill = group)) +
  geom_col() +
  scale_fill_manual(values = COLORS) +
  coord_flip() +
  labs(
    title = "Total Reads",
    x = "Sample",
    y = "Reads (millions)",
    fill = "Group"
  ) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

# 4B: Mapping rate
p4b <- ggplot(summary_stats, aes(x = group, y = mapping_rate * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Mapping Rate",
    subtitle = "Should be >90%",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 4C: Properly paired rate
p4c <- ggplot(summary_stats, aes(x = group, y = proper_pair_rate * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Properly Paired Rate",
    subtitle = "Both mates mapped correctly",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 4D: Duplicate rate
p4d <- ggplot(summary_stats, aes(x = group, y = duplicate_rate * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "orange") +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Duplicate Rate",
    subtitle = "PCR duplicates (should be <20%)",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 4E: Filter pass rate
p4e <- ggplot(summary_stats, aes(x = group, y = filter_pass_rate * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "QC Filter Pass Rate",
    subtitle = "Reads passing all filters",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 4F: MAPQ distribution
p4f <- ggplot(summary_stats, aes(x = group, y = mean_mapq, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Mean MAPQ",
    x = "Group",
    y = "MAPQ Score"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# Combine mapping stats
p4_combined <- (p4a | (p4b | p4c) / (p4d | p4e) / p4f) +
  plot_annotation(
    title = "Mapping & Alignment Quality (Flagstat Metrics)",
    tag_levels = "A"
  ) +
  plot_layout(widths = c(1.2, 1))

ggsave(file.path(FIG_DIR, "fig4_mapping_stats.pdf"), 
       p4_combined, width = 14, height = 10, dpi = 300)
ggsave(file.path(FIG_DIR, "fig4_mapping_stats.png"), 
       p4_combined, width = 14, height = 10, dpi = 300)

# =============================================================================
# Figure 5: cfDNA-Specific Quality Metrics
# =============================================================================

# 5A: Short fragment ratio (important for cfDNA quality)
p5a <- ggplot(summary_stats, aes(x = group, y = short_frag_ratio * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Short Fragment Ratio",
    subtitle = "Fragments < 150bp",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 5B: Mono-nucleosome ratio
p5b <- ggplot(summary_stats, aes(x = group, y = mono_nuc_ratio * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Mono-nucleosome Ratio",
    subtitle = "Fragments 140-200bp (expected peak)",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# 5C: Di-nucleosome ratio
p5c <- ggplot(summary_stats, aes(x = group, y = di_nuc_ratio * 100, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Di-nucleosome Ratio",
    subtitle = "Fragments 300-400bp",
    x = "Group",
    y = "Percentage (%)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# Combine cfDNA metrics
p5_combined <- (p5a | p5b | p5c) +
  plot_annotation(
    title = "cfDNA-Specific Quality Metrics",
    tag_levels = "A"
  )

ggsave(file.path(FIG_DIR, "fig5_cfDNA_metrics.pdf"), 
       p5_combined, width = 12, height = 4, dpi = 300)
ggsave(file.path(FIG_DIR, "fig5_cfDNA_metrics.png"), 
       p5_combined, width = 12, height = 4, dpi = 300)

# =============================================================================
# Figure 6: Statistical Comparisons (ALS vs Control)
# =============================================================================

# Perform t-tests for key metrics
run_ttest <- function(data, metric) {
  if (metric %in% names(data)) {
    vals <- data[[metric]]
    groups <- data$group
    
    if (any(!is.na(vals))) {
      tt <- tryCatch(
        t.test(vals[groups == "ALS"], vals[groups == "Ctrl"]),
        error = function(e) NULL
      )
      if (!is.null(tt)) {
        return(data.frame(
          metric = metric,
          als_mean = mean(vals[groups == "ALS"], na.rm = TRUE),
          ctrl_mean = mean(vals[groups == "Ctrl"], na.rm = TRUE),
          p_value = tt$p.value,
          significant = tt$p.value < 0.05
        ))
      }
    }
  }
  NULL
}

metrics_to_test <- c(
  # Mapping metrics
  "mapping_rate", "proper_pair_rate", "duplicate_rate",
  # Fragment metrics
  "mean_frag_length", "median_frag_length", "short_frag_ratio", "mono_nuc_ratio",
  # Methylation metrics
  "cpg_methylation", "bisulfite_conversion"
)

stat_results <- map_dfr(metrics_to_test, ~run_ttest(summary_stats, .x))

if (nrow(stat_results) > 0) {
  stat_results <- stat_results %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  readr::write_csv(stat_results, file.path(TABLE_DIR, "statistical_comparisons.csv"))
  
  # Plot statistical summary
  p6 <- stat_results %>%
    mutate(metric = factor(metric, levels = rev(metrics_to_test))) %>%
    ggplot(aes(x = metric, y = -log10(p_adj))) +
    geom_col(aes(fill = significant), width = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(aes(label = label), hjust = -0.2, size = 4) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray70")) +
    coord_flip() +
    labs(
      title = "Statistical Comparison: ALS vs Control",
      x = "Metric",
      y = "-log10(adjusted p-value)"
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  ggsave(file.path(FIG_DIR, "fig6_statistical_comparison.pdf"), 
         p6, width = 8, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, "fig6_statistical_comparison.png"), 
         p6, width = 8, height = 6, dpi = 300)
}

# =============================================================================
# Save Combined Metrics for Downstream Analysis
# =============================================================================

# Save fragment lengths for downstream use
frag_summary <- frag_df %>%
  group_by(sample_id, group) %>%
  summarise(
    fragment_lengths = list(fragment_length),
    .groups = "drop"
  )

saveRDS(frag_summary, file.path(DATA_DIR, "processed", "fragment_lengths.rds"))

# Save motif frequencies
saveRDS(motif_freq, file.path(DATA_DIR, "processed", "motif_frequencies.rds"))

# =============================================================================
# Final Report
# =============================================================================

message("\n")
message(strrep("=", 60))
message("QC Processing Complete")
message(strrep("=", 60))
message(sprintf("Samples processed: %d", length(all_metrics)))
message(sprintf("Figures saved to: %s", FIG_DIR))
message(sprintf("Tables saved to: %s", TABLE_DIR))
message("\nFigures generated:")
message("  - fig1_fragment_length.pdf/png    (Fragment size distribution)")
message("  - fig2_end_motifs.pdf/png         (End motif analysis)")
message("  - fig3_methylation_qc.pdf/png     (Methylation & bisulfite QC)")
message("  - fig4_mapping_stats.pdf/png      (Flagstat: mapping, dup rate)")
message("  - fig5_cfDNA_metrics.pdf/png      (cfDNA-specific metrics)")
message("  - fig6_statistical_comparison.pdf/png (ALS vs Ctrl stats)")
message("\nTables generated:")
message("  - qc_summary_stats.csv            (All QC metrics per sample)")
message("  - statistical_comparisons.csv    (t-test results)")
message(strrep("=", 60))
