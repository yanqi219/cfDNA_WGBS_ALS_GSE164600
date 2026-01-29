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
  "ggplot2", "dplyr", "tidyr", "readr", "purrr", "stringr", "tibble",
  "caret", "pROC", "patchwork", "here", "parallel", "ggpubr"
)
bioc_pkgs <- c(
  "Rsamtools", "cigarillo", "GenomicRanges", "IRanges",
  "Biostrings", "BSgenome.Hsapiens.UCSC.hg38"
)
invisible(lapply(cran_pkgs, install_and_load, bioc = FALSE))
invisible(lapply(bioc_pkgs, install_and_load, bioc = TRUE))

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
    Group = factor(Group, levels = c("Ctrl", "ALS"))
  )

# Load genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Figure theme
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
{
  bam_path <- file.path(RAW_DIR, metadata$Bam_file[1])
  sample_id = metadata$sample_id[1]
  frag_dir = FRAG_DIR
  chunk_size = 1e6
}
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
                       isDuplicate = FALSE,
                       isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE,
                       isNotPassingQualityControls = FALSE),
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "isize", "seq", "qwidth"),
    tag = c("XM")  # Bismark methylation tag
  )

  # Open file
  bf <- BamFile(bam_path, yieldSize = chunk_size)
  open(bf)

  frag_lengths <- integer()
  meth_calls <- list()
  total_reads <- 0
  mapq_values <- integer()
  gc_values <- numeric()

  # Output files: fragment BED.
  # Fragment BED consumed by downstream scripts in this repo:
  # chrom, start0, end0, name, mapq, strand, length_bp, gc
  # where start0 is 0-based, end0 is compatible with 1-based inclusive coordinate for BSgenome
  # conversion via: start1 = start0 + 1; end1 = end0; and length_bp == end0 - start0.
  fragment_bed_path <- file.path(frag_dir, paste0(sample_id, ".fragments.bed.gz"))
  frag_con <- gzfile(fragment_bed_path, open = "wt")
  on.exit(try(close(frag_con), silent = TRUE), add = TRUE)

  # Carry-over buffer for mate pairing across chunk boundaries
  pending_aln <- tibble::tibble(
    qname = character(),
    flag = integer(),
    rname = character(),
    pos = integer(),
    cigar = character(),
    mapq = integer()
  )

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

    # BEDPE -> fragment BED (bedtools-style pairing; robust to chunk boundaries) ----
    # We pair alignments by QNAME using both mates' alignment records.
    # Fragment coordinates come from the union span of the two mate alignments.
    if (length(chunk$qname) > 0) {
      aln <- tibble::tibble(
        qname = as.character(chunk$qname),
        flag = as.integer(chunk$flag),
        rname = as.character(chunk$rname),
        pos = as.integer(chunk$pos),
        cigar = as.character(chunk$cigar),
        mapq = as.integer(chunk$mapq)
      ) %>%
        dplyr::filter(
          !is.na(qname),
          !is.na(rname),
          !is.na(pos),
          !is.na(cigar),
          is.finite(pos),
          pos > 0L
        )

      if (nrow(aln) > 0) {
        aln <- aln %>%
          dplyr::mutate(
            start0 = pos - 1L,
            end0_excl = start0 + as.integer(cigarillo::cigar_extent_along_ref(cigar)),
            strand = dplyr::if_else(bitwAnd(flag, 16L) > 0L, "-", "+")
          ) %>%
          dplyr::filter(is.finite(start0), is.finite(end0_excl), end0_excl > start0)

        aln_all <- dplyr::bind_rows(pending_aln, aln)
        tab <- table(aln_all$qname)
        have_pair <- names(tab[tab >= 2L])

        if (length(have_pair) > 0) {
          pair_rows <- aln_all %>%
            dplyr::filter(qname %in% have_pair) %>%
            dplyr::group_by(qname) %>%
            dplyr::slice(1:2) %>%
            dplyr::ungroup()

          pending_aln <- aln_all %>%
            dplyr::filter(!(qname %in% have_pair))

          # Order mates to mimic bedtools BEDPE ordering:
          # - if same chrom: smaller start0 first
          # - else: lexicographically lower chrom first
          # Fast mate ordering + wide reshape (assumes exactly 2 rows per qname).
          # This replaces group_by()/mutate()/pivot_wider(), which are slow in tight loops.
          pair_rows <- pair_rows[order(pair_rows$qname), , drop = FALSE]
          mate_raw <- stats::ave(pair_rows$qname, pair_rows$qname, FUN = seq_along) # 1,2 within each qname

          idx1 <- which(mate_raw == 1L)
          idx2 <- which(mate_raw == 2L)

          # Determine swap per pair (bedtools-style):
          # - same chrom: smaller start first
          # - else: lexicographically lower chrom first
          same_chr <- pair_rows$rname[idx1] == pair_rows$rname[idx2]
          swap_per_pair <- ifelse(
            same_chr,
            pair_rows$start0[idx1] > pair_rows$start0[idx2],
            pair_rows$rname[idx1] > pair_rows$rname[idx2]
          )

          i1 <- ifelse(swap_per_pair, idx2, idx1)
          i2 <- ifelse(swap_per_pair, idx1, idx2)

          pair_w <- tibble::tibble(
            qname = pair_rows$qname[idx1],
            rname_1 = pair_rows$rname[i1],
            start0_1 = pair_rows$start0[i1],
            end0_excl_1 = pair_rows$end0_excl[i1],
            mapq_1 = pair_rows$mapq[i1],
            strand_1 = pair_rows$strand[i1],
            rname_2 = pair_rows$rname[i2],
            start0_2 = pair_rows$start0[i2],
            end0_excl_2 = pair_rows$end0_excl[i2],
            mapq_2 = pair_rows$mapq[i2],
            strand_2 = pair_rows$strand[i2]
          ) %>%
            dplyr::filter(!is.na(.data$rname_1) & !is.na(.data$rname_2))

          pair_w <- pair_w %>% dplyr::arrange(rname_1, start0_1)

          if (nrow(pair_w) > 0) {
            score <- pmin(pair_w$mapq_1, pair_w$mapq_2)

            # Derive fragment BED intervals
            same_chr <- pair_w$rname_1 == pair_w$rname_2
            if (any(same_chr)) {
              # Protocol-faithful "BEDPE -> fragment" definition:
              # fragment start = chrom1/start1 (BEDPE col 2), fragment end = chrom2/end2 (BEDPE col 6)
              # (this matches the common awk pattern: frag_size = $6 - $2)
              frag_start0 <- pair_w$start0_1[same_chr]
              frag_end0 <- pair_w$end0_excl_2[same_chr]
              length_bp <- frag_end0 - frag_start0

              # Keep "reasonable" fragment sizes for downstream plots/metrics
              keep_len <- is.finite(length_bp) & length_bp > 0L & length_bp <= 1000L
              if (any(keep_len)) {
                # GC content (bedtools nuc convention): (G+C) / (A+C+G+T), excluding Ns.
                gc <- rep(NA_real_, sum(keep_len))
                gr_frag <- GenomicRanges::GRanges(
                  seqnames = pair_w$rname_1[same_chr][keep_len],
                  ranges = IRanges::IRanges(
                    start = as.integer(frag_start0[keep_len] + 1L),
                    end = as.integer(frag_end0[keep_len])
                  ),
                  strand = "*"
                )

                gc <- tryCatch({
                  seqs <- BSgenome::getSeq(genome, gr_frag)
                  acgt <- Biostrings::letterFrequency(seqs, letters = c("A", "C", "G", "T"), as.prob = FALSE)
                  denom <- rowSums(acgt)
                  num_gc <- acgt[, "C"] + acgt[, "G"]
                  ifelse(denom > 0, num_gc / denom, NA_real_)
                }, error = function(e) {
                  rep(NA_real_, length(gr_frag))
                })

                frag_lines <- paste(
                  pair_w$rname_1[same_chr][keep_len],
                  as.integer(frag_start0[keep_len]),
                  as.integer(frag_end0[keep_len]),
                  pair_w$qname[same_chr][keep_len],
                  as.integer(score[same_chr][keep_len]),
                  pair_w$strand_1[same_chr][keep_len],
                  as.integer(length_bp[keep_len]),
                  formatC(gc, format = "f", digits = 4),
                  sep = "\t"
                )
                writeLines(frag_lines, con = frag_con, useBytes = TRUE)
                frag_lengths <- c(frag_lengths, as.integer(length_bp[keep_len]))
                gc_values <- c(gc_values, gc)
              }
            }
          }
        } else {
          pending_aln <- aln_all
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
    gc_content = gc_values,
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
  gc <- metrics$gc_content
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
    gc_content = mean(gc, na.rm = TRUE),
    n_fragments = length(fl),
    mean_frag_length = mean(fl),
    median_frag_length = median(fl),
    sd_frag_length = sd(fl),
    frag_mode = as.numeric(names(sort(table(fl), decreasing = TRUE))[1]),
    short_frag_ratio = sum(fl < 150) / length(fl),  # <150bp
    mono_nuc_ratio = sum(fl >= 150 & fl <= 220) / length(fl),  # mono-nucleosome
    di_nuc_ratio = sum(fl >= 300 & fl <= 400) / length(fl),  # di-nucleosome
    mean_mapq = mean(metrics$mapq),
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

saveRDS(all_metrics, file.path(DATA_DIR, "processed", "all_metrics.rds"))

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
    bisulfite_conversion,
    gc_content
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
    gc_content_mean = mean(gc_content, na.rm = TRUE),
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
    dplyr::select(sample_id, Group, total_reads, filtered_reads, median_frag_length, mean_mapq, cpg_methylation, bisulfite_conversion, gc_content)
)

message("\n================ QC key metrics (by group) ================\n")
print(key_metrics_by_group)

## plot ----
key_metrics_long_sample <- key_metrics %>%
  dplyr::select(sample_id, Group, total_reads_m, filtered_reads_m, median_frag_length, mean_mapq, cpg_methylation_pct, bisulfite_conversion_pct, gc_content) %>%
  tidyr::pivot_longer(
    cols = c(total_reads_m, filtered_reads_m, median_frag_length, mean_mapq, cpg_methylation_pct, bisulfite_conversion_pct, gc_content),
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
      levels = c("total_reads_m", "filtered_reads_m", "mean_mapq", "median_frag_length", "gc_content", "cpg_methylation_pct", "bisulfite_conversion_pct"),
      labels = c("Total reads (M)", "Filtered reads (M)", "Mean MAPQ", "Median fragment length (bp)", "GC content (%)", "CpG methylation (%)", "Bisulfite conversion (%)")
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
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "QCs by group",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(method = "t.test")

p_qc_combined <- p_qc_group / p_qc_sample +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.4))

ggsave(file.path(FIG_DIR, "fig1_qc_key_metrics.png"),
       p_qc_combined, width = 14, height = 18, dpi = 300)

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
