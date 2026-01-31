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

# Figure 0: Age distribution by group ----
p_age <- ggplot(metadata, aes(x = Group, y = Age, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
  scale_fill_manual(values = COLORS) +
  labs(title = "Age distribution by group", x = "Group", y = "Age") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.2),
                              margin = ggplot2::margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, margin = ggplot2::margin(b = 6)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3, linetype = "dotted"),
    plot.margin = ggplot2::margin(12, 12, 12, 12)
  )

ggsave(file.path(FIG_DIR, "fig0_age_distribution_by_group.png"), p_age, width = 6, height = 5, dpi = 450)

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