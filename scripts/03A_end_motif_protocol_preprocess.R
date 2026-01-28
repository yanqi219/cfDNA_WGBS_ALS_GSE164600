#!/usr/bin/env Rscript
#
# 03A_end_motif_protocol_preprocess.R
#
# Protocol-faithful end-motif preprocessing up to (but not including)
# "End motif analysis and visualization".
#
# Based on STAR Protocols (Liu et al., 2024) end-motif processing:
# - Generate 5' n-mer end motifs from each fragment end
# - Compute fragment-level GC content (reference genome)
# - Filter by fragment size and perform GC-bias correction via LOESS
# - Consolidate per-sample end-motif tables into a single RDS for downstream analysis
#
# Source protocol:
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC11489062/
#
# Notes for this project:
# - We start from gzipped fragment BEDs produced by `scripts/01_qc_processing.R`:
#     chrom, start0, end0, name, mapq, strand, length_bp
#   where start0 = fragment start (0-based), end0 = fragment end (1-based inclusive)
#   so that length_bp == end0 - start0.
# - We do not require external `bedtools`/FASTA files: we use BSgenome hg38 (already used elsewhere).
# - We implement the protocol algorithms directly (same binning, filtering, LOESS GC correction),
#   but avoid creating intermediate BED files on disk for efficiency.
#
# Output:
# - data/processed/end_motif_protocol/20_motif_gc/<sample>_<nmer>bpmotif_gc.rds
# - data/processed/end_motif_protocol/20_motif_gc_stats/<sample>_<nmer>bpmotif_stat.rds
# - data/processed/end_motif_protocol/20_gcbias_plots/<sample>/<sample>_<frac>_<nmer>bpmotif_gc.png
# - data/processed/end_motif_protocol/21_combine_motif/endmotif_<nmer>bp_gc_summary.rds
#
# Usage (from project root, e.g. Take_home_task/):
#   Rscript cfDNA_WGBS_ALS_GSE164600/scripts/03A_end_motif_protocol_preprocess.R \
#     --nmer 4 --chunk_size 200000 --cores 4
#
#

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(parallel))

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

## --- Minimal CLI parsing (protocol-style; no new dependencies) ---
parse_args <- function(argv) {
  # Accepts: --key value
  if (length(argv) == 0) return(list())
  if (any(grepl("^--[^=]+=.*", argv))) {
    stop("Please pass args as '--key value' (no '=') for consistency.")
  }
  keys <- argv[grepl("^--", argv)]
  vals <- argv[!grepl("^--", argv)]
  if (length(keys) != length(vals)) {
    stop("Argument parsing error. Use: --key value (repeat).")
  }
  out <- as.list(vals)
  names(out) <- sub("^--", "", keys)
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

nmer <- as.integer(args$nmer %||% 4L)
chunk_size <- as.integer(args$chunk_size %||% 200000L)
cores <- as.integer(args$cores %||% max(1L, parallel::detectCores(logical = TRUE) - 1L))

if (!is.finite(nmer) || nmer <= 0L) stop("nmer must be a positive integer.")
if (!is.finite(chunk_size) || chunk_size <= 0L) stop("chunk_size must be a positive integer.")
if (!is.finite(cores) || cores <= 0L) stop("cores must be a positive integer.")

## --- Paths ---
PROJECT_DIR <- here::here()
DATA_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600", "data")
FRAG_DIR <- file.path(DATA_DIR, "processed", "fragments")

OUT_BASE <- file.path(DATA_DIR, "processed", "end_motif_protocol")
OUT_GC_DIR <- file.path(OUT_BASE, "20_motif_gc")
OUT_PLOT_DIR <- file.path(OUT_BASE, "20_gcbias_plots")
OUT_STAT_DIR <- file.path(OUT_BASE, "20_motif_gc_stats")
OUT_COMBINE_DIR <- file.path(OUT_BASE, "21_combine_motif")

dir.create(OUT_GC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_STAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_COMBINE_DIR, recursive = TRUE, showWarnings = FALSE)

## --- Protocol functions (verbatim logic) ---
seqlast <- function(from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if (tail(vec, 1) != to) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

gc.correct <- function(coverage, bias) {
  gc_grid <- seqlast(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE), by = 0.01)
  coverage.trend <- loess(coverage ~ bias)

  # Two-stage smoothing as in the protocol, written explicitly to avoid
  # ambiguity about formula environments.
  df_model <- data.frame(
    gc = gc_grid,
    y = predict(coverage.trend, gc_grid)
  )
  coverage.model <- loess(y ~ gc, data = df_model)
  coverage.pred <- predict(coverage.model, data.frame(gc = bias))
  coverage - coverage.pred + median(coverage)
}

## --- Fragment reader (chunked; works for .gz and plain files) ---
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

read_fragments_chunked <- function(path, chunk_size, FUN) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
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
      comment.char = ""
    )
    if (nrow(x) == 0) break
    FUN(x)
  }
  invisible(TRUE)
}

## --- Core: per-sample protocol preprocessing ---
process_one_sample <- function(fragment_bed_path, sample_id, genome, nmer, chunk_size,
                               out_gc_dir, out_plot_dir, out_stat_dir) {
  message(sprintf("[end-motif protocol] sample=%s", sample_id))
  message(sprintf("  fragments: %s", fragment_bed_path))

  # Accumulate counts by (end_motif, frac, gc) without storing per-fragment rows.
  # Key: paste(end_motif, frac, gc, sep = "\t")
  counts_env <- new.env(parent = emptyenv(), hash = TRUE)

  # Running stats
  n_total <- 0L
  n_len_kept <- 0L
  n_gc_tail_kept <- NA_integer_  # computed after tail filtering (from binned totals)

  add_counts <- function(keys, vals) {
    # keys character, vals numeric
    for (i in seq_along(keys)) {
      k <- keys[[i]]
      v <- vals[[i]]
      if (!is.finite(v) || v <= 0) next
      cur <- counts_env[[k]]
      if (is.null(cur)) counts_env[[k]] <- v else counts_env[[k]] <- cur + v
    }
  }

  chunk_fun <- function(df) {
    n_total <<- n_total + nrow(df)

    # Enforce expected column types and basic validity
    df <- df %>%
      mutate(
        chrom = as.character(.data$chrom),
        start0 = as.integer(.data$start0),
        end0 = as.integer(.data$end0),
        length_bp = as.integer(.data$length_bp)
      ) %>%
      filter(
        is.finite(.data$start0),
        is.finite(.data$end0),
        is.finite(.data$length_bp),
        .data$length_bp > 0L,
        (.data$end0 - .data$start0) == .data$length_bp
      )

    # Protocol step 7d: filter by fragment size (100–650 bp), then label fractions
    df <- df %>%
      filter(.data$length_bp >= 100L & .data$length_bp <= 650L) %>%
      mutate(
        frac = case_when(
          .data$length_bp >= 100L & .data$length_bp <= 250L ~ "mono",
          .data$length_bp >= 251L & .data$length_bp <= 450L ~ "di",
          .data$length_bp >= 451L & .data$length_bp <= 650L ~ "tri",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(.data$frac))

    if (nrow(df) == 0) return(invisible(NULL))
    n_len_kept <<- n_len_kept + nrow(df)

    # Convert to 1-based inclusive coordinates used by BSgenome
    frag_start <- df$start0 + 1L
    frag_end <- df$end0

    # Whole-fragment GRanges for GC
    gr_frag <- GRanges(
      seqnames = df$chrom,
      ranges = IRanges(start = frag_start, end = frag_end),
      strand = "*"
    )

    # 5' end motifs (protocol step 5–6 logic, but derived from fragment interval)
    gr_left <- GRanges(
      seqnames = df$chrom,
      ranges = IRanges(start = frag_start, width = nmer),
      strand = "+"
    )
    gr_right <- GRanges(
      seqnames = df$chrom,
      ranges = IRanges(start = frag_end - (nmer - 1L), width = nmer),
      strand = "-"
    )

    # Validate coordinates (avoid negative/zero starts)
    valid_left <- start(gr_left) > 0L
    valid_right <- start(gr_right) > 0L
    valid_frag <- start(gr_frag) > 0L & end(gr_frag) >= start(gr_frag)
    keep_all <- valid_left & valid_right & valid_frag
    if (!any(keep_all)) return(invisible(NULL))

    df <- df[keep_all, , drop = FALSE]
    gr_frag <- gr_frag[keep_all]
    gr_left <- gr_left[keep_all]
    gr_right <- gr_right[keep_all]

    # Extract sequences from reference genome
    seq_frag <- getSeq(genome, gr_frag)
    seq_left <- getSeq(genome, gr_left)
    seq_right <- getSeq(genome, gr_right)

    # Fragment GC content, rounded to 0.01 as in protocol
    gc_mat <- Biostrings::letterFrequency(seq_frag, letters = c("G", "C"), as.prob = TRUE)
    gc <- rowSums(gc_mat)
    gc <- round(gc, 2)

    # Motifs uppercased (protocol uses toupper on getfasta output)
    motif_left <- toupper(as.character(seq_left))
    motif_right <- toupper(as.character(seq_right))

    # Aggregate counts for both ends (each fragment contributes two motifs)
    frac <- df$frac

    keys_left <- paste(motif_left, frac, gc, sep = "\t")
    keys_right <- paste(motif_right, frac, gc, sep = "\t")

    add_counts(keys_left, rep.int(1, length(keys_left)))
    add_counts(keys_right, rep.int(1, length(keys_right)))

    invisible(NULL)
  }

  read_fragments_chunked(fragment_bed_path, chunk_size = chunk_size, FUN = chunk_fun)

  # Materialize aggregated counts
  keys <- ls(envir = counts_env, all.names = TRUE)
  if (length(keys) == 0) {
    warning("No motifs found for sample: ", sample_id)
    df_empty <- tibble(
      end_motif = character(),
      frac = character(),
      gc = numeric(),
      count = numeric()
    )
    out_file <- file.path(out_gc_dir, sprintf("%s_%dbpmotif_gc.rds", sample_id, nmer))
    stat_file <- file.path(out_stat_dir, sprintf("%s_%dbpmotif_stat.rds", sample_id, nmer))
    saveRDS(df_empty, out_file)
    saveRDS(list(n_total = n_total, n_len_kept = n_len_kept, n_gc_tail_kept = 0L), stat_file)
    return(invisible(list(out_file = out_file, stat_file = stat_file)))
  }

  vals <- vapply(keys, function(k) counts_env[[k]], numeric(1))
  parts <- strsplit(keys, "\t", fixed = TRUE)
  df_motif <- tibble(
    end_motif = vapply(parts, `[[`, character(1), 1),
    frac = vapply(parts, `[[`, character(1), 2),
    gc = as.numeric(vapply(parts, `[[`, character(1), 3)),
    count = as.numeric(vals)
  ) %>%
    filter(is.finite(.data$gc), .data$gc >= 0, .data$gc <= 1) %>%
    mutate(frac = factor(.data$frac, levels = c("mono", "di", "tri")))

  # Protocol step 7e: GC-bias correction
  df_gc <- df_motif %>%
    group_by(.data$frac, .data$gc) %>%
    summarise(frac.count = sum(.data$count), .groups = "drop") %>%
    arrange(.data$frac, .data$gc) %>%
    group_by(.data$frac) %>%
    mutate(
      frac.total = sum(.data$frac.count),
      cum_count = cumsum(.data$frac.count),
      percentile = (.data$cum_count / .data$frac.total) * 100
    ) %>%
    ungroup() %>%
    filter(.data$percentile >= 5.00 & .data$percentile <= 95.00) %>%
    group_by(.data$frac) %>%
    mutate(
      frac.count.corrected = gc.correct(.data$frac.count, .data$gc),
      total.corrected = sum(.data$frac.count.corrected),
      gc.scale = .data$frac.count.corrected / .data$frac.count
    ) %>%
    ungroup() %>%
    select(dplyr::all_of(c("frac", "gc", "frac.count", "frac.count.corrected", "gc.scale")))

  # Apply GC correction scalar; tail bins are filtered out by inner_join (protocol behavior)
  df_out <- inner_join(df_motif, df_gc, by = c("frac", "gc")) %>%
    mutate(count.corrected = .data$gc.scale * .data$count) %>%
    group_by(.data$frac) %>%
    mutate(
      frac.total = sum(.data$count),
      frac.total.corrected = sum(.data$count.corrected)
    ) %>%
    ungroup() %>%
    mutate(
      total.count = sum(.data$count),
      total.corrected = sum(.data$count.corrected)
    ) %>%
    arrange(.data$frac, .data$gc, .data$end_motif)

  # Stats after GC tail filtering (counts table; each row is a motif-bin count)
  n_gc_tail_kept <- as.integer(sum(df_out$count) / 2) # approximate unique fragments after tail filter

  df_stats <- list(
    sample_id = sample_id,
    n_total_fragments_in_file = n_total,
    n_fragments_len_filtered = n_len_kept,
    n_fragments_after_gc_tail_filter = n_gc_tail_kept,
    nmer = nmer,
    chunk_size = chunk_size
  )

  # GC-bias plots per fraction (as in protocol)
  frac_levels <- c("mono", "di", "tri")
  for (f in frac_levels) {
    dfp <- df_gc %>% filter(as.character(.data$frac) == f)
    if (nrow(dfp) == 0) next

    p <- ggplot(dfp, aes(x = .data$gc, y = .data$frac.count)) +
      geom_line(linewidth = 0.6) +
      geom_line(aes(y = .data$frac.count.corrected, color = "corrected.counts"), linewidth = 0.6) +
      scale_color_manual(values = c("corrected.counts" = "red"), name = NULL) +
      labs(
        title = sprintf("%s: GC bias (n=%s-mer)", sample_id, nmer),
        subtitle = sprintf("Fraction: %s", f),
        x = "GC (binned, 0.01)",
        y = "Fragment count"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
      )

    out_dir_sample <- file.path(out_plot_dir, sample_id)
    dir.create(out_dir_sample, recursive = TRUE, showWarnings = FALSE)
    out_png <- file.path(out_dir_sample, sprintf("%s_%s_%dbpmotif_gc.png", sample_id, f, nmer))
    ggsave(filename = out_png, plot = p, width = 7, height = 5, dpi = 300)
  }

  out_file <- file.path(out_gc_dir, sprintf("%s_%dbpmotif_gc.rds", sample_id, nmer))
  stat_file <- file.path(out_stat_dir, sprintf("%s_%dbpmotif_stat.rds", sample_id, nmer))
  saveRDS(df_out, out_file)
  saveRDS(df_stats, stat_file)

  invisible(list(out_file = out_file, stat_file = stat_file))
}

## --- Discover inputs ---
if (!dir.exists(FRAG_DIR)) stop("Fragments directory not found: ", FRAG_DIR)
frag_files <- list.files(
  FRAG_DIR,
  pattern = "\\.fragments\\.bed(\\.gz)?$",
  full.names = TRUE
)

if (length(frag_files) == 0) {
  stop("No fragment BED files found in: ", FRAG_DIR, "\nExpected: *.fragments.bed.gz")
}

sample_ids <- basename(frag_files) %>%
  str_replace("\\.fragments\\.bed\\.gz$", "") %>%
  str_replace("\\.fragments\\.bed$", "")

names(frag_files) <- sample_ids

message(sprintf("Found %d fragment files.", length(frag_files)))
message(sprintf("nmer=%d, chunk_size=%d, cores=%d", nmer, chunk_size, cores))

## --- Run per-sample preprocessing (parallel where safe) ---
worker <- function(sid) {
  process_one_sample(
    fragment_bed_path = frag_files[[sid]],
    sample_id = sid,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    nmer = nmer,
    chunk_size = chunk_size,
    out_gc_dir = OUT_GC_DIR,
    out_plot_dir = OUT_PLOT_DIR,
    out_stat_dir = OUT_STAT_DIR
  )
}

if (.Platform$OS.type == "unix" && cores > 1L) {
  invisible(mclapply(sample_ids, worker, mc.cores = cores))
} else {
  invisible(lapply(sample_ids, worker))
}

## --- Protocol step 8 (combine per-sample RDS to one tibble) ---
gc_files <- list.files(OUT_GC_DIR, pattern = sprintf("_%dbpmotif_gc\\.rds$", nmer), full.names = TRUE)
if (length(gc_files) == 0) stop("No per-sample GC-corrected motif files produced in: ", OUT_GC_DIR)

files_list <- lapply(gc_files, readRDS)
ids <- basename(gc_files) %>%
  str_replace(sprintf("_%dbpmotif_gc\\.rds$", nmer), "")

tib_list <- map2(files_list, ids, ~ mutate(as_tibble(.x), id = .y)) %>%
  bind_rows() %>%
  select(id, everything())

out_summary <- file.path(OUT_COMBINE_DIR, sprintf("endmotif_%dbp_gc_summary.rds", nmer))
saveRDS(tib_list, out_summary)

message("Done.")
message("Wrote:")
message(sprintf("  - per-sample GC-corrected motifs: %s", OUT_GC_DIR))
message(sprintf("  - GC-bias plots: %s", OUT_PLOT_DIR))
message(sprintf("  - per-sample stats: %s", OUT_STAT_DIR))
message(sprintf("  - combined summary RDS: %s", out_summary))

