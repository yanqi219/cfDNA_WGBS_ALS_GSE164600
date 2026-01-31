# cfDNA_analysis_functions.R
#
# Reusable functions for cfDNA WGBS analysis
# Extracted from scripts 01-06 for use in single-sample analysis notebook
#
# This module ensures numerical reproducibility by using identical algorithms
# as the original cohort analysis scripts.

# Suppress R CMD check notes for NSE variables used in dplyr pipelines
utils::globalVariables(c(

  "qname", "rname", "pos", "cigar", "flag", "start0", "end0_excl", "rname_1", 
  "start0_1", "annotation", "annotation_simple", ".data", "chrom", "end0",
  "length_bp", "mapq", "strand", "frequency", "count", "motif"
))

# =============================================================================
# PUBLICATION THEME
# =============================================================================

#' Publication-quality ggplot2 theme
#' @param base_size Base font size
#' @param base_family Font family
#' @return ggplot2 theme object
theme_pub <- function(base_size = 14, base_family = "Helvetica") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, size = ggplot2::rel(1.25),
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, size = ggplot2::rel(1.07),
        margin = ggplot2::margin(b = 6)
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.08),
        margin = ggplot2::margin(t = 9)
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.08),
        margin = ggplot2::margin(r = 9)
      ),
      axis.text = ggplot2::element_text(color = "black", size = ggplot2::rel(1.00)),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.45),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.7),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.08)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.00)),
      legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      legend.justification = "center",
      legend.box.background = ggplot2::element_rect(color = NA, fill = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "gray85", linewidth = 0.30, linetype = "dotted"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
      strip.text = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.05),
        margin = ggplot2::margin(t = 3, b = 3)
      ),
      plot.margin = ggplot2::margin(16, 16, 16, 16)
    )
}

# Color palette
COLORS <- c("ALS" = "#E64B35", "Ctrl" = "#4DBBD5")

# =============================================================================
# QC FUNCTIONS (from 01_qc_processing.R)
# =============================================================================

#' Get BAM flag statistics (similar to samtools flagstat)
#' @param bam_path Path to BAM file
#' @return data.frame with flag statistics
get_flagstat <- function(bam_path) {
  param_all <- Rsamtools::ScanBamParam(what = "flag")
  bf <- Rsamtools::BamFile(bam_path, yieldSize = 5e6)
  open(bf)
  
  total <- 0
  mapped <- 0
  paired <- 0
  proper_pair <- 0
  duplicates <- 0
  secondary <- 0
  supplementary <- 0
  
  repeat {
    flags <- Rsamtools::scanBam(bf, param = param_all)[[1]]$flag
    if (length(flags) == 0) break
    
    total <- total + length(flags)
    mapped <- mapped + sum(bitwAnd(flags, 4) == 0)
    paired <- paired + sum(bitwAnd(flags, 1) > 0)
    proper_pair <- proper_pair + sum(bitwAnd(flags, 2) > 0)
    duplicates <- duplicates + sum(bitwAnd(flags, 1024) > 0)
    secondary <- secondary + sum(bitwAnd(flags, 256) > 0)
    supplementary <- supplementary + sum(bitwAnd(flags, 2048) > 0)
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

#' Parse Bismark XM methylation tags
#' @param xm_tags Vector of XM tag strings
#' @return data.frame with methylation counts by context
parse_xm_tags <- function(xm_tags) {
  all_chars <- paste(xm_tags, collapse = "")
  
  data.frame(
    Z = stringr::str_count(all_chars, "Z"),
    z = stringr::str_count(all_chars, "z"),
    X = stringr::str_count(all_chars, "X"),
    x = stringr::str_count(all_chars, "x"),
    H = stringr::str_count(all_chars, "H"),
    h = stringr::str_count(all_chars, "h")
  )
}

#' Extract BAM metrics and generate fragment BED
#' @param bam_path Path to BAM file
#' @param sample_id Sample identifier
#' @param genome BSgenome object
#' @param chunk_size Number of reads per chunk
#' @param frag_dir Output directory for fragment BED
#' @return List with QC metrics
extract_bam_metrics <- function(bam_path, sample_id, genome, chunk_size = 1e6, frag_dir = ".") {
  
  message(sprintf("Processing: %s", sample_id))
  message(sprintf("BAM file: %s", bam_path))
  
  if (!file.exists(bam_path)) {
    warning(sprintf("BAM file not found: %s", bam_path))
    return(NULL)
  }
  
  flagstat <- get_flagstat(bam_path)
  
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isPaired = TRUE, isProperPair = TRUE,
      isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
      isDuplicate = FALSE,
      isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE
    ),
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "isize", "seq", "qwidth"),
    tag = c("XM")
  )
  
  bf <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)
  open(bf)
  
  frag_lengths <- integer()
  meth_calls <- list()
  total_reads <- 0
  mapq_values <- integer()
  gc_values <- numeric()
  
  fragment_bed_path <- file.path(frag_dir, paste0(sample_id, ".fragments.bed.gz"))
  frag_con <- gzfile(fragment_bed_path, open = "wt")
  on.exit(try(close(frag_con), silent = TRUE), add = TRUE)
  
  pending_aln <- tibble::tibble(
    qname = character(), flag = integer(), rname = character(),
    pos = integer(), cigar = character(), mapq = integer()
  )
  
  repeat {
    chunk_raw <- Rsamtools::scanBam(bf, param = param)[[1]]
    if (length(chunk_raw$qname) == 0) break
    
    mapq_values <- c(mapq_values, chunk_raw$mapq)
    keep_idx <- which(chunk_raw$mapq >= 30)
    chunk <- lapply(chunk_raw, function(x) {
      if (length(x) == length(chunk_raw$mapq)) return(x[keep_idx])
      if (is.list(x)) return(lapply(x, function(tag_vector) tag_vector[keep_idx]))
      return(x)
    })
    
    total_reads <- total_reads + length(chunk$qname)
    
    if (length(chunk$qname) > 0) {
      aln <- tibble::tibble(
        qname = as.character(chunk$qname),
        flag = as.integer(chunk$flag),
        rname = as.character(chunk$rname),
        pos = as.integer(chunk$pos),
        cigar = as.character(chunk$cigar),
        mapq = as.integer(chunk$mapq)
      ) %>%
        dplyr::filter(!is.na(qname), !is.na(rname), !is.na(pos), !is.na(cigar),
                      is.finite(pos), pos > 0L)
      
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
          
          pending_aln <- aln_all %>% dplyr::filter(!(qname %in% have_pair))
          
          pair_rows <- pair_rows[order(pair_rows$qname), , drop = FALSE]
          mate_raw <- stats::ave(pair_rows$qname, pair_rows$qname, FUN = seq_along)
          
          idx1 <- which(mate_raw == 1L)
          idx2 <- which(mate_raw == 2L)
          
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
            same_chr <- pair_w$rname_1 == pair_w$rname_2
            
            if (any(same_chr)) {
              frag_start0 <- pair_w$start0_1[same_chr]
              frag_end0 <- pair_w$end0_excl_2[same_chr]
              length_bp <- frag_end0 - frag_start0
              
              keep_len <- is.finite(length_bp) & length_bp > 0L & length_bp <= 1000L
              if (any(keep_len)) {
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
                }, error = function(e) rep(NA_real_, length(gr_frag)))
                
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
    
    if (!is.null(chunk$tag$XM)) {
      xm_tags <- chunk$tag$XM[!is.na(chunk$tag$XM)]
      if (length(xm_tags) > 0) {
        meth_calls[[length(meth_calls) + 1]] <- parse_xm_tags(xm_tags)
      }
    }
    
    rm(chunk)
    gc(verbose = FALSE)
  }
  
  close(bf)
  
  meth_summary <- if (length(meth_calls) > 0) {
    do.call(rbind, meth_calls) %>% dplyr::summarise(dplyr::across(dplyr::everything(), sum))
  } else {
    data.frame(Z = 0, z = 0, X = 0, x = 0, H = 0, h = 0)
  }
  
  list(
    sample_id = sample_id,
    fragment_bed = fragment_bed_path,
    flagstat = flagstat,
    filtered_reads = total_reads,
    fragment_lengths = frag_lengths,
    gc_content = gc_values,
    mapq = mapq_values,
    methylation = meth_summary
  )
}

#' Calculate summary statistics from BAM metrics
#' @param metrics List from extract_bam_metrics
#' @return data.frame with summary statistics
calculate_summary_stats <- function(metrics) {
  if (is.null(metrics)) return(data.frame(sample_id = NA))
  
  fl <- metrics$fragment_lengths
  gc <- metrics$gc_content
  meth <- metrics$methylation
  fs <- metrics$flagstat
  
  cpg_total <- meth$Z + meth$z
  cpg_meth_rate <- if (cpg_total > 0) meth$Z / cpg_total else NA
  
  chh_total <- meth$H + meth$h
  chh_meth_rate <- if (chh_total > 0) meth$H / chh_total else NA
  bisulfite_conv <- if (!is.na(chh_meth_rate)) 1 - chh_meth_rate else NA
  
  data.frame(
    sample_id = metrics$sample_id,
    total_reads = fs$total_reads,
    mapped_reads = fs$mapped,
    unmapped_reads = fs$unmapped,
    properly_paired = fs$properly_paired,
    duplicates = fs$duplicates,
    mapping_rate = fs$mapping_rate,
    proper_pair_rate = fs$proper_pair_rate,
    duplicate_rate = fs$duplicate_rate,
    filtered_reads = metrics$filtered_reads,
    filter_pass_rate = metrics$filtered_reads / fs$total_reads,
    gc_content = mean(gc, na.rm = TRUE),
    n_fragments = length(fl),
    mean_frag_length = mean(fl),
    median_frag_length = median(fl),
    sd_frag_length = sd(fl),
    frag_mode = as.numeric(names(sort(table(fl), decreasing = TRUE))[1]),
    short_frag_ratio = sum(fl < 150) / length(fl),
    mono_nuc_ratio = sum(fl >= 150 & fl <= 220) / length(fl),
    di_nuc_ratio = sum(fl >= 300 & fl <= 400) / length(fl),
    mean_mapq = mean(metrics$mapq),
    cpg_methylation = cpg_meth_rate,
    chh_methylation = chh_meth_rate,
    bisulfite_conversion = bisulfite_conv
  )
}

# =============================================================================
# FRAGMENTATION FUNCTIONS (from 02_fragmentation_analysis.R)
# =============================================================================

# Fragment BED column definitions
FRAG_COLS <- c("chrom", "start0", "end0", "name", "mapq", "strand", "length_bp", "gc")
FRAG_COL_CLASSES <- c(
  chrom = "character", start0 = "integer", end0 = "integer", name = "character",
  mapq = "integer", strand = "character", length_bp = "integer", gc = "numeric"
)

#' Read fragment BED file in chunks
#' @param path Path to gzipped fragment BED
#' @param chunk_size Number of rows per chunk
#' @param FUN Function to apply to each chunk
read_fragments_chunked <- function(path, chunk_size = 200000L, FUN) {
  chunk_size <- as.integer(chunk_size)
  con <- gzfile(path, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  
  repeat {
    x <- utils::read.delim(
      file = con, header = FALSE, sep = "\t", nrows = chunk_size,
      col.names = FRAG_COLS, colClasses = FRAG_COL_CLASSES,
      stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE
    )
    if (nrow(x) == 0) break
    FUN(x)
  }
  invisible(TRUE)
}

#' Convert fragment BED chunk to midpoint GRanges
#' @param df Fragment BED data.frame chunk
#' @return GRanges object with fragment midpoints
chunk_to_midpoint_gr <- function(df) {
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

#' Convert fragment BED chunk to 5' start GRanges
#' @param df Fragment BED data.frame chunk
#' @return GRanges object with 5' start positions
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

#' Moving average smoothing
#' @param x Numeric vector
#' @param k Window size (odd number)
#' @return Smoothed vector
moving_average <- function(x, k = 21L) {
  k <- as.integer(k)
  if (k < 3L) return(x)
  if (k %% 2L == 0L) k <- k + 1L
  w <- rep(1 / k, k)
  as.numeric(stats::filter(x, w, sides = 2, method = "convolution"))
}

#' GC-bias correction using LOESS
#' @param counts Numeric vector of counts
#' @param gc_values GC content values per bin
#' @param fit_mask Logical mask for fitting
#' @param pseudocount Pseudocount for log transform
#' @param span LOESS span parameter
#' @return Corrected counts
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
  
  if (sum(fit_mask) < 20L) return(counts)
  
  fit <- stats::loess(
    y[fit_mask] ~ gc_values[fit_mask],
    span = span, degree = 1, family = "symmetric",
    control = stats::loess.control(surface = "direct")
  )
  
  pred <- stats::predict(fit, newdata = gc_values)
  med <- stats::median(counts[fit_mask] + pseudocount, na.rm = TRUE)
  corrected <- (counts + pseudocount) / exp(pred) * med
  
  corrected[!fit_mask_in] <- NA_real_
  corrected[!is.finite(pred)] <- NA_real_
  corrected
}

#' Compute bin-level short/long fragment counts for one sample
#' @param path Path to fragment BED
#' @param bins_gr GRanges of genomic bins
#' @param bins_df Data.frame with bin metadata including gc
#' @param bins_keep_for_gc Logical vector for GC fitting
#' @param short_range Short fragment length range
#' @param long_range Long fragment length range
#' @return tibble with bin-level counts and ratios
bin_counts_one_sample <- function(path, bins_gr, bins_df, bins_keep_for_gc,
                                   short_range = c(90L, 150L), long_range = c(151L, 220L)) {
  short_counts <- integer(length(bins_gr))
  long_counts <- integer(length(bins_gr))
  
  cb <- function(x) {
    l <- x$length_bp
    keep <- !is.na(l) & l >= 30L & l <= 500L
    if (!any(keep)) return(invisible())
    x <- x[keep, , drop = FALSE]
    l <- x$length_bp
    
    mid_gr <- chunk_to_midpoint_gr(x)
    
    idx_s <- which(l >= short_range[1] & l <= short_range[2])
    if (length(idx_s) > 0) {
      hits <- GenomicRanges::findOverlaps(mid_gr[idx_s], bins_gr, ignore.strand = TRUE)
      if (length(hits) > 0) {
        tab <- tabulate(S4Vectors::subjectHits(hits), nbins = length(bins_gr))
        short_counts <<- short_counts + tab
      }
    }
    
    idx_l <- which(l >= long_range[1] & l <= long_range[2])
    if (length(idx_l) > 0) {
      hits <- GenomicRanges::findOverlaps(mid_gr[idx_l], bins_gr, ignore.strand = TRUE)
      if (length(hits) > 0) {
        tab <- tabulate(S4Vectors::subjectHits(hits), nbins = length(bins_gr))
        long_counts <<- long_counts + tab
      }
    }
    invisible()
  }
  
  read_fragments_chunked(path = path, chunk_size = 200000L, FUN = cb)
  
  short_corrected <- correct_gc_bias(short_counts, bins_df$gc, fit_mask = bins_keep_for_gc)
  long_corrected <- correct_gc_bias(long_counts, bins_df$gc, fit_mask = bins_keep_for_gc)
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

#' Annotate fragment 5' starts with genomic features
#' @param fragment_bed_gz Path to fragment BED
#' @param sample_id Sample identifier
#' @param txdb TxDb object for annotation
#' @param chunk_size Chunk size for reading
#' @return tibble with annotation counts
annotate_fragment_starts_one_sample <- function(fragment_bed_gz, sample_id, txdb, chunk_size = 200000L) {
  counts <- c(Promoter = 0L, Exon = 0L, Intron = 0L, Distal_Intergenic = 0L,
              Three_UTR = 0L, Five_UTR = 0L, Downstream = 0L)
  total <- 0L
  
  read_fragments_chunked(
    path = fragment_bed_gz,
    chunk_size = chunk_size,
    FUN = function(x) {
      gr <- chunk_to_5p_start_gr(x)
      if (length(gr) == 0) return(invisible())
      total <<- total + length(gr)
      
      gr <- ChIPseeker::annotatePeak(gr, TxDb = txdb, level = "gene", addFlankGeneInfo = TRUE)
      gr_anno <- data.frame(gr@anno)
      gr_anno <- gr_anno %>%
        dplyr::mutate(annotation_simple = gsub(" \\(.*", "", annotation)) %>%
        dplyr::mutate(annotation_simple = dplyr::case_when(
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

#' Compute TSS enrichment profile from fragment 5' starts
#' @param fragment_bed_gz Path to fragment BED
#' @param tss_windows GRanges of TSS windows
#' @param tss_site TSS positions
#' @param tss_strand TSS strands
#' @param flank Flanking distance in bp
#' @param chunk_size Chunk size for reading
#' @param edge_width Width of edge region for normalization
#' @return tibble with enrichment profile
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
      pos <- GenomicRanges::start(gr[q])
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

# =============================================================================
# END MOTIF FUNCTIONS (from 03_end_motif_analysis.R)
# =============================================================================

#' Generate all 256 4-mer combinations
#' @return Character vector of 4-mers
all_4mers <- function() {
  bases <- c("A", "C", "G", "T")
  grid <- expand.grid(b1 = bases, b2 = bases, b3 = bases, b4 = bases, stringsAsFactors = FALSE)
  apply(grid, 1, paste0, collapse = "")
}

#' Calculate Shannon entropy in bits
#' @param p Probability vector
#' @return Entropy in bits
shannon_entropy_bits <- function(p) {
  p <- p[is.finite(p) & p > 0]
  -sum(p * log2(p))
}

#' Extract 4-mer end motif counts from fragment BED
#' @param frag_bed_gz Path to fragment BED
#' @param sample_id Sample identifier
#' @param genome BSgenome object
#' @param valid_seqnames Valid chromosome names
#' @param motifs_all Vector of all 4-mers
#' @param chunk_size Chunk size for reading
#' @return tibble with motif counts
extract_4mer_counts_from_frag_bed <- function(frag_bed_gz, sample_id, genome, valid_seqnames,
                                               motifs_all, chunk_size = 200000L) {
  counts <- setNames(integer(length(motifs_all)), motifs_all)
  con <- gzfile(frag_bed_gz, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  
  total_extracted <- 0L
  
  repeat {
    x <- tryCatch(
      utils::read.delim(
        file = con, header = FALSE, sep = "\t", nrows = as.integer(chunk_size),
        col.names = FRAG_COLS, colClasses = FRAG_COL_CLASSES,
        stringsAsFactors = FALSE, quote = ""
      ),
      error = function(e) NULL
    )
    if (is.null(x) || nrow(x) == 0) break
    
    x <- x %>%
      dplyr::filter(
        .data$chrom %in% valid_seqnames,
        is.finite(.data$start0), is.finite(.data$end0),
        .data$start0 >= 0L, .data$end0 > .data$start0
      )
    if (nrow(x) == 0) next
    
    gr_l <- GenomicRanges::GRanges(
      seqnames = x$chrom,
      ranges = IRanges::IRanges(start = x$start0 + 1L, width = 4L),
      strand = "+"
    )
    
    r_start0 <- x$end0 - 4L
    keep_r <- is.finite(r_start0) & r_start0 >= 0L
    gr_r <- GenomicRanges::GRanges(
      seqnames = x$chrom[keep_r],
      ranges = IRanges::IRanges(start = r_start0[keep_r] + 1L, width = 4L),
      strand = "-"
    )
    
    s_l <- BSgenome::getSeq(genome, gr_l)
    s_r <- if (length(gr_r) > 0) BSgenome::getSeq(genome, gr_r) else Biostrings::DNAStringSet()
    
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
    stop(sprintf("No motifs extracted for %s. Check contig names and genome.", sample_id), call. = FALSE)
  }
  
  message(sprintf("  %s: extracted %s motifs", sample_id, format(total_extracted, big.mark = ",")))
  tibble::tibble(sample_id = sample_id, motif = names(counts), count = as.integer(counts))
}

# =============================================================================
# CLASSIFICATION FUNCTIONS (from 06_classification_analysis.R)
# =============================================================================

#' Filter features by variance
#' @param mat Feature matrix (samples x features)
#' @param top_pct Top percentage to keep
#' @param min_features Minimum features to keep
#' @return Filtered matrix
filter_by_variance <- function(mat, top_pct = 0.5, min_features = 10) {
  vars <- apply(mat, 2, var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  n_keep <- min(max(min_features, floor(ncol(mat) * top_pct)), sum(vars > 0))
  mat[, order(vars, decreasing = TRUE)[seq_len(n_keep)], drop = FALSE]
}

#' Impute missing values with column medians
#' @param mat Feature matrix
#' @return Matrix with imputed values
impute_median <- function(mat) {
  apply(mat, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
}

#' Preprocess feature matrix for classification
#' @param mat Feature matrix (samples x features)
#' @param top_pct Top percentage of features by variance
#' @param min_feat Minimum features to keep
#' @param max_missing Maximum fraction of missing values per feature
#' @return Preprocessed and scaled matrix
preprocess_features <- function(mat, top_pct, min_feat, max_missing = 0.2) {
  miss_frac <- colMeans(is.na(mat))
  keep <- miss_frac <= max_missing
  
  if (!any(keep)) {
    o <- order(miss_frac, decreasing = FALSE)
    keep_idx <- o[seq_len(min(min_feat, length(o)))]
    mat <- mat[, keep_idx, drop = FALSE]
  } else {
    mat <- mat[, keep, drop = FALSE]
  }
  
  mat <- impute_median(mat)
  mat <- filter_by_variance(mat, top_pct, min_feat)
  mat <- scale(mat)
  mat[is.nan(mat)] <- 0
  mat
}

#' Predict class for new sample using pre-trained model
#' @param new_features Named numeric vector or single-row matrix of features
#' @param model_result Result from nested CV training
#' @return List with predicted class and probabilities
predict_single_sample <- function(new_features, model_result) {
  fit <- model_result$model
  selected_features <- model_result$selected_features
  
  if (is.vector(new_features)) {
    new_features <- matrix(new_features, nrow = 1, dimnames = list(NULL, names(new_features)))
  }
  
  missing_features <- setdiff(selected_features, colnames(new_features))
  if (length(missing_features) > 0) {
    warning(sprintf("Missing %d features: %s", length(missing_features),
                    paste(head(missing_features, 5), collapse = ", ")))
    fill_mat <- matrix(0, nrow = 1, ncol = length(missing_features),
                       dimnames = list(NULL, missing_features))
    new_features <- cbind(new_features, fill_mat)
  }
  
  new_features <- new_features[, selected_features, drop = FALSE]
  
  pred_class <- predict(fit$final_fit, newdata = as.data.frame(new_features))
  pred_prob <- predict(fit$final_fit, newdata = as.data.frame(new_features), type = "prob")
  
  list(
    predicted_class = as.character(pred_class),
    probabilities = pred_prob
  )
}
