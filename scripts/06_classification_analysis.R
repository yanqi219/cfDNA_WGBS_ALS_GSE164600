# 06_classification_analysis.R
#
# cfDNA WGBS Classification Analysis: ALS vs Control Diagnostic Predictor
#
# Goals:
# - Build diagnostic classifiers for ALS vs Control using Random Forest
# - Compare three feature sets: stackHMM methylation, end motifs, methylation tiles
# - Apply variance-based pre-filtering and Boruta feature selection
# - Use nested cross-validation (LOOCV outer, 5-fold inner) for small sample size
#
# References:
# - nestedcv: Lewis et al., Bioinformatics (2023)
# - Boruta: Kursa & Rudnicki, J Stat Software (2010)

# Setup ----
install_and_load <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "tibble", "purrr",
               "caret", "randomForest", "nestedcv", "Boruta",
               "pROC", "patchwork", "circlize", "here")
bioc_pkgs <- c("bsseq", "GenomicRanges", "IRanges", "GenomeInfoDb", 
               "BSgenome.Hsapiens.UCSC.hg38", "ComplexHeatmap")

invisible(lapply(cran_pkgs, install_and_load, bioc = FALSE))
invisible(lapply(bioc_pkgs, install_and_load, bioc = TRUE))

# Theme and colors
theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
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
}

COLORS <- c("ALS" = "#E64B35", "Ctrl" = "#4DBBD5")
FEATURE_COLORS <- c("stackHMM" = "#8491B4", "Motif" = "#91D1C2", "MethTile" = "#F39B7F")

set.seed(389)

# Paths
PROJECT_DIR <- here::here()
DATA_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/data")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
RESULTS_DIR <- file.path(PROJECT_DIR, "cfDNA_WGBS_ALS_GSE164600/results")
FIG_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
RDS_DIR <- file.path(DATA_DIR, "processed/methylation/rds")

for (d in c(FIG_DIR, TABLE_DIR, RDS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Helper functions ----
filter_by_variance <- function(mat, top_pct = 0.5, min_features = 10) {
  vars <- apply(mat, 2, var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  n_keep <- min(max(min_features, floor(ncol(mat) * top_pct)), sum(vars > 0))
  mat[, order(vars, decreasing = TRUE)[seq_len(n_keep)], drop = FALSE]
}

impute_median <- function(mat) {
  apply(mat, 2, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x })
}

# Load data ----
message("Loading data...")

metadata <- readr::read_csv(file.path(DATA_DIR, "sample_metadata.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Ctrl", "ALS")))

bs <- readRDS(file.path(RDS_DIR, "bsseq_object.rds"))

stackhmm_df <- readr::read_tsv(file.path(PROCESSED_DIR, "hg38_genome_100_segments.bed.gz"),
                               col_names = c("chr", "start", "end", "state"),
                               col_types = "ciic", progress = FALSE) %>%
  dplyr::filter(chr == "chr21")

stackhmm_gr <- GenomicRanges::GRanges(
  seqnames = stackhmm_df$chr,
  ranges = IRanges::IRanges(start = stackhmm_df$start + 1L, end = stackhmm_df$end),
  state = stackhmm_df$state
)
states <- unique(stackhmm_df$state)

motif_long <- readr::read_csv(file.path(TABLE_DIR, "end_motif_4mer_frequencies_long.csv"), show_col_types = FALSE)
tile_mat_raw <- readRDS(file.path(RDS_DIR, "dmrseq_tile_methylation_matrix.rds"))

message(sprintf("  %d samples, %d CpGs, %d stackHMM states, %d tiles",
                nrow(metadata), nrow(bs), length(states), nrow(tile_mat_raw)))

# Feature Engineering ----
message("Engineering features...")

# 1. stackHMM methylation (median per state)
cpg_gr <- GenomicRanges::granges(bs)
meth_mat <- bsseq::getMeth(bs, type = "raw", what = "perBase")
colnames(meth_mat) <- pData(bs)$sample_id

stackhmm_meth <- matrix(NA_real_, nrow = ncol(bs), ncol = length(states),
                        dimnames = list(colnames(meth_mat), states))
for (i in seq_along(states)) {
  cpg_idx <- S4Vectors::queryHits(GenomicRanges::findOverlaps(cpg_gr, stackhmm_gr[stackhmm_gr$state == states[i]]))
  if (length(cpg_idx) > 0) {
    stackhmm_meth[, states[i]] <- apply(meth_mat[cpg_idx, , drop = FALSE], 2, median, na.rm = TRUE)
  }
}

# 2. End motifs (log-transformed frequencies)
motif_wide <- motif_long %>%
  dplyr::select(sample_id, motif, frequency) %>%
  tidyr::pivot_wider(names_from = motif, values_from = frequency, values_fill = 0) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix()
motif_log <- log10(motif_wide[metadata$sample_id, ] + 1e-8)

# 3. Methylation tiles
tile_mat <- t(tile_mat_raw)[metadata$sample_id, ]

# Preprocessing: impute -> filter -> scale
message("Preprocessing features...")

preprocess <- function(mat, top_pct, min_feat, max_missing = 0.2) {
  # 1) Remove features missing in > max_missing fraction of samples
  miss_frac <- colMeans(is.na(mat))
  keep <- miss_frac <= max_missing
  if (!any(keep)) {
    # If everything fails missingness, keep the least-missing features (up to min_feat)
    o <- order(miss_frac, decreasing = FALSE)
    keep_idx <- o[seq_len(min(min_feat, length(o)))]
    mat <- mat[, keep_idx, drop = FALSE]
  } else {
    mat <- mat[, keep, drop = FALSE]
  }
  # 2) Impute remaining missing values (median per feature)
  mat <- impute_median(mat)
  # 3) Filter by variance, then scale
  mat <- filter_by_variance(mat, top_pct, min_feat)
  mat <- scale(mat)
  mat[is.nan(mat)] <- 0
  mat
}

stackhmm_scaled <- preprocess(stackhmm_meth[metadata$sample_id, , drop = FALSE], 0.5, 10)
motif_scaled <- preprocess(motif_log, 0.5, 50)
tile_scaled <- preprocess(tile_mat, 0.05, 100)

message(sprintf("  stackHMM: %d features, Motif: %d features, MethTile: %d features",
                ncol(stackhmm_scaled), ncol(motif_scaled), ncol(tile_scaled)))

feature_sets <- list(stackHMM = stackhmm_scaled, Motif = motif_scaled, MethTile = tile_scaled)
y <- setNames(as.factor(metadata$Group), metadata$sample_id)

# Model Training ----
message("Training Random Forest models with Boruta feature selection...")

N_INNER_FOLDS <- 5
all_results <- list()

for (feat_name in names(feature_sets)) {
  message(sprintf("  %s...", feat_name))
  X <- feature_sets[[feat_name]]
  
  # Boruta feature selection
  set.seed(389)
  boruta_result <- Boruta(x = X, y = y, doTrace = 0, maxRuns = 100)
  selected <- getSelectedAttributes(boruta_result, withTentative = TRUE)
  
  # Fallback if too few features selected
  if (length(selected) < 3) {
    imp <- attStats(boruta_result)
    imp <- imp[rownames(imp) %in% colnames(X), , drop = FALSE]
    selected <- rownames(imp)[order(imp$meanImp, decreasing = TRUE)[seq_len(min(10, nrow(imp)))]]
  }
  selected <- intersect(selected, colnames(X))
  if (length(selected) == 0) selected <- colnames(X)
  
  X_sel <- X[, selected, drop = FALSE]
  
  # Nested CV with Random Forest
  mtry_vals <- unique(c(2, floor(sqrt(ncol(X_sel))), floor(ncol(X_sel) / 3)))
  mtry_vals <- mtry_vals[mtry_vals > 0 & mtry_vals <= ncol(X_sel)]
  
  fit <- nestcv.train(y = y, x = X_sel, method = "rf", outer_method = "LOOCV",
                      n_inner_folds = N_INNER_FOLDS, tuneGrid = data.frame(mtry = mtry_vals),
                      importance = TRUE, cv.cores = 1, verbose = FALSE)
  
  all_results[[feat_name]] <- list(model = fit, boruta = boruta_result, selected_features = selected)
}

saveRDS(all_results, file.path(RDS_DIR, "nested_cv_results.rds"))

# Save tile_scaled matrix with scaling attributes for single-sample prediction
# This ensures consistent scaling when applying the model to new samples
saveRDS(tile_scaled, file.path(RDS_DIR, "tile_scaled_matrix.rds"))

# Save explicit scaling parameters for MethTile (used by single-sample notebook)
meththile_scaling_params <- list(
  center = attr(tile_scaled, "scaled:center"),
  scale = attr(tile_scaled, "scaled:scale"),
  feature_names = colnames(tile_scaled),
  selected_features = all_results[["MethTile"]]$selected_features
)
saveRDS(meththile_scaling_params, file.path(RDS_DIR, "meththile_scaling_params.rds"))

message(sprintf("  Saved tile_scaled matrix (%d samples x %d features) and scaling parameters", 
                nrow(tile_scaled), ncol(tile_scaled)))

# Performance Metrics ----
message("Computing performance metrics...")

extract_metrics <- function(result, feature_set) {
  fit <- result$model
  cm <- fit$summary$table
  
  TP <- if ("ALS" %in% rownames(cm) && "ALS" %in% colnames(cm)) cm["ALS", "ALS"] else 0
  TN <- if ("Ctrl" %in% rownames(cm) && "Ctrl" %in% colnames(cm)) cm["Ctrl", "Ctrl"] else 0
  FP <- if ("ALS" %in% rownames(cm) && "Ctrl" %in% colnames(cm)) cm["ALS", "Ctrl"] else 0
  FN <- if ("Ctrl" %in% rownames(cm) && "ALS" %in% colnames(cm)) cm["Ctrl", "ALS"] else 0
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  f1 <- if ((prec + sens) > 0) 2 * prec * sens / (prec + sens) else 0
  
  auc_val <- auc_lower <- auc_upper <- NA
  if (!is.null(fit$roc)) {
    auc_val <- as.numeric(pROC::auc(fit$roc))
    tryCatch({
      ci <- pROC::ci.auc(fit$roc, method = "delong")
      auc_lower <- ci[1]; auc_upper <- ci[3]
    }, error = function(e) NULL)
  }
  
  tibble::tibble(feature_set = feature_set, n_features = length(result$selected_features),
                 AUC = auc_val, AUC_lower = auc_lower, AUC_upper = auc_upper,
                 Sensitivity = sens, Specificity = spec, Precision = prec, F1 = f1,
                 TP = TP, TN = TN, FP = FP, FN = FN)
}

performance_df <- purrr::map_dfr(names(all_results), ~extract_metrics(all_results[[.x]], .x)) %>%
  dplyr::arrange(desc(AUC))

print(performance_df %>% dplyr::select(feature_set, n_features, AUC, Sensitivity, Specificity, F1))
readr::write_csv(performance_df, file.path(TABLE_DIR, "classification_performance_summary.csv"))

boruta_features_df <- purrr::map_dfr(names(all_results), ~tibble::tibble(
  feature_set = .x, feature = all_results[[.x]]$selected_features))
readr::write_csv(boruta_features_df, file.path(TABLE_DIR, "boruta_selected_features.csv"))

# Visualizations ----
message("Generating visualizations...")

# ROC Curves - use geom_path() and sort data properly for smooth curves
roc_data <- purrr::map_dfr(names(all_results), function(feat_name) {
  roc_obj <- all_results[[feat_name]]$model$roc
  if (is.null(roc_obj)) return(NULL)
  
  df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    feature_set = feat_name,
    auc = as.numeric(pROC::auc(roc_obj))
  )
  df[order(df$fpr, df$tpr), ]
})

## ROC Curves ----
auc_lookup <- roc_data %>% 
  dplyr::distinct(feature_set, auc)
feature_labels <- sapply(names(FEATURE_COLORS), function(fs) {
  auc_val <- auc_lookup$auc[auc_lookup$feature_set == fs]
  sprintf("%s (AUC=%.2f)", fs, auc_val)
})

p_roc <- ggplot(roc_data, aes(x = fpr, y = tpr, color = feature_set)) +
  geom_path(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = FEATURE_COLORS, labels = feature_labels) +
  labs(title = "ROC Curves: ALS vs Control Classification",
       subtitle = "Random Forest with Boruta feature selection (nested LOOCV)",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = "Feature Set") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_pub(base_size = 12) +
  theme(legend.position = c(0.7, 0.25))

ggsave(file.path(FIG_DIR, "fig7A_roc_curves_comparison.png"), p_roc, width = 7, height = 6, dpi = 450)
 
## Performance Barplot ----
perf_long <- performance_df %>%
  dplyr::select(feature_set, AUC, Sensitivity, Specificity, F1) %>%
  tidyr::pivot_longer(-feature_set, names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("AUC", "Sensitivity", "Specificity", "F1")))

p_perf <- ggplot(perf_long, aes(x = feature_set, y = value, fill = feature_set)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 3.5) +
  facet_wrap(~metric, nrow = 1) +
  scale_fill_manual(values = FEATURE_COLORS) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Classification Performance by Feature Set", x = NULL, y = "Score") +
  theme_pub(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_DIR, "fig7B_performance_barplot.png"), p_perf, width = 10, height = 5, dpi = 450)

## Feature Importance ----
# Uses importance from the final RF model trained on all selected features
imp_plots <- list()

for (feat_name in names(all_results)) {
  fit <- all_results[[feat_name]]$model
  
  # Extract RF importance from final model (MeanDecreaseGini for classification)
  rf_imp <- fit$final_fit$finalModel$importance
  
  imp_df <- data.frame(
    feature = rownames(rf_imp),
    importance = rf_imp[, "MeanDecreaseGini"]
  ) %>%
    dplyr::arrange(desc(importance))
  
  n_feat <- nrow(imp_df)
  
  imp_plots[[feat_name]] <- ggplot(imp_df, aes(x = reorder(feature, importance), y = importance)) +
    geom_col(fill = FEATURE_COLORS[feat_name], alpha = 0.85) +
    coord_flip() +
    labs(title = feat_name,
         subtitle = sprintf("n = %d features", n_feat),
         x = NULL, y = "Mean Decrease Gini") +
    theme_pub(base_size = 10) +
    theme(plot.title = element_text(size = 12))
}

# Combine plots: one row, three columns
p_imp_combined <- patchwork::wrap_plots(imp_plots, nrow = 1) +
  patchwork::plot_annotation(
    title = "Random Forest Feature Importance by Feature Set",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(FIG_DIR, "fig7C_rf_importance_combined.png"), 
       p_imp_combined, width = 14, height = 6, dpi = 450)

## stackHMM Heatmap ----
sample_order <- metadata %>% dplyr::arrange(Group) %>% dplyr::pull(sample_id)
anno_df <- data.frame(Group = metadata$Group[match(sample_order, metadata$sample_id)],
                      row.names = sample_order)

ha <- ComplexHeatmap::HeatmapAnnotation(
  df = anno_df, col = list(Group = COLORS),
  annotation_name_gp = grid::gpar(fontface = "bold", fontsize = 10)
)

ht <- ComplexHeatmap::Heatmap(
  t(stackhmm_scaled[sample_order, ]), name = "Z-score", top_annotation = ha,
  cluster_rows = TRUE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Samples", row_title = sprintf("stackHMM States (n=%d)", ncol(stackhmm_scaled)),
  col = circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426"))
)

png(file.path(FIG_DIR, "fig7D_stackhmm_methylation_heatmap.png"), width = 10, height = 10, units = "in", res = 450)
ComplexHeatmap::draw(ht)
invisible(dev.off())

## EhnA3 plot ----
plot_matrix <- stackhmm_meth %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(metadata, by = "sample_id")

p_enha3 <- ggplot(plot_matrix, aes(x = Group, y = `45_EnhA3`, fill = Group)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 2, alpha = 0.8) +
  scale_fill_manual(values = COLORS) +
  labs(x = NULL, y = "EhnA3 methylation level") +
  theme_pub(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig7E_enha3_methylation_boxplot.png"), p_enha3, width = 6, height = 5, dpi = 450)

# Combied
ht_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(ht, heatmap_legend_side = "left",
                        annotation_legend_side = "left")
)
p_ht <- patchwork::wrap_elements(full = ht_grob)

p_classification_combined <- p_perf / p_imp_combined / c(p_roc | p_ht | p_enha3) +
  patchwork::plot_layout(heights = c(1, 1, 1.2)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "fig7_classification_analysis_combined.png"), p_classification_combined, width = 14, height = 12, dpi = 450)

# Sample Predictions using MethTile RF model ----
all_results <- readRDS(file.path(RDS_DIR, "nested_cv_results.rds"))

message("\nGenerating sample predictions using MethTile model...")

# Get the MethTile model and selected features
methtile_result <- all_results[["MethTile"]]
methtile_fit <- methtile_result$model
methtile_features <- methtile_result$selected_features

# Get the final trained RF model
final_rf <- methtile_fit$final_fit

# Prepare feature matrix for prediction (use the scaled tile data with selected features)
X_pred <- tile_scaled[, methtile_features, drop = FALSE]

# Predict class and probabilities
pred_class <- predict(final_rf, newdata = as.data.frame(X_pred))
pred_probs <- predict(final_rf, newdata = as.data.frame(X_pred), type = "prob")

# Create prediction summary table
prediction_df <- data.frame(
  sample_id = metadata$sample_id,
  true_group = metadata$Group,
  predicted_group = pred_class,
  prob_Ctrl = round(pred_probs[, "Ctrl"], 3),
  prob_ALS = round(pred_probs[, "ALS"], 3),
  correct = ifelse(as.character(metadata$Group) == as.character(pred_class), "Yes", "No")
)

# Sort by ALS probability (descending)
prediction_df <- prediction_df %>%
  dplyr::arrange(desc(prob_ALS))

# Print summary
message("\nMethTile RF Model - Sample Predictions:")
message(sprintf("  Features used: %d (%s)", 
                length(methtile_features), 
                paste(methtile_features, collapse = ", ")))
message(sprintf("  Accuracy: %.1f%% (%d/%d correct)",
                100 * sum(prediction_df$correct == "Yes") / nrow(prediction_df),
                sum(prediction_df$correct == "Yes"),
                nrow(prediction_df)))

print(prediction_df)