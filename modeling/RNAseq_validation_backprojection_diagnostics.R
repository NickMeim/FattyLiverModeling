## Diagnostic script only.
## This does not source or modify RNAseq_validation_focused.R.
## It checks whether RNAval-derived directions actually span the LIV2TRANS
## phenotype-predictive and extra-LV directions used in backprojection.

suppressPackageStartupMessages({
  library(ropls)
})

results_dir <- "../results"
rnaval_dir <- "../RNAseq_validation"
diagnostics_dir <- "rnaval_backprojection_diagnostics"
dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)

contrast_conditions <- c("Beta", "BetaDMSO", "BetaM", "AlphaBeta")
condition_levels <- c("T2D", contrast_conditions)

sample_key <- data.frame(
  sampleID = c("BetaDMSO1", "BetaDMSO2", "BetaDMSO3",
               "BetaM1", "BetaM2", "BetaM3",
               "T2D1", "T2D2", "T2D3",
               "Beta1", "Beta2", "Beta3",
               "AlphaBeta1", "AlphaBeta2", "AlphaBeta3"),
  source_id = c("2C72ND_10", "2C72ND_11", "2C72ND_12",
                "2C72ND_13", "2C72ND_14", "2C72ND_15",
                "2C72ND_1", "2C72ND_2", "2C72ND_3",
                "2C72ND_4", "2C72ND_5", "2C72ND_6",
                "2C72ND_7", "2C72ND_8", "2C72ND_9"),
  condition = c("BetaDMSO", "BetaDMSO", "BetaDMSO",
                "BetaM", "BetaM", "BetaM",
                "T2D", "T2D", "T2D",
                "Beta", "Beta", "Beta",
                "AlphaBeta", "AlphaBeta", "AlphaBeta"),
  stringsAsFactors = FALSE
)

metadata_template <- data.frame(
  sampleID = sample_key$sampleID,
  condition = factor(sample_key$condition, levels = condition_levels),
  stringsAsFactors = FALSE
)

log2cpm <- function(m) {
  apply(m, 2, function(x) log2(1 + x / sum(x) * 1e6))
}

row_center <- function(m) {
  sweep(m, 1, rowMeans(m), "-")
}

read_rnaval_counts <- function() {
  expr_file <- file.path(rnaval_dir, "2C72ND-expression-matrix.tsv")
  gene_expression <- read.delim(expr_file, check.names = FALSE)
  measure_cols <- paste0(sample_key$source_id, "_count")

  gene_name <- as.character(gene_expression$gene_name)
  gene_name[gene_name == "" | is.na(gene_name)] <-
    as.character(gene_expression$gene_id[gene_name == "" | is.na(gene_name)])

  expr_mat <- as.matrix(gene_expression[, measure_cols, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  colnames(expr_mat) <- sample_key$sampleID

  mat <- rowsum(expr_mat, group = gene_name, reorder = FALSE)
  mat <- round(mat)
  storage.mode(mat) <- "numeric"
  mat
}

build_rnaval_liv2trans_basis <- function(count_matrix, metadata, genes_use) {
  counts <- count_matrix[genes_use, , drop = FALSE]
  lcpm <- log2cpm(counts)
  centered <- row_center(lcpm)

  meta <- data.frame(sampleName = metadata$sampleID,
                     condition = as.character(metadata$condition),
                     stringsAsFactors = FALSE)
  meta <- meta[match(colnames(centered), meta$sampleName), , drop = FALSE]
  stopifnot(identical(meta$sampleName, colnames(centered)))

  groups <- factor(meta$condition, levels = condition_levels)
  groups_use <- levels(groups)[levels(groups) %in% as.character(groups)]
  data_grouped <- sapply(groups_use, function(g) {
    rowMeans(centered[, groups == g, drop = FALSE])
  })
  rownames(data_grouped) <- rownames(centered)
  colnames(data_grouped) <- groups_use

  data_grouped_c <- row_center(data_grouped)
  pca_group <- prcomp(t(data_grouped_c), center = FALSE, scale. = FALSE)
  keep_pcs <- seq_len(ncol(pca_group$rotation) - 1)
  Wm_rnaval <- pca_group$rotation[, keep_pcs, drop = FALSE]
  colnames(Wm_rnaval) <- paste0("RNAval_PC", seq_len(ncol(Wm_rnaval)))

  list(Wm = Wm_rnaval,
       Xm = t(centered),
       grouped_centered = data_grouped_c,
       pca = pca_group)
}

orthonormalize_basis <- function(W, prefix = "basis", tol = 1e-10) {
  W <- as.matrix(W)
  keep <- colSums(!is.finite(W)) == 0
  W <- W[, keep, drop = FALSE]
  if (ncol(W) == 0) {
    out <- matrix(0, nrow = nrow(W), ncol = 0)
    rownames(out) <- rownames(W)
    return(out)
  }
  qr_w <- qr(W, tol = tol)
  if (qr_w$rank == 0) {
    out <- matrix(0, nrow = nrow(W), ncol = 0)
    rownames(out) <- rownames(W)
    return(out)
  }
  Q <- qr.Q(qr_w, complete = FALSE)[, seq_len(qr_w$rank), drop = FALSE]
  rownames(Q) <- rownames(W)
  colnames(Q) <- paste0(prefix, seq_len(ncol(Q)))
  Q
}

cbind_aligned <- function(..., prefixes = NULL) {
  mats <- list(...)
  mats <- lapply(mats, as.matrix)
  common <- Reduce(intersect, lapply(mats, rownames))
  out <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
  rownames(out) <- common
  if (!is.null(prefixes)) {
    colnames(out) <- unlist(mapply(function(prefix, m) {
      paste0(prefix, seq_len(ncol(m)))
    }, prefixes, mats, SIMPLIFY = FALSE))
  }
  out
}

projection_capture <- function(target, basis, basis_name) {
  common <- Reduce(intersect, list(rownames(target), rownames(basis)))
  target <- as.matrix(target[common, , drop = FALSE])
  target_names <- colnames(target)
  if (is.null(target_names) || length(target_names) != ncol(target)) {
    target_names <- paste0("target", seq_len(ncol(target)))
  }
  basis <- orthonormalize_basis(basis[common, , drop = FALSE], prefix = "Q")
  values <- sapply(seq_len(ncol(target)), function(i) {
    denom <- sum(target[, i]^2)
    if (denom == 0 || ncol(basis) == 0) return(NA_real_)
    sum(crossprod(basis, target[, i, drop = FALSE])^2) / denom
  })
  data.frame(
    basis = basis_name,
    target = target_names,
    capture = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

cosine_matrix <- function(A, B) {
  common <- Reduce(intersect, list(rownames(A), rownames(B)))
  A <- as.matrix(A[common, , drop = FALSE])
  B <- as.matrix(B[common, , drop = FALSE])
  A <- sweep(A, 2, sqrt(colSums(A^2)), "/")
  B <- sweep(B, 2, sqrt(colSums(B^2)), "/")
  crossprod(A, B)
}

predict_backprojected_plsr <- function(Xh, Yh, Wh, Bh, W_basis) {
  common <- Reduce(intersect, list(colnames(Xh), rownames(Wh), rownames(W_basis)))
  X <- as.matrix(Xh[, common, drop = FALSE])
  W <- orthonormalize_basis(W_basis[common, , drop = FALSE], prefix = "Q")
  Wh_common <- as.matrix(Wh[common, , drop = FALSE])

  projected_scores <- X %*% W %*% t(W) %*% Wh_common
  yhat <- cbind(1, projected_scores) %*%
    rbind(colMeans(Yh, na.rm = TRUE), Bh)
  colnames(yhat) <- colnames(Yh)
  rownames(yhat) <- rownames(Xh)
  yhat
}

prediction_metrics <- function(actual, predicted, input) {
  do.call(rbind, lapply(colnames(actual), function(pheno) {
    ok <- is.finite(actual[, pheno]) & is.finite(predicted[, pheno])
    data.frame(
      input = input,
      phenotype = pheno,
      r = cor(actual[ok, pheno], predicted[ok, pheno], method = "pearson"),
      rho = cor(actual[ok, pheno], predicted[ok, pheno], method = "spearman"),
      MAE = mean(abs(actual[ok, pheno] - predicted[ok, pheno])),
      RMSE = sqrt(mean((actual[ok, pheno] - predicted[ok, pheno])^2)),
      stringsAsFactors = FALSE
    )
  }))
}

cat("Loading saved LIV2TRANS objects and RNAval counts...\n")
tt <- readRDS(file.path(results_dir, "liv2trans_results_Govaere_Kostrzewski.rds"))
processed <- readRDS(file.path(results_dir, "processed_data_list_govaere_kostrzewski.rds"))
count_matrix <- read_rnaval_counts()

Xh <- as.matrix(processed$Xh)
Yh <- as.matrix(processed$Yh)
model_trained <- tt$model
Wh <- as.matrix(tt$Wh)
Bh <- t(model_trained@weightMN) %*% model_trained@coefficientMN
phi <- Wh %*% Bh
colnames(phi) <- colnames(Yh)

Wm_mps_tt <- as.matrix(tt$W_invitro)
Wm_mps_processed <- as.matrix(processed$Wm)
Wm_extra_saved <- as.matrix(tt$W_opt)
Wm_total_saved <- as.matrix(tt$W_tot)
Wm_tc_saved <- as.matrix(tt$W_translatable)

if (is.null(colnames(Wm_extra_saved))) {
  colnames(Wm_extra_saved) <- paste0("extraLV", seq_len(ncol(Wm_extra_saved)))
}
if (is.null(colnames(Wm_tc_saved))) {
  colnames(Wm_tc_saved) <- paste0("TC", seq_len(ncol(Wm_tc_saved)))
}

genes_for_rnaval_basis <- Reduce(intersect, list(
  colnames(Xh),
  rownames(Wh),
  rownames(Wm_mps_tt),
  rownames(count_matrix)
))

rnaval_liv2trans <- build_rnaval_liv2trans_basis(
  count_matrix = count_matrix,
  metadata = metadata_template,
  genes_use = genes_for_rnaval_basis
)
Wm_rnaval <- rnaval_liv2trans$Wm

condition_diffs <- sapply(contrast_conditions, function(condition_name) {
  rnaval_liv2trans$grouped_centered[, condition_name] -
    rnaval_liv2trans$grouped_centered[, "T2D"]
})
rownames(condition_diffs) <- rownames(rnaval_liv2trans$grouped_centered)
colnames(condition_diffs) <- paste0(contrast_conditions, "_minus_T2D")

basis_list <- list(
  MPS_tt_W_invitro = Wm_mps_tt,
  MPS_processed_Wm = Wm_mps_processed,
  TC_saved = Wm_tc_saved,
  Extra_saved_only = Wm_extra_saved,
  MPS_plus_saved_extra = cbind_aligned(Wm_mps_tt, Wm_extra_saved),
  RNAval_PCs = Wm_rnaval,
  MPS_plus_RNAval_PCs = cbind_aligned(Wm_mps_tt, Wm_rnaval),
  RNAval_condition_diffs = condition_diffs,
  MPS_plus_condition_diffs = cbind_aligned(Wm_mps_tt, condition_diffs)
)

for (condition_name in contrast_conditions) {
  single_diff <- condition_diffs[, paste0(condition_name, "_minus_T2D"), drop = FALSE]
  basis_list[[paste0("MPS_plus_", condition_name, "_diff")]] <-
    cbind_aligned(Wm_mps_tt, single_diff)
}

cat("Computing subspace capture diagnostics...\n")
phi_capture <- do.call(rbind, lapply(names(basis_list), function(nm) {
  projection_capture(phi, basis_list[[nm]], nm)
}))

extra_capture <- do.call(rbind, lapply(names(basis_list), function(nm) {
  projection_capture(Wm_extra_saved, basis_list[[nm]], nm)
}))

condition_diff_capture <- do.call(rbind, lapply(
  list(MPS_tt_W_invitro = Wm_mps_tt,
       TC_saved = Wm_tc_saved,
       Extra_saved_only = Wm_extra_saved,
       MPS_plus_saved_extra = cbind_aligned(Wm_mps_tt, Wm_extra_saved)),
  function(W) projection_capture(condition_diffs, W, deparse(substitute(W)))
))
condition_diff_capture$basis <- rep(
  c("MPS_tt_W_invitro", "TC_saved", "Extra_saved_only",
    "MPS_plus_saved_extra"),
  each = ncol(condition_diffs)
)

cos_rnaval_extra <- cosine_matrix(Wm_rnaval, Wm_extra_saved)
cos_condition_extra <- cosine_matrix(condition_diffs, Wm_extra_saved)

Wm_mps_q <- orthonormalize_basis(Wm_mps_tt, prefix = "MPS")
common_phi <- Reduce(intersect, list(rownames(phi), rownames(Wm_mps_q)))
phi_residual <- phi[common_phi, , drop = FALSE] -
  Wm_mps_q[common_phi, , drop = FALSE] %*%
  crossprod(Wm_mps_q[common_phi, , drop = FALSE], phi[common_phi, , drop = FALSE])
cos_rnaval_phi_residual <- cosine_matrix(Wm_rnaval, phi_residual)
cos_condition_phi_residual <- cosine_matrix(condition_diffs, phi_residual)

cat("Computing backprojected performance diagnostics...\n")
yhat_full <- as.matrix(predict(model_trained, Xh))
metrics <- prediction_metrics(Yh, yhat_full, "Full trained PLSR")
for (nm in names(basis_list)) {
  yhat <- predict_backprojected_plsr(Xh, Yh, Wh, Bh, basis_list[[nm]])
  metrics <- rbind(metrics, prediction_metrics(Yh, yhat, nm))
}

axis_basis <- cbind_aligned(Wm_tc_saved, Wm_extra_saved)
axis_common <- intersect(rownames(axis_basis), rownames(rnaval_liv2trans$grouped_centered))
axis_scores <- t(rnaval_liv2trans$grouped_centered[axis_common, , drop = FALSE]) %*%
  axis_basis[axis_common, , drop = FALSE]
axis_scores <- data.frame(condition = rownames(axis_scores),
                          axis_scores,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)

condition_diff_axis_scores <- t(condition_diffs[axis_common, , drop = FALSE]) %*%
  axis_basis[axis_common, , drop = FALSE]
condition_diff_norm <- sqrt(colSums(condition_diffs[axis_common, , drop = FALSE]^2))
condition_diff_axis_scores <- data.frame(
  contrast = rownames(condition_diff_axis_scores),
  norm = condition_diff_norm[rownames(condition_diff_axis_scores)],
  condition_diff_axis_scores,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pca_variance <- data.frame(
  PC = paste0("RNAval_PC", seq_along(rnaval_liv2trans$pca$sdev)),
  variance = rnaval_liv2trans$pca$sdev^2,
  pct_variance = 100 * rnaval_liv2trans$pca$sdev^2 /
    sum(rnaval_liv2trans$pca$sdev^2),
  cumulative_pct = cumsum(100 * rnaval_liv2trans$pca$sdev^2 /
                            sum(rnaval_liv2trans$pca$sdev^2)),
  stringsAsFactors = FALSE
)

dimensions <- data.frame(
  object = c("Xh genes", "W_invitro PCs", "processed Wm PCs",
             "RNAval PCs", "condition diffs", "saved extra LVs",
             "saved TCs", "saved total basis"),
  n_rows = c(ncol(Xh), nrow(Wm_mps_tt), nrow(Wm_mps_processed),
             nrow(Wm_rnaval), nrow(condition_diffs), nrow(Wm_extra_saved),
             nrow(Wm_tc_saved), nrow(Wm_total_saved)),
  n_cols = c(nrow(Xh), ncol(Wm_mps_tt), ncol(Wm_mps_processed),
             ncol(Wm_rnaval), ncol(condition_diffs), ncol(Wm_extra_saved),
             ncol(Wm_tc_saved), ncol(Wm_total_saved)),
  stringsAsFactors = FALSE
)

orthogonality <- data.frame(
  basis = names(basis_list),
  rank_after_qr = sapply(basis_list, function(W) {
    common <- Reduce(intersect, list(rownames(W), colnames(Xh), rownames(Wh)))
    ncol(orthonormalize_basis(W[common, , drop = FALSE], prefix = "Q"))
  }),
  max_crossprod_error_after_qr = sapply(basis_list, function(W) {
    common <- Reduce(intersect, list(rownames(W), colnames(Xh), rownames(Wh)))
    Q <- orthonormalize_basis(W[common, , drop = FALSE], prefix = "Q")
    if (ncol(Q) == 0) return(NA_real_)
    max(abs(crossprod(Q) - diag(ncol(Q))))
  }),
  stringsAsFactors = FALSE
)

ggplot(metrics_df, aes(x = input, y = rho, fill = input)) +
    geom_col(color = "black", width = 0.72) +
    geom_text(aes(label = sprintf("%.2f", rho)), vjust = -0.35, size = 3.5) +
    facet_wrap(~ phenotype, nrow = 1) +
    coord_cartesian(ylim = c(min(0, min(metrics_df$rho, na.rm = TRUE) - 0.08), 1)) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Backprojected PLSR Spearman performance",
         x = NULL, y = "Spearman rho") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
write.csv(metrics,
          file.path(diagnostics_dir, "metrics.csv"),
          row.names = FALSE)
write.csv(phi_capture,
          file.path(diagnostics_dir, "phi_capture.csv"),
          row.names = FALSE)
write.csv(extra_capture,
          file.path(diagnostics_dir, "extraLV_capture.csv"),
          row.names = FALSE)
write.csv(axis_scores,
          file.path(diagnostics_dir, "condition_axis_scores.csv"),
          row.names = FALSE)
write.csv(condition_diff_axis_scores,
          file.path(diagnostics_dir, "condition_diff_axis_scores.csv"),
          row.names = FALSE)
write.csv(condition_diff_capture,
          file.path(diagnostics_dir, "condition_diff_capture.csv"),
          row.names = FALSE)
write.csv(pca_variance,
          file.path(diagnostics_dir, "rnaval_grouped_pca_variance.csv"),
          row.names = FALSE)
write.csv(dimensions,
          file.path(diagnostics_dir, "dimensions.csv"),
          row.names = FALSE)
write.csv(orthogonality,
          file.path(diagnostics_dir, "orthogonality.csv"),
          row.names = FALSE)
write.csv(cos_rnaval_extra,
          file.path(diagnostics_dir, "cos_RNAvalPC_extraLV.csv"))
write.csv(cos_condition_extra,
          file.path(diagnostics_dir, "cos_conditionDiff_extraLV.csv"))
write.csv(cos_rnaval_phi_residual,
          file.path(diagnostics_dir, "cos_RNAvalPC_phiResidual.csv"))
write.csv(cos_condition_phi_residual,
          file.path(diagnostics_dir, "cos_conditionDiff_phiResidual.csv"))

saveRDS(
  list(
    metrics = metrics,
    phi_capture = phi_capture,
    extra_capture = extra_capture,
    condition_diff_capture = condition_diff_capture,
    axis_scores = axis_scores,
    condition_diff_axis_scores = condition_diff_axis_scores,
    pca_variance = pca_variance,
    dimensions = dimensions,
    orthogonality = orthogonality,
    cos_rnaval_extra = cos_rnaval_extra,
    cos_condition_extra = cos_condition_extra,
    cos_rnaval_phi_residual = cos_rnaval_phi_residual,
    cos_condition_phi_residual = cos_condition_phi_residual,
    Wm_rnaval = Wm_rnaval,
    condition_diffs = condition_diffs
  ),
  file.path(diagnostics_dir, "diagnostics.rds")
)

cat("\nDimensions:\n")
print(dimensions)

cat("\nBackprojection metrics:\n")
print(metrics[order(metrics$phenotype, metrics$input), ], row.names = FALSE)

cat("\nPhi capture:\n")
print(phi_capture[order(phi_capture$target, phi_capture$basis), ],
      row.names = FALSE)

cat("\nSaved extra-LV capture:\n")
print(extra_capture[order(extra_capture$target, extra_capture$basis), ],
      row.names = FALSE)

cat("\nRNAval condition mean scores on saved TC/extra axes:\n")
print(axis_scores, row.names = FALSE)

cat("\nCondition-difference capture by saved spaces:\n")
print(condition_diff_capture[order(condition_diff_capture$target,
                                   condition_diff_capture$basis), ],
      row.names = FALSE)

cat("\nRNAval grouped PCA variance:\n")
print(pca_variance, row.names = FALSE)

cat("\nWrote diagnostic outputs to rnaval_backprojection_diagnostics/.\n")
