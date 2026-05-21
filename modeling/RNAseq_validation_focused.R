library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(matrixStats)
library(ropls)
library(LIV2Trans)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)

source("enrichment_calculations.R")

fig_dir <- "../figures"
results_dir <- "../results"
rnaval_dir <- "../RNAseq_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

contrast_conditions <- c("Beta", "BetaDMSO", "BetaM", "AlphaBeta")
condition_levels <- c("T2D", contrast_conditions)
condition_plot_labels <- c(
  T2D = "T2D",
  Beta = expression(paste("TGF", beta)),
  BetaDMSO = expression(paste("TGF", beta," + DMSO")),
  BetaM = expression(paste("TGF", beta," + Mef. acid")),
  AlphaBeta = expression(paste("TGF", beta," + IFN",alpha))
)
contrast_plot_labels <- paste(condition_plot_labels[contrast_conditions], "vs T2D")
names(contrast_plot_labels) <- paste(contrast_conditions, "vs T2D")

sample_key <- tibble::tribble(
  ~sampleID,     ~source_id,  ~condition,
  "BetaDMSO1",  "2C72ND_10", "BetaDMSO",
  "BetaDMSO2",  "2C72ND_11", "BetaDMSO",
  "BetaDMSO3",  "2C72ND_12", "BetaDMSO",
  "BetaM1",     "2C72ND_13", "BetaM",
  "BetaM2",     "2C72ND_14", "BetaM",
  "BetaM3",     "2C72ND_15", "BetaM",
  "T2D1",       "2C72ND_1",  "T2D",
  "T2D2",       "2C72ND_2",  "T2D",
  "T2D3",       "2C72ND_3",  "T2D",
  "Beta1",      "2C72ND_4",  "Beta",
  "Beta2",      "2C72ND_5",  "Beta",
  "Beta3",      "2C72ND_6",  "Beta",
  "AlphaBeta1", "2C72ND_7",  "AlphaBeta",
  "AlphaBeta2", "2C72ND_8",  "AlphaBeta",
  "AlphaBeta3", "2C72ND_9",  "AlphaBeta"
)

metadata_template <- sample_key %>%
  transmute(sampleID, condition = factor(condition, levels = condition_levels))

save_plot <- function(filename, plot, width, height, dpi = 600) {
  stem <- file.path(fig_dir, tools::file_path_sans_ext(filename))
  ggsave(paste0(stem, ".png"), plot,
         width = width, height = height, dpi = dpi, units = "in")
  ggsave(paste0(stem, ".pdf"), plot,
         width = width, height = height, units = "in")
}

sig_label <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p <= 0.0001, "****",
                ifelse(p <= 0.001, "***",
                       ifelse(p <= 0.01, "**",
                              ifelse(p <= 0.05, "*", "")))))
}

log2cpm <- function(m) {
  apply(m, 2, function(x) log2(1 + x / sum(x) * 1e6))
}

row_center <- function(m) {
  sweep(m, 1, rowMeans(m), "-")
}

read_rnaval_expression <- function(measure = c("cpm", "count")) {
  measure <- match.arg(measure)
  expr_file <- file.path(rnaval_dir, "2C72ND-expression-matrix.tsv")
  gene_expression <- read.delim(expr_file, check.names = FALSE)
  measure_cols <- paste0(sample_key$source_id, "_", measure)

  gene_expression <- gene_expression %>%
    dplyr::select(gene_id, gene_name, gene_biotype, all_of(measure_cols))
  colnames(gene_expression) <- c("gene_id", "gene_name", "gene_biotype",
                                 sample_key$sampleID)

  gene_expression <- gene_expression %>%
    mutate(gene_name = ifelse(gene_name == "", gene_id, gene_name))

  gene_annot <- gene_expression %>%
    group_by(gene_name) %>%
    summarise(
      gene_ids = paste(sort(unique(gene_id)), collapse = ";"),
      gene_biotypes = paste(sort(unique(gene_biotype)), collapse = ";"),
      n_ids = n_distinct(gene_id),
      .groups = "drop"
    )

  expr_by_name <- gene_expression %>%
    dplyr::select(-gene_id, -gene_biotype) %>%
    group_by(gene_name) %>%
    summarise(across(all_of(sample_key$sampleID), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")

  mat <- expr_by_name %>%
    column_to_rownames("gene_name") %>%
    as.matrix()

  if (measure == "count") {
    mat <- round(mat)
    storage.mode(mat) <- "integer"
  } else {
    storage.mode(mat) <- "numeric"
  }

  list(matrix = mat, gene_annot = gene_annot)
}

name_basis <- function(W, prefix) {
  if (is.null(colnames(W)) || any(is.na(colnames(W))) || any(colnames(W) == "")) {
    colnames(W) <- paste0(prefix, seq_len(ncol(W)))
  }
  W
}

align_basis_and_matrix <- function(X, W) {
  shared <- intersect(colnames(X), rownames(W))
  if (length(shared) == 0) {
    stop("No shared genes between expression matrix and basis.")
  }
  list(X = X[, shared, drop = FALSE], W = W[shared, , drop = FALSE])
}

plot_projection <- function(df, x, y, title, xlab, ylab) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = condition)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 3.5, shape = 21, color = "black") +
    ggrepel::geom_text_repel(aes(label = sampleID), size = 3,
                             color = "black", max.overlaps = Inf) +
    scale_fill_brewer(palette = "Set2", drop = FALSE,
                      breaks = names(condition_plot_labels),
                      labels = condition_plot_labels) +
    labs(title = title, x = xlab, y = ylab, fill = "Condition") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
}

make_stat_matrix <- function(contrasts_list) {
  common_genes <- Reduce(intersect, lapply(contrasts_list, `[[`, "gene_name"))
  stat_mat <- sapply(contrasts_list, function(df) {
    df$stat[match(common_genes, df$gene_name)]
  })
  rownames(stat_mat) <- common_genes
  stat_mat
}

format_hallmark_name <- function(x) {
  x <- sub("^FL1000_MSIG_H_HALLMARK_", "", x)
  x <- str_replace_all(x, "_", " ")
  paste0(toupper(substr(tolower(x), 1, 1)), substring(tolower(x), 2))
}

run_hallmark_gsea <- function(stat_mat, n_permutations = 10000) {
  ids <- mapIds(org.Hs.eg.db, keys = rownames(stat_mat),
                column = "ENTREZID", keytype = "SYMBOL",
                multiVals = "first")
  keep <- which(!is.na(ids))
  meas <- stat_mat[keep, , drop = FALSE]
  rownames(meas) <- ids[keep]

  hm <- fastenrichment(colnames(meas), rownames(meas), meas,
                       enrichment_space = "msig_db_h",
                       n_permutations = n_permutations,
                       order_columns = FALSE)

  nes_mat <- as.matrix(hm$NES$`NES MSIG Hallmark`)
  pval_mat <- as.matrix(hm$Pval$`Pval MSIG Hallmark`)
  colnames(nes_mat) <- colnames(meas)
  colnames(pval_mat) <- colnames(meas)

  left_join(
    as.data.frame(nes_mat) %>%
      rownames_to_column("Hallmark") %>%
      pivot_longer(-Hallmark, names_to = "contrast", values_to = "NES"),
    as.data.frame(pval_mat) %>%
      rownames_to_column("Hallmark") %>%
      pivot_longer(-Hallmark, names_to = "contrast", values_to = "padj"),
    by = c("Hallmark", "contrast")
  ) %>%
    mutate(Hallmark = format_hallmark_name(Hallmark),
           condition = factor(sub(" vs T2D$", "", contrast),
                              levels = contrast_conditions))
}

plot_hallmark_gsea <- function(df_hm) {
  df_hm_plot <- df_hm %>%
    group_by(contrast) %>%
    group_modify(~{
      sig <- .x %>% filter(padj <= 0.05)
      if (nrow(sig) == 0) .x %>% slice_max(abs(NES), n = 15) else sig
    }) %>%
    ungroup() %>%
    mutate(contrast = factor(contrast, levels = names(contrast_plot_labels),
                             labels = unname(contrast_plot_labels))) %>%
    arrange(contrast, NES) %>%
    mutate(HM_ord = paste(Hallmark, contrast, sep = "___"),
           HM_ord = factor(HM_ord, levels = unique(HM_ord)))

  offset <- max(abs(df_hm_plot$NES), na.rm = TRUE) * 0.04

  ggplot(df_hm_plot, aes(x = NES, y = HM_ord, fill = NES)) +
    geom_col() +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_gradient2(low = "darkblue", high = "indianred",
                         mid = "whitesmoke", midpoint = 0) +
    geom_text(aes(label = sig_label(padj),
                  x = ifelse(NES < 0, NES - offset, NES + offset)),
              size = 5, color = "black") +
    facet_wrap(~ contrast, ncol = 2, scales = "free_y") +
    labs(title = "Hallmark GSEA by condition vs T2D",
         x = "Normalized enrichment score", y = "Hallmark") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

build_rnaval_liv2trans_basis <- function(count_matrix, metadata, genes_use) {
  counts <- count_matrix[genes_use, , drop = FALSE]
  lcpm <- log2cpm(counts)
  centered <- row_center(lcpm)

  meta <- metadata %>%
    transmute(sampleName = sampleID, condition = as.character(condition))
  meta <- meta[match(colnames(centered), meta$sampleName), ]
  stopifnot(identical(meta$sampleName, colnames(centered)))

  groups <- factor(meta$condition, levels = condition_levels)
  groups_use <- levels(groups)[levels(groups) %in% groups]
  data_grouped <- sapply(groups_use, function(g) {
    rowMeans(centered[, groups == g, drop = FALSE])
  })
  rownames(data_grouped) <- rownames(centered)
  colnames(data_grouped) <- groups_use

  data_grouped_c <- row_center(data_grouped)
  pca_group <- prcomp(t(data_grouped_c), center = FALSE, scale. = FALSE)
  keep_pcs <- seq_len(ncol(pca_group$rotation)-1)
  Wm_rnaval <- pca_group$rotation[, keep_pcs, drop = FALSE]
  colnames(Wm_rnaval) <- paste0("RNAval_PC", seq_len(ncol(Wm_rnaval)))

  list(Wm = Wm_rnaval,
       Xm = t(centered),
       grouped_centered = data_grouped_c,
       pca = pca_group)
}

orthonormalize_basis <- function(W, prefix = "aug_basis") {
  W <- W[, colSums(!is.finite(W)) == 0, drop = FALSE]
  qr_w <- qr(W)
  Q <- qr.Q(qr_w, complete = FALSE)
  Q <- Q[, seq_len(qr_w$rank), drop = FALSE]
  rownames(Q) <- rownames(W)
  colnames(Q) <- paste0(prefix, seq_len(ncol(Q)))
  Q
}

combine_orthonormal <- function(W1, W2, tol = 1e-10, reorth = TRUE) {
  # Component of W2 orthogonal to span(W1)
  W2_orth <- W2 - W1 %*% crossprod(W1, W2)   # crossprod(W1, W2) == t(W1) %*% W2
  
  # Optional re-projection ("twice is enough") to kill floating-point leakage
  if (reorth) {
    W2_orth <- W2_orth - W1 %*% crossprod(W1, W2_orth)
  }
  
  # Orthonormalize the residual; rank detection drops directions W1 already covers
  qrW <- qr(W2_orth, tol = tol)
  r   <- qrW$rank
  if (r == 0) return(W1)               # W2 fully explained by W1
  
  Q <- qr.Q(qrW)[, seq_len(r), drop = FALSE]
  cbind(W1, Q)
}
augment_by_regression <- function(W1, W2, tol = 1e-8, prefix = "aug") {
  basis <- W1
  added <- vector("list", ncol(W2))
  norms <- numeric(ncol(W2))          # how "new" each W2 PC is
  k <- 0L
  for (j in seq_len(ncol(W2))) {
    r  <- lm.fit(x = basis, y = W2[, j])$residuals   # part of PC ⟂ basis
    nr <- sqrt(sum(r^2))
    norms[j] <- nr
    if (nr > tol) {                   # skip PCs already fully explained
      q <- r / nr
      basis <- cbind(basis, q)        # grow basis BEFORE next column
      k <- k + 1L
      added[[k]] <- q
    }
  }
  aug <- if (k > 0L) {
    A <- do.call(cbind, added[seq_len(k)])
    colnames(A) <- paste0(prefix, seq_len(k))
    A
  } else NULL
  list(combined        = if (is.null(aug)) W1 else cbind(W1, aug),
       augmentation    = aug,
       residual_norms  = norms)
}

predict_backprojected_plsr <- function(Xh, Yh, model, Wh, W_basis) {
  shared <- Reduce(intersect, list(colnames(Xh), rownames(Wh), rownames(W_basis)))
  X <- Xh[, shared, drop = FALSE]
  W <- W_basis[shared, , drop = FALSE]
  Wh <- Wh[shared, , drop = FALSE]
  Bh <- t(model@weightMN) %*% model@coefficientMN

  projected_scores <- X %*% W %*% t(W) %*% Wh
  yhat <- cbind(1, projected_scores) %*%
    rbind(colMeans(Yh, na.rm = TRUE), Bh)
  colnames(yhat) <- colnames(Yh)
  rownames(yhat) <- rownames(Xh)
  yhat
}

prediction_long <- function(actual, predicted, input) {
  bind_rows(lapply(colnames(actual), function(pheno) {
    tibble(
      sample = rownames(actual),
      phenotype = pheno,
      measured = actual[, pheno],
      predicted = predicted[, pheno],
      input = input
    )
  }))
}

prediction_metrics <- function(pred_df) {
  pred_df %>%
    group_by(input, phenotype) %>%
    summarise(
      r = cor(measured, predicted, method = "pearson", use = "complete.obs"),
      rho = cor(measured, predicted, method = "spearman", use = "complete.obs"),
      MAE = mean(abs(measured - predicted), na.rm = TRUE),
      RMSE = sqrt(mean((measured - predicted)^2, na.rm = TRUE)),
      .groups = "drop"
    )
}

plot_backprojection_scatter <- function(pred_df) {
  ggplot(pred_df, aes(measured, predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey45") +
    geom_point(size = 2.2, alpha = 0.72, color = "steelblue") +
    ggpubr::stat_cor(method = "spearman", label.x.npc = "left",
                     label.y.npc = "top", size = 3.4) +
    facet_grid(phenotype ~ input, scales = "free") +
    labs(title = "Trained PLSR performance after backprojection",
         x = "Measured", y = "Predicted") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text.x = element_text(size = 9))
}

plot_backprojection_metrics <- function(metrics_df) {
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
}

## Project raw RNAval into the already-loaded/previously identified MPS TCs and -----------------
## extra LVs. This uses the saved Kostrzewski TC and extra-LV bases and keeps the
## sample-level RNAval projection plots only.
Wm_extra_loaded <- readRDS(file.path(results_dir, "Wm_kostrzewski_extra.rds")) %>%
  as.matrix() %>%
  name_basis("extraLV")
Wm_TC_loaded <- readRDS(file.path(results_dir, "Wm_kostrzewski_combo.rds")) %>%
  as.matrix() %>%
  name_basis("TC")

projection_basis_loaded <- cbind(Wm_TC_loaded, Wm_extra_loaded)
projection_basis_loaded <- projection_basis_loaded[, c("TC1", "TC2",
                                                       "extraLV1", "extraLV2")]

rnaval_cpm <- read_rnaval_expression("cpm")
Xval_raw <- log2(rnaval_cpm$matrix + 1) %>%
  t()
Xval_raw <- sweep(Xval_raw, 2, colMeans(Xval_raw), "-")

aligned_raw <- align_basis_and_matrix(Xval_raw, projection_basis_loaded)
Z_raw <- aligned_raw$X %*% aligned_raw$W
df_raw_projection <- as.data.frame(Z_raw) %>%
  rownames_to_column("sampleID") %>%
  left_join(metadata_template, by = "sampleID")

p_raw_tc <- plot_projection(df_raw_projection, "TC1", "TC2",
                            "Raw RNAval projected into loaded TCs",
                            "TC1", "TC2")
p_raw_extra <- plot_projection(df_raw_projection, "extraLV1", "extraLV2",
                               "Raw RNAval projected into loaded extra LVs",
                               "Extra LV1", "Extra LV2")
p_raw_tc_extra <- plot_projection(df_raw_projection, "TC1", "extraLV1",
                                  "Raw RNAval projected into TC1 and extra LV1",
                                  "TC1", "Extra LV1")

save_plot("RNAval_raw_loaded_TC_extraLV_projection.png",
          (p_raw_tc | p_raw_extra) / p_raw_tc_extra,
          width = 13, height = 10)
save_plot("RNAval_raw_loaded_TCs.png", p_raw_tc,
          width = 7.5, height = 5.5)
save_plot("RNAval_raw_loaded_extraLVs.png", p_raw_extra,
          width = 7.5, height = 5.5)
save_plot("RNAval_raw_loaded_TC1_extraLV1.png", p_raw_tc_extra,
          width = 7.5, height = 5.5)
saveRDS(df_raw_projection,
        file.path(results_dir, "rnaval_raw_loaded_TC_extraLV_projection.rds"))

## Run differential expression, GSEA, and inferred pathway activity. Then plot-------------------------
## how selected Hallmark and PROGENy pathways change in each RNAval condition
## relative to the T2D control.

rnaval_counts <- read_rnaval_expression("count")
count_matrix <- rnaval_counts$matrix
gene_annot <- rnaval_counts$gene_annot

metadata_deseq <- metadata_template %>%
  as.data.frame()
rownames(metadata_deseq) <- metadata_deseq$sampleID
metadata_deseq <- metadata_deseq[colnames(count_matrix), , drop = FALSE]
metadata_deseq$condition <- factor(metadata_deseq$condition,
                                   levels = condition_levels)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata_deseq,
                              design = ~ condition)
dds <- DESeq(dds)

get_results <- function(dds, condition_name, annot, ref = "T2D", alpha = 0.05) {
  res <- results(dds,
                 contrast = c("condition", condition_name, ref),
                 alpha = alpha)
  as.data.frame(res) %>%
    rownames_to_column("gene_name") %>%
    left_join(annot, by = "gene_name") %>%
    arrange(padj)
}

contrasts_list <- setNames(
  lapply(contrast_conditions, function(condition_name) {
    get_results(dds, condition_name, gene_annot) %>%
      filter(!is.na(stat))
  }),
  paste(contrast_conditions, "vs T2D")
)

stat_mat <- make_stat_matrix(contrasts_list)

net_prog <- decoupleR::get_progeny(organism = "human", top = 500)
pwy_act <- decoupleR::run_viper(stat_mat, net_prog,
                                minsize = 1, verbose = TRUE) %>%
  dplyr::select(-any_of("statistic"))

df_pwy <- pwy_act %>%
  dplyr::rename(contrast = condition, Pathway = source) %>%
  mutate(condition = factor(sub(" vs T2D$", "", contrast),
                            levels = contrast_conditions))

pathway_order <- c("JAK-STAT", "Hypoxia", "EGFR", "WNT", "MAPK", "NFkB",
                   "Androgen", "VEGF", "TNFa", "Trail", "PI3K", "Estrogen",
                   "TGFb", "p53")

lim_pwy <- max(abs(df_pwy$score), na.rm = TRUE) * 1.15
p_pwy_conditions <- ggplot(
  df_pwy %>% mutate(Pathway = factor(Pathway, levels = rev(pathway_order))),
  aes(x = condition, y = Pathway, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", score)), size = 3) +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = 0,
                       limits = c(-lim_pwy, lim_pwy)) +
  scale_x_discrete(labels = condition_plot_labels) +
  labs(title = "PROGENy pathway activity: each condition vs T2D",
       x = "Condition", y = "Pathway", fill = "Activity") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_plot("rnaval_condition_vs_T2D_pathway_heatmap.png",
          p_pwy_conditions, width = 8, height = 6.8)

progeny_pathways_of_interest <- c("p53", "JAK-STAT", "TGFb", "NFkB", "Androgen")
df_pwy_selected <- df_pwy %>%
  filter(Pathway %in% progeny_pathways_of_interest) %>%
  mutate(Pathway = factor(Pathway, levels = progeny_pathways_of_interest),
         condition = factor(condition, levels = contrast_conditions))

p_pwy_selected <- ggplot(df_pwy_selected,
                         aes(x = condition, y = score, fill = score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65") +
  geom_col(color = "black", width = 0.72) +
  facet_wrap(~ Pathway, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = condition_plot_labels) +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = 0) +
  labs(title = "Selected PROGENy pathway activities by condition vs T2D",
       x = "Condition", y = "Pathway activity", fill = "Activity") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_plot("rnaval_selected_PROGENy_activity_vs_T2D.png",
          p_pwy_selected, width = 12, height = 4.5)
write.csv(df_pwy, file.path(results_dir,
                            "rnaval_condition_vs_T2D_pathway_activity.csv"),
          row.names = FALSE)

df_hm <- run_hallmark_gsea(stat_mat, n_permutations = 10000)
p_hm <- plot_hallmark_gsea(df_hm)
save_plot("rnaval_hallmark_condition_vs_T2D.png", p_hm,
          width = 14, height = 12)

hallmark_pathways_of_interest <- c(
  "P53 pathway",
  "Interferon alpha response",
  "Interferon gamma response",
  "Epithelial mesenchymal transition",
  "Cholesterol homeostasis"
)
df_hm_selected <- df_hm %>%
  filter(Hallmark %in% hallmark_pathways_of_interest) %>%
  mutate(Hallmark = factor(Hallmark, levels = hallmark_pathways_of_interest),
         condition = factor(condition, levels = contrast_conditions))

p_hm_selected <- ggplot(df_hm_selected,
                        aes(x = condition, y = NES, fill = NES)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65") +
  geom_col(color = "black", width = 0.72) +
  facet_wrap(~ Hallmark, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = condition_plot_labels) +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = 0) +
  labs(title = "Selected Hallmark NES by condition vs T2D",
       x = "Condition", y = "Normalized enrichment score", fill = "NES") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_plot("rnaval_selected_Hallmark_NES_vs_T2D.png",
          p_hm_selected, width = 14, height = 4.5)
write.csv(df_hm, file.path(results_dir,
                           "rnaval_hallmark_condition_vs_T2D.csv"),
          row.names = FALSE)

df_volcano <- bind_rows(lapply(names(contrasts_list), function(nm) {
  contrasts_list[[nm]] %>%
    dplyr::select(gene_name, log2FoldChange, padj) %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(contrast = nm)
})) %>%
  mutate(
    direction = case_when(
      padj <= 0.05 & log2FoldChange > 1 ~ "Up",
      padj <= 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "ns"
    ),
    contrast = factor(contrast, levels = names(contrast_plot_labels),
                      labels = unname(contrast_plot_labels))
  )

df_top_labels <- df_volcano %>%
  filter(direction != "ns") %>%
  group_by(contrast) %>%
  slice_min(padj, n = 10, with_ties = FALSE) %>%
  ungroup()

p_volcano <- ggplot(df_volcano, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = direction), alpha = 0.55, size = 1.2) +
  scale_color_manual(values = c(Up = "indianred",
                                Down = "steelblue",
                                ns = "grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey40") +
  ggrepel::geom_text_repel(data = df_top_labels,
                           aes(label = gene_name),
                           size = 3, max.overlaps = 20,
                           box.padding = 0.3, segment.size = 0.3) +
  facet_wrap(~ contrast, ncol = 2, scales = "free") +
  labs(title = "Differential expression by condition vs T2D",
       x = "Log2 fold change", y = "-Log10 adjusted p-value",
       color = "Direction") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

save_plot("rnaval_volcano_condition_vs_T2D.png", p_volcano,
          width = 12, height = 10)
saveRDS(contrasts_list, file.path(results_dir,
                                  "rnaval_deseq2_condition_vs_T2D_results.rds"))

## Rebuild RNAval Wm in LIV2Trans format by averaging replicates by condition,----------------
## augment the saved MPS W_invitro basis from the Govaere/Kostrzewski run, then
## backproject human data and evaluate the already-trained PLSR model.
tt <- readRDS(file.path(results_dir,
                        "liv2trans_results_Govaere_Kostrzewski.rds"))
processed <- readRDS(file.path(results_dir,
                               "processed_data_list_govaere_kostrzewski.rds"))

Xh <- processed$Xh
Yh <- processed$Yh
model_trained <- tt$model
Wh <- tt$Wh
Wm_mps <- processed$Wm

Bh <- t(tt$model@weightMN) %*% tt$model@coefficientMN
Wh <- tt$Wh
phi <- Wh %*% Bh


genes_for_augmented_basis <- Reduce(intersect, list(
  colnames(Xh),
  rownames(Wh),
  rownames(Wm_mps),
  rownames(count_matrix)
))

rnaval_liv2trans <- build_rnaval_liv2trans_basis(
  count_matrix = count_matrix,
  metadata = metadata_template,
  genes_use = genes_for_augmented_basis
)

Wm_rnaval <- rnaval_liv2trans$Wm
# WPC_all_rnaval <- prcomp(aligned_raw$X,center = TRUE,scale. = FALSE)$rotation
Wm_mps_common <- Wm_mps[genes_for_augmented_basis, , drop = FALSE]
Wm_rnaval_common <- Wm_rnaval[genes_for_augmented_basis, , drop = FALSE]
# Wm_mps_common <- prcomp(processed$Xm[, genes_for_augmented_basis], center = TRUE, scale. = FALSE)$rotation
phi_common <- phi[genes_for_augmented_basis,]

Wm_mps_orth <- orthonormalize_basis(Wm_mps_common, prefix = "MPS_PC")
# Wm_rnaval_orth <- orthonormalize_basis(Wm_rnaval_common, prefix = "RNAval_PC")
# Wm_augmented <- orthonormalize_basis(cbind(Wm_mps_common, Wm_rnaval_common),
#                                      prefix = "MPS_RNAval_PC")
# Wm_augmented <- combine_orthonormal(Wm_mps_orth, Wm_rnaval_common)
augmentation_results <- augment_by_regression(Wm_mps_orth, Wm_rnaval_common)
Wm_augmented <- augmentation_results$combined
dim(Wm_augmented)                      # 10 x 7  (3 from W1 + 4 genuinely new directions)
max(abs(crossprod(Wm_augmented) - diag(ncol(Wm_augmented))))

yhat_full <- predict(model_trained, Xh)
# yhat_mps <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh, Wm_mps_orth)
# yhat_rnaval <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh, Wm_rnaval_orth)
# yhat_augmented <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh,
#                                              Wm_augmented)
yhat_mps <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh, Wm_mps)
yhat_mps_ortho <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh, Wm_mps_orth)
yhat_rnaval <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh, Wm_rnaval_common)
yhat_augmented <- predict_backprojected_plsr(Xh, Yh, model_trained, Wh,
                                             Wm_augmented)

pred_df <- bind_rows(
  prediction_long(Yh, yhat_full, "Full trained PLSR"),
  prediction_long(Yh, yhat_mps, "MPS Wm"),
  prediction_long(Yh, yhat_mps_ortho, "MPS Wm re-othonormalized"),
  prediction_long(Yh, yhat_rnaval, "RNAval Wm"),
  prediction_long(Yh, yhat_augmented, "MPS + RNAval Wm")
) %>%
  mutate(input = factor(input, levels = c("Full trained PLSR", "MPS Wm","MPS Wm re-othonormalized",
                                          "RNAval Wm", "MPS + RNAval Wm")),
         phenotype = ifelse(phenotype == "fibrosis", "Fibrosis", phenotype))

metrics_df <- prediction_metrics(pred_df)
metrics_summary <- metrics_df %>%
  group_by(input) %>%
  summarise(across(c(r, rho, MAE, RMSE), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

p_backprojection_scatter <- plot_backprojection_scatter(pred_df)
p_backprojection_metrics <- plot_backprojection_metrics(metrics_df)
print(p_backprojection_metrics)

save_plot("rnaval_augmented_Wm_backprojection_scatter.png",
          p_backprojection_scatter, width = 13.5, height = 6.8)
save_plot("rnaval_augmented_Wm_backprojection_performance.png",
          p_backprojection_metrics, width = 9, height = 4.8)

write.csv(metrics_df, file.path(results_dir,
                                "rnaval_augmented_Wm_backprojection_metrics.csv"),
          row.names = FALSE)
write.csv(metrics_summary, file.path(results_dir,
                                     "rnaval_augmented_Wm_backprojection_summary.csv"),
          row.names = FALSE)
saveRDS(
  list(
    Wm_mps = Wm_mps_orth,
    Wm_rnaval = Wm_rnaval_orth,
    Wm_augmented = Wm_augmented,
    Xm_rnaval = rnaval_liv2trans$Xm,
    predictions = pred_df,
    metrics = metrics_df,
    metrics_summary = metrics_summary
  ),
  file.path(results_dir, "rnaval_augmented_liv2trans_Wm_backprojection.rds")
)

message("Focused RNAseq validation completed.")
