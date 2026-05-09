library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(LIV2Trans)
library(matrixStats)

load_datasets <- function(dataset_names, dir_data){
  files <- dir(dir_data)
  files <- files[grepl("dataset.RData", files)]
  data_list <- list(NULL)
  if (length(files) > 0){
    for (name in dataset_names){
      if (any(grepl(name, files, ignore.case = TRUE))){
        load(file = paste0(dir_data, files[grep(name, files, ignore.case = TRUE)]))
        print(paste0("Loading file: ", paste0(dir_data, files[grep(name, files, ignore.case = TRUE)])))
        data_list[[name]]$counts <- data
        data_list[[name]]$metadata <- metadata
        data_list[[name]]$genes <- rownames(data)
        if (exists("exp_factors")){
          data_list[[name]]$exp_factors <- exp_factors
        }
      }
    }
    data_list[[1]] <- NULL
  }
  return(data_list)
}

# ============================================================
# CONFIGURATION
# ============================================================
N_REPS        <- 10
RUN_GSEA      <- TRUE
FIGURES_DIR   <- "../figures/"
RESULTS_DIR   <- "../results/"
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)


dataset_names  <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset    <- "Govaere"
target_dataset <- "Kostrzewski"

# ============================================================
# PART 1: LOAD REAL DATA AND RUN LIV2TRANS
# ============================================================
data_list <- load_datasets(dataset_names, dir_data = '../data/')
tmp       <- process_datasets(data_list, filter_variance = FALSE)
data_list <- tmp$data_list

Yh <- as.matrix(data_list[[ref_dataset]]$metadata %>%
                  dplyr::select(nas_score, Fibrosis_stage))
colnames(Yh) <- c('NAS', 'fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

counts_m  <- data_list[[target_dataset]]$counts
gm0  <- rownames(counts_m)
ngm0 <- length(gm0)
sm0  <- colnames(counts_m)
nsm0 <- length(sm0)

liv2trans_results <- liv2trans_run(Xh, Yh, Xm, Wm)
Wm_opt <- liv2trans_results$W_opt
rownames(Wm_opt) <- rownames(Wm)
colnames(Wm_opt) <- c("LV_extra1", "LV_extra2")

# ---------- NEW: real translatable components ----------
W_TC <- as.matrix(liv2trans_results$W_translatable)
rownames(W_TC) <- rownames(Wm)
colnames(W_TC) <- paste0("TC", seq_len(ncol(W_TC)))
n_TC <- ncol(W_TC)
cat("Number of real TCs:", n_TC, "\n")
# -------------------------------------------------------

# Compute phi_PLS from the real PLSR model (for rho calculations)
model_real <- liv2trans_results$model
Wh_real    <- liv2trans_results$Wh
Bh_real    <- t(model_real@weightMN) %*% model_real@coefficientMN
phi_real   <- Wh_real %*% Bh_real  # p x q

# Real rho per phenotype: rho_i = ||Wm^T phi_i||^2 / ||phi_i||^2
rho_real <- sapply(1:ncol(phi_real), function(i) {
  phi_i <- phi_real[, i]
  genes_common <- intersect(names(phi_i), rownames(Wm))
  proj  <- t(Wm[genes_common, ]) %*% phi_i[genes_common]
  sum(proj^2) / sum(phi_i^2)
})
names(rho_real) <- colnames(Yh)
cat("Real rho (fraction of phi captured by MPS PC space):\n")
print(rho_real)

# ---------- NEW: real rho through TC space ----------
rho_TC_real <- sapply(1:ncol(phi_real), function(i) {
  phi_i <- phi_real[, i]
  genes_common <- intersect(names(phi_i), rownames(W_TC))
  proj  <- t(W_TC[genes_common, ]) %*% phi_i[genes_common]
  sum(proj^2) / sum(phi_i^2)
})
names(rho_TC_real) <- colnames(Yh)
cat("Real rho through TCs (fraction of phi captured by TC space):\n")
print(rho_TC_real)
# ----------------------------------------------------

# Theoretical random baselines: k/p
k_mps   <- ncol(Wm)
p_genes <- nrow(Wm)
rho_theoretical    <- k_mps / p_genes
rho_TC_theoretical <- n_TC  / p_genes
cat("Theoretical E[rho] = k/p =", rho_theoretical, "\n")
cat("Theoretical E[rho_TC] = n_TC/p =", rho_TC_theoretical, "\n")

# Store a copy of original data_list for re-use
data_list_orig <- load_datasets(dataset_names, dir_data = '../data/')

# ============================================================
# PART 2: HELPER FUNCTIONS  (unchanged)
# ============================================================

generate_single_baseline <- function(type, n, m, real_counts = NULL) {
  switch(type,
         "rnaseq_synthetic" = matrix(
           rnbinom(n * m,
                   mu   = sample(5:500, m, replace = TRUE),
                   size = sample(seq(0.1, 2, by = 0.1), m, replace = TRUE)),
           nrow = n, ncol = m),
         "nb_low_disp"  = matrix(rnbinom(n * m, mu = 100,  size = 10),  nrow = n, ncol = m),
         "nb_high_disp" = matrix(rnbinom(n * m, mu = 100,  size = 0.5), nrow = n, ncol = m),
         "nb_low_mean"  = matrix(rnbinom(n * m, mu = 10,   size = 1),   nrow = n, ncol = m),
         "nb_high_mean" = matrix(rnbinom(n * m, mu = 1000, size = 2),   nrow = n, ncol = m),
         "poisson" = matrix(rpois(n * m, lambda = sample(5:500, m, replace = TRUE)),
                            nrow = n, ncol = m),
         "zip" = matrix(
           rbinom(n * m, 1, prob = 0.7) * rpois(n * m, lambda = 50),
           nrow = n, ncol = m),
         "zinb" = matrix(
           rbinom(n * m, 1, prob = 0.6) * rnbinom(n * m, mu = 100, size = 1),
           nrow = n, ncol = m),
         "geometric" = matrix(
           rgeom(n * m, prob = sample(seq(0.01, 0.2, by = 0.01), m, replace = TRUE)),
           nrow = n, ncol = m),
         "uniform_int" = matrix(sample(0:1000, n * m, replace = TRUE), nrow = n, ncol = m),
         "from_real_params" = {
           mu_est   <- pmax(colMeans(real_counts, na.rm = TRUE), 0.01)
           gene_vars <- apply(real_counts, 2, var, na.rm = TRUE)
           size_est <- ifelse(gene_vars > mu_est,
                              mu_est^2 / (gene_vars - mu_est),
                              1e6)
           size_est <- pmax(size_est, 0.01)
           mu_est[!is.finite(mu_est)]     <- 0.01
           size_est[!is.finite(size_est)] <- 1
           sapply(1:m, function(j) rnbinom(n, mu = mu_est[j], size = size_est[j]))
         },
         "permuted" = matrix(sample(as.vector(real_counts)), nrow = n, ncol = m),
         stop(paste("Unknown baseline type:", type))
  )
}

tanimoto <- function(a, b) {
  a <- unique(a); b <- unique(b)
  union_size <- length(union(a, b))
  if (union_size == 0) return(NA_real_)
  length(intersect(a, b)) / union_size
}

get_top_features <- function(x, k, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  x <- sort(x, decreasing = direction == "positive")
  if (direction == "positive") x <- x[x > 0] else x <- x[x < 0]
  if (length(x) < k) return(character(0))
  names(x)[seq_len(k)]
}

top_k_values <- c(5, 10, 30, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

# Generalized: works for both LVs and TCs (any column count)
calc_tanimoto_for_column <- function(col_index, W_old, W_new,
                                     label_prefix = "Column") {
  old_scores <- W_old[, col_index]
  new_scores <- W_new[, col_index]
  names(old_scores) <- rownames(W_old)
  names(new_scores) <- rownames(W_new)
  
  bind_rows(
    lapply(top_k_values, function(k) {
      old_pos <- get_top_features(old_scores, k, "positive")
      new_pos <- get_top_features(new_scores, k, "positive")
      data.frame(column = paste0(label_prefix, " ", col_index), k = k,
                 direction = "Positive",
                 tanimoto = if (length(old_pos) == 0 || length(new_pos) == 0)
                   NA_real_ else tanimoto(old_pos, new_pos))
    }),
    lapply(top_k_values, function(k) {
      old_neg <- get_top_features(old_scores, k, "negative")
      new_neg <- get_top_features(new_scores, k, "negative")
      data.frame(column = paste0(label_prefix, " ", col_index), k = k,
                 direction = "Negative",
                 tanimoto = if (length(old_neg) == 0 || length(new_neg) == 0)
                   NA_real_ else tanimoto(old_neg, new_neg))
    })
  )
}

gsea_thresholds <- c(30, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

calc_gsea_distance_for_column <- function(col_index, W_old, W_new,
                                          thresholds = gsea_thresholds,
                                          label_prefix = "Column") {
  num_table <- cbind(
    W_opt     = W_old[, col_index],
    W_opt_new = W_new[, col_index]
  )
  rownames(num_table) <- rownames(W_old)
  
  bind_rows(lapply(thresholds, function(thres) {
    dist_mat <- distance_scores(
      num_table = num_table,
      threshold_count = thres,
      names = colnames(num_table)
    )
    data.frame(column = paste0(label_prefix, " ", col_index),
               threshold = thres,
               gsea_distance = dist_mat["W_opt", "W_opt_new"])
  }))
}

library(GeneExpressionSignature)
library(Biobase)
source("distance_scores.R")

# ============================================================
# PART 3: MAIN LOOP — 10 REPLICATES PER BASELINE TYPE
# ============================================================

baseline_types <- c("rnaseq_synthetic", "nb_low_disp", "nb_high_disp",
                    "nb_low_mean", "nb_high_mean", "poisson", "zip",
                    "zinb", "geometric", "uniform_int",
                    "from_real_params", "permuted")

# Containers for extra LVs
all_tanimoto   <- list()
all_gsea       <- list()
all_diag_cor   <- list()
all_spear_true <- list()
all_rho        <- list()

# ---------- NEW: Containers for TCs ----------
all_tanimoto_TC   <- list()
all_gsea_TC       <- list()
all_diag_cor_TC   <- list()
all_spear_true_TC <- list()
all_rho_TC        <- list()
# ---------------------------------------------

total_iters <- length(baseline_types) * N_REPS
iter <- 0

for (bl_name in baseline_types) {
  cat("\n========================================\n")
  cat("Baseline:", bl_name, "\n")
  cat("========================================\n")
  
  for (rep_i in 1:N_REPS) {
    iter <- iter + 1
    cat(sprintf("  [%d/%d] %s rep %d/%d ... ",
                iter, total_iters, bl_name, rep_i, N_REPS))
    
    tryCatch({
      # Generate fresh random counts
      random_counts <- generate_single_baseline(bl_name, nsm0, ngm0, counts_m)
      new_counts_m <- t(random_counts)
      colnames(new_counts_m) <- colnames(counts_m)
      rownames(new_counts_m) <- rownames(counts_m)
      
      # Re-process with random MPS counts
      data_list_new <- data_list_orig
      data_list_new[[target_dataset]]$counts <- new_counts_m
      tmp_new <- process_datasets(data_list_new, filter_variance = FALSE)
      data_list_new <- tmp_new$data_list
      
      Xm_new <- data_list_new[[target_dataset]]$data_center %>% t()
      Wm_new <- data_list_new[[target_dataset]]$Wm_group %>% as.matrix()
      
      # Run LIV2TRANS with random Wm
      liv2trans_new <- liv2trans_run(Xh, Yh, Xm_new, Wm_new)
      
      # --- Extra LVs (random) ---
      Wm_opt_new <- liv2trans_new$W_opt
      rownames(Wm_opt_new) <- rownames(Wm_new)
      colnames(Wm_opt_new) <- c("new_LV_extra1", "new_LV_extra2")
      
      # --- NEW: TCs (random) ---
      W_TC_new <- as.matrix(liv2trans_new$W_translatable)
      rownames(W_TC_new) <- rownames(Wm_new)
      colnames(W_TC_new) <- paste0("new_TC", seq_len(ncol(W_TC_new)))
      
      # ============================================================
      # METRICS FOR EXTRA LVs (existing)
      # ============================================================
      
      # Metric 1: Tanimoto
      tan_res <- bind_rows(
        calc_tanimoto_for_column(1, Wm_opt, Wm_opt_new, "Column"),
        calc_tanimoto_for_column(2, Wm_opt, Wm_opt_new, "Column")
      )
      tan_res$baseline  <- bl_name
      tan_res$replicate <- rep_i
      all_tanimoto[[length(all_tanimoto) + 1]] <- tan_res
      
      # Metric 2: GSEA distance
      if (RUN_GSEA) {
        gsea_res <- tryCatch({
          bind_rows(
            calc_gsea_distance_for_column(1, Wm_opt, Wm_opt_new, label_prefix = "Column"),
            calc_gsea_distance_for_column(2, Wm_opt, Wm_opt_new, label_prefix = "Column")
          )
        }, error = function(e) NULL)
        
        if (!is.null(gsea_res)) {
          gsea_res$baseline  <- bl_name
          gsea_res$replicate <- rep_i
          all_gsea[[length(all_gsea) + 1]] <- gsea_res
        }
      }
      
      # Metric 3: Diagonal cor(Yhat_real, Yhat_rand) using full Wm
      Yhat_real <- predict(model_real, Xh %*% Wm %*% t(Wm))
      Yhat_rand <- predict(model_real, Xh %*% Wm_new %*% t(Wm_new))
      cor_mat   <- cor(Yhat_real, Yhat_rand, method = 'spearman')
      diag_vals <- diag(cor_mat)
      
      all_diag_cor[[length(all_diag_cor) + 1]] <- data.frame(
        baseline  = bl_name,
        replicate = rep_i,
        phenotype = colnames(Yh),
        cor_value = diag_vals
      )
      
      # Metric 3b: cor(Yh, Yhat_rand)
      cor_true_rand <- cor(Yh, Yhat_rand, method = 'spearman')
      diag_true     <- diag(cor_true_rand)
      
      all_spear_true[[length(all_spear_true) + 1]] <- data.frame(
        baseline  = bl_name,
        replicate = rep_i,
        phenotype = colnames(Yh),
        cor_value = diag_true
      )
      
      # Metric 4: rho via Wm
      rho_rand <- sapply(1:ncol(phi_real), function(i) {
        phi_i <- phi_real[, i]
        genes_common <- intersect(names(phi_i), rownames(Wm_new))
        proj  <- t(Wm_new[genes_common, ]) %*% phi_i[genes_common]
        sum(proj^2) / sum(phi_i^2)
      })
      
      all_rho[[length(all_rho) + 1]] <- data.frame(
        baseline  = bl_name,
        replicate = rep_i,
        phenotype = colnames(Yh),
        rho_value = rho_rand,
        k_random  = ncol(Wm_new)
      )
      
      # ============================================================
      # NEW: METRICS FOR TRANSLATABLE COMPONENTS (TCs)
      # ============================================================
      
      # Match column count between real and random TCs (use min)
      n_TC_common <- min(ncol(W_TC), ncol(W_TC_new))
      
      if (n_TC_common > 0) {
        # Metric 1 (TC): Tanimoto similarity
        tan_res_TC <- bind_rows(lapply(seq_len(n_TC_common), function(j) {
          calc_tanimoto_for_column(j, W_TC, W_TC_new, label_prefix = "TC")
        }))
        tan_res_TC$baseline  <- bl_name
        tan_res_TC$replicate <- rep_i
        all_tanimoto_TC[[length(all_tanimoto_TC) + 1]] <- tan_res_TC
        
        # Metric 2 (TC): GSEA distance
        if (RUN_GSEA) {
          gsea_res_TC <- tryCatch({
            bind_rows(lapply(seq_len(n_TC_common), function(j) {
              calc_gsea_distance_for_column(j, W_TC, W_TC_new, label_prefix = "TC")
            }))
          }, error = function(e) NULL)
          
          if (!is.null(gsea_res_TC)) {
            gsea_res_TC$baseline  <- bl_name
            gsea_res_TC$replicate <- rep_i
            all_gsea_TC[[length(all_gsea_TC) + 1]] <- gsea_res_TC
          }
        }
      }
      
      # Metric 3 (TC): Diagonal cor(Yhat_realTC, Yhat_randTC)
      # Backproject Xh through TC space (W_TC is orthonormal in gene space)
      Yhat_real_TC <- predict(model_real, Xh %*% W_TC     %*% t(W_TC))
      Yhat_rand_TC <- predict(model_real, Xh %*% W_TC_new %*% t(W_TC_new))
      cor_mat_TC   <- cor(Yhat_real_TC, Yhat_rand_TC, method = 'spearman')
      diag_vals_TC <- diag(cor_mat_TC)
      
      all_diag_cor_TC[[length(all_diag_cor_TC) + 1]] <- data.frame(
        baseline  = bl_name,
        replicate = rep_i,
        phenotype = colnames(Yh),
        cor_value = diag_vals_TC
      )
      
      # Metric 3b (TC): cor(Yh, Yhat_randTC)
      cor_true_rand_TC <- cor(Yh, Yhat_rand_TC, method = 'spearman')
      diag_true_TC     <- diag(cor_true_rand_TC)
      
      all_spear_true_TC[[length(all_spear_true_TC) + 1]] <- data.frame(
        baseline  = bl_name,
        replicate = rep_i,
        phenotype = colnames(Yh),
        cor_value = diag_true_TC
      )
      
      # Metric 4 (TC): rho via random TCs
      rho_rand_TC <- sapply(1:ncol(phi_real), function(i) {
        phi_i <- phi_real[, i]
        genes_common <- intersect(names(phi_i), rownames(W_TC_new))
        proj  <- t(W_TC_new[genes_common, ]) %*% phi_i[genes_common]
        sum(proj^2) / sum(phi_i^2)
      })
      
      all_rho_TC[[length(all_rho_TC) + 1]] <- data.frame(
        baseline  = bl_name,
        replicate = rep_i,
        phenotype = colnames(Yh),
        rho_value = rho_rand_TC,
        k_random  = ncol(W_TC_new)
      )
      
      cat("done\n")
      
    }, error = function(e) {
      cat("FAILED:", e$message, "\n")
    })
  }
}

# Combine all results into data frames
df_tanimoto   <- bind_rows(all_tanimoto)
df_gsea       <- bind_rows(all_gsea)
df_diag_cor   <- bind_rows(all_diag_cor)
df_spear_true <- bind_rows(all_spear_true)
df_rho        <- bind_rows(all_rho)

# NEW: TC dataframes
df_tanimoto_TC   <- bind_rows(all_tanimoto_TC)
df_gsea_TC       <- bind_rows(all_gsea_TC)
df_diag_cor_TC   <- bind_rows(all_diag_cor_TC)
df_spear_true_TC <- bind_rows(all_spear_true_TC)
df_rho_TC        <- bind_rows(all_rho_TC)

# Save raw results
save(df_tanimoto, df_gsea, df_diag_cor, df_spear_true, df_rho,
     df_tanimoto_TC, df_gsea_TC, df_diag_cor_TC, df_spear_true_TC, df_rho_TC,
     rho_real, rho_theoretical, rho_TC_real, rho_TC_theoretical,
     Wm_opt, W_TC,
     file = file.path(RESULTS_DIR, "specificity_v3_results.RData"))

# ============================================================
# PART 4: PLOTS — ONE SET PER BASELINE TYPE (LVs + TCs)
# ============================================================

# Helper: reusable plotting routines (so we don't duplicate code for LVs and TCs)

# Some result tables can be empty if every replicate failed for a metric
# or if GSEA fails silently. These guards prevent downstream plotting errors.
has_plot_data <- function(df, required_cols = character(0)) {
  if (is.null(df)) return(FALSE)
  if (!is.data.frame(df)) return(FALSE)
  if (nrow(df) == 0) return(FALSE)
  if (!all(required_cols %in% colnames(df))) return(FALSE)
  TRUE
}

plot_tanimoto <- function(df, bl_name, label, fname_suffix) {
  if (!has_plot_data(df, c("baseline", "column", "k", "direction", "tanimoto"))) return(invisible(NULL))
  df_bl <- df %>% filter(baseline == bl_name)
  if (nrow(df_bl) == 0) return(invisible(NULL))
  s <- df_bl %>%
    group_by(column, k, direction) %>%
    summarise(
      mean_tan = mean(tanimoto, na.rm = TRUE),
      se_tan   = sd(tanimoto, na.rm = TRUE) / sqrt(sum(!is.na(tanimoto))),
      .groups  = "drop"
    )
  p <- ggplot(s, aes(x = k, y = mean_tan, color = direction, group = direction)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_tan - se_tan, ymax = mean_tan + se_tan),
                  width = 20, linewidth = 0.5) +
    facet_wrap(~ column) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Top k features",
      y = "Tanimoto similarity (mean +/- SEM)",
      color = "Feature set",
      title = paste0("Tanimoto: real vs random ", label, " [", bl_name, "]"),
      subtitle = paste0("N = ", N_REPS, " replicates")
    ) +
    theme_bw(base_size = 12)
  ggsave(paste0(FIGURES_DIR, "tanimoto_", fname_suffix, "_", bl_name, ".png"),
         plot = p, width = 10, height = 5, dpi = 300)
}

plot_gsea <- function(df, bl_name, label, fname_suffix) {
  if (!RUN_GSEA) return(invisible(NULL))
  if (!has_plot_data(df, c("baseline", "column", "threshold", "gsea_distance"))) return(invisible(NULL))
  df_bl <- df %>% filter(baseline == bl_name)
  if (nrow(df_bl) == 0) return(invisible(NULL))
  s <- df_bl %>%
    group_by(column, threshold) %>%
    summarise(
      mean_dist = mean(gsea_distance, na.rm = TRUE),
      se_dist   = sd(gsea_distance, na.rm = TRUE) / sqrt(sum(!is.na(gsea_distance))),
      .groups   = "drop"
    )
  p <- ggplot(s, aes(x = threshold, y = mean_dist, group = column, color = column)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_dist - se_dist, ymax = mean_dist + se_dist),
                  width = 20, linewidth = 0.5) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(limits = c(0, 2)) +
    labs(
      x = "Threshold count",
      y = "GSEA-based distance (mean +/- SEM)",
      color = "Component",
      title = paste0("GSEA distance: real vs random ", label, " [", bl_name, "]"),
      subtitle = paste0("N = ", N_REPS, " replicates")
    ) +
    theme_bw(base_size = 12)
  ggsave(paste0(FIGURES_DIR, "gsea_distance_", fname_suffix, "_", bl_name, ".png"),
         plot = p, width = 8, height = 5, dpi = 300)
}

plot_diag_cor <- function(df, bl_name, label, fname_suffix) {
  if (!has_plot_data(df, c("baseline", "phenotype", "cor_value"))) return(invisible(NULL))
  df_bl <- df %>% filter(baseline == bl_name)
  if (nrow(df_bl) == 0) return(invisible(NULL))
  p <- ggplot(df_bl, aes(x = cor_value, fill = phenotype)) +
    geom_histogram(bins = 15, alpha = 0.6, position = "identity",
                   color = "black", linewidth = 0.3) +
    facet_wrap(~ phenotype, scales = "free_x") +
    labs(
      x = "Spearman correlation (diagonal)",
      y = "Count",
      fill = "Phenotype",
      title = paste0("cor(Yhat_real, Yhat_random) using ", label, " [", bl_name, "]"),
      subtitle = paste0("N = ", N_REPS, " replicates")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave(paste0(FIGURES_DIR, "diag_cor_histogram_", fname_suffix, "_", bl_name, ".png"),
         plot = p, width = 8, height = 4, dpi = 300)
}

plot_spear_true <- function(df, bl_name, W_real, label, fname_suffix) {
  if (!has_plot_data(df, c("baseline", "phenotype", "cor_value"))) return(invisible(NULL))
  df_bl <- df %>% filter(baseline == bl_name)
  if (nrow(df_bl) == 0) return(invisible(NULL))
  Yhat_real_ref  <- predict(model_real, Xh %*% W_real %*% t(W_real))
  cor_real_ref   <- diag(cor(Yh, Yhat_real_ref, method = 'spearman'))
  ref_real_lines <- data.frame(
    phenotype = colnames(Yh),
    value     = cor_real_ref,
    label     = paste0("Real (", round(cor_real_ref, 3), ")")
  )
  p <- ggplot(df_bl, aes(x = cor_value, fill = phenotype)) +
    geom_histogram(bins = 15, alpha = 0.6, position = "identity",
                   color = "black", linewidth = 0.3) +
    geom_vline(data = ref_real_lines, aes(xintercept = value),
               linetype = "dashed", color = "red", linewidth = 1) +
    facet_wrap(~ phenotype, scales = "free_x") +
    labs(
      x = "Spearman correlation (diagonal)",
      y = "Count",
      fill = "Phenotype",
      title = paste0("cor(Yh, Yhat_random) using ", label, " [", bl_name, "]"),
      subtitle = paste0("N = ", N_REPS,
                        " replicates; dashed red = real ", label, " backprojection")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave(paste0(FIGURES_DIR, "spear_true_vs_rand_histogram_", fname_suffix, "_", bl_name, ".png"),
         plot = p, width = 8, height = 4, dpi = 300)
}

plot_rho <- function(df, bl_name, rho_real_ref, rho_theo_ref_fn,
                     label, fname_suffix) {
  if (!has_plot_data(df, c("baseline", "phenotype", "rho_value", "k_random"))) return(invisible(NULL))
  df_bl <- df %>% filter(baseline == bl_name)
  if (nrow(df_bl) == 0) return(invisible(NULL))
  
  k_rand_median <- median(df_bl$k_random)
  rho_theo_rand <- rho_theo_ref_fn(k_rand_median)
  
  ref_lines <- data.frame(
    phenotype = rep(colnames(Yh), 2),
    value     = c(rho_real_ref, rep(rho_theo_rand, ncol(Yh))),
    label     = c(paste0("Real (", round(rho_real_ref, 5), ")"),
                  rep(paste0("E[rho]=k/p (", round(rho_theo_rand, 5), ")"), ncol(Yh))),
    linetype  = c(rep("Real", ncol(Yh)),
                  rep("Theoretical k/p", ncol(Yh)))
  )
  
  p <- ggplot(df_bl, aes(x = rho_value)) +
    geom_histogram(bins = 15, fill = "steelblue", alpha = 0.7,
                   color = "black", linewidth = 0.3) +
    geom_vline(data = ref_lines,
               aes(xintercept = value, linetype = linetype, color = linetype),
               linewidth = 1) +
    scale_linetype_manual(values = c("Real" = "dashed",
                                     "Theoretical k/p" = "dotted")) +
    scale_color_manual(values = c("Real" = "red",
                                  "Theoretical k/p" = "darkgreen")) +
    facet_wrap(~ phenotype, scales = "free_x") +
    labs(
      x        = expression(rho ~ "(fraction of " ~ phi[PLS] ~ " captured)"),
      y        = "Count",
      linetype = "Reference",
      color    = "Reference",
      title    = paste0("Fraction of phi captured by random ", label,
                        " [", bl_name, "]"),
      subtitle = paste0("N = ", N_REPS, " replicates; k_random = ",
                        k_rand_median, ", p = ", p_genes)
    ) +
    theme_bw(base_size = 12)
  
  ggsave(paste0(FIGURES_DIR, "rho_histogram_", fname_suffix, "_", bl_name, ".png"),
         plot = p, width = 10, height = 5, dpi = 300)
}

# ============================================================
# Loop over baselines and produce both LV and TC plots
# ============================================================
for (bl_name in baseline_types) {
  cat("Plotting:", bl_name, "\n")
  
  # ---------- Extra LV plots ----------
  plot_tanimoto(df_tanimoto,   bl_name, "extra LVs", "LV")
  plot_gsea(    df_gsea,       bl_name, "extra LVs", "LV")
  plot_diag_cor(df_diag_cor,   bl_name, "extra LVs (Wm)", "LV")
  plot_spear_true(df_spear_true, bl_name, Wm, "extra LVs (Wm)", "LV")
  plot_rho(df_rho, bl_name,
           rho_real_ref     = rho_real,
           rho_theo_ref_fn  = function(k) k / p_genes,
           label            = "Wm",
           fname_suffix     = "LV")
  
  # ---------- TC plots ----------
  plot_tanimoto(df_tanimoto_TC,   bl_name, "TCs", "TC")
  plot_gsea(    df_gsea_TC,       bl_name, "TCs", "TC")
  plot_diag_cor(df_diag_cor_TC,   bl_name, "TCs", "TC")
  plot_spear_true(df_spear_true_TC, bl_name, W_TC, "TCs", "TC")
  plot_rho(df_rho_TC, bl_name,
           rho_real_ref     = rho_TC_real,
           rho_theo_ref_fn  = function(k) k / p_genes,
           label            = "TCs",
           fname_suffix     = "TC")
}

cat("\n============================================\n")
cat("All plots saved to:", FIGURES_DIR, "\n")
cat("Results saved to: ", file.path(RESULTS_DIR, "specificity_v3_results.RData"), "\n", sep = "")
cat("============================================\n")