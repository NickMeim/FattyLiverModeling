# Concentration-metric test for LIV2TRANS revision (Reviewer 1)
# For each Wm (real, generative-baseline, random orthonormal),
# compute per-PC alignment with phi_PLS and derive:
#   - rho_total      : fraction of phi_PLS captured by Wm
#   - effective rank : r_eff = (sum a)^2 / sum(a^2)
#   - concentration  : 1 - r_eff / k    (k/p-independent specificity)
#
# Two null distributions:
#   (a) Generative baselines (existing types)
#   (b) 1000 random orthonormal Wm of identical (p, k)

library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(LIV2Trans)
library(matrixStats)

# ---------------------------------------------------------------
# Reused loader from existing script
# ---------------------------------------------------------------
load_datasets <- function(dataset_names, dir_data){
  files <- dir(dir_data)
  files <- files[grepl("dataset.RData", files)]
  data_list <- list(NULL)
  if (length(files) > 0){
    for (name in dataset_names){
      if (any(grepl(name, files, ignore.case = TRUE))){
        load(file = paste0(dir_data, files[grep(name, files, ignore.case = TRUE)]))
        print(paste0("Loading file: ", paste0(dir_data, files[grep(name, files, ignore.case = TRUE)])))
        data_list[[name]]$counts   <- data
        data_list[[name]]$metadata <- metadata
        data_list[[name]]$genes    <- rownames(data)
        if (exists("exp_factors")){
          data_list[[name]]$exp_factors <- exp_factors
        }
      }
    }
    data_list[[1]] <- NULL
  }
  return(data_list)
}

# ===============================================================
# CONFIGURATION
# ===============================================================
N_RANDOM_ORTHO  <- 1000   # random orthonormal Wm
N_BASELINE_REPS <- 30     # repetitions per generative baseline

FIGURES_DIR <- "../figures/"
RESULTS_DIR <- "../results/"
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

dataset_names  <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset    <- "Govaere"
target_dataset <- "Kostrzewski"

clean_filename <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
dataset_tag <- paste0("ref-",     clean_filename(ref_dataset),
                      "__target-", clean_filename(target_dataset))

# ===============================================================
# PART 1: LOAD REAL DATA, RUN LIV2TRANS, COMPUTE phi_PLS
# ===============================================================
data_list <- load_datasets(dataset_names, dir_data = '../data/')
tmp       <- process_datasets(data_list, filter_variance = FALSE)
data_list <- tmp$data_list

Yh <- as.matrix(data_list[[ref_dataset]]$metadata %>%
                  dplyr::select(nas_score, Fibrosis_stage))
colnames(Yh) <- c('NAS', 'fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
Wm <- data_list[[target_dataset]]$Wm_group   %>% as.matrix()

counts_m   <- data_list[[target_dataset]]$counts
gm0  <- rownames(counts_m); ngm0 <- length(gm0)
sm0  <- colnames(counts_m); nsm0 <- length(sm0)

liv2trans_results <- liv2trans_run(Xh, Yh, Xm, Wm)
Wm_opt <- liv2trans_results$W_opt
rownames(Wm_opt) <- rownames(Wm)
colnames(Wm_opt) <- c("LV_extra1", "LV_extra2")

W_TC <- as.matrix(liv2trans_results$W_translatable)
rownames(W_TC) <- rownames(Wm)
colnames(W_TC) <- paste0("TC", seq_len(ncol(W_TC)))
n_TC <- ncol(W_TC)
cat("Number of real TCs:", n_TC, "\n")

# phi_PLS from real PLSR model
model_real <- liv2trans_results$model
Wh_real    <- liv2trans_results$Wh
Bh_real    <- t(model_real@weightMN) %*% model_real@coefficientMN
phi_real   <- Wh_real %*% Bh_real          # p x q
colnames(phi_real) <- colnames(Yh)

k_mps   <- ncol(Wm)
p_genes <- nrow(Wm)
cat("k =", k_mps, "  p =", p_genes, "  k/p =", round(k_mps/p_genes, 6), "\n")

# Cached copy of original data_list for re-use in baseline loop
data_list_orig <- load_datasets(dataset_names, dir_data = '../data/')

# ===============================================================
# HELPERS
# ===============================================================

# Per-PC alignment: alpha_{ij} = (w_i^T phi_j)^2 / ||phi_j||^2
# Returns a (k x q) matrix; columns sum to rho_j.
compute_alpha <- function(W, phi) {
  genes_common <- intersect(rownames(W), rownames(phi))
  if (length(genes_common) == 0) {
    stop("No common genes between W and phi.")
  }
  W_c   <- W[genes_common, , drop = FALSE]
  phi_c <- phi[genes_common, , drop = FALSE]
  alpha <- sapply(1:ncol(phi_c), function(j) {
    phi_j <- phi_c[, j]
    proj  <- as.vector(t(W_c) %*% phi_j)
    proj^2 / sum(phi_j^2)
  })
  if (is.null(dim(alpha))) alpha <- matrix(alpha, ncol = ncol(phi_c))
  rownames(alpha) <- if (is.null(colnames(W))) seq_len(ncol(W)) else colnames(W)
  colnames(alpha) <- colnames(phi_c)
  alpha
}

# Effective rank (participation ratio): r_eff = (sum a)^2 / sum(a^2)
# r_eff = 1 if all alignment is in one PC, r_eff = k if uniform across all PCs
effective_rank <- function(alpha_vec) {
  s1 <- sum(alpha_vec)
  s2 <- sum(alpha_vec^2)
  if (s2 < .Machine$double.eps) return(NA_real_)
  (s1^2) / s2
}

# Concentration metric: 1 - r_eff/k, ranges in (0, 1).
# Close to 1 = highly concentrated (biologically specific alignment)
# Close to 0 = uniform / dimensional artifact
concentration <- function(alpha_vec) {
  r <- effective_rank(alpha_vec)
  k <- length(alpha_vec)
  if (is.na(r) || k <= 1) return(NA_real_)
  1 - (r / k)
}

# Random orthonormal p x k matrix via QR of Gaussian random
random_orthonormal <- function(p, k, gene_names = NULL) {
  G <- matrix(rnorm(p * k), nrow = p, ncol = k)
  W <- qr.Q(qr(G))
  if (!is.null(gene_names)) rownames(W) <- gene_names
  W
}

# Empirical p-value (one-sided)
emp_pvalue <- function(null_vals, real_val, alternative = c("less", "greater")) {
  alternative <- match.arg(alternative)
  null_vals <- null_vals[is.finite(null_vals)]
  n_null <- length(null_vals)
  if (n_null == 0) return(NA_real_)
  if (alternative == "less") {
    (sum(null_vals <= real_val) + 1) / (n_null + 1)
  } else {
    (sum(null_vals >= real_val) + 1) / (n_null + 1)
  }
}

compute_pvals <- function(df_null, real_value, var_name, alternative) {
  ph <- unique(df_null$phenotype)
  bind_rows(lapply(ph, function(p) {
    null_vals <- df_null[df_null$phenotype == p, var_name, drop = TRUE]
    real_p    <- real_value[p]
    data.frame(phenotype = p,
               real      = unname(real_p),
               null_mean = mean(null_vals, na.rm = TRUE),
               null_sd   = sd(null_vals,   na.rm = TRUE),
               n_null    = sum(is.finite(null_vals)),
               p_value   = emp_pvalue(null_vals, real_p, alternative))
  }))
}

# ===============================================================
# REAL CONCENTRATION METRIC
# ===============================================================
alpha_real    <- compute_alpha(Wm, phi_real)            # k x q
r_eff_real    <- apply(alpha_real, 2, effective_rank)
conc_real     <- apply(alpha_real, 2, concentration)
rho_real_vec  <- colSums(alpha_real)
names(r_eff_real)   <- colnames(phi_real)
names(conc_real)    <- colnames(phi_real)
names(rho_real_vec) <- colnames(phi_real)

cat("\n--- REAL Wm ---\n")
cat("  rho_total      :", rho_real_vec, "\n")
cat("  effective rank :", r_eff_real, " (k =", k_mps, ")\n")
cat("  concentration  :", conc_real, "\n")

# ===============================================================
# RANDOM ORTHONORMAL Wm
# ===============================================================
cat("\n=== Random orthonormal Wm  (N =", N_RANDOM_ORTHO, ") ===\n")
results_ortho <- vector("list", N_RANDOM_ORTHO)
for (i in seq_len(N_RANDOM_ORTHO)) {
  if (i %% 100 == 0) cat("  ", i, "/", N_RANDOM_ORTHO, "\n")
  W_rand   <- random_orthonormal(p_genes, k_mps, rownames(Wm))
  alpha_r  <- compute_alpha(W_rand, phi_real)
  r_eff_r  <- apply(alpha_r, 2, effective_rank)
  conc_r   <- apply(alpha_r, 2, concentration)
  rho_r    <- colSums(alpha_r)
  results_ortho[[i]] <- data.frame(
    rep           = i,
    phenotype     = colnames(phi_real),
    r_eff         = r_eff_r,
    concentration = conc_r,
    rho_total     = rho_r,
    k_random      = k_mps,
    source        = "random_orthonormal",
    stringsAsFactors = FALSE
  )
}
df_ortho <- bind_rows(results_ortho)

# ===============================================================
# GENERATIVE BASELINES
# ===============================================================
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
         "zip" = matrix(rbinom(n * m, 1, prob = 0.7) * rpois(n * m, lambda = 50),
                        nrow = n, ncol = m),
         "zinb" = matrix(rbinom(n * m, 1, prob = 0.6) *
                           rnbinom(n * m, mu = 100, size = 1), nrow = n, ncol = m),
         "geometric" = matrix(
           rgeom(n * m, prob = sample(seq(0.01, 0.2, by = 0.01), m, replace = TRUE)),
           nrow = n, ncol = m),
         "uniform_int" = matrix(sample(0:1000, n * m, replace = TRUE),
                                nrow = n, ncol = m),
         "from_real_params" = {
           mu_est    <- pmax(colMeans(real_counts, na.rm = TRUE), 0.01)
           gene_vars <- apply(real_counts, 2, var, na.rm = TRUE)
           size_est  <- ifelse(gene_vars > mu_est,
                               mu_est^2 / (gene_vars - mu_est), 1e6)
           size_est  <- pmax(size_est, 0.01)
           mu_est[!is.finite(mu_est)]     <- 0.01
           size_est[!is.finite(size_est)] <- 1
           sapply(1:m, function(j) rnbinom(n, mu = mu_est[j], size = size_est[j]))
         },
         "permuted" = matrix(sample(as.vector(real_counts)), nrow = n, ncol = m),
         stop(paste("Unknown baseline type:", type))
  )
}

baseline_types <- c("rnaseq_synthetic", "nb_low_disp", "nb_high_disp",
                    "nb_low_mean", "nb_high_mean", "poisson", "zip",
                    "zinb", "geometric", "uniform_int",
                    "from_real_params", "permuted")

cat("\n=== Generative baselines  (N =", N_BASELINE_REPS, " per type) ===\n")
results_baseline <- list()
for (bl_name in baseline_types) {
  cat("Baseline:", bl_name, "\n")
  for (rep_i in seq_len(N_BASELINE_REPS)) {
    tryCatch({
      random_counts <- generate_single_baseline(bl_name, nsm0, ngm0, counts_m)
      new_counts_m  <- t(random_counts)
      colnames(new_counts_m) <- colnames(counts_m)
      rownames(new_counts_m) <- rownames(counts_m)
      
      data_list_new <- data_list_orig
      data_list_new[[target_dataset]]$counts <- new_counts_m
      tmp_new <- process_datasets(data_list_new, filter_variance = FALSE)
      data_list_new <- tmp_new$data_list
      Wm_new <- data_list_new[[target_dataset]]$Wm_group %>% as.matrix()
      
      alpha_b <- compute_alpha(Wm_new, phi_real)
      r_eff_b <- apply(alpha_b, 2, effective_rank)
      conc_b  <- apply(alpha_b, 2, concentration)
      rho_b   <- colSums(alpha_b)
      
      results_baseline[[length(results_baseline) + 1]] <- data.frame(
        baseline      = bl_name,
        rep           = rep_i,
        phenotype     = colnames(phi_real),
        r_eff         = r_eff_b,
        concentration = conc_b,
        rho_total     = rho_b,
        k_random      = ncol(Wm_new),
        source        = bl_name,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("  rep", rep_i, "FAILED:", e$message, "\n")
    })
  }
}
df_baseline <- bind_rows(results_baseline)

# ===============================================================
# REPORT P-VALUES
# ===============================================================
cat("\n========================================\n")
cat("vs random orthonormal Wm\n")
cat("========================================\n")
cat("\nP(r_eff_real <= r_eff_rand)\n")
print(compute_pvals(df_ortho, r_eff_real, "r_eff", "less"))
cat("\nP(concentration_real >= concentration_rand)\n")
print(compute_pvals(df_ortho, conc_real, "concentration", "greater"))
cat("\nP(rho_total_real >= rho_total_rand)\n")
print(compute_pvals(df_ortho, rho_real_vec, "rho_total", "greater"))

# ===============================================================
# SAVE
# ===============================================================
save(df_ortho, df_baseline,
     alpha_real, r_eff_real, conc_real, rho_real_vec,
     phi_real, Wm, Wm_opt, W_TC,
     k_mps, p_genes, n_TC,
     ref_dataset, target_dataset, dataset_tag,
     file = file.path(RESULTS_DIR,
                      paste0(dataset_tag,
                             "__concentration_results.RData")))

# ===============================================================
# PLOTS
# ===============================================================

# Combined comparison plots: real Wm (red dashed) vs all null models
df_combined <- bind_rows(
  df_ortho    %>% transmute(phenotype, r_eff, concentration, rho_total,
                            group = "Random orthonormal"),
  df_baseline %>% transmute(phenotype, r_eff, concentration, rho_total,
                            group = paste0("Baseline: ", baseline))
)

# Effective rank
ref_lines <- data.frame(phenotype = names(r_eff_real),
                        value     = unname(r_eff_real))
p_reff <- ggplot(df_combined, aes(x = group, y = r_eff)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  geom_hline(data = ref_lines, aes(yintercept = value),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ phenotype, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Effective rank of alignment",
       subtitle = "Red dashed = real Wm") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGURES_DIR,
                 paste0(dataset_tag, "__reff_comparison.png")),
       plot = p_reff, width = 12, height = 8, dpi = 300)

# Concentration
ref_lines2 <- data.frame(phenotype = names(conc_real),
                         value     = unname(conc_real))
p_conc <- ggplot(df_combined, aes(x = group, y = concentration)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  geom_hline(data = ref_lines2, aes(yintercept = value),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ phenotype, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Concentration  (1 - r_eff/k)",
       subtitle = "Red dashed = real Wm") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGURES_DIR,
                 paste0(dataset_tag, "__concentration_comparison.png")),
       plot = p_conc, width = 12, height = 8, dpi = 300)

cat("\n============================================\n")
cat("Plots:   ", FIGURES_DIR, "\n")
cat("Results: ",
    file.path(RESULTS_DIR,
              paste0(dataset_tag, "__concentration_results.RData")),
    "\n", sep = "")
cat("============================================\n")