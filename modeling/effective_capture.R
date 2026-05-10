# ===============================================================
#   Layer 1 (k/p-independent specificity):
#       Concentration metric (effective rank of alignment of Wm
#       columns with phi_PLS), tested against:
#         (a) all generative baselines from randomBaselines script
#         (b) 1000 random orthonormal Wm matrices
#
#   Layer 2 (permutation test on condition labels):
#       Permute the experimental-condition labels in MPS metadata,
#       re-aggregate replicates, recompute Wm, and compare the real
#       Wm against this null distribution on:
#         - rho_total (fraction of phi_PLS captured)
#         - effective rank
#         - concentration
#         - downstream predictive performance
#
#   Both layers use the same phi_PLS computed from the real PLSR
#   model, so the comparisons are apples-to-apples.
# ===============================================================

library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(LIV2Trans)
library(matrixStats)

# ---------------------------------------------------------------
# Data Loader
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
N_RANDOM_ORTHO  <- 1000   # Layer 1 (b): random orthonormal Wm
N_PERMUTATIONS  <- 1000   # Layer 2: condition-label permutations
N_BASELINE_REPS <- 30     # Layer 1 (a): generative baselines

FIGURES_DIR <- "../figures/"
RESULTS_DIR <- "../results/"
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

dataset_names  <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset    <- "Govaere"
target_dataset <- "Kostrzewski"

# ---------------------------------------------------------------
# IMPORTANT (Layer 2): which column in MPS metadata identifies
# the experimental condition. This is what gets shuffled in the
# permutation test. Adjust to match your dataset.
#
# The script tries the names below in order and uses the first
# one that exists in metadata.
# ---------------------------------------------------------------
candidate_condition_cols <- c("condition", "exp_condition", "Condition",
                              "treatment", "Treatment", "group", "Group")

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
metadata_m <- data_list[[target_dataset]]$metadata
gm0  <- rownames(counts_m); ngm0 <- length(gm0)
sm0  <- colnames(counts_m); nsm0 <- length(sm0)

# Identify the condition column in metadata
condition_col <- NULL
for (cn in candidate_condition_cols) {
  if (cn %in% colnames(metadata_m)) {
    condition_col <- cn; break
  }
}
if (is.null(condition_col)) {
  stop(sprintf(
    "Could not find a condition column in MPS metadata. Tried: %s.\nAvailable columns: %s.\nEdit `candidate_condition_cols` at the top of this script.",
    paste(candidate_condition_cols, collapse = ", "),
    paste(colnames(metadata_m), collapse = ", ")
  ))
}
cat("Using condition column for permutation test:", condition_col, "\n")

# Real LIV2TRANS run
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

# Cached copy of original data_list for re-use in baseline / permutation loops
data_list_orig <- load_datasets(dataset_names, dir_data = '../data/')

# ===============================================================
# HELPERS
# ===============================================================

# Per-PC alignment: alpha_{ij} = (w_i^T phi_j)^2 / ||phi_j||^2
# Returns a (k x q) matrix; columns sum to rho_j (total captured fraction
# for phenotype j).
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

# Aggregate (mean) gene expression Xm by group labels, then center,
# and compute its principal-component basis with k components.
# Used in the permutation test to bypass process_datasets and make
# the permutation logic fully transparent.
pca_after_aggregation <- function(X, group_labels, k) {
  if (length(group_labels) != nrow(X)) {
    stop("group_labels length must match nrow(X) (samples).")
  }
  groups <- unique(group_labels)
  agg <- t(sapply(groups, function(g) {
    idx <- which(group_labels == g)
    if (length(idx) == 1) X[idx, ] else colMeans(X[idx, , drop = FALSE])
  }))
  rownames(agg) <- groups
  # Center across groups
  agg_c <- scale(agg, center = TRUE, scale = FALSE)
  # PCA via SVD
  sv <- svd(agg_c, nu = 0, nv = min(k, ncol(agg_c)))
  W  <- sv$v
  rownames(W) <- colnames(X)
  colnames(W) <- paste0("PC", seq_len(ncol(W)))
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
r_eff_real    <- apply(alpha_real, 2, effective_rank)   # one per phenotype
conc_real     <- apply(alpha_real, 2, concentration)
rho_real_vec  <- colSums(alpha_real)
names(r_eff_real) <- colnames(phi_real)
names(conc_real)  <- colnames(phi_real)
names(rho_real_vec) <- colnames(phi_real)

cat("\n--- REAL Wm ---\n")
cat("  rho_total      :", rho_real_vec, "\n")
cat("  effective rank :", r_eff_real, " (k =", k_mps, ")\n")
cat("  concentration  :", conc_real, "\n")

# ===============================================================
# LAYER 1 (b): 1000 RANDOM ORTHONORMAL Wm
# ===============================================================
cat("\n=== Layer 1 (b): random orthonormal Wm  (N =", N_RANDOM_ORTHO, ") ===\n")
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
# LAYER 1 (a): GENERATIVE BASELINES (same types as existing script)
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

cat("\n=== Layer 1 (a): generative baselines  (N =", N_BASELINE_REPS, " per type) ===\n")
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
# LAYER 2: PERMUTATION TEST ON CONDITION LABELS
# ===============================================================
# Logic:
#   - Take the ORIGINAL MPS counts and metadata.
#   - Permute the values in the condition column (i.e., randomly
#     re-assign samples to conditions, preserving the number of
#     samples per condition).
#   - Aggregate (mean) samples within each *permuted* condition.
#   - Run PCA on the aggregated, centered matrix.
#   - This preserves: number of samples, number of conditions k,
#     gene-marginal distributions, gene-gene covariance.
#   - This destroys: the relationship between specific perturbations
#     and gene-expression responses.
#
# We use direct PCA (not process_datasets) for transparency.
# ===============================================================
cat("\n=== Layer 2: condition-label permutation test  (N =", N_PERMUTATIONS, ") ===\n")

# Sanity check: real Wm vs direct PCA on real metadata
# (this should produce something equivalent to data_list[[target]]$Wm_group)
real_groups <- as.character(metadata_m[[condition_col]])

# Sample-level centered MPS expression (rows = samples, cols = genes)
# Note: data_center in your pipeline is genes x samples, so we transpose
Xm_samples <- data_list[[target_dataset]]$data_center %>% as.matrix() %>% t()

# Reference Wm computed by the same direct path (used only for sanity logging)
Wm_check <- pca_after_aggregation(Xm_samples, real_groups, k = k_mps)
cos_sim_check <- diag(abs(t(Wm_check[rownames(Wm), 1:k_mps]) %*% Wm[, 1:k_mps]))
cat("Sanity check (direct-PCA reproducibility): mean |cos sim| of leading PC pairs =",
    round(mean(cos_sim_check, na.rm = TRUE), 3), "\n")
# The leading PCs should be very similar; small rotations are normal.

results_perm <- vector("list", N_PERMUTATIONS)
for (i in seq_len(N_PERMUTATIONS)) {
  if (i %% 100 == 0) cat("  permutation", i, "/", N_PERMUTATIONS, "\n")
  tryCatch({
    perm_groups <- sample(real_groups)   # break sample-condition link
    Wm_perm     <- pca_after_aggregation(Xm_samples, perm_groups, k = k_mps)
    
    alpha_p <- compute_alpha(Wm_perm, phi_real)
    r_eff_p <- apply(alpha_p, 2, effective_rank)
    conc_p  <- apply(alpha_p, 2, concentration)
    rho_p   <- colSums(alpha_p)
    
    # Predictive performance using permuted Wm
    Yhat_perm  <- predict(model_real, Xh %*% Wm_perm %*% t(Wm_perm))
    spear_perm <- diag(cor(Yh, Yhat_perm, method = 'spearman'))
    
    # Note: we do NOT run liv2trans_run on the permuted Wm here, to keep
    # the loop fast and well-defined. The relevant null comparison is
    # rho/r_eff/concentration on permuted Wm vs real Wm. If you also
    # want to test extra-LV pathway loadings under the null, see the
    # optional block below.
    
    results_perm[[i]] <- data.frame(
      rep           = i,
      phenotype     = colnames(phi_real),
      r_eff         = r_eff_p,
      concentration = conc_p,
      rho_total     = rho_p,
      spearman      = spear_perm,
      k_random      = ncol(Wm_perm),
      source        = "permuted_conditions",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("  perm", i, "FAILED:", e$message, "\n")
  })
}
df_perm <- bind_rows(results_perm)

# ===============================================================
# REPORT P-VALUES
# ===============================================================
cat("\n========================================\n")
cat("LAYER 1 (b): vs random orthonormal Wm\n")
cat("========================================\n")
cat("\nP(r_eff_real <= r_eff_rand)  -- real should be SMALLER (more concentrated)\n")
print(compute_pvals(df_ortho, r_eff_real, "r_eff", "less"))
cat("\nP(concentration_real >= concentration_rand)\n")
print(compute_pvals(df_ortho, conc_real, "concentration", "greater"))
cat("\nP(rho_total_real >= rho_total_rand)\n")
print(compute_pvals(df_ortho, rho_real_vec, "rho_total", "greater"))

cat("\n========================================\n")
cat("LAYER 2: vs permuted condition labels\n")
cat("========================================\n")
cat("\nP(r_eff_real <= r_eff_perm)\n")
print(compute_pvals(df_perm, r_eff_real, "r_eff", "less"))
cat("\nP(concentration_real >= concentration_perm)\n")
print(compute_pvals(df_perm, conc_real, "concentration", "greater"))
cat("\nP(rho_total_real >= rho_total_perm)\n")
print(compute_pvals(df_perm, rho_real_vec, "rho_total", "greater"))

# ===============================================================
# SAVE
# ===============================================================
save(df_ortho, df_baseline, df_perm,
     alpha_real, r_eff_real, conc_real, rho_real_vec,
     phi_real, Wm, Wm_opt, W_TC,
     k_mps, p_genes, n_TC,
     ref_dataset, target_dataset, dataset_tag,
     condition_col,
     file = file.path(RESULTS_DIR,
                      paste0(dataset_tag,
                             "__layer1_layer2_concentration_results.RData")))

# ===============================================================
# PLOTS
# ===============================================================

# Helper: histogram of null distribution with red dashed line for real value
plot_null_with_real <- function(df_null, real_value, var, xlab, title, fname,
                                facet_var = "phenotype") {
  if (nrow(df_null) == 0) return(invisible(NULL))
  ref_lines <- data.frame(phenotype = names(real_value),
                          value     = unname(real_value))
  p <- ggplot(df_null, aes(x = .data[[var]])) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7,
                   color = "black", linewidth = 0.3) +
    geom_vline(data = ref_lines, aes(xintercept = value),
               color = "red", linetype = "dashed", linewidth = 1) +
    facet_wrap(as.formula(paste0("~", facet_var)), scales = "free") +
    labs(x = xlab, y = "Count", title = title,
         subtitle = paste0("Red dashed = real Wm  (N = ", nrow(df_null) /
                             length(unique(df_null$phenotype)), " replicates)")) +
    theme_bw(base_size = 12)
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__", fname, ".png")),
         plot = p, width = 10, height = 5, dpi = 300)
}

# --- Layer 1 (b): random orthonormal ---
plot_null_with_real(df_ortho, r_eff_real, "r_eff",
                    "Effective rank of alignment",
                    "Layer 1: Effective rank vs random orthonormal Wm",
                    "layer1_reff_random_orthonormal")
plot_null_with_real(df_ortho, conc_real, "concentration",
                    "Concentration  (1 - r_eff/k)",
                    "Layer 1: Concentration vs random orthonormal Wm",
                    "layer1_concentration_random_orthonormal")
plot_null_with_real(df_ortho, rho_real_vec, "rho_total",
                    expression(rho ~ "(fraction of " * phi[PLS] * " captured)"),
                    "Layer 1: rho_total vs random orthonormal Wm",
                    "layer1_rho_random_orthonormal")

# --- Layer 1 (a): generative baselines, faceted by baseline x phenotype ---
if (nrow(df_baseline) > 0) {
  ref_lines <- data.frame(phenotype = names(r_eff_real),
                          value     = unname(r_eff_real))
  p_bl <- ggplot(df_baseline, aes(x = r_eff)) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7,
                   color = "black", linewidth = 0.3) +
    geom_vline(data = ref_lines, aes(xintercept = value),
               color = "red", linetype = "dashed", linewidth = 1) +
    facet_grid(baseline ~ phenotype, scales = "free") +
    labs(x = "Effective rank of alignment", y = "Count",
         title = "Layer 1 (a): Effective rank across generative baselines",
         subtitle = "Red dashed = real Wm") +
    theme_bw(base_size = 9)
  ggsave(file.path(FIGURES_DIR,
                   paste0(dataset_tag, "__layer1_reff_baselines_grid.png")),
         plot = p_bl, width = 10, height = 14, dpi = 300)
}

# --- Layer 2: permutation ---
plot_null_with_real(df_perm, r_eff_real, "r_eff",
                    "Effective rank of alignment",
                    "Layer 2: Effective rank vs permuted condition labels",
                    "layer2_reff_permutation")
plot_null_with_real(df_perm, conc_real, "concentration",
                    "Concentration  (1 - r_eff/k)",
                    "Layer 2: Concentration vs permuted condition labels",
                    "layer2_concentration_permutation")
plot_null_with_real(df_perm, rho_real_vec, "rho_total",
                    expression(rho ~ "(fraction of " * phi[PLS] * " captured)"),
                    "Layer 2: rho_total vs permuted condition labels",
                    "layer2_rho_permutation")

# Also: predictive performance under permutation
if ("spearman" %in% colnames(df_perm)) {
  Yhat_real_full <- predict(model_real, Xh %*% Wm %*% t(Wm))
  spear_real_full <- diag(cor(Yh, Yhat_real_full, method = 'spearman'))
  ref_lines <- data.frame(phenotype = names(spear_real_full),
                          value     = unname(spear_real_full))
  p_sp <- ggplot(df_perm, aes(x = spearman)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7,
                   color = "black", linewidth = 0.3) +
    geom_vline(data = ref_lines, aes(xintercept = value),
               color = "red", linetype = "dashed", linewidth = 1) +
    facet_wrap(~ phenotype, scales = "free") +
    labs(x = "Spearman correlation (true Y vs Yhat through permuted Wm)",
         y = "Count",
         title = "Layer 2: Predictive performance under permuted Wm",
         subtitle = "Red dashed = real Wm performance") +
    theme_bw(base_size = 12)
  ggsave(file.path(FIGURES_DIR,
                   paste0(dataset_tag, "__layer2_spearman_permutation.png")),
         plot = p_sp, width = 10, height = 5, dpi = 300)
}

# --- Combined comparison plot: r_eff across all null models ---
df_combined <- bind_rows(
  df_ortho    %>% transmute(phenotype, r_eff, concentration, rho_total,
                            group = "Random orthonormal"),
  df_baseline %>% transmute(phenotype, r_eff, concentration, rho_total,
                            group = paste0("Baseline: ", baseline)),
  df_perm     %>% transmute(phenotype, r_eff, concentration, rho_total,
                            group = "Permuted conditions")
)
ref_lines <- data.frame(phenotype = names(r_eff_real),
                        value     = unname(r_eff_real))
p_comb <- ggplot(df_combined, aes(x = group, y = r_eff)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  geom_hline(data = ref_lines, aes(yintercept = value),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ phenotype, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Effective rank of alignment",
       title = "Layer 1 + Layer 2: r_eff across all null models",
       subtitle = "Red dashed = real Wm") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGURES_DIR,
                 paste0(dataset_tag, "__combined_reff_comparison.png")),
       plot = p_comb, width = 12, height = 8, dpi = 300)

# Same comparison for concentration
ref_lines2 <- data.frame(phenotype = names(conc_real),
                         value     = unname(conc_real))
p_comb2 <- ggplot(df_combined, aes(x = group, y = concentration)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  geom_hline(data = ref_lines2, aes(yintercept = value),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ phenotype, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Concentration  (1 - r_eff/k)",
       title = "Layer 1 + Layer 2: concentration across all null models",
       subtitle = "Red dashed = real Wm") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGURES_DIR,
                 paste0(dataset_tag, "__combined_concentration_comparison.png")),
       plot = p_comb2, width = 12, height = 8, dpi = 300)

cat("\n============================================\n")
cat("All plots:   ", FIGURES_DIR, "\n")
cat("Results:     ",
    file.path(RESULTS_DIR,
              paste0(dataset_tag, "__layer1_layer2_concentration_results.RData")),
    "\n", sep = "")
cat("============================================\n")

