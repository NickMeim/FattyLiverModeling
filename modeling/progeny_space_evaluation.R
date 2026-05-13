# LIV2TRANS in PROGENy pathway-activity space (high k/p regime)
# This version adds Tanimoto similarity and GSEA-based distance at
# multiple top-N values (2, 3, 5) for both TCs and extra LVs.
# Cosine-similarity figures are unchanged; new figures are saved
# separately.

library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(LIV2Trans)
library(matrixStats)
library(decoupleR)
library(fgsea)

# ---------------------------------------------------------------
# Loader
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
N_RANDOM_ORTHO   <- 1000
N_RANDOM_REPS    <- 50
TOP_N_CONDITIONS <- 6
TOP_N_VALUES     <- c(2, 3, 5)        # for PROGENy (p_path ~= 14)
PROGENY_TOP      <- 500
FGSEA_NPERM      <- 100

FIGURES_DIR <- "../figures/"
RESULTS_DIR <- "../results/"
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

dataset_names  <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset    <- "Govaere"
target_dataset <- "Kostrzewski"

clean_filename <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
dataset_tag <- paste0("ref-",     clean_filename(ref_dataset),
                      "__target-", clean_filename(target_dataset),
                      "__PROGENYspace")

# ===============================================================
# PART 1: LOAD DATA
# ===============================================================
data_list <- load_datasets(dataset_names, dir_data = '../data/')
tmp       <- process_datasets(data_list, filter_variance = FALSE)
data_list <- tmp$data_list

Yh <- as.matrix(data_list[[ref_dataset]]$metadata %>%
                  dplyr::select(nas_score, Fibrosis_stage))
colnames(Yh) <- c('NAS', 'fibrosis')

Xh_gene <- data_list[[ref_dataset]]$data_center %>% as.matrix()
Xm_gene <- data_list[[target_dataset]]$data_center %>% as.matrix()

metadata_m <- data_list[[target_dataset]]$metadata

condition_col <- 'treatment'
cat("Using condition column:", condition_col, "\n")
real_groups <- as.character(metadata_m[[condition_col]])
cat("Number of unique MPS conditions:", length(unique(real_groups)), "\n")

# ===============================================================
# PART 2: PATHWAY ACTIVITY INFERENCE
# ===============================================================
cat("\n=== Pathway activity inference (PROGENy + VIPER) ===\n")
network <- tryCatch(
  decoupleR::get_progeny(organism = "human", top = PROGENY_TOP),
  error = function(e) {
    cat("  decoupleR::get_progeny failed; trying OmnipathR fallback\n")
    OmnipathR::progeny(organism = "human", top = PROGENY_TOP)
  }
)
cat("  PROGENy regulon loaded:", nrow(network), "interactions,",
    length(unique(network$source)), "pathways\n")

run_pathway_inference <- function(expr_mat) {
  res <- decoupleR::run_viper(
    mat       = expr_mat,
    network   = network,
    .source   = "source",
    .target   = "target",
    .mor      = "weight",
    minsize   = 5
  )
  pmat <- res %>%
    dplyr::filter(statistic == "viper") %>%
    tidyr::pivot_wider(id_cols    = condition,
                       names_from = source,
                       values_from = score) %>%
    tibble::column_to_rownames("condition") %>%
    as.matrix()
  pmat
}

P_h_raw <- run_pathway_inference(Xh_gene)
P_m_raw <- run_pathway_inference(Xm_gene)

common_paths <- intersect(colnames(P_h_raw), colnames(P_m_raw))
P_h <- P_h_raw[, common_paths, drop = FALSE]
P_m <- P_m_raw[, common_paths, drop = FALSE]
P_m <- P_m[rownames(metadata_m), , drop = FALSE]
ref_meta <- data_list[[ref_dataset]]$metadata
# P_h      <- P_h[rownames(ref_meta), , drop = FALSE]

p_path <- ncol(P_h)
cat("  Pathway activities computed:", p_path, "pathways;",
    nrow(P_h), "human samples;", nrow(P_m), "MPS samples\n")

# Ensure top-N values are sensible relative to p_path
TOP_N_VALUES <- TOP_N_VALUES[TOP_N_VALUES < p_path]
if (length(TOP_N_VALUES) == 0)
  stop("All TOP_N_VALUES >= p_path; nothing to compute.")

Xh_path <- scale(P_h, center = TRUE, scale = FALSE)
Xm_path <- scale(P_m, center = TRUE, scale = FALSE)

# ===============================================================
# PART 3: HELPERS
# ===============================================================

compute_Wm_from_matrix <- function(M_samples_by_features,
                                   condition_labels,
                                   k = NULL) {
  if (length(condition_labels) != nrow(M_samples_by_features)) {
    stop("condition_labels length must match nrow(M).")
  }
  groups <- unique(condition_labels)
  agg <- t(sapply(groups, function(g) {
    idx <- which(condition_labels == g)
    if (length(idx) == 1) M_samples_by_features[idx, ]
    else colMeans(M_samples_by_features[idx, , drop = FALSE])
  }))
  rownames(agg) <- groups
  agg_c <- scale(agg, center = TRUE, scale = FALSE)
  k_max <- min(nrow(agg_c) - 1, ncol(agg_c))
  if (is.null(k)) k <- k_max else k <- min(k, k_max)
  sv <- svd(agg_c, nu = 0, nv = k)
  W <- sv$v[, seq_len(k), drop = FALSE]
  rownames(W) <- colnames(M_samples_by_features)
  colnames(W) <- paste0("PC", seq_len(ncol(W)))
  W
}

compute_alpha <- function(W, phi) {
  feats <- intersect(rownames(W), rownames(phi))
  if (length(feats) == 0) stop("No common features between W and phi.")
  W_c   <- W[feats, , drop = FALSE]
  phi_c <- phi[feats, , drop = FALSE]
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

effective_rank <- function(a) {
  s1 <- sum(a); s2 <- sum(a^2)
  if (s2 < .Machine$double.eps) return(NA_real_)
  (s1^2) / s2
}
concentration <- function(a) {
  r <- effective_rank(a); k <- length(a)
  if (is.na(r) || k <= 1) return(NA_real_)
  1 - (r / k)
}

random_orthonormal <- function(p, k, feature_names = NULL) {
  G <- matrix(rnorm(p * k), nrow = p, ncol = k)
  W <- qr.Q(qr(G))
  if (!is.null(feature_names)) rownames(W) <- feature_names
  W
}

matched_cosine <- function(A, B) {
  if (ncol(A) == 0 || ncol(B) == 0) return(data.frame())
  feats <- intersect(rownames(A), rownames(B))
  if (length(feats) == 0) return(data.frame())
  A_c <- A[feats, , drop = FALSE]
  B_c <- B[feats, , drop = FALSE]
  norm_cols <- function(M) {
    nrm <- sqrt(colSums(M^2)); nrm[nrm == 0] <- 1
    sweep(M, 2, nrm, "/")
  }
  A_n <- norm_cols(A_c); B_n <- norm_cols(B_c)
  S <- t(A_n) %*% B_n
  out <- data.frame()
  used_B <- integer(0)
  for (i in seq_len(min(ncol(A_n), ncol(B_n)))) {
    cand <- abs(S[i, , drop = TRUE])
    cand[used_B] <- -Inf
    j <- which.max(cand)
    out <- rbind(out, data.frame(
      A_col = colnames(A_n)[i], B_col = colnames(B_n)[j],
      A_idx = i, B_idx = j,
      cosine = S[i, j], abs_cosine = abs(S[i, j])
    ))
    used_B <- c(used_B, j)
  }
  out
}

# Tanimoto/Jaccard on top-N feature sets (by |loading|)
tanimoto_topN <- function(vA, vB, top_n) {
  if (length(vA) < top_n || length(vB) < top_n) return(NA_real_)
  top_A <- names(sort(abs(vA), decreasing = TRUE))[seq_len(top_n)]
  top_B <- names(sort(abs(vB), decreasing = TRUE))[seq_len(top_n)]
  length(intersect(top_A, top_B)) / length(union(top_A, top_B))
}

# Symmetric GSEA: mean(|ES_AB|, |ES_BA|); distance = 1 - mean
gsea_topN <- function(vA, vB, top_n) {
  if (length(vA) < top_n || length(vB) < top_n) {
    return(list(es_mean = NA_real_, distance = NA_real_))
  }
  feats <- intersect(names(vA), names(vB))
  if (length(feats) < top_n + 1) {
    return(list(es_mean = NA_real_, distance = NA_real_))
  }
  vA <- vA[feats]; vB <- vB[feats]
  top_A <- names(sort(abs(vA), decreasing = TRUE))[seq_len(top_n)]
  top_B <- names(sort(abs(vB), decreasing = TRUE))[seq_len(top_n)]
  abs_vA <- abs(vA); abs_vB <- abs(vB)
  
  run_fgsea <- function(set_top, stats_vec) {
    s <- stats_vec + rnorm(length(stats_vec), sd = .Machine$double.eps * 10)
    fg <- tryCatch(
      suppressWarnings(fgsea::fgsea(
        pathways    = list(set = set_top),
        stats       = s,
        minSize     = 1,
        maxSize     = length(set_top) + 1,
        scoreType   = "pos",
        nPermSimple = FGSEA_NPERM
      )),
      error = function(e) NULL
    )
    if (is.null(fg) || nrow(fg) == 0) return(NA_real_)
    fg$ES[1]
  }
  es_AB <- run_fgsea(top_A, abs_vB)
  es_BA <- run_fgsea(top_B, abs_vA)
  es_mean <- mean(c(abs(es_AB), abs(es_BA)), na.rm = TRUE)
  list(es_mean = es_mean, distance = 1 - es_mean)
}

matched_topN_metrics <- function(A, B, top_ns) {
  cs <- matched_cosine(A, B)
  if (nrow(cs) == 0) return(data.frame())
  feats <- intersect(rownames(A), rownames(B))
  if (length(feats) == 0) return(data.frame())
  A_f <- A[feats, , drop = FALSE]
  B_f <- B[feats, , drop = FALSE]
  out <- list()
  for (i in seq_len(nrow(cs))) {
    a_name <- cs$A_col[i]; b_name <- cs$B_col[i]
    vA <- A_f[, a_name]; names(vA) <- feats
    vB <- B_f[, b_name]; names(vB) <- feats
    for (n in top_ns) {
      tan <- tanimoto_topN(vA, vB, n)
      gs  <- gsea_topN(vA, vB, n)
      out[[length(out) + 1]] <- data.frame(
        A_col          = a_name,
        B_col          = b_name,
        top_n          = n,
        tanimoto       = tan,
        gsea_es        = gs$es_mean,
        gsea_distance  = gs$distance,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

# ===============================================================
# PART 4: REAL LIV2TRANS RUN IN PATHWAY SPACE
# ===============================================================
cat("\n=== Real LIV2TRANS in PROGENy space ===\n")

Wm_path <- compute_Wm_from_matrix(P_m, real_groups)
k_path  <- ncol(Wm_path)
cat("  k_path =", k_path, "  p_path =", p_path,
    "  k/p =", round(k_path/p_path, 4), "\n")

liv2trans_real <- liv2trans_run(Xh_path, Yh, Xm_path, Wm_path, LV_number = 2)

Wm_opt_real <- liv2trans_real$W_opt
rownames(Wm_opt_real) <- rownames(Wm_path)
colnames(Wm_opt_real) <- paste0("LV_extra", seq_len(ncol(Wm_opt_real)))

W_TC_real <- as.matrix(liv2trans_real$W_translatable)
rownames(W_TC_real) <- rownames(Wm_path)
colnames(W_TC_real) <- paste0("TC", seq_len(ncol(W_TC_real)))
n_TC_real <- ncol(W_TC_real)
cat("  number of real TCs:", n_TC_real, "\n")

model_real <- liv2trans_real$model
Wh_real    <- liv2trans_real$Wh
Bh_real    <- t(model_real@weightMN) %*% model_real@coefficientMN
phi_real   <- Wh_real %*% Bh_real
colnames(phi_real) <- colnames(Yh)
rownames(phi_real) <- rownames(Wm_path)

alpha_real <- compute_alpha(Wm_path, phi_real)
r_eff_real <- apply(alpha_real, 2, effective_rank)
conc_real  <- apply(alpha_real, 2, concentration)
rho_real   <- colSums(alpha_real)
names(r_eff_real) <- colnames(phi_real)
names(conc_real)  <- colnames(phi_real)
names(rho_real)   <- colnames(phi_real)

cat("\n  Real Wm_path metrics:\n")
cat("    rho_total      :", rho_real, "\n")
cat("    effective rank :", r_eff_real, " (k =", k_path, ")\n")
cat("    concentration  :", conc_real, "\n")

# ===============================================================
# PART 5: NULL DISTRIBUTIONS
# ===============================================================
results_align   <- list()
results_TC      <- list()
results_TC_topN <- list()

cat("\n=== Null (a): random orthonormal Wm  (N =", N_RANDOM_ORTHO, ") ===\n")
for (i in seq_len(N_RANDOM_ORTHO)) {
  if (i %% 200 == 0) cat("  ", i, "/", N_RANDOM_ORTHO, "\n")
  W_rand <- random_orthonormal(p_path, k_path, rownames(Wm_path))
  a_r    <- compute_alpha(W_rand, phi_real)
  results_align[[length(results_align) + 1]] <- data.frame(
    source = "Random orthonormal", rep = i,
    phenotype     = colnames(phi_real),
    r_eff         = apply(a_r, 2, effective_rank),
    concentration = apply(a_r, 2, concentration),
    rho_total     = colSums(a_r),
    stringsAsFactors = FALSE
  )
}

cat("\n=== Null (b): random Gaussian pathway matrix  (N =", N_RANDOM_REPS, ") ===\n")
mu_path  <- mean(P_m); sig_path <- sd(P_m)
random_path_runs <- list()
for (i in seq_len(N_RANDOM_REPS)) {
  if (i %% 10 == 0) cat("  ", i, "/", N_RANDOM_REPS, "\n")
  tryCatch({
    P_m_rand <- matrix(rnorm(nrow(P_m) * p_path, mean = mu_path, sd = sig_path),
                       nrow = nrow(P_m), ncol = p_path,
                       dimnames = dimnames(P_m))
    Xm_path_rand <- scale(P_m_rand, center = TRUE, scale = FALSE)
    Wm_rand      <- compute_Wm_from_matrix(P_m_rand, real_groups, k = k_path)
    
    a_r <- compute_alpha(Wm_rand, phi_real)
    results_align[[length(results_align) + 1]] <- data.frame(
      source = "Random Gaussian", rep = i,
      phenotype     = colnames(phi_real),
      r_eff         = apply(a_r, 2, effective_rank),
      concentration = apply(a_r, 2, concentration),
      rho_total     = colSums(a_r),
      stringsAsFactors = FALSE
    )
    
    res_rand <- liv2trans_run(Xh_path, Yh, Xm_path_rand, Wm_rand, LV_number = 2)
    W_TC_rand <- as.matrix(res_rand$W_translatable)
    rownames(W_TC_rand) <- rownames(Wm_rand)
    colnames(W_TC_rand) <- paste0("TC", seq_len(ncol(W_TC_rand)))
    Wm_opt_rand <- res_rand$W_opt
    rownames(Wm_opt_rand) <- rownames(Wm_rand)
    colnames(Wm_opt_rand) <- paste0("LV_extra", seq_len(ncol(Wm_opt_rand)))
    
    if (ncol(W_TC_rand) > 0 && ncol(W_TC_real) > 0) {
      cs <- matched_cosine(W_TC_real, W_TC_rand)
      cs$source <- "Random Gaussian"; cs$rep <- i; cs$kind <- "TC"
      results_TC[[length(results_TC) + 1]] <- cs
      
      tn <- matched_topN_metrics(W_TC_real, W_TC_rand, TOP_N_VALUES)
      if (nrow(tn) > 0) {
        tn$source <- "Random Gaussian"; tn$rep <- i; tn$kind <- "TC"
        results_TC_topN[[length(results_TC_topN) + 1]] <- tn
      }
    }
    
    random_path_runs[[i]] <- list(
      P_m = P_m_rand, Xm = Xm_path_rand,
      Wm = Wm_rand, W_TC = W_TC_rand, Wm_opt = Wm_opt_rand,
      source = "Random Gaussian"
    )
  }, error = function(e) cat("    rep", i, "FAILED:", e$message, "\n"))
}

cat("\n=== Null (c): shuffled pathway table  (N =", N_RANDOM_REPS, ") ===\n")
shuffled_path_runs <- list()
for (i in seq_len(N_RANDOM_REPS)) {
  if (i %% 10 == 0) cat("  ", i, "/", N_RANDOM_REPS, "\n")
  tryCatch({
    P_m_sh <- matrix(sample(as.vector(P_m)),
                     nrow = nrow(P_m), ncol = p_path,
                     dimnames = dimnames(P_m))
    Xm_path_sh <- scale(P_m_sh, center = TRUE, scale = FALSE)
    Wm_sh      <- compute_Wm_from_matrix(P_m_sh, real_groups, k = k_path)
    
    a_r <- compute_alpha(Wm_sh, phi_real)
    results_align[[length(results_align) + 1]] <- data.frame(
      source = "Shuffled pathway table", rep = i,
      phenotype     = colnames(phi_real),
      r_eff         = apply(a_r, 2, effective_rank),
      concentration = apply(a_r, 2, concentration),
      rho_total     = colSums(a_r),
      stringsAsFactors = FALSE
    )
    
    res_sh <- liv2trans_run(Xh_path, Yh, Xm_path_sh, Wm_sh, LV_number = 2)
    W_TC_sh <- as.matrix(res_sh$W_translatable)
    rownames(W_TC_sh) <- rownames(Wm_sh)
    colnames(W_TC_sh) <- paste0("TC", seq_len(ncol(W_TC_sh)))
    Wm_opt_sh <- res_sh$W_opt
    rownames(Wm_opt_sh) <- rownames(Wm_sh)
    colnames(Wm_opt_sh) <- paste0("LV_extra", seq_len(ncol(Wm_opt_sh)))
    
    if (ncol(W_TC_sh) > 0 && ncol(W_TC_real) > 0) {
      cs <- matched_cosine(W_TC_real, W_TC_sh)
      cs$source <- "Shuffled pathway table"; cs$rep <- i; cs$kind <- "TC"
      results_TC[[length(results_TC) + 1]] <- cs
      
      tn <- matched_topN_metrics(W_TC_real, W_TC_sh, TOP_N_VALUES)
      if (nrow(tn) > 0) {
        tn$source <- "Shuffled pathway table"; tn$rep <- i; tn$kind <- "TC"
        results_TC_topN[[length(results_TC_topN) + 1]] <- tn
      }
    }
    
    shuffled_path_runs[[i]] <- list(
      P_m = P_m_sh, Xm = Xm_path_sh,
      Wm = Wm_sh, W_TC = W_TC_sh, Wm_opt = Wm_opt_sh,
      source = "Shuffled pathway table"
    )
  }, error = function(e) cat("    rep", i, "FAILED:", e$message, "\n"))
}

df_align    <- bind_rows(results_align)
df_TC       <- if (length(results_TC)      > 0) bind_rows(results_TC)      else data.frame()
df_TC_topN  <- if (length(results_TC_topN) > 0) bind_rows(results_TC_topN) else data.frame()

# ===============================================================
# PART 6: CONDITION REDUCTION + RE-RUN LIV2TRANS
# ===============================================================
cat("\n=== Condition reduction + re-run LIV2TRANS ===\n")

select_driving_conditions <- function(Xm_centered, condition_labels,
                                      W_TC, top_n = TOP_N_CONDITIONS) {
  if (ncol(W_TC) < 2) {
    tc_cols <- seq_len(ncol(W_TC))
  } else {
    tc_cols <- 1:2
  }
  feats <- intersect(rownames(W_TC), colnames(Xm_centered))
  proj  <- Xm_centered[, feats] %*% W_TC[feats, tc_cols, drop = FALSE]
  cond_proj <- t(sapply(unique(condition_labels), function(g) {
    idx <- which(condition_labels == g)
    if (length(idx) == 1) proj[idx, ] else colMeans(proj[idx, , drop = FALSE])
  }))
  rownames(cond_proj) <- unique(condition_labels)
  mag <- sqrt(rowSums(cond_proj^2))
  keep_conditions <- names(sort(mag, decreasing = TRUE))[seq_len(min(top_n, length(mag)))]
  list(keep_conditions = keep_conditions,
       keep_samples    = condition_labels %in% keep_conditions,
       cond_proj       = cond_proj,
       mag             = mag)
}

reduce_and_rerun <- function(P_mat, Xm_mat, W_TC, label) {
  sel <- select_driving_conditions(Xm_mat, real_groups, W_TC,
                                   top_n = TOP_N_CONDITIONS)
  kept_samples <- which(sel$keep_samples)
  if (length(unique(real_groups[kept_samples])) < 2) {
    cat("  [", label, "] only", length(unique(real_groups[kept_samples])),
        "condition(s) selected; skipping rerun\n")
    return(NULL)
  }
  Xm_sub  <- Xm_mat[kept_samples, , drop = FALSE]
  groups_sub <- real_groups[kept_samples]
  Wm_sub  <- compute_Wm_from_matrix(P_mat[kept_samples, , drop = FALSE],
                                    groups_sub)
  res_sub <- liv2trans_run(Xh_path, Yh, Xm_sub, Wm_sub, LV_number = 2)
  Wm_opt_sub <- res_sub$W_opt
  rownames(Wm_opt_sub) <- rownames(Wm_sub)
  colnames(Wm_opt_sub) <- paste0("LV_extra", seq_len(ncol(Wm_opt_sub)))
  list(
    Wm        = Wm_sub,
    Wm_opt    = Wm_opt_sub,
    W_TC      = as.matrix(res_sub$W_translatable),
    kept      = sel$keep_conditions,
    n_kept    = length(sel$keep_conditions),
    k_sub     = ncol(Wm_sub),
    p_sub     = nrow(Wm_sub)
  )
}

cat("  Reducing real MPS dataset...\n")
real_reduced <- reduce_and_rerun(P_m, Xm_path, W_TC_real, "real")
cat("    kept", real_reduced$n_kept, "conditions:",
    paste(real_reduced$kept, collapse = ", "), "\n")
cat("    new k =", real_reduced$k_sub, "  p =", real_reduced$p_sub,
    "  (k/p =", round(real_reduced$k_sub/real_reduced$p_sub, 3), ")\n")

results_LVopt      <- list()
results_LVopt_topN <- list()

cat("  Reducing random Gaussian runs...\n")
for (i in seq_along(random_path_runs)) {
  r <- random_path_runs[[i]]
  if (is.null(r) || is.null(r$W_TC) || ncol(r$W_TC) < 1) next
  red <- reduce_and_rerun(r$P_m, r$Xm, r$W_TC, paste0("rand_", i))
  if (is.null(red)) next
  cs <- matched_cosine(real_reduced$Wm_opt, red$Wm_opt)
  if (nrow(cs) > 0) {
    cs$source <- "Random Gaussian"; cs$rep <- i; cs$kind <- "LV_extra"
    results_LVopt[[length(results_LVopt) + 1]] <- cs
    
    tn <- matched_topN_metrics(real_reduced$Wm_opt, red$Wm_opt, TOP_N_VALUES)
    if (nrow(tn) > 0) {
      tn$source <- "Random Gaussian"; tn$rep <- i; tn$kind <- "LV_extra"
      results_LVopt_topN[[length(results_LVopt_topN) + 1]] <- tn
    }
  }
}

cat("  Reducing shuffled pathway runs...\n")
for (i in seq_along(shuffled_path_runs)) {
  r <- shuffled_path_runs[[i]]
  if (is.null(r) || is.null(r$W_TC) || ncol(r$W_TC) < 1) next
  red <- reduce_and_rerun(r$P_m, r$Xm, r$W_TC, paste0("shuf_", i))
  if (is.null(red)) next
  cs <- matched_cosine(real_reduced$Wm_opt, red$Wm_opt)
  if (nrow(cs) > 0) {
    cs$source <- "Shuffled pathway table"; cs$rep <- i; cs$kind <- "LV_extra"
    results_LVopt[[length(results_LVopt) + 1]] <- cs
    
    tn <- matched_topN_metrics(real_reduced$Wm_opt, red$Wm_opt, TOP_N_VALUES)
    if (nrow(tn) > 0) {
      tn$source <- "Shuffled pathway table"; tn$rep <- i; tn$kind <- "LV_extra"
      results_LVopt_topN[[length(results_LVopt_topN) + 1]] <- tn
    }
  }
}

df_LVopt      <- if (length(results_LVopt)      > 0) bind_rows(results_LVopt)      else data.frame()
df_LVopt_topN <- if (length(results_LVopt_topN) > 0) bind_rows(results_LVopt_topN) else data.frame()

self_cs <- matched_cosine(real_reduced$Wm_opt, real_reduced$Wm_opt)
cat("\n  Sanity: real-vs-real extra-LV cosine =",
    round(self_cs$cosine, 3), "\n")

# ===============================================================
# PART 7: PLOTS
# ===============================================================

# 7a: alignment metrics
ref_lines_reff <- data.frame(phenotype = names(r_eff_real),
                             value     = unname(r_eff_real))
p_reff <- ggplot(df_align, aes(x = source, y = r_eff)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  geom_hline(data = ref_lines_reff, aes(yintercept = value),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ phenotype, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Effective rank of alignment",
       subtitle = "Red dashed = real Wm (PROGENy space)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__reff_comparison.png")),
       plot = p_reff, width = 11, height = 6, dpi = 300)

ref_lines_conc <- data.frame(phenotype = names(conc_real),
                             value     = unname(conc_real))
p_conc <- ggplot(df_align, aes(x = source, y = concentration)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  geom_hline(data = ref_lines_conc, aes(yintercept = value),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ phenotype, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Concentration  (1 - r_eff/k)",
       subtitle = "Red dashed = real Wm (PROGENy space)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__concentration_comparison.png")),
       plot = p_conc, width = 11, height = 6, dpi = 300)

# 7b: cosine
if (nrow(df_TC) > 0) {
  p_TCcos <- ggplot(df_TC, aes(x = abs_cosine, fill = source)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "black", linewidth = 0.3,
                   position = "identity") +
    facet_grid(source ~ A_col, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "|cosine similarity| to real TC", y = "Count",
         subtitle = "Best-matched pairing per replicate") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__TC_cosine_similarity.png")),
         plot = p_TCcos, width = 11, height = 6, dpi = 300)
}

if (nrow(df_LVopt) > 0) {
  p_LVcos <- ggplot(df_LVopt, aes(x = abs_cosine, fill = source)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "black", linewidth = 0.3,
                   position = "identity") +
    facet_grid(source ~ A_col, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "|cosine similarity| to real extra LV (post condition-reduction)",
         y = "Count",
         subtitle = paste0("Top ", TOP_N_CONDITIONS,
                           " conditions kept per scenario; LIV2TRANS rerun")) +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__extraLV_cosine_similarity.png")),
         plot = p_LVcos, width = 11, height = 6, dpi = 300)
}

# 7c: NEW Tanimoto / GSEA for TC
if (nrow(df_TC_topN) > 0) {
  df_TC_topN$top_n_f <- factor(df_TC_topN$top_n, levels = TOP_N_VALUES)
  p_TCtan <- ggplot(df_TC_topN, aes(x = top_n_f, y = tanimoto, fill = source)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
    facet_wrap(~ A_col, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Top-N pathways", y = "Tanimoto similarity",
         fill = "",
         subtitle = "Top-N selected by |loading|; real vs null per replicate") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__TC_tanimoto_similarity.png")),
         plot = p_TCtan, width = 11, height = 6, dpi = 300)
  
  p_TCgsea <- ggplot(df_TC_topN, aes(x = top_n_f, y = gsea_distance, fill = source)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
    facet_wrap(~ A_col, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Top-N pathways", y = "GSEA-based distance  (1 - mean |ES|)",
         fill = "",
         subtitle = "Symmetric: top-N of one as set, |loading| of other as rank") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__TC_gsea_distance.png")),
         plot = p_TCgsea, width = 11, height = 6, dpi = 300)
}

# 7d: NEW Tanimoto / GSEA for extra LV
if (nrow(df_LVopt_topN) > 0) {
  df_LVopt_topN$top_n_f <- factor(df_LVopt_topN$top_n, levels = TOP_N_VALUES)
  p_LVtan <- ggplot(df_LVopt_topN, aes(x = top_n_f, y = tanimoto, fill = source)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
    facet_wrap(~ A_col, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Top-N pathways", y = "Tanimoto similarity",
         fill = "",
         subtitle = paste0("Post condition-reduction; top ",
                           TOP_N_CONDITIONS, " conditions kept")) +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__extraLV_tanimoto_similarity.png")),
         plot = p_LVtan, width = 11, height = 6, dpi = 300)
  
  p_LVgsea <- ggplot(df_LVopt_topN, aes(x = top_n_f, y = gsea_distance, fill = source)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
    facet_wrap(~ A_col, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Top-N pathways", y = "GSEA-based distance  (1 - mean |ES|)",
         fill = "",
         subtitle = paste0("Post condition-reduction; top ",
                           TOP_N_CONDITIONS, " conditions kept")) +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGURES_DIR, paste0(dataset_tag, "__extraLV_gsea_distance.png")),
         plot = p_LVgsea, width = 11, height = 6, dpi = 300)
}

# ===============================================================
# PART 8: SAVE
# ===============================================================
save(df_align, df_TC, df_TC_topN, df_LVopt, df_LVopt_topN,
     alpha_real, r_eff_real, conc_real, rho_real,
     Wm_path, Wm_opt_real, W_TC_real, phi_real,
     real_reduced,
     P_h, P_m, real_groups, condition_col,
     k_path, p_path, TOP_N_CONDITIONS, TOP_N_VALUES,
     ref_dataset, target_dataset, dataset_tag,
     file = file.path(RESULTS_DIR,
                      paste0(dataset_tag,
                             "__progeny_space_results_with_topN.RData")))

cat("\n============================================\n")
cat("Plots:   ", FIGURES_DIR, "\n")
cat("Results: ",
    file.path(RESULTS_DIR,
              paste0(dataset_tag, "__progeny_space_results_with_topN.RData")),
    "\n", sep = "")
cat("============================================\n")