library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(LIV2Trans)
# library(AnnotationDbi)
# library(org.Hs.eg.db)
library(matrixStats)
# Function for loading datasets and appending them to a list
load_datasets <- function(dataset_names, dir_data){
  # Look at files in data directory
  files <- dir(dir_data)
  # Trim to datasets
  files <- files[grepl("dataset.RData", files)]
  # Loop for dataset names given and load one by one, appending to a list.
  # dataset names don't have to match perfectly
  data_list <- list(NULL)
  if (length(files) > 0){ # Loop if any match was found
    for (name in dataset_names){
      if (any(grepl(name, files, ignore.case = TRUE))){ # Avoid loading bad matches
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
    # Remove redundant first element
    data_list[[1]] <- NULL
  }
  
  return(data_list)
  
} # End function

## Load my actual data-------------
# Load data and compute extra basis as in main_example_run.R
# Xh: human gene-expression matrix (samples × genes)
# Yh: vector of MASLD or fibrosis scores
# Wm: PC rotation matrix of MPS (genes × nPC)
# Wh: human PLSR loading matrix (genes × nComp)
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = FALSE)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% dplyr::select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
if (ref_dataset=='Hoang'){
  sex_inferred <- data_list[[ref_dataset]]$metadata$sex
}else if (ref_dataset=='Pantano') {
  sex_inferred <- data_list[[ref_dataset]]$metadata$Sex
}else{
  sex_inferred <- apply(as.matrix(Xh[,c('RPS4Y1')]),2,sign)
  sex_inferred <- 1*(sex_inferred>0)
  sex_inferred <- ifelse(sex_inferred==1,'male','female')
}
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# geneVars <- rowVars(t(Xm))
# ind <- order(-geneVars)[1:100]
# Compute extra basis (two extra latent variables)
liv2trans_results <- liv2trans_run(Xh,Yh,Xm,Wm)
Wm_opt <- liv2trans_results$W_opt
rownames(Wm_opt) <- rownames(Wm)
colnames(Wm_opt) <- c("LV_extra1", "LV_extra2")
## also keep counts of in vitro subsetted
counts_m <- data_list[[target_dataset]]$counts
gm0 <- rownames(counts_m)
ngm0 <- length(gm0)
sm0 <- colnames(counts_m)
nsm0 <- length(sm0)

### Generation of random data-------------------
generate_all_baselines <- function(n, m, real_counts = NULL) {
  
  baselines <- list()
  
  # 1) Random synthetic bulk RNA-seq counts
  baselines$rnaseq_synthetic <- matrix(
    rnbinom(n * m, mu = sample(5:500, m, replace = TRUE),
            size = sample(seq(0.1, 2, by = 0.1), m, replace = TRUE)),
    nrow = n, ncol = m
  )
  
  # 2) Negative binomial with different parameterizations
  baselines$nb_low_disp  <- matrix(rnbinom(n * m, mu = 100,  size = 10),  nrow = n, ncol = m)
  baselines$nb_high_disp <- matrix(rnbinom(n * m, mu = 100,  size = 0.5), nrow = n, ncol = m)
  baselines$nb_low_mean  <- matrix(rnbinom(n * m, mu = 10,   size = 1),   nrow = n, ncol = m)
  baselines$nb_high_mean <- matrix(rnbinom(n * m, mu = 1000, size = 2),   nrow = n, ncol = m)
  
  # 3) Other integer count distributions
  baselines$poisson <- matrix(rpois(n * m, lambda = sample(5:500, m, replace = TRUE)),
                              nrow = n, ncol = m)
  baselines$zip <- matrix(
    rbinom(n * m, 1, prob = 0.7) * rpois(n * m, lambda = 50),
    nrow = n, ncol = m
  )
  baselines$zinb <- matrix(
    rbinom(n * m, 1, prob = 0.6) * rnbinom(n * m, mu = 100, size = 1),
    nrow = n, ncol = m
  )
  baselines$geometric <- matrix(rgeom(n * m, prob = sample(seq(0.01, 0.2, by = 0.01), m, replace = TRUE)),
                                nrow = n, ncol = m)
  baselines$uniform_int <- matrix(sample(0:1000, n * m, replace = TRUE), nrow = n, ncol = m)
  
  # 4) NB from real data params (method of moments, handles non-integer data)
  if (!is.null(real_counts)) {
    mu_est   <- pmax(colMeans(real_counts, na.rm = TRUE), 0.01)
    gene_vars <- apply(real_counts, 2, var, na.rm = TRUE)
    # NB size = mu^2 / (var - mu); clamp when var <= mu (Poisson-like)
    size_est <- ifelse(gene_vars > mu_est,
                       mu_est^2 / (gene_vars - mu_est),
                       1e6) # large size ≈ Poisson limit
    size_est <- pmax(size_est, 0.01)
    # Replace any remaining NA/NaN/Inf
    mu_est[!is.finite(mu_est)]     <- 0.01
    size_est[!is.finite(size_est)] <- 1
    
    baselines$from_real_params <- sapply(1:m, function(j) {
      rnbinom(n, mu = mu_est[j], size = size_est[j])
    })
    
    # 5) Full matrix permutation (works with non-integer averaged data too)
    baselines$permuted <- matrix(sample(as.vector(real_counts)), nrow = n, ncol = m)
  }
  
  return(baselines)
}

# --- Usage
set.seed(42)
baselines <- generate_all_baselines(n = nsm0, m = ngm0, real_counts = counts_m)

## Run LIV2TRANS with random and save loss of backprojection and cosine similarity of extra LVs---------
data_list_new <- load_datasets(dataset_names, dir_data = '../data/')
new_counts_m <- t(baselines$uniform_int)
colnames(new_counts_m) <- colnames(counts_m)
rownames(new_counts_m) <- rownames(counts_m)
data_list_new[[target_dataset]]$counts <- new_counts_m
tmp <- process_datasets(data_list_new, filter_variance = F)
data_list_new <- tmp$data_list
plt_list_new <- tmp$plt_list
# Xh_new <- data_list_new[[ref_dataset]]$data_center %>% t()
Xm_new <- data_list_new[[target_dataset]]$data_center %>% t()
Xm_grouped_new <- data_list_new$Kostrzewski$Xm_grouped
Wm_new <- data_list_new[[target_dataset]]$Wm_group %>% as.matrix()


# Compute extra basis (two extra latent variables)
# W <- matrix(rnorm(26001*47),26001,47)
# W <- matrix(sample(as.vector(Wm)), nrow = 26001, ncol = 47)
W <- Wm_new
liv2trans_results_new <- liv2trans_run(Xh,Yh,Xm_new,W)
Wm_opt_new <- liv2trans_results_new$W_opt
rownames(Wm_opt_new) <- rownames(Wm_new)
colnames(Wm_opt_new) <- c("new_LV_extra1", "new_LV_extra2")
cor(Wm_opt_new,Wm_opt)

### Check what happens when looking at top weights
top_k <- c(5,10,30,50,100,200,300,400,500,600,700,800,900,1000)
## 1) Tanimoto similarity
tanimoto <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  
  union_size <- length(union(a, b))
  if (union_size == 0) return(NA_real_)
  
  length(intersect(a, b)) / union_size
}

get_top_features <- function(x, k, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  
  x <- sort(x, decreasing = direction == "positive")
  
  if (direction == "positive") {
    x <- x[x > 0]
  } else {
    x <- x[x < 0]
  }
  
  if (length(x) < k) return(character(0))
  
  names(x)[seq_len(k)]
}

calc_tanimoto_for_column <- function(col_index,W_opt,W_opt_new) {
  old_scores <- W_opt[, col_index]
  new_scores <- W_opt_new[, col_index]
  
  names(old_scores) <- rownames(W_opt)
  names(new_scores) <- rownames(W_opt_new)
  
  bind_rows(
    lapply(top_k, function(k) {
      old_pos <- get_top_features(old_scores, k, "positive")
      new_pos <- get_top_features(new_scores, k, "positive")
      
      data.frame(
        column = paste0("Column ", col_index),
        k = k,
        direction = "Positive",
        tanimoto = if (length(old_pos) == 0 || length(new_pos) == 0) NA_real_
        else tanimoto(old_pos, new_pos)
      )
    }),
    lapply(top_k, function(k) {
      old_neg <- get_top_features(old_scores, k, "negative")
      new_neg <- get_top_features(new_scores, k, "negative")
      
      data.frame(
        column = paste0("Column ", col_index),
        k = k,
        direction = "Negative",
        tanimoto = if (length(old_neg) == 0 || length(new_neg) == 0) NA_real_
        else tanimoto(old_neg, new_neg)
      )
    })
  )
}

tanimoto_results <- bind_rows(
  calc_tanimoto_for_column(1,Wm_opt,Wm_opt_new),
  calc_tanimoto_for_column(2,Wm_opt,Wm_opt_new)
)

## Plot both columns as facets
p_tanimoto <- ggplot(
  tanimoto_results,
  aes(x = k, y = tanimoto, color = direction, group = direction)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ column) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(limits = c(0,1))+
  labs(
    x = "Top k features",
    y = "Tanimoto similarity",
    color = "Feature set",
    title = "Tanimoto similarity between W_opt and W_opt_new"
  ) +
  theme_bw()
print(p_tanimoto)
ggsave(
  filename = "../figures/tanimoto_similarity_W_opt_vs_W_opt_from_random_data.png",
  plot = p_tanimoto,
  width = 8,
  height = 4
)

## 2) GSEA-based distance
library(GeneExpressionSignature)
library(Biobase)
source("distance_scores.R")
gsea_thresholds <- c(30,50,100,200,300,400,500,600,700,800,900,1000)

calc_gsea_distance_for_column <- function(col_index, W_opt,W_opt_new,thresholds = gsea_thresholds) {
  
  num_table <- cbind(
    W_opt     = W_opt[, col_index],
    W_opt_new = W_opt_new[, col_index]
  )
  
  rownames(num_table) <- rownames(W_opt)
  
  bind_rows(lapply(thresholds, function(thres) {
    
    dist_mat <- distance_scores(
      num_table = num_table,
      threshold_count = thres,
      names = colnames(num_table)
    )
    
    data.frame(
      column = paste0("Column ", col_index),
      threshold = thres,
      gsea_distance = dist_mat["W_opt", "W_opt_new"]
    )
  }))
}

gsea_distance_results <- bind_rows(
  calc_gsea_distance_for_column(1,Wm_opt,Wm_opt_new),
  calc_gsea_distance_for_column(2,Wm_opt,Wm_opt_new)
)

p_gsea <- ggplot(
  gsea_distance_results,
  aes(x = threshold, y = gsea_distance, group = column, color = column)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(limits = c(0,2))+
  labs(
    x = "Threshold count",
    y = "GSEA-based distance",
    color = "Matrix column",
    title = "GSEA-based distance between W_opt and W_opt_new"
  ) +
  theme_bw()

print(p_gsea)
ggsave(
  filename = "../figures/gsea_distance_W_opt_vs_W_opt_new.png",
  plot = p_gsea,
  width = 7,
  height = 4
)

Yhat1 <- predict(liv2trans_results$model,Xh %*% Wm %*% t(Wm))
cor(Yh,Yhat1,method='spearman')
Yhat2 <- predict(liv2trans_results_new$model,Xh %*% Wm_new %*% t(Wm_new))
cor(Yh,Yhat2,method='spearman')
cor(Yhat1,Yhat2,method='spearman')

# Train PLSR once to get phi_PLS
model <- liv2trans_results$model
phi <- coef(model)  # or however you extract W_h %*% B_h
# For each in vitro dataset:
rho_real <- sum((Wm[rownames(phi),] %*% (t(Wm[rownames(phi),]) %*% phi))^2) / sum(phi^2)
rho_random <- sum((Wm_new[rownames(phi),] %*% (t(Wm_new[rownames(phi),]) %*% phi))^2) / sum(phi^2)
cat("Fraction of phi captured by real MPS:", rho_real, "\n")
cat("Fraction of phi captured by random:", rho_random, "\n")

