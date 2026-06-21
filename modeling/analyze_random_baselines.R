# This script compares the results from LIV2TRANS with the expected values (and empirically distributed values) obtained with random orthonormal matrices
# These results are part of a response to reviewers and inform a supplemental note in the manuscript
library(dplyr)
library(ggplot2)
library(matrixStats)
source('modeling/functions_translation.R')
source("utils/plotting_functions.R")

set.seed(123)
## Load the matrices obtained with LIV2TRANS for the main example
processed_data <- readRDS("results/processed_data_list_govaere_kostrzewski.rds")
Xh <- processed_data$Xh
Xm <- processed_data$Xm
Yh <- processed_data$Yh
Wm <- processed_data$Wm
# Load PLSR model
plsr_model <- readRDS("results/PLSR_model_govaere.rds")
# Get Wh of PLSR
Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
rownames(Wh) <- colnames(Xh)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)){
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
}
# Get regression coefficients
Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
phi <- Wh %*% Bh
# Load the grouped CN Bio dataset
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = 'data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
Xm_group <- data_list$Kostrzewski$Xm_grouped %>% t()
# Load LIV2TRANS results - specifically rotation matrices
Wm_TC <- readRDS("results/Wm_kostrzewski_combo.rds")
Wm_extra <- readRDS("results/Wm_kostrzewski_extra.rds")

## For a different number of observations (up to 50), make random orthonormal matrices and get stats from projection/backprojection
# Do not try to run this for more than 5000 observations - the othonormalization might run out of memory

df_rand <- data.frame()
n_feature <- nrow(Wm)
n_reps <- 20

for (n_obs in seq(1,51,by=5)){
  print(n_obs)
  for (rr in 1:n_reps){
   # Directly make a random orthonormal basis
    random_matrix <- matrix(rnorm(n_feature * n_obs), nrow = n_feature, ncol = n_obs)
    q_matrix <- qr.Q(qr(random_matrix))
    
    # Project/bakproject to get phenotypic scores
    Yhat <- (Xh %*% q_matrix) %*% (t(q_matrix) %*% phi)
    # Project to get variance of the data on random basis
    Sfull <- Xh %*% q_matrix
    Strunc <- Xh %*% Wh %*% t(Wh) %*% q_matrix
    # Append
    df_rand <- rbind(df_rand,
                    data.frame(n_obs = n_obs,
                               var_Sfull = sum(colVars(Sfull)),
                               var_Strunc = sum(colVars(Strunc)),
                               var_NAS = var(Yhat[,1]),
                               var_Fib = var(Yhat[,2]),
                               cor_NAS = cor(Yh[,1], Yhat[,1]),
                               cor_Fib = cor(Yh[,2], Yhat[,2]),
                               cov_NAS = cov(Yh[,1], Yhat[,1]),
                               cov_Fib = cov(Yh[,2], Yhat[,2])))
  }
}

# Now pic subsets of the MPS dataset to get their stats when doing projection/backprojection using different number of samples picked
sample_sizes <- seq(2, nrow(Xm_group), by = 4)  
n_partitions <- 10

df_results <- data.frame()

for (size in sample_sizes) {
  print(paste("Testing sample size:", size))
  for (partition in 1:n_partitions) {
    # Random subset of Xm
    random_idx <- sample(1:nrow(Xm_group), size, replace=FALSE)
    Xm_subset <- Xm_group[random_idx, ]
    
    # PCA on subset
    W <- prcomp(Xm_subset, center=TRUE, scale.=FALSE)$rotation
    
    # Evaluate on held-out data
    Ypred <- Xh %*% W %*% t(W) %*% phi
    
    # Append
    df_results <- rbind(df_results, data.frame(
      n_obs = size - 1, #subtract one since the number of PCs is one less than the number of observations
      partition = partition,
      cov_NAS = cov(Yh[,1], Ypred[,1]),
      cov_Fib = cov(Yh[,2], Ypred[,2]),
      var_Sfull = colVars(Xh %*% W) %>% sum(),
      var_Strunc = colVars((Xh %*% Wh) %*% t(Wh) %*% W) %>% sum()
    ))
  }
}

# Make plots


# Plot covariance between Yh and Yh pred for NAS score
plt_cov_NAS <- df_rand %>% ggplot(aes(x = n_obs, y = abs(cov_NAS))) + 
  geom_point(alpha = 0.5, color = "black", size = size_dot*0.75) + # Add empirical data for random matrices
  geom_line(data = data.frame(x = 1:51, y = 1:51/n_feature*cov(Yh[,1], Xh %*% phi[,1])), mapping = aes(x = x,y=y), color = "black", linewidth = size_line) + # Add analytical line for expected value
  geom_point(data = df_results, color = "steelblue", size = size_dot*0.75) + # Add empirical results from subsets of Xm data 
  stat_summary(data = df_results, color = "steelblue", geom = "line", fun = mean, size = size_line) +  #Plot means
  geom_hline(yintercept =  cov(Yh[,1], (Xh %*% Wm) %*% (t(Wm) %*% phi[,1])), color = "indianred", linetype = 2, linewidth = size_line) + # Covariance between prediction with all PCs and true data
  geom_hline(yintercept =  cov(Yh[,1], Xh %*% phi[,1]), color = "darkgreen", linetype = 2, linewidth = size_line) + # Covariance between prediction with without backprojection 
  scale_y_log10()  + scale_x_log10() + labs(x = "Number of principal components", y = "cov(Yh, Yh predicted)") + theme_bw()

# Plot covariance between Yh and Yh pred for Fib score
plt_cov_Fib <- df_rand %>% ggplot(aes(x = n_obs, y = abs(cov_Fib))) + 
  geom_point(alpha = 0.5, color = "black", size = size_dot*0.75) + # Add empirical data for random matrices
  geom_line(data = data.frame(x = 1:51, y = 1:51/n_feature*cov(Yh[,2], Xh %*% phi[,2])), mapping = aes(x = x,y=y), color = "black", linewidth = size_line) + # Add analytical line for expected value
  geom_point(data = df_results, color = "steelblue", size = size_dot*0.75) + # Add empirical results from subsets of Xm data 
  stat_summary(data = df_results, color = "steelblue", geom = "line", fun = mean, size = size_line) +  #Plot means
  geom_hline(yintercept =  cov(Yh[,2], (Xh %*% Wm) %*% (t(Wm) %*% phi[,2])), color = "indianred", linetype = 2, linewidth = size_line) + # Covariance between prediction with all PCs and true data
  geom_hline(yintercept =  cov(Yh[,2], Xh %*% phi[,2]), color = "darkgreen", linetype = 2, linewidth = size_line) + # Covariance between prediction with without backprojection 
  scale_y_log10()  + scale_x_log10() + labs(x = "Number of principal components", y = "cov(Yh, Yh predicted)") + theme_bw()

# Plot captured data variance for Xh %*% Wh
plt_var_rand <- df_rand %>% ggplot(aes(x = n_obs, y = var_Strunc)) + 
  geom_point(alpha = 0.5, color = "black", size = size_dot*0.75) + # Add empirical data for random matrices
  geom_line(data = data.frame(x = 1:51, y = 1:51/n_feature*sum(colVars(Xh %*% Wh))), mapping = aes(x = x,y=y), color = "black", linewidth = size_line) + # Add analytical line for expected value
  geom_point(data = df_results, color = "steelblue", size = size_dot*0.75) + # Add empirical results from subsets of Xm data 
  stat_summary(data = df_results, color = "steelblue", geom = "line", fun = mean, size = size_line) +  #Plot means
  geom_hline(yintercept =  sum(colVars(Xh %*% Wh)), color = "darkgreen", linetype = 2, linewidth = size_line) + # Covariance between prediction with without backprojection 
  scale_y_log10()  + scale_x_log10() + labs(x = "Number of principal components", y = "Total variance (Xh*W)") + theme_bw()


# Now make random matrices of the size of Wm to calculate their TCs and extraLVs
n_obs <- ncol(Wm)
cos_dist_TC <- rep(0,100)
cos_dist_extra <- cos_dist_TC
rand_TC <- NULL
rand_extra <- NULL
for (ii in 1:100){
  random_matrix <- matrix(rnorm(n_feature * n_obs), nrow = n_feature, ncol = n_obs)
  q_matrix <- qr.Q(qr(random_matrix))
  # Get TCs
  q_TC <- get_translatable_LV_2phenotype(Xh, Yh, Wh, q_matrix, Bh)$Wm_TC
  # Get cosine distance between actual TC1 and random TC1
  cos_dist_TC[ii] <- t(q_TC[,1]) %*% Wm_TC[,1]
  
  # Get extraLV
  q_extra <- analytical_solution_opt(y=Yh, W_invitro = q_matrix, phi = phi)
  cos_dist_extra[ii] <- t(q_extra[,1]) %*% Wm_extra[,1]
  
  # Append first vectors just for fun to compare their similarity with each other
  rand_TC <- cbind(rand_TC, q_TC[,1])
  rand_extra <- cbind(rand_extra, q_extra[,1])
}

# Plot
pheatmap::pheatmap(t(rand_extra) %*% rand_extra)
pheatmap::pheatmap(t(rand_extra) %*% rand_TC)
data.frame(x = cos_dist_TC) %>% ggplot(aes(x = x)) + geom_histogram(fill = "steelblue") + labs(x = "cosine similarity WmTC - randTC") + theme_bw()
data.frame(x = cos_dist_extra) %>% ggplot(aes(x = x)) + geom_histogram(fill = "steelblue") + labs(x = "cosine similarity Wm extra - rand extra") + theme_bw()


