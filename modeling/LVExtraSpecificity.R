library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
source('CrossValidationUtilFunctions.R')
source('functions_translation.R')
source("../utils/plotting_functions.R")
source("vector_space_interpretation.R")

## Initial PLSR and LIV2TRANS run-------------------------------
# Load data and compute extra basis as in main_example_run.R
# Xh: human gene-expression matrix (samples × genes)
# Yh: vector of MASLD or fibrosis scores
# Wm: PC rotation matrix of MPS (genes × nPC)
# Wh: human PLSR loading matrix (genes × nComp)
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
data_list <- load_datasets(dataset_names, dir_data = '../data/')
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
Xm_grouped <- data_list$Kostrzewski$Xm_grouped
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()
# Get Wh of PLSR
plsr_model <- opls(x = Xh, 
                   y = Yh,
                   predI = 8,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
rownames(Wh) <- colnames(Xh)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)){
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
}
# Get regression coefficients
Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
phi <- Wh %*% Bh

# Compute extra basis (two extra latent variables)
Wm_opt <- analytical_solution_opt(y=Yh,
                                  W_invitro = Wm,
                                  phi = phi)
rownames(Wm_opt) <- rownames(Wm)
colnames(Wm_opt) <- c("LV_extra1", "LV_extra2")

# Append extra basis to MPS basis and evaluate baseline performance
Wm_tot  <- cbind(Wm, Wm_opt)
# Predicted scores using PLSR after backprojection
Ypred   <- cbind(1, Xh %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
avg_spearman <- mean(diag(cor(Yh,Ypred,method='spearman')))

# Randomization test: generate random orthonormal bases of the same size-----------------
k <- ncol(Wm_opt)
nrep <- 1000  # number of random draws
spear_random <- numeric(nrep)
max_cors <- NULL
pb <- txtProgressBar(min = 0, max = nrep, style = 3)
for (i in 1:nrep) {
  Z <- matrix(rnorm(nrow(Wm) * k), nrow = nrow(Wm))
  # remove component in span(Wm)
  Z_orth <- Z - Wm %*% (t(Wm) %*% Z)
  # orthonormalize what's left
  rand_basis <- qr.Q(qr(Z_orth))[, 1:k]
  Wm_rand_tot <- cbind(Wm, rand_basis)
  max_cor <- cor(Wm, rand_basis)
  max_cors[i] <- max(abs(max_cor[upper.tri(max_cor)]))
  # Evaluate the model
  Ypred_rand    <- cbind(1, Xh %*% Wm_rand_tot %*% t(Wm_rand_tot) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
  spear_random[i]  <- mean(diag(cor(Yh,Ypred_rand,method='spearman')))
  setTxtProgressBar(pb, i)
}
close(pb)
hist(max_cors)

# Compare baseline to random distribution
p_value <- mean(spear_random >= avg_spearman)  # empirical false-positive rate
cat("Baseline Spearman`s correlation with extra LVs:", avg_spearman, "\n")
cat("Mean Spearman`s correlation for random bases:", mean(spear_random), "\n")
cat("Fraction of random draws exceeding baseline (p‑value):", p_value, "\n")

# Random gene perturbations------------
# fast score (no huge projector)
Bfull <- rbind(apply(Yh, 2, mean), Bh)
score_spearman_basis <- function(Wm_tot, Xh, Wh, Yh, Bfull) {
  # Xh %*% Wm_tot %*% t(Wm_tot) %*% Wh computed efficiently
  Z <- Xh %*% Wm_tot            # samples x r
  M <- t(Wm_tot) %*% Wh         # r x q
  F <- Z %*% M                  # samples x q
  Ypred <- cbind(1, F) %*% Bfull
  mean(diag(cor(Yh, Ypred, method = "spearman")))
}
augment_with_two_conditions <- function(Xm_sg, genes1, genes2, infl) {
  mu <- colMeans(Xm_sg)
  sdg <- apply(Xm_sg, 2, sd)
  
  new1 <- mu
  new2 <- mu
  
  # perturb around baseline mean using SD scale (mean preserved)
  new1[genes1] <- mu[genes1] + infl * sdg[genes1]
  new2[genes2] <- mu[genes2] + infl * sdg[genes2]
  
  rbind(Xm_sg, new1, new2)  # adds 2 new "conditions"
}
recompute_Wm_from_Xm <- function(Xm_sg, nPC) {
  # center genes
  Xc <- scale(Xm_sg, center = TRUE, scale = FALSE)
  pca <- prcomp(Xc, scale. = FALSE)
  pca$rotation[, 1:nPC, drop = FALSE]  # genes x nPC
}


## Parameters
k_extra  <- ncol(Wm_opt)
n_size   <- 40
n_infl   <- 3
n_pools  <- 10
infl_min <- 1.1
infl_max <- 5
nPC_use <- ncol(Wm) + k_extra

## Align Xm_grouped to samples x genes (conditions x genes)
G <- nrow(Wm)
if (ncol(Xh) != G) stop("Expected Xh samples×genes with ncol(Xh)==nrow(Wm).")

if (ncol(Xm_grouped) == G) {
  Xm_grouped_sg <- Xm_grouped
} else if (nrow(Xm_grouped) == G) {
  Xm_grouped_sg <- t(Xm_grouped)
} else {
  stop("Cannot align Xm_grouped to genes in Wm.")
}

max_genes <- ncol(Xm_grouped_sg)

# Ensure gene-set size at least 2 so we can split into two pools
subset_sizes <- pmax(2L, sample.int(max_genes, n_size, replace = TRUE))

# storage (include pool_idx)
out <- expand.grid(size_idx = 1:n_size, pool_idx = 1:n_pools, infl_idx = 1:n_infl)
out$ngenes   <- NA_integer_
out$inflate  <- NA_real_
out$spearman <- NA_real_

pb <- txtProgressBar(min=0, max=nrow(out), style=3)
rowi <- 0

for (s in 1:n_size) {
  ng <- subset_sizes[s]
  
  for (p in 1:n_pools) {
    # sample ONE pool of ng genes, then split into two sets => two independent perturbation directions
    genes_s <- sample.int(max_genes, ng)
    half <- floor(ng / 2)
    genes1 <- genes_s[1:half]
    genes2 <- genes_s[(half + 1):ng]
    # if ng==2, each gets 1 gene; if ng is odd, genes2 has one more gene.
    
    for (j in 1:n_infl) {
      rowi <- rowi + 1
      infl <- runif(1, infl_min, infl_max)
      
      # Add two new conditions
      Xm_aug <- augment_with_two_conditions(Xm_grouped_sg, genes1, genes2, infl)
      
      # Recompute PCA basis: this Wm_aug is now your "Wm_tot"
      # Cap to feasible PCs:
      maxPC <- min(nrow(Xm_aug) - 1L, ncol(Xm_aug))
      nPC_eff <- min(nPC_use, maxPC)
      
      Wm_tot <- recompute_Wm_from_Xm(Xm_aug, nPC = nPC_eff)
      
      out$ngenes[rowi]   <- ng
      out$inflate[rowi]  <- infl
      out$spearman[rowi] <- score_spearman_basis(Wm_tot, Xh, Wh, Yh, Bfull)
      
      setTxtProgressBar(pb, rowi)
    }
  }
}
close(pb)

summary(out$spearman)

## Visualize all---------
# Combine into one long data frame
plot_df <- bind_rows(
  tibble(approach = "Random orthonormal basis",      spearman = spear_random),
  tibble(approach = "Gene perturbation", spearman = out$spearman)
) %>%
  mutate(
    approach = factor(approach, levels = c("Random orthonormal basis", "Gene perturbation"))
  )

# Compute "p-values" relative to the same baseline (avg_spearman)
pvals_df <- plot_df %>%
  group_by(approach) %>%
  summarise(
    p_value = mean(spearman >= avg_spearman, na.rm = TRUE),
    y_max   = max(spearman, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # place label a bit above the top of each violin
  mutate(
    y_lab = y_max + 0.05 * diff(range(plot_df$spearman, na.rm = TRUE)),
    label = sprintf("p-value = %.3g", p_value),
    size=5
  )

# Plot
p <- ggplot(plot_df, aes(x = approach, y = spearman)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.25, size = 1.0) +
  geom_hline(yintercept = avg_spearman, linetype = "dashed", linewidth = 1.0) +
  geom_text(
    data = pvals_df,
    aes(x = approach, y = y_lab, label = label),
    inherit.aes = FALSE,
    vjust = 0
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Performance after backprojection with random perturbations",
    subtitle = sprintf("Dashed line = average Spearman correlation using LV extra = %.3f)", avg_spearman),
    x = NULL,
    y = "Spearman correlation"
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(5.5, 25, 5.5, 5.5),
    axis.text = element_text(size=12),
    axis.title = element_text(size=16)
  )

print(p)
ggsave("../figures/supplementary_randomization_basis_test_performance.png",
       plot = p, 
       width = 15,
       height = 9,
       units = 'cm',
       dpi = 600)
ggsave("../figures/supplementary_randomization_basis_test_performance.pdf",
       plot = p, 
       width = 15,
       height = 9,
       units = 'cm',
       dpi = 600)
