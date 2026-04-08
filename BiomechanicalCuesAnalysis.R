library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak)
library(patchwork)
library(ggforce)
library(ggradar)
library(ggiraphExtra)
library(grid)
library(ggExtra)
library(matrixStats)
source('modeling/CrossValidationUtilFunctions.R')
source('modeling/functions_translation.R')
source("utils/plotting_functions.R")
source("modeling/vector_space_interpretation.R")

## Load data ---------------------------------------------------------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset    <- "Govaere"
target_dataset <- "Kostrzewski"   # keep Kostrzewski as the system LIV2TRANS was built on
wang_dataset   <- "Wang"          # Wang = the dataset we want to project

data_list <- load_datasets(dataset_names, dir_data = 'data/')
tmp       <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list

# Reference (human) matrices
Yh <- as.matrix(data_list[[ref_dataset]]$metadata %>%
                  select(nas_score, Fibrosis_stage))
colnames(Yh) <- c('NAS', 'fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()

# Kostrzewski — the system LIV2TRANS learns the basis from
Xm  <- data_list[[target_dataset]]$data_center %>% t()
Wm  <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# Wang — the system we want to PROJECT
Xw  <- data_list[[wang_dataset]]$data_center %>% t()
Ww  <- data_list[[wang_dataset]]$Wm_group %>% as.matrix()

# Wang metadata for colouring plots later
meta_wang <- data_list[[wang_dataset]]$metadata

## PLSR model ---------------------------------------------------------------
plsr_model <- opls(x = Xh,
                   y = Yh,
                   predI = 8,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")

# Build Wh (genes × LVs)
Wh <- matrix(0, nrow = ncol(Xh), ncol = ncol(plsr_model@weightMN))
rownames(Wh) <- colnames(Xh)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)) {
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii, ]
}

# Regression coefficients
Bh     <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
Wh_red <- plsr_model@weightMN

## LIV2TRANS — Extra LVs (analytical solution) ------------------------------
phi    <- Wh %*% Bh
Wm_opt <- analytical_solution_opt(y = Yh,
                                  W_invitro = Wm,
                                  phi = phi)
rownames(Wm_opt) <- rownames(Wm)
colnames(Wm_opt) <- c('V1', 'V2')

## LIV2TRANS — Translatable Components (evolutionary algorithm) -------------
Wm_combo           <- get_translatable_LV_2phenotype(Xh, Yh, Wh, Wm, Bh)
Wm_combo           <- Wm_combo$Wm_TC
rownames(Wm_combo) <- rownames(Wm)
colnames(Wm_combo) <- c('V1', 'V2')

## Project Wang into Extra LVs and TCs -------------------------------------
# Align Wang genes to the basis (only shared genes)
shared_genes_opt   <- intersect(colnames(Xw), rownames(Wm_opt))
shared_genes_combo <- intersect(colnames(Xw), rownames(Wm_combo))

# Project: samples × 2 score matrices
Zw_opt   <- Xw[, shared_genes_opt]   %*% Wm_opt[shared_genes_opt, ]
colnames(Zw_opt) <- c('extraLV1','extraLV2')
Zw_combo <- Xw[, shared_genes_combo] %*% Wm_combo[shared_genes_combo, ]
colnames(Zw_combo) <- c('TC1','TC2')

# Tidy: attach Wang metadata
df_wang_opt <- as.data.frame(Zw_opt) %>%
  rownames_to_column('sample') %>%
  left_join(meta_wang %>% mutate(sample = paste0('X',filename)), by = 'sample')

df_wang_combo <- as.data.frame(Zw_combo) %>%
  rownames_to_column('sample') %>%
  left_join(meta_wang %>% mutate(sample = paste0('X',filename)), by = 'sample')

## Variance captured in Wang by each space ---------------------------------
var_opt   <- round(sum(colVars(Zw_opt))   / sum(colVars(Xw[, shared_genes_opt]))   * 100, 2)
var_combo <- round(sum(colVars(Zw_combo)) / sum(colVars(Xw[, shared_genes_combo])) * 100, 2)
cat("Variance in extra LVs (Wang): ", var_opt,   "%\n")
cat("Variance in TCs       (Wang): ", var_combo, "%\n")

## Visualise Wang in Extra LV space ----------------------------------------
# Adjust colour_by to a column that exists in your Wang metadata
# e.g. condition, treatment, dose, time_point etc.
p_wang_opt <- ggplot(df_wang_opt,
                     aes(x = extraLV1, y = extraLV2,
                         fill  = scaffold)) +
  geom_point(size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_viridis_d() +
  labs(title = paste0("Wang projected into Extra LVs(",
                      var_opt, "% variance captured)"),
       x    = "Extra LV1",
       y    = "Extra LV2") +
  theme_bw()
print(p_wang_opt)

## Visualise Wang in TC space -----------------------------------------------
p_wang_combo <- ggplot(df_wang_combo,
                       aes(x = TC1, y = TC2,
                           fill  = scaffold)) +
  geom_point(size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_viridis_d() +
  labs(title = paste0("Wang projected into TCs(",
                      var_combo, "% variance captured)"),
       x    = "TC1",
       y    = "TC2") +
  theme_bw()
print(p_wang_combo)

## Side-by-side comparison --------------------------------------------------
p_wang_opt + p_wang_combo +
  plot_annotation(title = "Wang et al. in vitro — LIV2TRANS projection")

## Radial analysis: which direction in extra LV space predicts NAS/Fibrosis?
# (same logic as the example's df_correlation_radial section)
Zh_wang <- as.data.frame(Zw_opt)
Zh_wang$NAS      <- NA    # Wang has no Yh — leave NA or fill if available
Zh_wang$fibrosis <- NA

thetas   <- seq(0, 360, 5)
df_radial_wang <- data.frame()

for (theta in thetas) {
  u     <- as.matrix(c(cos(theta * pi / 180), sin(theta * pi / 180)))
  Wproj <- Wm_opt %*% u     # project basis to 1D direction
  
  # Score each Wang sample along this direction
  proj_scores <- Xw[, shared_genes_opt] %*% Wproj
  
  # Variance captured in this direction
  var_dir <- var(proj_scores) / sum(colVars(Xw[, shared_genes_opt]))
  
  df_radial_wang <- rbind(df_radial_wang,
                          data.frame(theta    = theta,
                                     var_dir  = var_dir,
                                     proj_sd  = sd(proj_scores)))
}

# Plot: how much variance does each angular direction capture in Wang?
p_radial_wang <- ggplot(df_radial_wang,
                        aes(x = theta, y = var_dir * 100)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_area(fill  = "steelblue", alpha = 0.2) +
  scale_x_continuous(breaks = seq(0, 360, 45)) +
  labs(title = "Wang variance by direction in Extra LV space",
       x     = "Angle θ (degrees)",
       y     = "% variance captured") +
  theme_bw()
print(p_radial_wang)