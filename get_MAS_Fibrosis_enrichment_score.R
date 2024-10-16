# Supplementary script for gene signature enrichment using VIPER

### Load data and packages
library(decoupleR)
library(dplyr)
library(ggpubr)

root_dir <- "E:/Jose Luis/Documents/GitHub/FattyLiverModeling"
setwd(root_dir)
source("./utils/plotting_functions.R")
source('modeling/functions_translation.R')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(c(ref_dataset, target_dataset), dir_data = 'data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Extract gene expression
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Extract phenotype
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS

### Naive comparison - correlate all MPS samples to all human samples based on normalized and centered RNAseq
cor_Xh_Xm <- cor(t(Xh), t(Xm))

cor_Xh_Xm <- cor_Xh_Xm %>% 
  as.data.frame() %>% 
  mutate(human_sample = rownames(.)) %>%
  pivot_longer(cols = 1:nrow(Xm), names_to = "mps_sample", values_to = "cor")


cor_Xh_Xm <- merge(cor_Xh_Xm, 
                   data_frame(human_sample = rownames(Xh), 
                              NAS = Yh[,1], fibrosis = Yh[,2]), 
                   by = "human_sample")

# Group by MPS sample and get the median NAS and fibrosis for the top 10 correlations
cor_Xh_Xm_top <- cor_Xh_Xm %>% 
  group_by(mps_sample) %>% 
  top_n(10, cor) %>% 
  summarize(NAS = median(NAS), fibrosis = median(fibrosis)) %>%
  mutate(X = mps_sample)

cor_Xh_Xm_top <- merge(cor_Xh_Xm_top, data_list[[target_dataset]]$metadata, by = "X")

# Make plots
plt_cor_NAS <- cor_Xh_Xm_top %>%
  ggplot(aes(x = NAS, y = reorder(Combination_name, NAS))) +
  geom_boxplot(size = size_line) +
  labs(x = "Top correlated human MAS", y = NULL) +
  xlim(0,8)

plt_cor_NAS <- add_theme(plt_cor_NAS)

plt_cor_fib <- cor_Xh_Xm_top %>%
  ggplot(aes(x = fibrosis, y = reorder(Combination_name, fibrosis))) +
  geom_boxplot(size = size_line) +
  labs(x = "Top correlated human Fibrosis", y = NULL) +
  xlim(0,4)

plt_cor_fib <- add_theme(plt_cor_fib)

# Export
ggsave('./figures/plt_cor_NAS.pdf', 
       plot = add_theme(plt_cor_NAS),
       device = "pdf",
       height = 6,
       width=7,
       units = 'cm')

ggsave('./figures/plt_cor_fib.pdf', 
       plot = add_theme(plt_cor_fib),
       device = "pdf",
       height = 6,
       width=7,
       units = 'cm')

### Look for enrichment using genes from Hoang et al and VIPER

# Signatures
gNAS <- c("ARPC5", "YWHAH", "ARF4", "TNFRSF12A", "ADHFE1",
          "USP33", "CD52", "ACVR2B", "ING5", "ASB3", "IFI30",
          "ERVW-1", "YWHAZ", "ERBB3", "KPNA2", "COQ10B", "MAGI1",
          "MAPRE1", "ABCA6")

gFib <- c("TIMP1", "MYLI2B", "LUM", "ZNF395", "AKAP9",
          "ACTR2", "LGALS3", "MAPRE1", "FRK", "ANKRD28",
          "IGFBP7", "YWHAZ", "USP33", "CD59", "TAX1BP3",
          "FAM221A", "ADHFE1", "TNFRSF12A")

# Build regulon-like object to do non-parametric enrichment with VIPER
sig_regulon <- rbind(
  data.frame(source = "MAS", target = gNAS, mor = c(1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1)),
  data.frame(source = "Fibrosis stage", target = gFib, mor = c(1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,1,-1,-1,-1)))


# Run on human data
sig_viper <- decoupleR::run_viper(t(Xh), net = sig_regulon)
colnames(sig_viper)[3] <- "GEO_Accession..exp."

sig_viper <- merge(sig_viper, data_list[[ref_dataset]]$metadata[, c("GEO_Accession..exp.", "Fibrosis_stage", "nas_score")], by = "GEO_Accession..exp.")

sig_viper <- sig_viper %>% 
  pivot_longer(cols = c("Fibrosis_stage", "nas_score"), 
               names_to = "Phenotype", 
               values_to = "val") %>%
  mutate(Phenotype = factor(ifelse(Phenotype == "nas_score", "MAS", "Fibrosis stage")))

# Run on MPS data
sig_viper_MPS <- decoupleR::run_viper(t(Xm), net = sig_regulon)
colnames(sig_viper_MPS)[3] <- "X"
sig_viper_MPS <- merge(sig_viper_MPS, data_list[[target_dataset]]$metadata, by = "X")

# Plot
plt_sig_viper <- sig_viper %>%
  ggplot(aes(x = val, y = score)) +
  geom_jitter(width = 0.1, color = "steelblue", size = 0.8*size_dot) +
  facet_wrap(~Phenotype, scales = "free") +
  geom_smooth(method = "lm", se = T, color = "indianred", linewidth = size_line, linetype = 2) +
  labs(x = "Histological score", y = "Gene signature enrichment score") +
  stat_cor(size = size_annotation*0.8)

plt_sig_viper_MPS <- sig_viper_MPS %>%
  ggplot(aes(x = score, y = reorder(Combination_name, score))) +
  geom_boxplot(size = size_line) +
  facet_wrap(~source) +
  labs(x = "Gene signature enrichment score", y = NULL)

plt_sig_viper <- add_theme(plt_sig_viper)
plt_sig_viper_MPS <- add_theme(plt_sig_viper_MPS)


ggsave('./figures/plt_sig_viper.pdf', 
       plot = add_theme(plt_sig_viper),
       device = "pdf",
       height = 6,
       width=8,
       units = 'cm')

ggsave('./figures/plt_sig_viper_MPS.pdf', 
       plot = add_theme(plt_sig_viper_MPS),
       device = "pdf",
       height = 6,
       width=8,
       units = 'cm')
