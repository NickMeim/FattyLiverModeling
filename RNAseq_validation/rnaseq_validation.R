# Script for processing RNAseq data obtained from liver spheroids
library(tidyverse)
library(matrixStats)
source("utils/plotting_functions.R")
source("modeling/vector_space_interpretation.R")

# Read count matrix
data_counts <- read.table("RNAseq_validation/2C72ND-expression-matrix.tsv", header = T, sep = "\t")
sample_key <- read.table("RNAseq_validation/sampleIDs.txt", header = T, sep = "\t")

# Load processed data from main paper
processed_data <- readRDS("results/processed_data_list_govaere_kostrzewski.rds")
Xh <- processed_data$Xh
Yh <- processed_data$Yh
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

# Load LIV2TRANS results - specifically rotation matrices
Wm_TC <- readRDS("results/Wm_kostrzewski_combo.rds")
Wm_extra <- readRDS("results/Wm_kostrzewski_extra.rds")
colnames(Wm_extra) <- c("extraLV1", "extraLV2")
genes_Wm <- rownames(Wm_TC)

# Process RNAseq data:

# Use genenames; aggregate counts per gene
data_counts <- data_counts[,-c(1,3)]
data_counts <- data_counts %>% 
                group_by(gene_name) %>%
                summarize_all(sum) %>%
                as.data.frame()
# Subset to genes used in the LIV2TRANS part
data_counts <- data_counts %>% filter(gene_name %in% genes_Wm)
genes <- data_counts$gene_name
# Keep only counts, not cpms
data_counts <- data_counts[,grep("count", colnames(data_counts))]
rownames(data_counts) <- genes
colnames(data_counts) <- gsub("_count","",colnames(data_counts))
# Some genes not found; add them with zero value so that we can do projection/backprojection
data_counts <- data_counts[genes_Wm,]
rownames(data_counts) <- genes_Wm
data_counts[is.na(data_counts)] <- 0

# Convert to log2 cpm + 1
data_cpm <- apply(data_counts, MARGIN = 2, FUN = function(x){log2(1 + x/sum(x)*1e6)})

# PCA
Xnew <- data_cpm %>% t() %>% scale(., center = T, scale = F)
pca_rnaseq <- prcomp(Xnew)
Wnew <- pca_rnaseq$rotation
per_var <- round(100*pca_rnaseq$sdev^2/sum(pca_rnaseq$sdev^2), digits = 2)
# Plot PCA
plt_PCA_rnaseq <- data.frame(pca_rnaseq$x, 
                             sample_key) %>%
                  ggplot(aes(x = PC1, y = PC2, fill=condition)) +
                  geom_point(size = size_dot, shape = 21, stroke = size_stroke, color = "black")+
                  xlab(paste0("PC1 (", per_var[1],"%)")) + ylab(paste0("PC2 (", per_var[2],"%)"))

plt_PCA_rnaseq <- add_theme(plt_PCA_rnaseq) + scale_fill_brewer(palette = "Dark2")
ggsave("figures/spheroid_RNAseq_PCA.pdf", plot = plt_PCA_rnaseq, units = "cm", height = 5, width = 6)

# Project new data onto TC and LV
Xnew_project <-  Xnew %*% cbind(Wm_TC, Wm_extra)
per_var_project <- round(100*colVars(Xnew_project)/sum(pca_rnaseq$sdev^2), digits = 2)

plt_project <- data.frame(x = Xnew_project[,1], 
                          y = Xnew_project[,3], 
                          sample_key) %>%
                ggplot(aes(x = x, y = y, fill=condition)) +
                geom_point(size = size_dot, shape = 21, stroke = size_stroke, color = "black")+
                xlab(paste0("TC1 (", per_var_project[1],"%)")) + ylab(paste0("LV extra 1 (", per_var_project[3],"%)"))

plt_project <- add_theme(plt_project) + scale_fill_brewer(palette = "Dark2")
ggsave("figures/spheroid_RNAseq_projection.pdf", plot = plt_project, units = "cm", height = 5, width = 6)


# Interpret PCs to compare with LIV2TRANS latent variables

# Get pathway activity from PC loadings
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
pathway_acts <- decoupleR::run_viper(cbind(Wm_TC, Wm_extra, Wnew), net_prog, minsize = 1,verbose = TRUE)
pathway_acts_wide <- pathway_acts %>% pivot_wider(names_from = source, values_from = score, id_cols = condition) %>% as.data.frame()
rownames(pathway_acts_wide) <- pathway_acts_wide$condition
pathway_acts_wide <- pathway_acts_wide[,-1] %>% t() %>% as.data.frame()
pathway_acts_wide$pwy <- rownames(pathway_acts_wide)
# Plot
plt_pwy_PC1 <- pathway_acts %>% filter(condition %in% c("PC1")) %>%
  ggplot(aes(x = score,y = reorder(source, score),fill = score)) + 
  geom_bar(stat = 'identity', show.legend = F, color = "black", linewidth = size_col) +
  scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                       midpoint = 0)+
  labs(y = "Pathway", x = "Activity score", fill = "Score")

plt_pwy_PC2 <- pathway_acts %>% filter(condition %in% c("PC2")) %>%
  ggplot(aes(x = score,y = reorder(source, score),fill = score)) + 
  geom_bar(stat = 'identity', show.legend = F, color = "black", linewidth = size_col) +
  scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                       midpoint = 0)+
  labs(y = "Pathway", x = "Activity score", fill = "Score")

plt_pwy_PC1 <- add_theme(plt_pwy_PC1)
plt_pwy_PC2 <- add_theme(plt_pwy_PC2)
ggsave("figures/spheroid_RNAseq_pwyPC1.pdf", plot = plt_pwy_PC1, units = "cm", height = 4.5, width = 4)
ggsave("figures/spheroid_RNAseq_pwyPC2.pdf", plot = plt_pwy_PC2, units = "cm", height = 4.5, width = 4)
