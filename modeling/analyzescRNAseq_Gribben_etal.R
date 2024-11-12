################################################################################
#### Scripts for processing single cell data from GSE202379 and map it to 
#### translatable and extra latent variables

library(tidyverse)
library(corrplot)
source("./utils/plotting_functions.R")
target_dataset <- "Kostrzewski"
### Part 1: Open "GSE202379_SeuratObject_AllCells.rds" downloaded from GEO and 
### pseudobulk based on the given cell type annotation and patient ID.
### Only run if you have the .rds file. We are including the pseudobulked file 
### in the repo since this part is slow

library(Seurat)
seur <- readRDS(file = "data/GSE202379_SeuratObject_AllCells.rds")

features_meta <- c("orig.ident", "Patient.ID", "Disease.status",
                   "Lobe", "SAF.Score",
                   "Steatosis", "Ballooning", "Inflammation",
                   "Fibrosis.score..F0.4.",
                   "Gender", "cell.annotation")

metadata <- seur@meta.data[, features_meta] %>%
  mutate(Group = paste0(Patient.ID,"_",cell.annotation))

# Pseudo-bulk by patient and cell type
count_mat <- seur@assays$RNA$counts

# Create groupings
groups <- metadata$Group %>% unique()

count_bulk <- matrix(0, nrow = nrow(count_mat), ncol = length(groups))
for (ii in 1:ncol(count_bulk)){
  idx <- which(metadata$Group == groups[ii])
  if (length(idx) > 1){
    count_bulk[,ii] <- rowSums(count_mat[, idx])
  } else {
    count_bulk[,ii] <- count_mat[, idx]
  }
}

colnames(count_bulk) <- groups
rownames(count_bulk) <- rownames(count_mat)

# Summarize metadata for groupings and count how many cells were added
metadata_group <- metadata %>%
  group_by(Patient.ID, cell.annotation, Group,
           Disease.status, SAF.Score, Steatosis, 
           Ballooning, Inflammation, Fibrosis.score..F0.4.) %>%
  summarize(nCells = n())

# Rename fibrosis score column
colnames(metadata_group)[colnames(metadata_group) == "Fibrosis.score..F0.4."] <- "Fibrosis.score"

# Add NAS score. Assuming that healthy control means zero for all scores
metadata_group <- metadata_group %>%
  mutate(NAS = ifelse(Steatosis %in% c("healthy control", "end stage"),
                      0,
                      as.numeric(Steatosis) + as.numeric(Ballooning) + as.numeric(Inflammation)))

metadata_group$NAS[metadata_group$Steatosis == "end stage"] <- "end stage"

# Group by patient only
metadata_patient <- metadata_group %>% 
  ungroup() %>%
  select(Patient.ID,
         Disease.status, SAF.Score, Steatosis, 
         Ballooning, Inflammation, Fibrosis.score, NAS) %>%
  unique()

# Organize columns to match metadata order
count_bulk <- count_bulk[,metadata_group$Group]

# Save
save(count_bulk, metadata_group, metadata_patient, file = "data/GSE202379_pseudobulk_data.RData")

################################################################################
### Part 2: Load pseudo-bulk data and process. Normalize as log2 + 1 cpm

load("data/GSE202379_pseudobulk_data.RData")
cpm_pseudo <- apply(count_bulk, MARGIN = 2, FUN = function(x){log2(1 + x/sum(x)*1e6)})
# Load extra latent variables in target dataset
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
colnames(Wm_opt) <- c("LVextra1", "LVextra2")

# Organize cpm in the same order and content as the LVs
cpm_center <- matrix(0, nrow = nrow(Wm_opt), ncol = ncol(cpm_pseudo))
rownames(cpm_center) <- rownames(Wm_opt)
for (ii in 1:ncol(cpm_pseudo)){
  cpm_center[rownames(cpm_pseudo)[rownames(cpm_pseudo) %in% rownames(Wm_opt)], ii] <- cpm_pseudo[rownames(cpm_pseudo)[rownames(cpm_pseudo) %in% rownames(Wm_opt)], ii]
}
# Center rows by cell type since we'll project by cell types
for (ii in 1:nrow(cpm_center)){
  for (jj in unique(metadata_group$cell.annotation)){
    idx <- metadata_group$cell.annotation == jj
    cpm_center[ii,idx] <- cpm_center[ii,idx] - mean(cpm_center[ii,idx])
  }
}

### Part 3: Project data onto LVs and plot
pseudobulk_project <- t(cpm_center) %*% Wm_opt
# Add metadata and make long table to plot - remove end stage because it's not covered
# in our original dataset and remove cell types with very few cells (<10) per patient
pseudobulk_project <- data.frame(pseudobulk_project, metadata_group) %>% 
  filter(Disease.status != "end stage" &
           !(cell.annotation %in% c("unknown", "B-cell 1", "B-cell 2", "Neutrophils")) &
           nCells > 10)

pseudobulk_project_long <- pseudobulk_project %>%
                            pivot_longer(cols = all_of(colnames(Wm_opt)), 
                                         names_to = "LV", values_to = "score") 


# For each cell type, get correlation between projected score and phenotypic score
pseudobulk_cor_MAS <- pseudobulk_project_long %>%
                      group_by(cell.annotation, LV) %>%
                      summarize(cor = cor.test(score, as.numeric(NAS))[['estimate']], 
                                cor_pval = cor.test(score, as.numeric(NAS))[['p.value']]) %>%
                      mutate(Pheno = "MAS")

pseudobulk_cor_Fib <- pseudobulk_project_long %>%
                      group_by(cell.annotation, LV) %>%
                      summarize(cor = cor.test(score, as.numeric(Fibrosis.score))[['estimate']], 
                                cor_pval = cor.test(score, as.numeric(Fibrosis.score))[['p.value']]) %>%
                      mutate(Pheno = "Fibrosis stage")

# Append and adjust pvalues within each LV and phenotype
pseudobulk_cor <- rbind(pseudobulk_cor_MAS, pseudobulk_cor_Fib) %>%
                  group_by(Pheno, LV) %>%
                  mutate(cor_pval_adj = p.adjust(cor_pval, method = 'hochberg')) %>%
                  ungroup() %>%
                  mutate(annot = ifelse(cor_pval_adj <= 0.0001, "****",
                                        ifelse(cor_pval_adj <= 0.001,"***",
                                               ifelse(cor_pval_adj<=0.01,"**",
                                                      ifelse(cor_pval_adj<=0.05,'*',
                                                             ifelse(cor_pval_adj<=0.1,'\u2022', #\u2219
                                                                    'ns'))))))


# Plot
plt_pseudobulk_cor <- pseudobulk_cor %>% 
                      ggplot(aes(y = cell.annotation, x = cor, fill = LV, label = annot)) +
                        geom_col(position = position_dodge(width = 0.9), size = size_col, color = 'black') +
                        geom_text(aes(x = ifelse(cor < 0, cor - 0.06, cor + 0.06)),
                                  color = 'black', size = size_annotation*0.7,
                                  angle = 0, show.legend = FALSE, hjust = "center", vjust = "center",
                                  position = position_dodge(width = 0.9)) +
                        facet_wrap(~Pheno) +
                        labs(x = 'Pearson correlation', y = NULL, fill = NULL) +
                        scale_fill_brewer(palette = "Dark2") +
                        coord_cartesian(xlim = c(-0.3, 0.75))
    
plt_pseudobulk_cor <- add_theme(plt_pseudobulk_cor)

ggsave(filename = "./Figures/figure5/plt_pseudobulk_cor.pdf", 
       plot = add_theme(plt_pseudobulk_cor), units = "cm", width = 10, height = 4.5)  
