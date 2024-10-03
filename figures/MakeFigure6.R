### Functions for making elements of figure 6 - General optimization of variance
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets
library(tidyverse)
library(matrixStats)
library(ggExtra)
root_dir <- "C:/Users/nmeim/Documents/LiverModeling/FattyLiverModeling"
setwd(root_dir)
source("utils/plotting_functions.R")
source("modeling/functions_translation.R")
source("modeling/vector_space_interpretation.R")
target_dataset <- "Kostrzewski"
ref_dataset <-  "Govaere"
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
processed_data_list <- readRDS(paste0('results/processed_data_list_',
                                      tolower(ref_dataset),'_',
                                      tolower(target_dataset),'.rds'))
Xh <- processed_data_list$Xh
Yh <- processed_data_list$Yh
sex_inferred <- processed_data_list$sex_inferred
Xm <- processed_data_list$Xm
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
colnames(Wm_opt) <- c("LV extra 1", "LV extra 2")
#################################################################################
### Panel - Human data on LV optimized through perturbation
# Load perturbation
Wm <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
dx_lean <- data.table::fread('./results/optimized_mps/dx_lean_govaere_kostrzewski_corrected.csv') %>%
            dplyr::select(-V1) %>%
            t()

# Plot pathway activities
dx_pwy_activity <- pathway_activity_interpretation(dx_lean, Wm)

plt_pwy_dx_lean <- plot_pwy_activity(dx_pwy_activity %>% filter(condition == "V1"),
                                 plt_lim = 8.5, show_fill_legend = T, offset_annot = 1.2, n.breaks = 6)

plt_pwy_dx_lean <- add_theme(plt_pwy_dx_lean)

# Plot hallmarks
df_msig <- readRDS(paste0("results/hallmark_enrichment_", tolower(target_dataset),"_dx_lean.rds"))


plt_hallmarks_dx_lean <- ggplot(df_msig %>% arrange(NES) %>%
                              filter(padj<=0.1), aes(x=NES,y=reorder(Hallmark,NES),fill=NES))+
  geom_bar(stat = 'identity', size = size_col, color = "black", show.legend = F) +
  # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
  scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0) +
  labs(x = 'Normalized Enrichment Score', y = 'Hallmark')

plt_hallmarks_dx_lean <- add_theme(plt_hallmarks_dx_lean)

#################################################################################
### Panel - Human PLSR labelled by (inferred) sex
sex_inferred <- readRDS(paste0("results/", tolower(ref_dataset),"_sex_inferred.rds"))
# plsr_model <- readRDS(paste0('results/PLSR_model_',tolower(ref_dataset),'.rds'))
# Zh_plsr <- plsr_model@scoreMN
Z <- Xh %*% Wm_opt
lvs <- c(1,2)
lv_vars <- round(100 * colVars(Z)/sum(colVars(Xh)),2)


plt_PLSR_sex <- ggplot(cbind(as.data.frame(Z[,c(lvs[1],lvs[2])]),
                             as.data.frame(sex_inferred) %>% dplyr::select(c('sex'='V1'))),
                       aes(x=`LV extra 1`,y=`LV extra 2`,color=sex)) +
                geom_point(size = size_dot)+
                scale_color_viridis_d()+
                stat_ellipse(size = size_line, show.legend = F, linetype = 2) +
                labs(color = 'Sex')
plt_PLSR_sex <- add_theme(plt_PLSR_sex)

### Save panels as figures
ggsave(filename = "./Figures/figure6/plt_pwy_dx_lean.pdf", plot = plt_pwy_dx_lean, units = "cm", width = 5.5, height = 5)
ggsave(filename = "./Figures/figure6/plt_hallmarks_dx_lean.pdf", plot = plt_hallmarks_dx_lean, units = "cm", width = 7, height = 5)
ggsave(filename = "./Figures/figure6/plt_PLSR_sex.pdf", plot = plt_PLSR_sex, units = "cm", width = 6, height = 5)



