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

### Load all data------------------------------------------------------------------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = 'data/')
# Manually load also the other clinical datasets I have from ARCHS4
geo <- 'GSE162694' # only this one has multiple NAS and fibrosis scores
meta_data <- read.delim('data/ARCHS4/FattyLiver_meta_data.tsv',row.names = 1)
meta_data <- meta_data %>% filter(series_id==geo)
old_cols <- colnames(meta_data)
meta_data <- meta_data %>% separate_rows(characteristics_ch1, sep = ",") %>%
  separate(characteristics_ch1, into = c("key", "value"), sep = ":") %>%
  mutate(key = str_trim(key), value = str_trim(value)) %>%
  spread(key, value)
new_cols <- colnames(meta_data)
new_cols <- new_cols[which(!(new_cols %in% old_cols))]
meta_data <- meta_data %>%
  mutate_at(vars(new_cols), ~ ifelse(grepl("\\d", .), as.numeric(.), .))
meta_data <- meta_data %>% filter(!(is.na(`nas score`) & is.na(`fibrosis stage`))) %>%
  mutate(`nas score`=ifelse(`nas score`=='NA',NA,`nas score`)) %>% 
  mutate(`fibrosis stage`=ifelse(`fibrosis stage`=='normal liver histology',0,`fibrosis stage`))
meta_data <- meta_data %>% filter(!is.na(`nas score`)) %>% filter(!is.na(`fibrosis stage`)) %>%
  filter(`nas score`!='NA') %>% filter(`fibrosis stage`!='NA')
expression_matrix <- readRDS('data/ARCHS4/FattyLiver_expression_matrix.rds')
expression_matrix <- expression_matrix[,meta_data$sample]
data_list[['Pantano']] <- list(counts = expression_matrix,
                               metadata = meta_data,
                               genes = rownames(expression_matrix))
### Pre-process all
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list

### Iterate through all possible combinations of in-vitro / in-vivo datasets-------------------
in_vivos <-  c("Govaere","Hoang","Pantano")
in_vitros <- c("Kostrzewski","Wang", "Feaver")
all_results <- data.frame()
for (dataset in in_vivos){
  message(paste0('Begun ',dataset ,' dataset'))
  Xh <- data_list[[dataset]]$data_center %>% t()
  if (dataset == 'Hoang'){
    Yh <- as.matrix(data_list[[dataset]]$metadata  %>% select(nafld_activity_score,Fibrosis_stage))
    
  }else if (dataset == 'Pantano'){
    Yh <- as.matrix(data_list[[dataset]]$metadata  %>% select(`nas score`,`fibrosis stage`))
    Yh <- apply(Yh,c(1,2),as.numeric)
  }else{
    Yh <- as.matrix(data_list[[dataset]]$metadata  %>% select(nas_score,Fibrosis_stage))
  }
  colnames(Yh) <- c('NAS','fibrosis')
  plsr_model <- opls(x = Xh, 
                     y = Yh,
                     predI = 8,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
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
  for (experimental_model in in_vitros){
    Wm <- data_list[[experimental_model]]$Wm_group %>% as.matrix()
    Wm_opt <- analytical_solution_opt(y=Yh,
                                      W_invitro = Wm,
                                      phi = phi)
    colnames(Wm_opt) <- c("LV_extra1", "LV_extra2")
    rownames(Wm_opt) <- rownames(Wm)
    
    translatable_components_progenies <- pathway_activity_interpretation(Wm_opt, Wm)
    
    if (experimental_model=='Kostrzewski' & dataset=='Govaere'){
      # get figure of pathway activities
      plt_pwy_LV1 <- plot_pwy_activity(translatable_components_progenies %>% filter(condition == "LV_extra1"),
                                       plt_lim = 9, show_fill_legend = T, offset_annot = 0.8, n.breaks = 6)
      
      plt_pwy_LV1 <- add_theme(plt_pwy_LV1)
      pathways_order_1 <- plt_pwy_LV1$data$Pathway[order(plt_pwy_LV1$data$score)]
      
      plt_pwy_LV2 <- plot_pwy_activity(translatable_components_progenies %>% filter(condition == "LV_extra2"),
                                       plt_lim = 9, show_fill_legend = T, offset_annot = 0.8, n.breaks = 6)
      
      plt_pwy_LV2 <- add_theme(plt_pwy_LV2)
      pathways_order_2 <- plt_pwy_LV2$data$Pathway[order(plt_pwy_LV2$data$score)]
    }else{
      all_results <- rbind(all_results,
                           translatable_components_progenies %>% 
                             mutate(clinical=dataset) %>%
                             mutate(invitro=experimental_model))
    }
    print(paste0('Finished in vitro model: ',experimental_model))
  }
}
saveRDS(all_results,'results/all_invitro_invivo_combos_extraLVs_pathways.rds')

##### Visualize---------------------------
## plot LV1
df_LV1 <- all_results %>% filter(condition=='LV_extra1')
df_LV1$Pathway <- factor(df_LV1$Pathway,levels = pathways_order_1)
# Set plot limits
plt_lim = max(abs(df_LV1$score)) + 1.5

# Annotate pvals with asterisks
df_LV1 <- df_LV1 %>%
  mutate(annot = ifelse(p_value <= 0.0001, "****",
                        ifelse(p_value <= 0.001,"***",
                               ifelse(p_value<=0.01,"**",
                                      ifelse(p_value<=0.05,'*',
                                             ifelse(p_value<=0.1,'\u2022', #\u2219
                                                    'ns')))))) %>%
  mutate(offset = ifelse(annot == 'ns', 2, 0.5*2))
# Make plot
plt_lv1 <- ggplot(df_LV1,aes(x = score,y = Pathway,fill = score, label = annot)) + 
  geom_bar(stat = 'identity', show.legend = FALSE, color = "black", size = size_col) +
  scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                       midpoint = 0,limits = c(-plt_lim,plt_lim))+
  scale_x_continuous(n.breaks = 8,limits = c(-plt_lim,plt_lim))+
  labs(y = "Pathway", x = "Activity score", fill = "Score") +
  geom_text(aes(x = ifelse(score < 0, score - offset, score + offset)),
            color = 'black', size = size_annotation*0.7,
            angle = 90, show.legend = FALSE, hjust = "center", vjust = "center")  +
  facet_wrap(vars(clinical,invitro))
plt_lv1 <- add_theme(plt_lv1)
print(plt_lv1)
ggsave(filename = "figures/supplementary_all_invitro_invivo_combos_LV1_pathways.pdf",
       plot = plt_lv1, units = "cm", width = 12, height = 12)

## plot LV2
df_LV2 <- all_results %>% filter(condition=='LV_extra2')
df_LV2$Pathway <- factor(df_LV2$Pathway,levels = pathways_order_2)
# Set plot limits
plt_lim = max(abs(df_LV2$score)) + 1.5
# Annotate pvals with asterisks
df_LV2 <- df_LV2 %>%
  mutate(annot = ifelse(p_value <= 0.0001, "****",
                        ifelse(p_value <= 0.001,"***",
                               ifelse(p_value<=0.01,"**",
                                      ifelse(p_value<=0.05,'*',
                                             ifelse(p_value<=0.1,'\u2022', #\u2219
                                                    'ns')))))) %>%
  mutate(offset = ifelse(annot == 'ns', 2, 0.5*2))
# Make plot
plt_lv2 <- ggplot(df_LV2,aes(x = score,y = Pathway,fill = score, label = annot)) + 
  geom_bar(stat = 'identity', show.legend = FALSE, color = "black", size = size_col) +
  scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                       midpoint = 0,limits = c(-plt_lim,plt_lim))+
  scale_x_continuous(n.breaks = 8,limits = c(-plt_lim,plt_lim))+
  labs(y = "Pathway", x = "Activity score", fill = "Score") +
  geom_text(aes(x = ifelse(score < 0, score - offset, score + offset)),
            color = 'black', size = size_annotation*0.7,
            angle = 90, show.legend = FALSE, hjust = "center", vjust = "center")  +
  facet_wrap(vars(clinical,invitro))
plt_lv2 <- add_theme(plt_lv2)
print(plt_lv2)
ggsave(filename = "figures/supplementary_all_invitro_invivo_combos_LV2_pathways.pdf",
       plot = plt_lv2, units = "cm", width = 12, height = 12)
