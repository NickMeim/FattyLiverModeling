library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(caret)
library(ropls)
source("../utils/plotting_functions.R")
source("../modeling/functions_translation.R")

### Load data---------------------------------------------------------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')

### Run PCA----------------------------------------------------------
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Govaere$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# Get grouped data and center them 
Xm_grouped <- data_list[[target_dataset]]$Xm_grouped %>% as.matrix()
Xm_grouped <- Xm_grouped  - rowMeans(Xm_grouped)
Xm_grouped <- t(Xm_grouped)
Xm_grouped <- as.data.frame(Xm_grouped)
data.table::fwrite(Xm_grouped,paste0('X_',target_dataset,'_grouped.csv'),row.names = TRUE)
data.table::fwrite(as.data.frame(Xh),paste0('X_',ref_dataset,'.csv'),row.names = TRUE)
data.table::fwrite(as.data.frame(Yh),paste0('Y_',ref_dataset,'.csv'),row.names = TRUE)

### Calculate percentage variance of each pair of conditions-------

# First calculate total human variance
var_human <- sum(apply(Xh, 2, var))

combinations <- combn(rownames(Xm_grouped), 2)
# combinations_list <- apply(combinations, 2, function(pair) list(pair))
perc_pair_var_explained <- NULL 
pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
for (j in 1:ncol(combinations)){
  X <- Xm_grouped[combinations[,j],]
  pca_selected <- prcomp(X,scale. = FALSE,center = TRUE)
  X_proj <- Xh %*% pca_selected$rotation
  var_proj <- sum(apply(X_proj, 2, var))
  perc_pair_var_explained[j] <- var_proj/var_human
  # Update the progress bar
  setTxtProgressBar(pb, j)
}
df_res <- as.data.frame(t(combinations))
df_res$pair_id <- seq(1,nrow(df_res))
df_res$var_explained <- perc_pair_var_explained
df_res <- df_res %>% arrange(-var_explained)
## calculate cumulative
cumulative_var_explained <- NULL
for (j in 2:nrow(df_res)){
  X <- Xm_grouped[unique(c(as.matrix(df_res[1:j,1:2]))),]
  pca_selected <- prcomp(X,scale. = FALSE,center = TRUE)
  X_proj <- Xh %*% pca_selected$rotation
  var_proj <- sum(apply(X_proj, 2, var))
  cumulative_var_explained[j] <- var_proj/var_human
  # Update the progress bar
  setTxtProgressBar(pb, j)
}
cumulative_var_explained[1] <- df_res$var_explained[1]
df_res$cumulative_var_explained <- cumulative_var_explained
df_res <- df_res %>% arrange(cumulative_var_explained)
df_res$pair_id <- factor(df_res$pair_id,levels = df_res$pair_id)
pathviewr::find_curve_elbow(df_res %>% select(pair_id,cumulative_var_explained) %>% mutate(pair_id = as.numeric(pair_id)),plot_curve = TRUE)
ggplot(df_res,aes(x=pair_id,y=100 * cumulative_var_explained, group = 1)) + 
  geom_point(color = 'blue' , size=1.5) + geom_line(color='blue',lwd=1)+
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  geom_hline(yintercept = 100*df_res$cumulative_var_explained[42])+
  geom_vline(xintercept = 42)+
  xlab('') + ylab('cumulative explained variance (%)')+
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_blank(),
        axis.ticks.x =  element_blank())
ggsave('cumulative_human_variance_explained_by_invitro.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)  

kept_conditions <- df_res[1:42,]
df_res1 <- kept_conditions %>% separate(V1,into = data_list[[target_dataset]]$exp_factors,sep = '_',remove = FALSE) %>%
  gather('stimuli','value',-pair_id,-V1,-V2,-var_explained,-cumulative_var_explained) %>% mutate(group = 'V1')
df_res2 <- kept_conditions %>% separate(V2,into = data_list[[target_dataset]]$exp_factors,sep = '_',remove = FALSE) %>%
  gather('stimuli','value',-pair_id,-V1,-V2,-var_explained,-cumulative_var_explained) %>% mutate(group = 'V2')
df_res_gathered <- rbind(df_res1,df_res2)
df_res_gathered <- df_res_gathered %>% arrange(-cumulative_var_explained)
df_res_gathered <- df_res_gathered %>% mutate(stimuli = ifelse(stimuli=='NPC',paste0('NPC ',value),
                                                               ifelse(stimuli=='Background',value,
                                                                      stimuli))) %>%
  mutate(value = ifelse(value %in% c('TRUE','FALSE'),value,'TRUE'))
df_res_gathered <- df_res_gathered %>% filter(value == 'TRUE') %>% 
  select(pair_id,stimuli,var_explained,cumulative_var_explained) %>% unique()
ggplot(df_res_gathered %>% select(pair_id,stimuli) %>% unique() %>% 
         group_by(stimuli) %>%
         count() %>%
         mutate(percentage = n / length(unique(df_res_gathered$pair_id)) * 100),
       aes(x=reorder(stimuli,-percentage),y=percentage)) + geom_bar(stat = "identity", position = "dodge") +
  ylab('counts (%)') + xlab('stimuli') +
  theme_pubr(base_size = 15,base_family = 'Arial')+
  theme(panel.grid.major = element_line(),
        axis.line = element_blank())
ggsave('percentage_invitro_stimuli_in_high_explainable_samples.png',
       width = 9,
       height = 6,
       units = 'in',
       dpi=600) 
## see only controls
df_res_controls <- rbind(df_res %>% filter(pair_id %in% df_res_gathered$pair_id) %>%
  mutate(control = ifelse(grepl('FALSE_FALSE_FALSE_FALSE',V1),V1,
                          ifelse(grepl('FALSE_FALSE_FALSE_FALSE',V2),V2,NA))) %>%
    filter(!is.na(control)) %>% select(pair_id,cumulative_var_explained,var_explained,control) %>% unique(),
  df_res %>% filter(pair_id %in% df_res_gathered$pair_id) %>%
    mutate(control = ifelse(grepl('FALSE_FALSE_FALSE_FALSE',V2),V2,
                            ifelse(grepl('FALSE_FALSE_FALSE_FALSE',V1),V1,NA))) %>%
    filter(!is.na(control)) %>% select(pair_id,cumulative_var_explained,var_explained,control) %>% unique()) %>%
  unique()
df_res_controls <- df_res_controls %>% mutate(control = sub('_FALSE_FALSE_FALSE_FALSE','',control))
ggplot(df_res_controls %>% select(pair_id,control) %>% unique() %>% 
         group_by(control) %>%
         count() %>%
         mutate(percentage = n / length(unique(df_res_gathered$pair_id)) * 100),
       aes(x=reorder(control,-percentage),y=percentage)) + geom_bar(stat = "identity", position = "dodge") +
  ylab('counts (%)') + xlab('control condition') +
  theme_pubr(base_size = 15,base_family = 'Arial')+
  theme(panel.grid.major = element_line(),
        axis.line = element_blank())
ggsave('controls_percentage_invitro_stimuli_in_high_explainable_samples.png',
       width = 9,
       height = 6,
       units = 'in',
       dpi=600) 
