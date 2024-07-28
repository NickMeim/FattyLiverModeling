library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(ggforce)
library(caret)
library(ropls)
library(factoextra)
library(Rtsne)
# library(GA)
source("../utils/plotting_functions.R")
source("functions_translation.R")
source("CrossValidationUtilFunctions.R")

### Load in-vitro and in-vivo datasets-----------------------------------------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Manually load also the other clinical datasets I have from ARCHS4
geo <- 'GSE162694' # only this one has multiple NAS and fibrosis scores
meta_data <- read.delim('../data/ARCHS4/FattyLiver_meta_data.tsv',row.names = 1)
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
expression_matrix <- readRDS('../data/ARCHS4/FattyLiver_expression_matrix.rds')
expression_matrix <- expression_matrix[,meta_data$sample]
data_list[['Pantano']] <- list(counts = expression_matrix,
                               metadata = meta_data,
                               genes = rownames(expression_matrix))


# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list

## X_all and PCA and t-SNE visualize
clinical <- c('Govaere','Hoang','Pantano')
meata_data_all <- data.frame()
for (i in 1:length(clinical)){
  set <- clinical[i]
  if (i==1){
    X_all <- data_list[[set]]$data_center %>% t() %>% scale()
  }else{
    X_all <- rbind(X_all,data_list[[set]]$data_center %>% t() %>% scale())
  }
  if (set=='Govaere'){
    meata_data_all <- rbind(meata_data_all,
                            data_list[[set]]$metadata  %>% 
                              select(c('NAS'='nas_score'),c('fibrosis'='Fibrosis_stage')) %>%
                              mutate(dataset = set))
  }else if (set=='Hoang'){
    meata_data_all <- rbind(meata_data_all,
                            data_list[[set]]$metadata  %>% 
                              select(c('NAS'='nafld_activity_score'),c('fibrosis'='Fibrosis_stage')) %>%
                              mutate(dataset = set))
  }else{
    meata_data_all <- rbind(meata_data_all,
                            data_list[[set]]$metadata  %>% 
                              select(c('NAS'='nas score'),c('fibrosis'='fibrosis stage')) %>%
                              mutate(dataset = set))
  }
}
inds <- apply(X_all,2,function(x)(return(sum(is.na(x)))))
inds <- which(inds==0)
X_all <- X_all[,inds]
## visualize to see batch effects
pca_all_humans <- prcomp(X_all,center = TRUE,scale. = FALSE)
fviz_screeplot(pca_all_humans)
df_pca <- as.data.frame(pca_all_humans$x[,1:10])
df_pca <- cbind(df_pca,meata_data_all %>% mutate(NAS=as.numeric(NAS)) %>% mutate(fibrosis=as.numeric(fibrosis)))
ggplot(df_pca,aes(x=PC1,y=PC2,fill=dataset))+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  # scale_fill_viridis_c()+
  theme_minimal(base_size = 20,base_family = 'Arial')
ggplot(df_pca,aes(x=PC1,y=PC3,fill=dataset))+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  theme_minimal(base_size = 20,base_family = 'Arial')
ggplot(df_pca,aes(x=PC2,y=PC3,fill=dataset))+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  theme_minimal(base_size = 20,base_family = 'Arial')
ggplot(df_pca,aes(x=PC4,y=PC7,fill=dataset))+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  theme_minimal(base_size = 20,base_family = 'Arial')
anov <- aov(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)~NAS+fibrosis+dataset,data=df_pca %>% 
              mutate(dataset=as.numeric(as.factor(dataset))))
stats_results <- df_pca %>%  mutate(dataset=as.numeric(as.factor(dataset))) %>%
  gather('PC','value',-NAS,-fibrosis,-dataset) %>% group_by(PC) %>%
  rstatix::anova_test(value~NAS+fibrosis+dataset) %>% 
  rstatix::adjust_pvalue(method = 'BH') %>%
  ungroup()
coeffs <- coef(anov)
coeffs <- as.data.frame(coeffs)[2:nrow(coeffs),] %>% rownames_to_column('Effect') %>%
  gather('PC','coefficient',-Effect)
stats_results <- as.data.frame(stats_results)
stats_results <- stats_results %>% 
  mutate(significance = ifelse(p.adj<0.001,'***',
                               ifelse(p.adj<0.01,'**',
                                      ifelse(p.adj<0.05,'*',
                                             ifelse(p.adj<0.1,'.','')))))
stats_results <- left_join(stats_results,coeffs)
p_F <- ggplot(stats_results,aes(x=PC,y=Effect,fill=`F`)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient(low = 'white',high='red')+
  labs(fill = "F-statistic") + 
  ggtitle('ANOVA results')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
p_coeff <- ggplot(stats_results,aes(x=PC,y=Effect,fill=coefficient)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient2(low = 'blue',high='red',mid = 'white',midpoint = 0)+
  ggtitle('Linear model coefficients')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
p <- p_F + p_coeff
print(p)

### run t-SNE
perpl = DescTools::RoundTo(sqrt(nrow(X_all)), multiple = 5, FUN = round)
emb_size = ncol(X_all)
set.seed(42)
tsne_all_humans <- Rtsne(X_all, 
                  dims = 2, perplexity=perpl, 
                  verbose=TRUE, max_iter = 1000,initial_dims = 30,check_duplicates = T,
                  pca_scale = FALSE,
                  num_threads = 15)
df_tsne <- as.data.frame(tsne_all_humans$Y[,1:2])
df_tsne <- cbind(df_tsne,meata_data_all %>% mutate(NAS=as.numeric(NAS)) %>% mutate(fibrosis=as.numeric(fibrosis)))
ggplot(df_tsne,aes(x=V1,y=V2,fill=dataset))+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  theme_minimal(base_size = 20,base_family = 'Arial')+
  theme(legend.position = 'top')
(ggplot(df_tsne,aes(x=V1,y=V2,fill=dataset))+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  theme_minimal(base_size = 20,base_family = 'Arial')+
    theme(legend.position = 'top')) +
  (ggplot(df_tsne,aes(x=V1,y=V2,fill=NAS))+
     scale_fill_viridis_c()+
     geom_point(size=2.8,shape=21,stroke=1.2)+
     theme_minimal(base_size = 20,base_family = 'Arial')+
     theme(legend.position = 'top'))+
  ggplot(df_tsne,aes(x=V1,y=V2,fill=fibrosis))+
  scale_fill_viridis_c()+
  geom_point(size=2.8,shape=21,stroke=1.2)+
  theme_minimal(base_size = 20,base_family = 'Arial')+
  theme(legend.position = 'top')

# Define matrices of interest
# Yh <- as.matrix(meata_data_all %>% select(NAS,fibrosis)) #keep both Fibrosis and NAS
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% select(nas_score,Fibrosis_stage))
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
# Xh <- X_all
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()
mps_meta_data <- data_list[[target_dataset]]$metadata
rownames(mps_meta_data) <- NULL
mps_meta_data <- mps_meta_data %>% select(sampleName,treatment,condition) %>% unique()
print(all(rownames(Xm)==mps_meta_data$sampleName))
# saveRDS(mps_meta_data,paste0('../preprocessing/',tolower(target_dataset),'_meta_data.rds'))
# saveRDS(Xm,paste0('../preprocessing/',tolower(target_dataset),'_centered_data.rds'))
# saveRDS(Xh,paste0('../preprocessing/',tolower(ref_dataset),'_centered_data.rds'))
# saveRDS(Yh,paste0('../preprocessing/',tolower(ref_dataset),'_phenotype_data.rds'))


#### Train PLSR model and score MPS samples----------------------------------------------------
plsr_model <- opls(x = Xh, 
                   y = Yh,
                   predI = 8,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
total_plsr_var <- sum(apply(plsr_model@scoreMN,2,var))/sum(apply(Xh,2,var))
print(paste0('PLSR with 8 LVs accounts for ',round(100*total_plsr_var,2),'% of total variance'))
# Get Wh of PLSR
Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
rownames(Wh) <- colnames(Xh)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)){
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
}
# Get regression coefficients
Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN

Yhat_MPS <- predict(plsr_model,Xm)
df_mps <- rbind(as.data.frame(Yhat_MPS) %>% rownames_to_column('sample') %>% 
                  gather('phenotype','score',-sample) %>% mutate(system = 'MPS'),
                    as.data.frame(Yh) %>% rownames_to_column('sample') %>% 
                  gather('phenotype','score',-sample) %>% mutate(system = 'human')) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))

custom_breaks <- list('Fibrosis stage'=seq(0,4,1),
                      NAS=seq(0,8,1))
# Create a custom scale function to set breaks based on phenotype
custom_scale_x <- function(phenotype) {
  scale_x_continuous(breaks = custom_breaks[[phenotype]], limits = c(0, NA))
}
p <- ggplot(df_mps, aes(x = score, fill = system)) +
  geom_histogram(color = 'black', alpha = 0.4, bins = 20, position = 'identity') +
  facet_wrap_paginate(~phenotype, scales = 'free', nrow = 2, ncol = 1, page = 1) +
  theme_pubr(base_family = 'Arial', base_size = 20) +
  theme(text = element_text(size = 20, family = 'Arial'))
# Print the plot for each phenotype with custom x-axis breaks
for (i in seq_along(custom_breaks)) {
  print(p + custom_scale_x(names(custom_breaks)[i]))
}
ggsave('../figures/plsr_scored_mps_samples.eps',
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)


## Move sample in LV extra 1---------------------------------
#### edo rikse to Z new kai to dX kai des gia to dX ti pathways vgainoun kai gia to Z pou exei paei sto latent space
perturbations <- data.table::fread('../results/FindPerturbation4TargetNAS_Kostrzewski_for_govaere.csv') %>% select(-V1)
pert_info <- perturbations %>% select(treatment,NAS,mean_distance) %>% mutate(id = paste0(treatment,':',NAS)) %>%
  column_to_rownames('id')
perturbations <- perturbations %>% mutate(id = paste0(treatment,':',NAS)) %>% select(-treatment,-NAS,-mean_distance) %>%
  column_to_rownames('id')
perturbations <- t(perturbations)
gc()

### Infer pathway activity
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
path_activity  <- decoupleR::run_viper(perturbations, net_prog,minsize = 1,verbose = TRUE) %>% select(-statistic)
path_activity_mat <- as.matrix(path_activity %>% select(source,condition,score) %>% spread('source','score') %>% 
  column_to_rownames('condition'))
dend1 <- hclust(dist(path_activity_mat))
dend2 <- hclust(dist(t(path_activity_mat)))
ordered_matrix <- path_activity_mat[dend1$order, dend2$order]
cond_groups <- ifelse(grepl('Fat',rownames(ordered_matrix)),'Fat','Lean')
nas_groups <- data.frame(group = rownames(ordered_matrix))
nas_groups <- nas_groups %>% mutate(NAS = str_split_fixed(group,':',2)[,2])
nas_groups <- nas_groups %>%mutate(NAS=as.numeric(NAS))
nas_groups <- nas_groups$NAS
groups_row <- data.frame(Condition = cond_groups,NAS=nas_groups)
rownames(groups_row) <- rownames(ordered_matrix)
# png('../results/similarity_of_optimal_analytical_lvs_for_many_partitions.png',
#     height = 14,width = 16,units = 'in',res=600)
pheatmap::pheatmap(ordered_matrix, 
                   color = colorRampPalette(c("blue", "white", "red"))(100), 
                   breaks = seq(-4, 4, length.out = 100),
                   annotation_row = groups_row,
                   fontsize = 18,
                   show_colnames = TRUE,  
                   show_rownames = FALSE,
                   filename = '../figures/perturbations_clustermap.png',
                   height = 12,
                   width = 14,
                   units = 'in',
                   res=600
)


## Project onto exrta basis space
Wm_opt <- readRDS('../results/Wm_kostrzewski_extra.rds')
Zpert <- data.table::fread('../results/Z_pertNAS_Kostrzewski_for_govaere.csv') 
Zh <- as.data.frame(Xh %*% Wm_opt)
Zh$NAS <- Yh[,'NAS']
colnames(Zh)[1:2] <- c('Z0','Z1')
Z_all <- rbind(Zh %>% mutate(type='human'),Zpert %>% select(Z0,Z1,NAS) %>% mutate(type='MPS perturbed'))
# (ggplot(Z_all,aes(x=Z0,y=Z1,fill=type))+
#   geom_point(size=2.8,shape=21,stroke=1.2)+
#   theme_minimal(base_size = 20,base_family = 'Arial') )/ 
#   (ggplot(Z_all %>% filter(type=='MPS perturbed'),aes(x=Z0,y=Z1,fill=NAS))+
#      geom_point(size=2.8,shape=21,stroke=1.2)+
#      scale_fill_viridis_c()+
#      theme_minimal(base_size = 20,base_family = 'Arial'))/
#   (ggplot(Zh,aes(x=Z0,y=Z1,fill=NAS))+
#   geom_point(size=2.8,shape=21,stroke=1.2)+
#   scale_fill_viridis_c()+
#   theme_minimal(base_size = 20,base_family = 'Arial'))

Z_all$type <- factor(Z_all$type,levels = c('MPS perturbed','human')) 
ggplot(Z_all,aes(x=Z0,y=Z1,fill=NAS,size=type))+
  geom_point(shape=21,stroke=1.2)+
  scale_fill_viridis_c()+
  scale_size_manual(values = c(1.5,3))+
  theme_minimal(base_size = 20,base_family = 'Arial')
  
