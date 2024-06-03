library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(ggradar)
library(ggiraphExtra)
source('modeling/CrossValidationUtilFunctions.R')
source('modeling/functions_translation.R')
source("utils/plotting_functions.R")
source("modeling/vector_space_interpretation.R")

### Load the all the data to be used-----------------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = 'data/')
# Run PCA
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





### Check current TF and pathway activity in the data------------------------------
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)

path_acitivity <- decoupleR::run_viper(t(Xm), net_prog,minsize = 1,verbose = FALSE)
path_acitivity <- path_acitivity %>% select(c('Pathway'='source'),condition,score,p_value)
path_acitivity <- left_join(path_acitivity,
                            data_list$Kostrzewski$metadata %>% select(sampleName,treatment,TGF,LPS,Cholesterol,Fructose),
                            by=c('condition'='sampleName'))
path_acitivity <- left_join(path_acitivity,
                            as.data.frame(Xm %*% Wm) %>% select(PC1,PC2) %>% rownames_to_column('condition'))
# path_acitivity <- path_acitivity %>% filter(p_value<0.01)
ggscatter(path_acitivity %>% mutate(p=ifelse(p_value<0.0001,'<0.0001',
                                                     ifelse(p_value<0.001,'<0.001',
                                                           ifelse(p_value<0.01,'<0.01',
                                                                  ifelse(p_value<0.05,'<0.05',
                                                                        'ns'))))) %>%
            mutate(p = factor(p,levels = c('ns','<0.05','<0.01','<0.001','<0.0001'))) ,
          x='PC1','PC2',fill = 'score',size ='p',shape=21) +
  scale_fill_gradient2(low='blue',high='green',midpoint = 0,mid='white') +
  theme(text = element_text(size=20,family='Arial'))+
  facet_wrap(~Pathway)
ggsave(paste0('results/',tolower(target_dataset),'_pathway_activity_in_pc_space.png'),
       height = 12,
       width = 16,
       units = 'in',
       dpi=600)

ggscatter(left_join(path_acitivity,
                    as.data.frame(Xm %*% Wm) %>% rownames_to_column('condition')) %>%
            mutate(p=ifelse(p_value<0.0001,'<0.0001',
                                              ifelse(p_value<0.001,'<0.001',
                                                     ifelse(p_value<0.01,'<0.01',
                                                            ifelse(p_value<0.05,'<0.05',
                                                                   'ns'))))) %>%
            mutate(p = factor(p,levels = c('ns','<0.05','<0.01','<0.001','<0.0001'))) %>%
            gather('PC','value',-p,-condition,-score,-treatment,-TGF,-LPS,-Cholesterol,-Fructose,-Pathway,-p_value) %>%
            filter(PC %in% c('PC1','PC2','PC3','PC4','PC5')) %>% filter(Pathway=='JAK-STAT'),
          x='value',y='score',cor.coef = TRUE) +
  # scale_fill_gradient2(low='blue',high='green',midpoint = 0,mid='white') +
  theme(text = element_text(size=20,family='Arial'))+
  facet_wrap(~PC)
ggscatter(left_join(path_acitivity,
                    as.data.frame(Xm %*% Wm) %>% rownames_to_column('condition')) %>%
            mutate(p=ifelse(p_value<0.0001,'<0.0001',
                                              ifelse(p_value<0.001,'<0.001',
                                                     ifelse(p_value<0.01,'<0.01',
                                                            ifelse(p_value<0.05,'<0.05',
                                                                   'ns'))))) %>%
            mutate(p = factor(p,levels = c('ns','<0.05','<0.01','<0.001','<0.0001'))) %>%
            gather('PC','value',-p,-condition,-score,-treatment,-TGF,-LPS,-Cholesterol,-Fructose,-Pathway,-p_value) %>%
            filter(PC %in% c('PC1','PC2','PC3','PC4','PC5')) %>% filter(Pathway=='TGFb'),
          x='value',y='score',cor.coef = TRUE) +
  # scale_fill_gradient2(low='blue',high='green',midpoint = 0,mid='white') +
  theme(text = element_text(size=20,family='Arial'))+
  facet_wrap(~PC)


### Run PLSR and find extra basis--------------------------
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

# Find extra basis
phi <- Wh %*% Bh
Wm_opt <- analytical_solution_opt(y=Yh,
                                  W_invitro = Wm,
                                  phi = phi)
Wm_tot <- cbind(Wm, Wm_opt)
# predict and evaluate briefly
y_hat <- cbind(1, Xh %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
plot(Yh[,1],y_hat[,1],xlab='True NAS',ylab='Predicted NAS',main=paste0('r = ',cor(Yh[,1],y_hat[,1])))
plot(Yh[,2],y_hat[,2],xlab='True fibrosis',ylab='Predicted fibrosis',main=paste0('r = ',cor(Yh[,2],y_hat[,2])))

### Plot the directions of correlation with NAS ans Fibrosis in this space
Zh <- as.data.frame(Xh %*% Wm_opt)
Zh$NAS <- Yh[,1]
Zh$fibrosis <- Yh[,2]
Zh <- Zh %>% rownames_to_column('sample')
thetas <- seq(0,360,5)
df_proj <- data.frame()
df <- data.frame()
for (theta in thetas){
  u <- c(cos(theta * pi / 180),sin(theta * pi / 180))
  u <- as.matrix(u)
  proj <- u %*% t(u) # since it is unit norm vector
  tmp <- as.matrix(Zh%>% select(V1,V2))
  Z_proj <- tmp %*% proj
  if (theta %in% c(90,270)){
    corr_nas <- cor(Z_proj[,2],Zh$NAS)
    corr_fib <- cor(Z_proj[,2],Zh$fibrosis)
  }else{
    corr_nas <- cor(Z_proj[,1],Zh$NAS)
    corr_fib <- cor(Z_proj[,1],Zh$fibrosis)
  }
  df <- rbind(df,
              data.frame(theta = theta,phenotype = 'NAS',corr=corr_nas),
              data.frame(theta = theta,phenotype = 'fibrosis',corr=corr_fib))
  df_proj <- rbind(df_proj,
                   data.frame(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1],phenotype = 'NAS',corr=corr_nas),
                   data.frame(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1],phenotype = 'fibrosis',corr=corr_fib))
}
p1 <- ggplot(df %>% spread('phenotype','corr') %>% group_by(theta) %>% 
               mutate(`absolute average`=0.5*(abs(NAS)+abs(fibrosis))) %>% 
               ungroup() %>% gather('phenotype','corr',-theta),
             aes(x=theta,y=corr,colour=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = seq(-1,1,0.25))+
  scale_x_continuous(breaks = seq(0,360,30))+
  xlab(expression(theta*" (\u00B0)"))+
  theme_minimal(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family='Arial'),
        legend.position = 'top')

p2 <- ggplot(df_proj,aes(x=LV_opt_1,y=LV_opt_2,fill=corr,colour=corr))+
  geom_segment(aes(x=0,y=0,xend =LV_opt_1,yend=LV_opt_2),
               arrow = arrow(length = unit(0.03, "npc")))+
  scale_fill_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-1,1))+
  scale_color_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-1,1))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~phenotype)+
  theme_minimal(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family='Arial'),
        legend.position = 'right')
p1 / p2
ggsave(paste0('results/pc_loadings_scores_analysis/',
              tolower(ref_dataset),
              'optimal_space_phenotypes_correlations_directions',
              tolower(target_dataset),
              '.png'),
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)

(ggplot(Zh ,aes(x=V1,y=V2,fill=NAS))+
    geom_point(size=3,shape=21,stroke=1.3)+
    scale_fill_viridis_c()+
    theme_minimal(base_size=20,base_family = 'Arial')+
    theme(text= element_text(size=20,family = 'Arial'),
          legend.position = 'right')) +
  (ggplot(Zh ,aes(x=V1,y=V2,fill=fibrosis))+
     geom_point(size=3,shape=21,stroke=1.3)+
     scale_fill_viridis_c()+
     theme_minimal(base_size=20,base_family = 'Arial')+
     theme(text= element_text(size=20,family = 'Arial'),
           legend.position = 'right'))
ggsave(paste0('results/projected_',
              tolower(ref_dataset),
              '_samples_on_extra_basis',
              tolower(target_dataset),
              '.png'),
       height = 12,
       width = 16,
       units = 'in',
       dpi=600)

### Find translatable LV of the in vitro system
### Run evolutionary algorithm
Wm_combo <- get_translatable_LV(Xh, Yh, Wh, Wm,
                              rbind(apply(Yh,2,mean),Bh),
                              find_extra = FALSE,
                              verbose = TRUE)
Wm_combo <- Wm_combo$Wm_new
rownames(Wm_combo) <- rownames(Wm)
Zh <- as.data.frame(Xh %*% Wm_combo)
Zh$NAS <- Yh[,1]
Zh$fibrosis <- Yh[,2]
Zh <- Zh %>% rownames_to_column('sample')
cat(paste0('Correlation with NAS:',
             '\nLV1 = ',cor(Zh$LV_data1,Zh$NAS),' , LV2 = ',cor(Zh$LV_data2,Zh$NAS),
             '\nCorrelation with fibrosis:',
             '\nLV1 = ',cor(Zh$LV_data1,Zh$fibrosis),' , LV2 = ',cor(Zh$LV_data2,Zh$fibrosis),
             '\n'))
y_hat <- cbind(1, Xh %*% cbind(Wm, Wm_combo) %*% t(cbind(Wm, Wm_combo)) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
plot(Yh[,1],y_hat[,1],xlab='True NAS',ylab='Predicted NAS',main=paste0('r = ',cor(Yh[,1],y_hat[,1])))
plot(Yh[,2],y_hat[,2],xlab='True fibrosis',ylab='Predicted fibrosis',main=paste0('r = ',cor(Yh[,2],y_hat[,2])))
(ggplot(Zh ,aes(x=LV_data1,y=LV_data2,colour=NAS))+
    geom_point()+
    scale_color_viridis_c()+
    theme_minimal(base_size=20,base_family = 'Arial')+
    theme(text= element_text(size=20,family = 'Arial'),
          legend.position = 'right')) +
  (ggplot(Zh ,aes(x=LV_data1,y=LV_data2,colour=fibrosis))+
     geom_point()+
     scale_color_viridis_c()+
     theme_minimal(base_size=20,base_family = 'Arial')+
     theme(text= element_text(size=20,family = 'Arial'),
           legend.position = 'right'))
ggsave(paste0('results/',tolower(target_dataset),'_projected_',tolower(ref_dataset),'_samples_on_translatable_lvs.png'),
       height = 12,
       width = 16,
       units = 'in',
       dpi=600)
## Check directions for translatable LVs
thetas <- seq(0,360,5)
df_proj <- data.frame()
df <- data.frame()
for (theta in thetas){
  u <- c(cos(theta * pi / 180),sin(theta * pi / 180))
  u <- as.matrix(u)
  proj <- u %*% t(u) # since it is unit norm vector
  tmp <- as.matrix(Zh%>% select(LV_data1,LV_data2))
  Z_proj <- tmp %*% proj
  if (theta %in% c(90,270)){
    corr_nas <- cor(Z_proj[,2],Zh$NAS)
    corr_fib <- cor(Z_proj[,2],Zh$fibrosis)
  }else{
    corr_nas <- cor(Z_proj[,1],Zh$NAS)
    corr_fib <- cor(Z_proj[,1],Zh$fibrosis)
  }
  df <- rbind(df,
              data.frame(theta = theta,phenotype = 'NAS',corr=corr_nas),
              data.frame(theta = theta,phenotype = 'fibrosis',corr=corr_fib))
  df_proj <- rbind(df_proj,
                   data.frame(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1],phenotype = 'NAS',corr=corr_nas),
                   data.frame(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1],phenotype = 'fibrosis',corr=corr_fib))
}
p1 <- ggplot(df %>% spread('phenotype','corr') %>% group_by(theta) %>% 
               mutate(`absolute average`=0.5*(abs(NAS)+abs(fibrosis))) %>% 
               ungroup() %>% gather('phenotype','corr',-theta),
             aes(x=theta,y=corr,colour=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = seq(-1,1,0.25))+
  scale_x_continuous(breaks = seq(0,360,15))+
  xlab(expression(theta*" (\u00B0)"))+
  theme_minimal(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family='Arial'),
        legend.position = 'top')
p2 <- ggplot(df_proj,aes(x=LV_opt_1,y=LV_opt_2,fill=corr,colour=corr))+
  # geom_point()+
  geom_segment(aes(x=0,y=0,xend =LV_opt_1,yend=LV_opt_2),
               arrow = arrow(length = unit(0.03, "npc")))+
  scale_fill_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-1,1))+
  scale_color_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-1,1))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~phenotype)+
  theme_minimal(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family='Arial'),
        legend.position = 'right')
p1 / p2
ggsave(paste0('results/pc_loadings_scores_analysis/',tolower(target_dataset),'_translatable_lvs_space_phenotypes_correlations_directions.png'),
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)


### Save found result
rownames(Wm_tot) <- rownames(Wm)
rownames(Wm_opt) <- rownames(Wm)
rownames(Wm_combo) <- rownames(Wm)
rownames(Wh) <- rownames(Wm)

# saveRDS(Wm_tot,paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
# saveRDS(Wm_opt,paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
# saveRDS(Wm_combo,paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
# data.table::fwrite(as.data.frame(Wh),paste0('results/Wh_',tolower(ref_dataset),'.csv'),row.names = TRUE)
# data.table::fwrite(as.data.frame(Wm_opt),paste0('results/Wm_',tolower(target_dataset),'_extra.csv'),row.names = TRUE)


### Visualize activites in the new dimension
path_acitivity <- decoupleR::run_viper(t(Xm), net_prog,minsize = 1,verbose = FALSE)
path_acitivity <- path_acitivity %>% select(c('Pathway'='source'),condition,score,p_value)
path_acitivity <- left_join(path_acitivity,
                            data_list$Kostrzewski$metadata %>% select(sampleName,treatment,TGF,LPS,Cholesterol,Fructose),
                            by=c('condition'='sampleName'))
path_acitivity <- left_join(path_acitivity,
                            as.data.frame(Xm %*% Wm_opt) %>% rownames_to_column('condition'))
# path_acitivity <- path_acitivity %>% filter(p_value<0.01)
ggscatter(path_acitivity %>% mutate(p=ifelse(p_value<0.0001,'<0.0001',
                                             ifelse(p_value<0.001,'<0.001',
                                                    ifelse(p_value<0.01,'<0.01',
                                                           ifelse(p_value<0.05,'<0.05',
                                                                  'ns'))))) %>%
            mutate(p = factor(p,levels = c('ns','<0.05','<0.01','<0.001','<0.0001'))) ,
          x='V1','V2',fill = 'score',size ='p',shape=21) +
  scale_fill_gradient2(low='blue',high='green',midpoint = 0,mid='white') +
  theme(text = element_text(size=20,family='Arial'))+
  facet_wrap(~Pathway)
ggsave(paste0('results/',tolower(target_dataset),'pathway_activity_in_extra_space.png'),
       height = 12,
       width = 16,
       units = 'in',
       dpi=600)

ggscatter(left_join(path_acitivity,
                    as.data.frame(Xm %*% Wm_opt) %>% rownames_to_column('condition')) %>%
            mutate(p=ifelse(p_value<0.0001,'<0.0001',
                            ifelse(p_value<0.001,'<0.001',
                                   ifelse(p_value<0.01,'<0.01',
                                          ifelse(p_value<0.05,'<0.05',
                                                 'ns'))))) %>%
            mutate(p = factor(p,levels = c('ns','<0.05','<0.01','<0.001','<0.0001'))) %>%
            gather('LV','value',-p,-condition,-score,-treatment,-TGF,-LPS,-Cholesterol,-Fructose,-Pathway,-p_value) %>%
            filter(LV=='V1'),
          x='value',y='score',cor.coef = TRUE) +
  # scale_fill_gradient2(low='blue',high='green',midpoint = 0,mid='white') +
  theme(text = element_text(size=20,family='Arial'))+
  facet_wrap(~Pathway)



### Analyze gene loadings----------------------------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
plot_extra_gene_loadings_lv1 <- plot_gene_loadings(loadings = Wm_opt,
                                               selection='V1',
                                               y_lab = 'weight in extra LV1')
ggsave(paste0('results/pc_loadings_scores_analysis/gene_woptimal_LV1_',
              tolower(target_dataset),
              '_loadings.png'),
       plot = plot_extra_gene_loadings_lv1,
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

plot_extra_gene_loadings_lv2 <- plot_gene_loadings(Wm_opt,
                                               selection='V2',
                                               y_lab = 'weight in extra LV2')
ggsave(paste0('results/pc_loadings_scores_analysis/gene_woptimal_LV2_',
              tolower(target_dataset),
              '_loadings.png'),
       plot = plot_extra_gene_loadings_lv2,
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

plot_translatable_gene_loadings_lv1 <- plot_gene_loadings(Wm_combo,
                                                      colnames(Wm_combo)[1],
                                                      'translatable LV1')
ggsave(paste0('results/pc_loadings_scores_analysis/gene_wcombo_LV1',
              tolower(target_dataset),
              '_loadings.png'),
       plot = plot_translatable_gene_loadings_lv1,
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

plot_translatable_gene_loadings_lv2 <- plot_gene_loadings(Wm_combo,
                                                      colnames(Wm_combo)[2],
                                                      'translatable LV2')
ggsave(paste0('results/pc_loadings_scores_analysis/gene_wcombo_LV2',
              tolower(target_dataset),
              '_loadings.png'),
       plot = plot_translatable_gene_loadings_lv2,
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)


### Analyze loadings at the TF activity level--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
dorotheaData = read.table('data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
extra_basis_TF_activity <- TF_activity_interpretation(Wm_opt,
                                                      Wm,
                                                      dorotheaData)

ggsave(paste0('results/pc_loadings_scores_analysis/tfs_only_optimal_loadings_',
              tolower(target_dataset),
              '_barplot.png'),
       plot = extra_basis_TF_activity$figure,
       width = 16,
       height = 12,
       dpi = 600)
### Analyze loadings at the Pathway activity level--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
extra_basis_pathway_activity <- pathway_activity_interpretation(Wm_opt,
                                                                Wm)
ggsave(paste0('results/pc_loadings_scores_analysis/progenies_only_optimal_loadings_',
              tolower(target_dataset),
              '_barplot.png'),
       plot = extra_basis_pathway_activity$figure,
       width = 16,
       height = 12,
       dpi = 600)

### Visualize pathways in various directions of the new space-----------------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
thetas <- seq(0,360,5)
df <- data.frame()
for (theta in thetas){
  u <- c(cos(theta * pi / 180),sin(theta * pi / 180))
  u <- as.matrix(u)
  proj <- u %*% t(u) # since it is unit norm vector
  W_proj <- Wm_opt %*% proj
  extra_basis_pathway_activity <- pathway_activity_interpretation(W_proj,
                                                                  Wm,
                                                                  plotting = FALSE)
  
  tmp <- extra_basis_pathway_activity[[2]]
  tmp <- tmp %>% filter(condition == 'V1') %>% ungroup() %>% 
    select(Pathway,score) %>% mutate(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1]) %>%
    mutate(theta = theta)
  df <- rbind(df,
              tmp)
}
df_plot <- df %>% select(Pathway,theta,score) %>% mutate(score = abs(score)) %>% filter(theta<=90) %>%
  spread('theta','score')
df_labels <- data.frame(score = rep(seq(0,6,2),length(unique(df_plot$Pathway))))
df_paths <- data.frame(Pathway = rep(unique(df_plot$Pathway),4))
df_paths <- df_paths %>% arrange(Pathway)
df_labels <- cbind(df_labels,df_paths)
p1 <- ggRadar(df_plot, aes(group = Pathway), rescale = FALSE, use.label = TRUE,size = 1) +
  # scale_x_discrete(breaks = seq(0,180,30))+
  geom_text(data=df_labels,aes(x=17,y=score,label = score),color='black',alpha=0.7)+
  theme_light() +
  labs(title = "Pathway activity contribution in each direction of the extra basis")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank()) +
  facet_wrap(~Pathway)
# p + geom_text(aes(label = score), size = 5)
print(p1)
ggsave(paste0('results/pc_loadings_scores_analysis/pathways_in_extra_basis_',tolower(target_dataset),'.png'),
       plot = p1,
       width = 9,
       height = 9,
       units = 'in',
       dpi =600)
p2 <- ggplot(df,aes(x=LV_opt_1,y=LV_opt_2,fill=score,colour=score))+
  # geom_point()+
  geom_segment(aes(x=0,y=0,xend =LV_opt_1,yend=LV_opt_2),
               arrow = arrow(length = unit(0.03, "npc")))+
  scale_fill_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-8,8))+
  scale_color_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-8,8))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Pathway)+
  theme_light(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family='Arial'),
        legend.position = 'right')
print(p2)
ggsave(paste0('results/pc_loadings_scores_analysis/pathways_in_extra_basis_arrowplot_',tolower(target_dataset),'.png'),
       plot = p2,
       width = 9,
       height = 9,
       units = 'in',
       dpi =600)
### Analyze loadings with GSEA on MSIG Hallmarks genesets--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Wm_opt), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Wm_opt[inds,])
rownames(meas) <- entrez_ids
msig <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes <- as.data.frame(msig$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark')  #%>% gather('PC','NES',-Hallmark)
msig_nes <- msig_nes %>% gather('LV','NES',-Hallmark)
msig_pval <- as.data.frame(msig$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')#%>% gather('PC','padj',-Hallmark)
msig_pval <- msig_pval %>% gather('LV','padj',-Hallmark)
df_msig <- left_join(msig_nes,msig_pval)
df_msig <- df_msig %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))

p1 <- (ggplot(df_msig %>% filter(LV=='V1') %>% arrange(NES),aes(x=NES,y=reorder(Hallmark,-NES),fill=padj))+ 
         geom_bar(stat = 'identity',color='black',size=1.5) +
         scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in extra basis 1')+
         theme_pubr(base_family = 'Arial',base_size = 15)+
         theme(text = element_text(family = 'Arial',size=15),
               axis.text.y = element_text(size=10),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'right',
               legend.justification = "center"))+
  (ggplot(df_msig %>% filter(LV=='V2') %>% arrange(NES),aes(x=NES,y=reorder(Hallmark,-NES),fill=padj))+ 
     geom_bar(stat = 'identity',color='black',size=1.5) +
     scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
     xlab('Normalized Enrichment Score') +ylab('Hallmark')+
     ggtitle('Hallmarks enriched in extra basis 2')+
     theme_pubr(base_family = 'Arial',base_size = 15)+
     theme(text = element_text(family = 'Arial',size=15),
           axis.text.y = element_text(size=10),
           plot.title = element_text(hjust = 0.5),
           legend.key.size = unit(1.5, "lines"),
           legend.position = 'right',
           legend.justification = "center",
           axis.title.y = element_blank()))
print(p1)
ggsave(paste0('results/pc_loadings_scores_analysis/hallmark_',
              tolower(target_dataset),
              'on_optimal_loadings.png'),
       plot=p1,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)
### Analyze loadings with GSEA on KEGG Pathways genesets--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Wm_opt), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Wm_opt[inds,])
rownames(meas) <- entrez_ids
keggs <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'kegg',
                       n_permutations = 10000,
                       order_columns=F)
kegg_nes <- as.data.frame(keggs$NES$`NES KEGG`) %>% rownames_to_column('pathway')
kegg_nes <- kegg_nes %>% gather('LV','NES',-pathway)
kegg_pval <- as.data.frame(keggs$Pval$`Pval KEGG`) %>% rownames_to_column('pathway')
kegg_pval <- kegg_pval %>% gather('LV','padj',-pathway)
df_keggs <- left_join(kegg_nes,kegg_pval)
df_keggs <- df_keggs %>% mutate(pathway=strsplit(pathway,"_"))
df_keggs <- df_keggs %>% unnest(pathway) %>% filter(!(pathway %in% c("KEGG","FL1000")))
df_keggs <- df_keggs %>% mutate(pathway=as.character(pathway))
df_keggs <- df_keggs %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))

p1 <- (ggplot(df_keggs%>% filter(abs(NES)>1.30) %>% #filter(padj<0.5) %>% 
                filter(LV=='V1') %>% arrange(NES),aes(x=NES,y=reorder(pathway,-NES),fill=padj))+ 
         geom_bar(stat = 'identity',color='black',size=1.5) +
         scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_keggs$padj),1)) +
         xlab('Normalized Enrichment Score') + ylab('KEGG pathway')+
         ggtitle('KEGG Pathways enriched in extra basis 1')+
         theme_pubr(base_family = 'Arial',base_size = 15)+
         theme(text = element_text(family = 'Arial',size=15),
               axis.text.y = element_text(size=10),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'right',
               legend.justification = "center"))+
  (ggplot(df_keggs %>% filter(abs(NES)>1.30) %>% #filter(padj<0.5) %>% 
            filter(LV=='V2') %>% arrange(NES),aes(x=NES,y=reorder(pathway,-NES),fill=padj))+ 
     geom_bar(stat = 'identity',color='black',size=1.5) +
     scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_keggs$padj),1)) +
     xlab('Normalized Enrichment Score') +ylab('KEGG pathway')+
     ggtitle('KEGG Pathways enriched in extra basis 2')+
     theme_pubr(base_family = 'Arial',base_size = 15)+
     theme(text = element_text(family = 'Arial',size=15),
           axis.text.y = element_text(size=10),
           plot.title = element_text(hjust = 0.5),
           legend.key.size = unit(1.5, "lines"),
           legend.position = 'right',
           legend.justification = "center",
           axis.title.y = element_blank()))
print(p1)
ggsave(paste0('results/pc_loadings_scores_analysis/keggs_',
              tolower(target_dataset),
              'on_optimal_loadings.png'),
       plot=p1,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)
### Analyze loadings with GSEA on GO Terms BP genesets--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
colnames(Wm_opt) <- c('V1','V2')
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
gos <- fastenrichment(colnames(Wm_opt),
                        rownames(Wm_opt),
                        Wm_opt,
                        enrichment_space = 'go_bp',
                        gene_id_type = 'symbol',
                        n_permutations = 10000,
                        order_columns=F)
go_nes <- as.data.frame(gos$NES$`NES GO BP`) %>% rownames_to_column('GO')
go_nes <- go_nes %>% gather('LV','NES',-GO)
go_pval <- as.data.frame(gos$Pval$`Pval GO BP`) %>% rownames_to_column('GO')
go_pval <- go_pval %>% gather('LV','padj',-GO)
df_gos <- left_join(go_nes,go_pval)
df_gos <- df_gos %>% mutate(GO=strsplit(GO,"_"))
df_gos <- df_gos %>% unnest(GO) %>% filter(!(GO %in% c("GO","FL1000","BP")))
df_gos <- df_gos %>% mutate(GO=as.character(GO))
df_gos <- left_join(go_annotations %>% select(GO,`GO Terms`) %>% unique(),df_gos)
df_gos <- df_gos %>% filter(!is.na(NES))

p1 <- (ggplot(df_gos%>% filter(abs(NES)>1.75) %>% filter(padj<0.5) %>% 
                filter(LV=='V1') %>% arrange(NES),aes(x=NES,y=reorder(`GO Terms`,-NES),fill=padj))+ 
         geom_bar(stat = 'identity',color='black',size=1.5) +
         scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_gos$padj),1)) +
         xlab('Normalized Enrichment Score') + ylab('GO Terms')+
         ggtitle('GO Terms enriched in extra basis 1')+
         theme_pubr(base_family = 'Arial',base_size = 15)+
         theme(text = element_text(family = 'Arial',size=15),
               axis.text.y = element_text(size=10),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'right',
               legend.justification = "center"))+
  (ggplot(df_gos %>% filter(abs(NES)>1.75) %>% filter(padj<0.5) %>% 
            filter(LV=='V2') %>% arrange(NES),aes(x=NES,y=reorder(`GO Terms`,-NES),fill=padj))+ 
     geom_bar(stat = 'identity',color='black',size=1.5) +
     scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_gos$padj),1)) +
     xlab('Normalized Enrichment Score') +ylab('GO Terms')+
     ggtitle('GO Terms enriched in extra basis 2')+
     theme_pubr(base_family = 'Arial',base_size = 15)+
     theme(text = element_text(family = 'Arial',size=15),
           axis.text.y = element_text(size=10),
           plot.title = element_text(hjust = 0.5),
           legend.key.size = unit(1.5, "lines"),
           legend.position = 'right',
           legend.justification = "center",
           axis.title.y = element_blank()))
print(p1)
ggsave(paste0('results/pc_loadings_scores_analysis/gos_',
              tolower(target_dataset),
              'on_optimal_loadings.png'),
       plot=p1,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)


### Identify external perturbations with ChemPert "regulon"--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
dorotheaData = read.table('data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
### Load relevent TF perturbations
tf_responses <- readRDS("data/ChemPert/Transcriptional_responses.rds")
colnames(tf_responses)[1] <- "Response_ID"
metadata <- read.csv("data/ChemPert/Information_for_transcriptional_responses.csv")
metadata_human <- metadata %>% filter(Species == "Human" & grepl("Hepatocy", Cell_Source))
tf_responses_hep <- tf_responses %>% filter(Response_ID %in% metadata_human$Response_ID)
resp_net <- NULL
for (ii in 1:nrow(tf_responses_hep)){
  if (tf_responses_hep$`Up-regulation number`[ii]>0){
    dummy <- data.frame(source = tf_responses_hep$Response_ID[ii],
                        target = strsplit(tf_responses_hep$`Up-regulation`[ii], split = "; ", fixed = T)
                        %>% unlist(),
                        mor = 1)
  }
  if (tf_responses_hep$`Down-regulation number`[ii]>0){
    dummy <- rbind(dummy,
                   data.frame(source = tf_responses_hep$Response_ID[ii],
                              target = strsplit(tf_responses_hep$`Down-regulation`[ii], split = "; ", fixed = T)
                              %>% unlist(),
                              mor = -1))
  }
  resp_net <- rbind(resp_net, dummy)
}
resp_net <- resp_net %>% mutate(Response_ID = source)
resp_net <- merge(resp_net, metadata_human, by = "Response_ID")

### Run analysis with ChemPert
extra_basis_inferred_perts <- perturnation_activity_inference(Wm_opt,metadata_human,dorotheaData,resp_net)

ggsave(paste0('results/pc_loadings_scores_analysis/optimal_direction_',
              tolower(target_dataset),
              '_perturbation_activity.png'),
       plot = extra_basis_inferred_perts$figure,
       width = 16,
       height = 12,
       dpi = 600)


