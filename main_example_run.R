library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
source('modeling/functions_translation_jose.R')
source('modeling/CrossValidationUtilFunctions.R')
source("utils/plotting_functions.R")
### Load the all the data to be used-----------------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver","Govaere")
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

### Check current TF and pathway activity
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
ggsave('results/pathway_activity_in_pc_space.png',
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
ggsave('results/projected_govaere_samples_on_translatable_lvs.png',
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
ggsave('results/pc_loadings_scores_analysis/translatable_lvs_space_phenotypes_correlations_directions.png',
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)


### Save found result
saveRDS(Wm_tot,'results/Wm_total.rds')
saveRDS(Wm_opt,'results/Wm_extra.rds')
saveRDS(Wm_combo,'results/Wm_combo.rds')

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
ggsave('results/pathway_activity_in_extra_space.png',
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

### Analyze loadings at the TF activity level--------------

### Analyze loadings at the Pathway activity level--------------

### Analyze loadings with GSEA on MSIG Hallmarks genesets--------------

### Analyze loadings with GSEA on KEGG Pathways genesets--------------

### Analyze loadings with GSEA on GO Terms BP genesets--------------