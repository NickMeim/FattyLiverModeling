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

################################################################################
## Load the all the data to be used
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
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
if (ref_dataset=='Hoang'){
  sex_inferred <- data_list[[ref_dataset]]$metadata$sex
}else if (ref_dataset=='Pantano') {
  sex_inferred <- data_list[[ref_dataset]]$metadata$Sex
}else{
  sex_inferred <- apply(as.matrix(Xh[,c('RPS4Y1')]),2,sign) 
  sex_inferred <- 1*(sex_inferred>0)
  sex_inferred <- ifelse(sex_inferred==1,'male','female')
}
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# For supplementary data: Print some statistics about NAS and fibrosis distribution between sexes
print(all(rownames(Xh)==rownames(sex_inferred)))
pheno_stats <- as.data.frame(cbind(Yh,sex_inferred))
colnames(pheno_stats) <- c('NAS','fibrosis','sex')
pheno_stats$NAS <- as.numeric(pheno_stats$NAS)
pheno_stats$fibrosis <- as.numeric(pheno_stats$fibrosis)

plt_pheno_stats <- ggboxplot(pheno_stats %>% gather('phenotype','score',-sex) %>% 
                             mutate(phenotype= ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype)),
                              x='sex',y='score',add='jitter', add.params = list(size = size_dot), size = size_line)+
                      stat_compare_means(comparisons = list(c('male','female')),
                      method='wilcox',
                      size = size_annotation)+
                      facet_wrap(~phenotype,scales="free_y")
saveRDS(plt_pheno_stats, "Figures/plt_pheno_stats.rds")


### Project in-vitro data to their group-derived PC space
Zm <- Xm %*% Wm
per_var <- round(100 * colVars(Zm)/sum(colVars(Xm)), 2)
plt_PCA_MPS <- ggplot(data.frame(Zm), aes(x = PC1,y = PC2)) +
                geom_point(color = '#FD0714',size = size_dot)+
                xlab(paste0('MPS PC1 (',per_var[1],'%)')) + ylab(paste0('MPS PC2 (',per_var[2],'%)')) 
saveRDS(plt_PCA_MPS, "Figures/plt_PCA_MPS.rds")


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

### Repeat for Wm_combo
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
colnames(Wm_combo) <- c('TC1','TC2')
path_acitivity <- decoupleR::run_viper(t(Xm), net_prog,minsize = 1,verbose = FALSE)
path_acitivity <- path_acitivity %>% select(c('Pathway'='source'),condition,score,p_value)
path_acitivity <- left_join(path_acitivity,
                            data_list$Kostrzewski$metadata %>% select(sampleName,treatment,TGF,LPS,Cholesterol,Fructose),
                            by=c('condition'='sampleName'))
path_acitivity <- left_join(path_acitivity,
                            as.data.frame(Xm %*% Wm_combo) %>% select(TC1,TC2) %>% rownames_to_column('condition'))
ggscatter(path_acitivity %>% mutate(p=ifelse(p_value<0.0001,'<0.0001',
                                             ifelse(p_value<0.001,'<0.001',
                                                    ifelse(p_value<0.01,'<0.01',
                                                           ifelse(p_value<0.05,'<0.05',
                                                                  'ns'))))) %>%
            mutate(p = factor(p,levels = c('ns','<0.05','<0.01','<0.001','<0.0001'))) ,
          x='TC1','TC2',fill = 'score',size ='p',shape=21) +
  scale_fill_gradient2(low='blue',high='green',midpoint = 0,mid='white') +
  theme(text = element_text(size=20,family='Arial'))+
  facet_wrap(~Pathway)
ggsave(paste0('results/',tolower(target_dataset),'_pathway_activity_in_translatable_components.png'),
       height = 12,
       width = 16,
       units = 'in',
       dpi=600)


### Run PLSR and find extra basis--------------------------
plsr_model <- opls(x = Xh, 
                   y = Yh,
                   predI = 8,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
total_plsr_var <- sum(colVars(plsr_model@scoreMN))/sum(colVars(Xh))
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

### Plot human PLSR space and projected truncated data
Zh_plsr <- plsr_model@scoreMN
mcor <- cor(cbind(Zh_plsr,Yh))
mcor[upper.tri(mcor,diag = TRUE)] <- 100
mcor <- reshape2::melt(mcor)
mcor <- mcor %>% filter(value != 100)
mcor <- mcor %>% mutate(keep = ifelse(Var1=='fibrosis' | Var2=='fibrosis' | Var1=='NAS' | Var2=='NAS',TRUE,FALSE)) %>%
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(keep = ifelse((Var1=='fibrosis' & Var2=='NAS') | (Var1=='NAS' & Var2=='fibrosis'),FALSE,TRUE))%>%
  filter(keep==TRUE) %>% select(-keep) 
mcor <- mcor %>% group_by(Var2) %>% mutate(avg_abs_cor = mean(abs(value))) %>% select(c('LV'='Var2'),avg_abs_cor) %>%
  unique()
lvs <- mcor$LV[order(-mcor$avg_abs_cor)]
lvs <- as.character(lvs[1:2])
lv_vars <- round(100 * colVars(Zh_plsr)/sum(colVars(Xh)),2) 
lv1_var <- round(100 * var(Zh_plsr[,c(lvs[1])])/sum(apply(Xm,2,var)),2) 
lv2_var <- round(100 * var(Zh_plsr[,c(lvs[2])])/sum(apply(Xm,2,var)),2) 
invivo_plsr <- as.data.frame(Zh_plsr[,c(lvs[1],lvs[2])])
invivo_plsr <- cbind(invivo_plsr,as.data.frame(Yh))
invivo_plsr <- invivo_plsr %>% gather('phenotype','Score',-all_of(lvs))
colnames(invivo_plsr)[1:2] <- c('V1','V2')
invivo_plsr <- invivo_plsr %>% 
                mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype)) %>%
                group_by(phenotype) %>% 
                mutate(normed_score=Score/max(Score))


plt_PLSR_human <- ggplot(invivo_plsr, aes(x=V1,y=V2,color=normed_score)) +
                    geom_point(size = size_dot)+
                    scale_color_viridis_c()+
                    labs(color = 'Score')+
                    xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv_vars[1],'%)')) + ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv_vars[2],'%)')) +
                    facet_wrap(~phenotype)
saveRDS(plt_PCA_MPS, "Figures/plt_PLSR_human.rds")

### visualize sex separation in PLSR
plt_PLSR_sex <- ggplot(cbind(as.data.frame(Zh_plsr[,c(lvs[1],lvs[2])]), #c(lvs[1],lvs[2])
                             as.data.frame(sex_inferred) %>% select(c('sex'='V1'))),
                         aes(x=p1,y=p2,color=sex)) +
  geom_point(size = size_dot)+
  scale_color_viridis_d()+
  stat_ellipse(size = size_line, show.legend = F) +
  labs(color = 'Sex')+
  xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv_vars[1],'%)')) + ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv_vars[2],'%)'))

saveRDS(plt_PLSR_sex, "Figures/plt_PLSR_sex.rds")


## Visualize plsr for backprojected data
# Extract the necessary matrices from the opls object
P <- plsr_model@loadingMN  # Loadings matrix (P)
# plsr_scores <- plsr_model@scoreMN    # Scores matrix (T)
W <- plsr_model@weightMN   # Weights matrix (W)
# C <- plsr_model@loadingYMN # Y-loadings matrix (C)
# E <- plsr_model@residualsMN # Residuals matrix (E)
# U <- plsr_model@orthoScoreMN # Orthogonal scores matrix (U)
# Manually calculate the scores (T)
# T = X * W * (P' * W)^-1
Zh_plsr_backprj <- Xh %*% Wm %*% t(Wm)  %*% Wh %*% solve(t(P) %*% W)
invivo_plsr_backprj <- as.data.frame(Zh_plsr_backprj[,c(lvs[1],lvs[2])])
invivo_plsr_backprj <- cbind(invivo_plsr_backprj,as.data.frame(Yh))
invivo_plsr_backprj <- invivo_plsr_backprj %>% gather('phenotype','Score',-all_of(lvs))
colnames(invivo_plsr_backprj)[1:2] <- c('V1','V2')
invivo_plsr_backprj <- invivo_plsr_backprj %>%
                        mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype)) %>%
                        group_by(phenotype) %>% 
                        mutate(normed_score=Score/max(Score))

plt_PLSR_backproject <- ggplot(invivo_plsr_backprj,aes(x=V1,y=V2,color=normed_score)) +
                        geom_point(size = size_dot)+
                        scale_color_viridis_c()+
                        labs(color = 'Score')+
                        xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv_vars[1],'%)')) + ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv_vars[2],'%)')) +
                        facet_wrap(~phenotype)
saveRDS(plt_PLSR_backproject, "Figures/plt_PLSR_backproject.rds")


################################################################################
## Find extra latent variables
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
  # proj <- u %*% t(u) # since it is unit norm vector
  # tmp <- as.matrix(Zh %>% select(V1,V2))
  # Z_proj <- tmp %*% proj
  Wproj <- Wm_opt %*% u
  # if (theta %in% c(90,270)){
  #   corr_nas <- cor(Z_proj[,2],Zh$NAS)
  #   corr_fib <- cor(Z_proj[,2],Zh$fibrosis)
  # }else{
  #   corr_nas <- cor(Z_proj[,1],Zh$NAS)
  #   corr_fib <- cor(Z_proj[,1],Zh$fibrosis)
  # }
  Wtmp <- cbind(Wm,Wproj)
  Yhat <- cbind(1, Xh %*% Wtmp %*% t(Wtmp) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
  corr_nas <- cor(Yhat[,1],Zh$NAS)
  corr_fib <- cor(Yhat[,2],Zh$fibrosis)
  df <- rbind(df,
              data.frame(theta = theta,phenotype = 'NAS',corr=corr_nas),
              data.frame(theta = theta,phenotype = 'fibrosis',corr=corr_fib))
  # df_proj <- rbind(df_proj,
  #                  data.frame(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1],phenotype = 'NAS',corr=corr_nas),
  #                  data.frame(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1],phenotype = 'fibrosis',corr=corr_fib))
}
# p1 <- ggplot(df %>% spread('phenotype','corr') %>% group_by(theta) %>% 
#                mutate(`absolute average`=0.5*(abs(NAS)+abs(fibrosis))) %>% 
#                ungroup() %>% gather('phenotype','corr',-theta),
#              aes(x=theta,y=corr,colour=phenotype)) +
#   geom_line()+
#   geom_point()+
#   geom_hline(yintercept = 0) +
#   scale_y_continuous(breaks = seq(-1,1,0.25))+
#   scale_x_continuous(breaks = seq(0,360,30))+
#   xlab(expression(theta*" (\u00B0)"))+
#   ylab('correlation') +
#   theme_minimal(base_size = 24,base_family = 'Arial')+
#   theme(text = element_text(size=24,family='Arial'),
#         legend.position = 'top',
#         legend.text =  element_text(size=26))
# 
# p2 <- ggplot(df_proj,aes(x=LV_opt_1,y=LV_opt_2,fill=corr,colour=corr))+
#   geom_segment(aes(x=0,y=0,xend =LV_opt_1,yend=LV_opt_2),
#                arrow = arrow(length = unit(0.03, "npc")),size=1)+
#   scale_fill_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-1,1))+
#   scale_color_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-1,1))+
#   xlab('extra basis 1') + ylab('extra basis 2') +
#   labs(fill = 'correlation',colour='correlation')+
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   facet_wrap(~phenotype)+
#   theme_minimal(base_size = 24,base_family = 'Arial')+
#   theme(text = element_text(size=24,family='Arial'),
#         strip.text=element_text(size=24),
#         legend.position = 'right')
# p1 / p2
# ggsave(paste0('results/pc_loadings_scores_analysis/',
#               tolower(ref_dataset),
#               'optimal_space_phenotypes_correlations_directions',
#               tolower(target_dataset),
#               '.png'),
#        height = 12,
#        width = 12,
#        units = 'in',
#        dpi=600)

scatter_box_plot <- function(df,legend_title,
                             font_size =20,font_family = 'Arial',
                             point_shape = 21,point_size=2.8,point_stroke=1.2,
                             x_axis='LV extra 1',y_axis='LV extra 2',
                             box_y_width = 0.2,jitter_y_width=0.1,
                             jitter_x_height = 0.2,
                             theme_use = 'minimal',
                             plotting=TRUE){
  if (theme_use=='bw'){
    th <- ggplot2::theme_bw(base_size = font_size, base_family = font_family)
  }else{
    th <- ggplot2::theme_minimal(base_size = font_size, base_family = font_family)
  }
  scatter_plot <- ggplot(df, aes(x = V1, y = V2, fill = pheno)) +
    geom_point(size = point_size, shape = point_shape, stroke = point_stroke) +
    scale_fill_viridis_c() +
    xlab(x_axis) +
    ylab(y_axis) +
    labs(fill=legend_title)+
    th +
    theme(
      text = element_text(size = font_size, family = font_family),
      legend.position = 'top',
    )
  
  # Boxplot for V1 (x-axis)
  boxplot_x <- ggplot(df, aes(x = V1, y = "", color = pheno)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(height = jitter_x_height)+
    scale_color_viridis_c() +
    theme_minimal(base_size = font_size, base_family = font_family) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    )
  
  # Boxplot for V2 (y-axis)
  boxplot_y <- ggplot(df, aes(x = "", y = V2, color = pheno)) +
    geom_boxplot(outliers = FALSE,width = box_y_width) +
    geom_jitter(width=jitter_y_width)+
    scale_color_viridis_c() +
    theme_minimal(base_size = font_size, base_family = font_family) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    )
  # Combine plots for NAS
  combo1 <- scatter_plot + boxplot_y + plot_layout(widths = c(3, 1))
  combo2 <- boxplot_x + plot_spacer() + plot_layout(widths = c(3, 1))
  combined_plot <- (combo1) / (combo2) +
    plot_layout(heights = c(4, 1))
  if (plotting==TRUE){
    print(combined_plot)
  }
  return(combined_plot)
}
# ggMarginal(scatter_plot,type = 'box', groupColour = TRUE, groupFill = TRUE)
combined_plot_nas <- scatter_box_plot(Zh %>% mutate(pheno=NAS),'NAS')
combined_plot_fibrosis <- scatter_box_plot(Zh %>% mutate(pheno=fibrosis),'Fibrosis stage')
# p <- (ggplot(Zh ,aes(x=V1,y=V2,fill=NAS))+
#     geom_point(size=2.8,shape=21,stroke=1.2)+
#     scale_fill_viridis_c()+
#     xlab('LV extra 1')+ylab('LV extra 2')+  
#     theme_minimal(base_size=20,base_family = 'Arial')+
#     theme(text= element_text(size=20,family = 'Arial'),
#           legend.position = 'top')) +
#   (ggplot(Zh ,aes(x=V1,y=V2,fill=fibrosis))+
#      geom_point(size=2.8,shape=21,stroke=1.2)+
#      scale_fill_viridis_c()+
#      xlab('LV extra 1')+ylab('LV extra 2')+ 
#      labs(fill='Fibrosis stage') +
#      theme_minimal(base_size=20,base_family = 'Arial')+
#      theme(text= element_text(size=20,family = 'Arial'),
#            legend.position = 'top'))
# print(p)
ggsave(paste0('figures/projected_',
              tolower(ref_dataset),
              '_samples_on_extra_basis_NAS',
              tolower(target_dataset),
              '.png'),
       plot = combined_plot_nas,
       height = 5,
       width = 12,
       units = 'in',
       dpi=600)
ggsave(paste0('figures/projected_',
              tolower(ref_dataset),
              '_samples_on_extra_basis_NAS',
              tolower(target_dataset),
              '.eps'),
       plot = combined_plot_nas,
       device = cairo_ps,
       height = 5,
       width = 12,
       units = 'in',
       dpi=600)
ggsave(paste0('figures/projected_',
              tolower(ref_dataset),
              '_samples_on_extra_basis_fibrosis',
              tolower(target_dataset),
              '.png'),
       plot = combined_plot_fibrosis,
       height = 5,
       width = 12,
       units = 'in',
       dpi=600)
ggsave(paste0('figures/projected_',
              tolower(ref_dataset),
              '_samples_on_extra_basis_fibrosis',
              tolower(target_dataset),
              '.eps'),
       plot = combined_plot_fibrosis,
       device = cairo_ps,
       height = 5,
       width = 12,
       units = 'in',
       dpi=600)

## separation of sex in the extra basis
ggplot(left_join(Zh , as.data.frame(sex_inferred) %>% select(c('sex'='V1')) %>%rownames_to_column('sample')),
       aes(x=V1,y=V2,fill=sex,color=sex)) +
  geom_point(size=2.8,shape=21,stroke=1.2,color='black')+
  scale_fill_viridis_d()+
  xlab('LV extra 1') + ylab('LV extra 2') +
  theme_bw(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family = 'Arial'),
        panel.background = element_rect(fill='white',linetype = 0 ),
        panel.border = element_blank(),
        legend.position = 'right')+
  stat_ellipse(aes(color = sex),size=1.2) +
  scale_color_viridis_d()+
  scale_fill_viridis_d()
ggsave(paste0('figures/',target_dataset,'_to_',ref_dataset,'_extra_basis_sex_separete.eps'),
       device = cairo_ps,
       height = 6,
       width = 9,
       units = 'in',
       dpi = 600)
ggsave(paste0('figures/',target_dataset,'_to_',ref_dataset,'_extra_basis_sex_separete.png'),
       height = 6,
       width = 9,
       units = 'in',
       dpi = 600)
### Can you even predict sex extra LVs ?
sex_model <- opls(x = as.matrix(left_join(Zh , 
                                as.data.frame(sex_inferred) %>% 
                                  select(c('sex'='V1')) %>%
                                  rownames_to_column('sample')) %>%
                    select(V1,V2)), 
                   y = as.matrix(left_join(Zh , 
                                 as.data.frame(sex_inferred) %>% 
                                   select(c('sex'='V1')) %>%
                                   rownames_to_column('sample')) %>%
                    select(sex) %>% mutate(sex=factor(sex))),
                   predI = 8,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
sex_predicted <- predict(sex_model,
                         as.matrix(left_join(Zh , 
                                             as.data.frame(sex_inferred) %>% 
                                               select(c('sex'='V1')) %>%
                                               rownames_to_column('sample')) %>%
                                     select(V1,V2)))
sex_predicted <- factor(sex_predicted,levels = c('male','female'))
sex_inferred2 <- factor(sex_inferred,levels = levels(sex_predicted))
print(confusionMatrix(data=sex_predicted,
                reference = sex_inferred2,
                positive = 'male'))
# try also with RF
rf_fit <- train(sex ~ ., 
                data = left_join(Zh , 
                                 as.data.frame(sex_inferred) %>%
                                   select(c('sex'='V1')) %>%
                                   rownames_to_column('sample')) %>%
                  select(V1,V2,sex), 
                method = "rf")
print(rf_fit$results)

df_mu_radial <- df %>% group_by(theta) %>% mutate(mu = mean(corr)) %>% ungroup() %>% 
  select(-corr,-phenotype) %>% unique() %>%
  mutate(phenotype = 'average') %>% select(theta,phenotype,c('corr'='mu')) %>%
  mutate(x = corr*cos(theta), y = corr*sin(theta))
df_radial <- df %>% 
  mutate(x = corr*cos(theta), y = corr*sin(theta)) %>%
  mutate(phenotype = ifelse(phenotype == "fibrosis", "Fibrosis stage",phenotype)) #%>% mutate(theta = pi*theta/180)
plt_cor_radial <-   ggplot(df_radial %>% filter(theta<=90), 
                           aes(x = theta, y = corr, color = phenotype, group = phenotype)) + 
  geom_line(linewidth = 1.5) +
  geom_line(data = df_mu_radial %>% filter(theta<=90),
             aes(x = theta, y = corr),
             color='black',lwd=1.5,linetype='dashed',
            inherit.aes = FALSE)+
  scale_x_continuous(breaks = seq(0,90,15))+
  geom_hline(yintercept = 0.4,color='black',lwd=1)+
  geom_vline(xintercept = 90,color='black',lwd=1)+
  coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F, direction = 1) +
  labs(y = "Pearson correlation", x = 'LV extra 2') +
  ggtitle('LV extra 1')+
  theme_bw() +
  theme(text = element_text(size = 20),
        plot.title = element_text(vjust = -15,hjust = 0.1,size=18,face = 'bold'),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(hjust = 0.66),
        axis.text.x = element_blank(),
        axis.title.x = element_text(vjust = 9,hjust = 1.05,size=18,face='bold'),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -1.5, unit = "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 1.2),
        panel.grid.minor = element_line(linewidth = 1.2)) +
  scale_color_brewer(palette = "Dark2")
print(plt_cor_radial)
ggsave('figures/quarter_radial_correlation.png',
       plot = plt_cor_radial,
       height = 6,
       width = 9,
       units = 'in',
       dpi = 600)
ggsave('figures/quarter_radial_correlation.eps',
       plot = plt_cor_radial,
       device = cairo_ps,
       height = 6,
       width = 9,
       units = 'in',
       dpi = 600)
# plt_cor_radial2 <-   ggplot(df_radial, aes(x = theta/(360), y = corr, color = phenotype, group = phenotype)) + 
#   geom_line(linewidth = 1.5) +
#   theme(text = element_text(size = 20),
#         axis.text.y = element_text(color = "black"),
#         axis.text.x = element_text(color = "black"),
#         legend.key.size = unit(size_legend, "cm")) +
#   scale_color_brewer(palette = "Dark2") +
#   labs(x = "Mixing proportion", y = "Pearson correlation", color = "Phenotype")
# plt_cor_radial2 <-  add_theme(plt_cor_radial2,size_text=20)
# print(plt_cor_radial2)

### Find translatable LV of the in vitro system
### Run evolutionary algorithm
Wm_combo <- get_translatable_LV(Xh, Yh, Wh, Wm,
                              rbind(apply(Yh,2,mean),Bh),
                              find_extra = FALSE,
                              verbose = TRUE)
Wm_combo <- Wm_combo$Wm_new
rownames(Wm_combo) <- rownames(Wm)

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
## Project into translatable latent variables
Wm_combo <- readRDS("results/Wm_kostrzewski_combo.rds")
# Project human data on these components
df_tc <- rbind(data.frame(x = Xh %*% Wm_combo[,1], y = Xh %*% Wm_combo[,2], phenotype = "Fibrosis stage", Score = Yh[,2]),
               data.frame(x = Xh %*% Wm_combo[,1], y = Xh %*% Wm_combo[,2], phenotype = "NAS", Score = Yh[,1]))
tc_fibrosis_scatter_boxplot <- scatter_box_plot(df_tc %>% filter(phenotype!='NAS') %>%
                                             mutate(pheno=Score) %>%
                                             mutate(V1=x) %>%
                                             mutate(V2=y) %>%
                                             select(-phenotype,-x,-y,-Score) %>%
                                             unique(),
                                           legend_title = 'Fibrosis stage',
                                           x_axis = 'MPS TC1',
                                           y_axis = 'MPS TC2',
                                           theme_use = 'bw')
ggsave('figures/AllData_TCs_Scatterplot_MPS_fibrosis.eps',
       plot = tc_fibrosis_scatter_boxplot,
       device = cairo_ps,
       height = 6,
       width=9,
       units = 'in',
       dpi=600)
tc_nas_scatter_boxplot <- scatter_box_plot(df_tc %>% filter(phenotype=='NAS') %>%
                                             mutate(pheno=Score) %>%
                                             mutate(V1=x) %>%
                                             mutate(V2=y) %>%
                                             select(-phenotype,-x,-y,-Score) %>%
                                             unique(),
                                           legend_title = 'NAS',
                                           x_axis = 'MPS TC1',
                                           y_axis = 'MPS TC2',
                                           theme_use = 'bw')
ggsave('figures/AllData_TCs_Scatterplot_MPS_nas.eps',
       plot = tc_nas_scatter_boxplot,
       device = cairo_ps,
       height = 6,
       width=9,
       units = 'in',
       dpi=600)
# ptc <- (ggplot(df_tc %>% filter(phenotype!='NAS'),aes(x=x,y=y,fill=Score))+
#           geom_point(size=2.8,shape=21,stroke=1.2)+
#           scale_fill_viridis_c()+
#           xlab('MPS TC1')+ylab('MPS TC2')+  
#           labs(fill='Fibrosis stage')+
#           theme_bw(base_size=28,base_family = 'Arial')+
#           theme(text= element_text(size=28,family = 'Arial'),
#                 axis.title.x = element_blank(),
#                 legend.position = 'top')) +
#   (ggplot(df_tc%>% filter(phenotype =='NAS') ,aes(x=x,y=y,fill=Score))+
#      geom_point(size=2.8,shape=21,stroke=1.2)+
#      scale_fill_viridis_c()+
#      labs(fill='NAS') +
#      xlab('MPS TC1')+ylab('MPS TC2')+  
#      theme_bw(base_size=28,base_family = 'Arial')+
#      theme(text= element_text(size=28,family = 'Arial'),
#            axis.title.y = element_blank(),
#            axis.title.x = element_text(hjust = -0.8),
#            legend.position = 'top'))
# print(ptc)
# ggsave('figures/AllData_TCs_Scatterplot_human.eps',
#        plot = ptc,
#        device = cairo_ps,
#        height = 6,
#        width=12,
#        units = 'in',
#        dpi=600)
# take a look also at invitro separation in the TCs
df_tc_mps <- data.frame(x = Xm %*% Wm_combo[,1], y = Xm %*% Wm_combo[,2], TGF = data_list[[target_dataset]]$metadata$TGF)
ptc_mps <- ggplot(df_tc_mps,aes(x=x,y=y,fill=TGF))+
          geom_point(size=2.8,shape=21,stroke=1.2)+
          scale_fill_manual(values = c('#66C2A5','#FC8D62'))+
          xlab('MPS TC1')+ylab('MPS TC2')+  
          theme_bw(base_size=28,base_family = 'Arial')+
          theme(text= element_text(size=28,family = 'Arial'),
                plot.title = element_text(hjust = 0.5,face='bold'),
                legend.position = 'top')
print(ptc_mps)
ggsave('figures/AllData_TCs_Scatterplot_MPS.eps',
       plot = ptc_mps,
       device = cairo_ps,
       height = 6,
       width=6,
       units = 'in',
       dpi=600)

## take a look also at performance of back-projection through the TCs
Th <- Xh %*% Wm_combo %*% t(Wm_combo) %*% Wh
Th0 <- Xh %*% Wm %*% t(Wm) %*% Wh
Y_pred_backproj <- cbind(1, Th)  %*% rbind(apply(Yh,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
Y_pred_truncated <-  cbind(1, Th0)  %*% rbind(apply(Yh,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
all_scatter_plot_backproj <- left_join(data.frame(Y_pred_truncated) %>% 
                                         mutate(id = seq(1,nrow(Y_pred_truncated))) %>% 
                                         gather('phenotype','true',-id),
                                       data.frame(Y_pred_backproj) %>% 
                                         mutate(id = seq(1,nrow(Y_pred_backproj))) %>% 
                                         gather('phenotype','prediction',-id)) %>%
  select(-id) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
all_cor_results_backproj <- all_scatter_plot_backproj %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
(ggplot(all_scatter_plot_backproj %>% filter(phenotype=='Fibrosis stage'),aes(x = true,y=prediction)) +
    geom_jitter(width = 0.05,color='#4682B4') + 
    geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
    geom_text(data = all_cor_results_backproj%>% filter(phenotype=='Fibrosis stage'), 
              aes(x = 1, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
              hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
    ylab('prediction using TCs') + xlab('prediction using all PCs')+
    ylim(c(1,3))+ xlim(c(1,3))+
    facet_wrap(~phenotype,scales = 'free')+
    theme_pubr(base_family = 'Arial',base_size=25)+
    theme(text = element_text(family = 'Arial',size=25),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face='bold'),
          panel.grid.major = element_line()))+
  (ggplot(all_scatter_plot_backproj %>% filter(phenotype!='Fibrosis stage'),aes(x = true,y=prediction)) +
     geom_jitter(width = 0.05,color='#4682B4') + 
     geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
     geom_text(data = all_cor_results_backproj%>% filter(phenotype!='Fibrosis stage'), 
               aes(x = 1, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
               hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
     ylab('prediction using TCs') + xlab('prediction using all PCs')+
     ylim(c(1,7))+xlim(c(1,7))+
     facet_wrap(~phenotype,scales = 'free')+
     theme_pubr(base_family = 'Arial',base_size=25)+
     theme(text = element_text(family = 'Arial',size=25),
           axis.title.x = element_text(hjust = -0.6,face='bold'),
           axis.title.y = element_blank(),
           panel.grid.major = element_line())) + 
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(size = 25, family = "Arial", hjust = 0.5,face='bold'))
  )
ggsave('figures/AllData_TCs_VS_PCs_Scatterplot_human_plsr.png',
       height = 6,
       width=9,
       units = 'in',
       dpi=600)
(ggplot(all_scatter_plot_backproj %>% filter(phenotype=='Fibrosis stage'),aes(x = true,y=prediction)) +
    geom_jitter(width = 0.05,color='#4682B4') + 
    geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
    geom_text(data = all_cor_results_backproj%>% filter(phenotype=='Fibrosis stage'), 
              aes(x = 1, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
              hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
    ylab('prediction using TCs') + xlab('prediction using all PCs')+
    ylim(c(1,3))+ xlim(c(1,3))+
    facet_wrap(~phenotype,scales = 'free')+
    theme_pubr(base_family = 'Arial',base_size=25)+
    theme(text = element_text(family = 'Arial',size=25),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face='bold'),
          panel.grid.major = element_line()))+
  (ggplot(all_scatter_plot_backproj %>% filter(phenotype!='Fibrosis stage'),aes(x = true,y=prediction)) +
     geom_jitter(width = 0.05,color='#4682B4') + 
     geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
     geom_text(data = all_cor_results_backproj%>% filter(phenotype!='Fibrosis stage'), 
               aes(x = 1, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
               hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
     ylab('prediction using TCs') + xlab('prediction using all PCs')+
     ylim(c(1,7))+xlim(c(1,7))+
     facet_wrap(~phenotype,scales = 'free')+
     theme_pubr(base_family = 'Arial',base_size=25)+
     theme(text = element_text(family = 'Arial',size=25),
           axis.title.x = element_text(hjust = -3,face='bold'),
           axis.title.y = element_blank(),
           panel.grid.major = element_line())) + 
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(size = 25, family = "Arial", hjust = 0.5,face='bold'))
  )
ggsave('figures/AllData_TCs_VS_PCs_Scatterplot_human_plsr.eps',
       device = cairo_ps,
       height = 6,
       width=12,
       units = 'in',
       dpi=600)

### Analyze gene loadings----------------------------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
plot_extra_gene_loadings_lv1 <- plot_gene_loadings(loadings = Wm_opt,
                                               selection='V1',
                                               y_lab = 'weight in extra LV1',
                                               top=20)
ggsave(paste0('figures/gene_woptimal_LV1_',
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
ggsave(paste0('figures/gene_woptimal_LV2_',
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
ggsave(paste0('figures/gene_wcombo_LV1',
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
ggsave(paste0('figures/gene_wcombo_LV2',
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

ggsave(paste0('figures/tfs_only_optimal_loadings_',
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
ggsave(paste0('figures/progenies_only_optimal_loadings_',
              tolower(target_dataset),
              '_barplot_LV1.png'),
       plot = extra_basis_pathway_activity$figure[[1]] + theme(plot.title = element_blank()),
       width = 9,
       height = 9,
       dpi = 600)
ggsave(paste0('figures/progenies_only_optimal_loadings_',
              tolower(target_dataset),
              '_barplot_LV2.png'),
       plot = extra_basis_pathway_activity$figure[[2]]+ theme(plot.title = element_blank()),
       width = 9,
       height = 9,
       dpi = 600)
setEPS()
postscript(paste0('figures/progenies_only_optimal_loadings_',
                  tolower(target_dataset),
                  '_barplot_LV1.eps'),
           height = 9, width = 9)
print(extra_basis_pathway_activity$figure[[1]]+ theme(plot.title = element_blank()))
dev.off()
setEPS()
postscript(paste0('figures/progenies_only_optimal_loadings_',
                  tolower(target_dataset),
                  '_barplot_LV2.eps'),
           height = 9, width = 9)
print(extra_basis_pathway_activity$figure[[2]]+ theme(plot.title = element_blank()))
dev.off()

colnames(Wm_combo) <- c('V1','V2')
translatable_components_progenies <- pathway_activity_interpretation(Wm_combo,
                                                                       Wm,
                                                                     lim=15)
p1 <- translatable_components_progenies$figure[[1]]
p1 <- p1 + ggtitle('Translatable Component 1')
p2 <- translatable_components_progenies$figure[[2]] + ggtitle('Translatable Component 2')
print(p1+p2)
p <- p1+p2

ggsave(paste0('figures/progenies_only_translatable_components_',
              tolower(target_dataset),
              '.png'),
       plot = p,
       width = 16,
       height = 9,
       dpi = 600)

setEPS()
postscript(paste0('figures/progenies_only_translatable_components_',
                  tolower(target_dataset),
                  '.eps'),
           height = 9, width = 16)
print(p)
dev.off()

### Visualize pathways in various directions of the new space-----------------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
thetas <- seq(0,360,5)
df <- data.frame()
for (theta in thetas){
  u <- c(cos(theta * pi / 180),sin(theta * pi / 180))
  u <- as.matrix(u)
  # proj <- u %*% t(u) # since it is unit norm vector
  # W_proj <- Wm_opt %*% proj
  W_proj <- Wm_opt %*% u
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
# df_plot <- df %>% select(Pathway,theta,score) %>% mutate(score = abs(score)) %>% 
#   mutate(score = ifelse(theta<=90,score,NA)) %>%
#   spread('theta','score')
df_labels <- data.frame(score = rep(seq(0,6,2),length(unique(df_plot$Pathway))))
df_paths <- data.frame(Pathway = rep(unique(df_plot$Pathway),4))
df_paths <- df_paths %>% arrange(Pathway)
df_labels <- cbind(df_labels,df_paths)
lvs <- unique(df_plot$Pathway)
plt_pw_radial <-  ggplot(df_plot  %>% gather('theta','activity',-Pathway) %>% 
                           mutate(theta=as.numeric(theta)) ,
       aes(x = theta, y = activity, color = Pathway,fill=Pathway, group = Pathway)) + 
  geom_line(linewidth = 1.5) +
  geom_area(alpha=0.2)+
  scale_x_continuous(breaks = seq(0,90,15))+
  scale_y_continuous(limits = c(0,8))+
  geom_hline(yintercept = 0,color='black',lwd=1)+
  geom_vline(xintercept = 90,color='black',lwd=1)+
  labs(y = "absolute activity", x = NULL) +
  # ggtitle('LV extra 1')+
  coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F, direction = 1)+
  theme_bw() +
  theme(text = element_text(size = 16),
        # plot.title = element_text(vjust = -15,hjust = 0.1,size=18,face = 'bold'),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(vjust = 4),
        axis.text.x = element_blank(),
        # axis.title.x = element_text(vjust = 9,hjust = 1.05,size=18,face='bold'),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 1.2),
        panel.grid.minor = element_line(linewidth = 1.2))+
  facet_wrap(~Pathway)
print(plt_pw_radial)
ggsave('figures/quarter_radial_pathway_act.png',
       plot = plt_pw_radial,
       height = 6,
       width = 9,
       units = 'in',
       dpi = 600)
ggsave('figures/quarter_radial_pathway_act.eps',
       plot = plt_pw_radial,
       device = cairo_ps,
       height = 6,
       width = 9,
       units = 'in',
       dpi = 600)
plt_pw_radial2 <- ggplot(df_plot  %>% gather('theta','activity',-Pathway) %>% 
                           mutate(theta=as.numeric(theta))%>% 
                           filter(Pathway %in% c('JAK-STAT','p53','NFkB')),
                         aes(x = theta, y = activity, color = Pathway,fill=Pathway, group = Pathway)) + 
  geom_line(linewidth = 1.5) +
  geom_area(alpha=0.2)+
  scale_x_continuous(breaks = seq(0,90,15))+
  scale_y_continuous(limits = c(0,8))+
  scale_colour_manual(values = c('#53B400','#00C094','#00BFC4'))+
  scale_fill_manual(values = c('#53B400','#00C094','#00BFC4'))+
  geom_hline(yintercept = 0,color='black',lwd=1)+
  geom_vline(xintercept = 90,color='black',lwd=1)+
  labs(y = "absolute activity", x = NULL) +
  coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F, direction = 1)+
  # ggtitle('LV extra 1')+
  theme_bw() +
  theme(text = element_text(size = 20),
        # plot.title = element_text(vjust = -15,hjust = 0.1,size=18,face = 'bold'),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(hjust = 0.66),
        axis.text.x = element_blank(),
        # axis.title.x = element_text(vjust = 9,hjust = 1.05,size=18,face='bold'),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        # legend.key.size = unit(0.3, "cm"),
        # legend.margin = margin(t = 0, r = 0, b = 0, l = -1.5, unit = "cm"),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 1.2),
        panel.grid.minor = element_line(linewidth = 1.2))+facet_wrap(~Pathway)
print(plt_pw_radial2)
# p1 <- ggRadar(df_plot, aes(group = Pathway), rescale = FALSE, use.label = TRUE,size = 1) +
#   scale_x_discrete(breaks = seq(0,90,15))+
#   geom_text(data=df_labels,aes(x=17,y=score,label = score),color='black',alpha=0.7,size=6)+
#   theme_light() +
#   labs(title = "Pathway activity along each direction of the extra basis space")+
#   theme(text = element_text(family = 'Arial',size=20),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_blank()) +
#   facet_wrap(~Pathway)
ggsave(paste0('figures/pathways_in_extra_basis_',tolower(target_dataset),'.png'),
       plot = plt_pw_radial2,
       width = 9,
       height = 6,
       units = 'in',
       dpi =600)
ggsave(paste0('figures/pathways_in_extra_basis_',tolower(target_dataset),'.eps'),
       plot = plt_pw_radial2,
       device = cairo_ps,
       width = 9,
       height = 6,
       units = 'in',
       dpi =600)
# p2 <- ggplot(df,aes(x=LV_opt_1,y=LV_opt_2,fill=score,colour=score))+
#   geom_segment(aes(x=0,y=0,xend =LV_opt_1,yend=LV_opt_2),
#                arrow = arrow(length = unit(0.03, "npc")),size=1)+
#   scale_fill_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-8,8))+
#   scale_color_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0,limits = c(-8,8))+
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   xlab('extra basis 1') + ylab('extra basis 2') +
#   labs(fill = 'activity',colour='activity')+
#   facet_wrap(~Pathway)+
#   theme_light(base_size = 24,base_family = 'Arial')+
#   labs(title = "Pathway activity along each direction of the extra basis space")+
#   theme(text = element_text(size=24,family='Arial'),
#         plot.title = element_text(hjust = 0.5),
#         strip.text = element_text(size=24),
#         legend.position = 'right')
# print(p2)
# ggsave(paste0('results/pc_loadings_scores_analysis/pathways_in_extra_basis_arrowplot_',tolower(target_dataset),'.png'),
#        plot = p2,
#        width = 12,
#        height = 12,
#        units = 'in',
#        dpi =600)
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
df_msig <- df_msig %>% mutate(Hallmark=str_replace_all(Hallmark,'_',' '))
df_msig <- df_msig %>% mutate(Hallmark = tolower(Hallmark)) %>% 
  mutate(Hallmark = paste0(toupper(substr(Hallmark, 1, 1)), tolower(substr(Hallmark, 2, nchar(Hallmark)))))
p1 <- (ggplot(df_msig %>% filter(LV=='V1') %>% arrange(NES) %>%
                filter(padj<=0.1),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV extra 1')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p1)
ggsave(paste0('figures/hallmark_',
              tolower(target_dataset),
              'on_optimal_lv1.png'),
       plot=p1,
       width=9,
       height=9,
       units = 'in',
       dpi = 600)
ggsave(paste0('figures/hallmark_',
              tolower(target_dataset),
              'on_optimal_lv1.eps'),
       device = cairo_ps,
       plot=p1,
       width=9,
       height=9,
       units = 'in',
       dpi = 600)
p2 <- (ggplot(df_msig %>% filter(LV=='V2') %>% arrange(NES)%>%
                filter(padj<=0.1),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') +ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV extra 2')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p2)
ggsave(paste0('figures/hallmark_',
              tolower(target_dataset),
              'on_optimal_lv2.png'),
       plot=p2,
       width=9,
       height=9,
       units = 'in',
       dpi = 600)

## repeat for TCs
colnames(Wm_combo) <- c('V1','V2')
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Wm_combo), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Wm_combo[inds,])
rownames(meas) <- entrez_ids
msig_tcs <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes_tcs <- as.data.frame(msig_tcs$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark')  #%>% gather('PC','NES',-Hallmark)
msig_nes_tcs <- msig_nes_tcs %>% gather('LV','NES',-Hallmark)
msig_pval_tcs <- as.data.frame(msig_tcs$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')#%>% gather('PC','padj',-Hallmark)
msig_pval_tcs <- msig_pval_tcs %>% gather('LV','padj',-Hallmark)
df_msig_tcs <- left_join(msig_nes_tcs,msig_pval_tcs)
df_msig_tcs <- df_msig_tcs %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))
df_msig_tcs <- df_msig_tcs %>% mutate(Hallmark=str_replace_all(Hallmark,'_',' '))
df_msig_tcs <- df_msig_tcs %>% mutate(Hallmark = tolower(Hallmark)) %>% 
  mutate(Hallmark = paste0(toupper(substr(Hallmark, 1, 1)), tolower(substr(Hallmark, 2, nchar(Hallmark)))))
p1 <- (ggplot(df_msig_tcs %>% filter(LV=='V1') %>% arrange(NES) %>%
                filter(padj<=0.1),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in TC 1')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p1)
p2 <- (ggplot(df_msig_tcs %>% filter(LV=='V2') %>% arrange(NES)%>%
                filter(padj<=0.1),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') +ylab('Hallmark')+
         ggtitle('Hallmarks enriched in TC 2')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p2)
p <- p1 + p2
ggsave(paste0('figures/hallmark_',
              tolower(target_dataset),
              'on_translatable_components.png'),
       plot=p,
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

ggsave(paste0('figures/optimal_LV_1_',
              tolower(target_dataset),
              '_perturbation_activity.png'),
       plot = extra_basis_inferred_perts$figure[[1]]+ theme(plot.title = element_blank()),
       width = 9,
       height = 9,
       dpi = 600)
ggsave(paste0('figures/optimal_LV_1_',
              tolower(target_dataset),
              '_perturbation_activity.eps'),
       plot = extra_basis_inferred_perts$figure[[1]]+ theme(plot.title = element_blank()),
       device = cairo_ps,
       width = 9,
       height = 9,
       dpi = 600)

ggsave(paste0('figures/optimal_LV_2_',
              tolower(target_dataset),
              '_perturbation_activity.png'),
       plot = extra_basis_inferred_perts$figure[[2]]+ theme(plot.title = element_blank()),
       width = 9,
       height = 9,
       dpi = 600)
ggsave(paste0('figures/optimal_LV_2_',
              tolower(target_dataset),
              '_perturbation_activity.eps'),
       plot = extra_basis_inferred_perts$figure[[2]]+ theme(plot.title = element_blank()),
       device = cairo_ps,
       width = 9,
       height = 9,
       dpi = 600)


