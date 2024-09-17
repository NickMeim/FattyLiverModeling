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
# Save sex inferred
saveRDS(sex_inferred, file = paste0("results/", tolower(ref_dataset),"_sex_inferred.rds"))
# Get matrices of target dataset
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
# Save
saveRDS(file = paste0('results/Wh_',tolower(target_dataset),'.rds'), object = Wh)
saveRDS(file = paste0('results/PLSR_model_',tolower(target_dataset),'.rds'), object = plsr_model)

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
  Wproj <- Wm_opt %*% u
  Wtmp <- cbind(Wm,Wproj)
  Yhat <- cbind(1, Xh %*% Wtmp %*% t(Wtmp) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
  corr_nas <- cor(Yhat[,1],Zh$NAS)
  corr_fib <- cor(Yhat[,2],Zh$fibrosis)
  df <- rbind(df,
              data.frame(theta = theta,phenotype = 'NAS',corr=corr_nas),
              data.frame(theta = theta,phenotype = 'fibrosis',corr=corr_fib))
}
# Save
saveRDS(df, paste0("results/df_correlation_radial_", tolower(target_dataset), ".rds"))


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

### Find translatable LV of the in vitro system
### Run evolutionary algorithm
Wm_combo <- get_translatable_LV(Xh, Yh, Wh, Wm,
                              rbind(apply(Yh,2,mean),Bh),
                              find_extra = FALSE,
                              verbose = TRUE)
Wm_combo <- Wm_combo$Wm_new
rownames(Wm_combo) <- rownames(Wm)

# Print percentage of variance captured in TC
print(paste0("Total variance captured in TCs: ", round(sum(colVars(Xm %*% Wm_combo))/sum(colVars(Xm))*100,2), " (%)"))

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

### Analyze gene loadings----------------------------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
plot_extra_gene_loadings_lv1 <- plot_gene_loadings(loadings = Wm_opt,
                                               selection='V1',
                                               y_lab = 'weight in extra LV1',
                                               top=20)


plot_extra_gene_loadings_lv2 <- plot_gene_loadings(Wm_opt,
                                               selection='V2',
                                               y_lab = 'weight in extra LV2')


plot_translatable_gene_loadings_lv1 <- plot_gene_loadings(Wm_combo,
                                                      colnames(Wm_combo)[1],
                                                      'translatable LV1')


plot_translatable_gene_loadings_lv2 <- plot_gene_loadings(Wm_combo,
                                                      colnames(Wm_combo)[2],
                                                      'translatable LV2')



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


### Analyze loadings at the Pathway activity level--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
extra_basis_pathway_activity <- pathway_activity_interpretation(Wm_opt,
                                                                Wm)


colnames(Wm_combo) <- c('V1','V2')
translatable_components_progenies <- pathway_activity_interpretation(Wm_combo,
                                                                       Wm,
                                                                     lim=15)
p1 <- translatable_components_progenies$figure[[1]]
p1 <- p1 + ggtitle('Translatable Component 1')
p2 <- translatable_components_progenies$figure[[2]] + ggtitle('Translatable Component 2')
print(p1+p2)

### Visualize pathways in various directions of the new space-----------------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
thetas <- seq(0,360,5)
df <- data.frame()
for (theta in thetas){
  u <- c(cos(theta * pi / 180),sin(theta * pi / 180))
  u <- as.matrix(u)
  W_proj <- Wm_opt %*% u
  extra_basis_pathway_activity <- pathway_activity_interpretation(W_proj,
                                                                  Wm)
  
  
  tmp <- extra_basis_pathway_activity %>% filter(condition == 'V1') %>% ungroup() %>% 
    select(Pathway,score) %>% mutate(LV_opt_1 = u[1,1],LV_opt_2 =u[2,1]) %>%
    mutate(theta = theta)
  df <- rbind(df,
              tmp)
}
df_plot <- df %>% select(Pathway,theta,score) %>% mutate(score = abs(score)) %>% filter(theta<=90) %>%
  spread('theta','score')
# Save
saveRDS(df_plot, paste0("results/df_pwy_radial_", tolower(target_dataset),"_.rds"))

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
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(hjust = 0.66),
        axis.text.x = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 1.2),
        panel.grid.minor = element_line(linewidth = 1.2))+facet_wrap(~Pathway)
print(plt_pw_radial2)


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
# Save
saveRDS(df_msig, paste0("results/hallmark_enrichment_", tolower(target_dataset),"_LVopt.rds" ))


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

# Save
saveRDS(df_msig_tcs, paste0("results/hallmark_enrichment_", tolower(target_dataset),"_TC.rds" ))


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
saveRDS(extra_basis_inferred_perts, paste0("results/extra_basis_", tolower(target_dataset),"_inferred_perts.rds"))

extra_TC_inferred_perts <- perturnation_activity_inference(Wm_combo,metadata_human,dorotheaData,resp_net)