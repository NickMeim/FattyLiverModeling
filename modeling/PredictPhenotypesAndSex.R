library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(caret)
library(ropls)
source("../utils/plotting_functions.R")
source("functions_translation.R")
source("CrossValidationUtilFunctions.R")
source('vector_space_interpretation.R')

### Load data--------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver","Hoang")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% 
                  select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
Xh <- data_list[[ref_dataset]]$data_center %>% t()
sex_inferred <- apply(as.matrix(Xh[,c('RPS4Y1')]),2,sign) 
sex_inferred <- 1*(sex_inferred>0)
Yh <- cbind(Yh,sex_inferred)
colnames(Yh) <- c('NAS','fibrosis','sex')
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

#### Train PLSR model----------------
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

### Plot human PLSR space and projected truncated data
Zh_plsr <- plsr_model@scoreMN
mcor <- cor(cbind(Zh_plsr,Yh))
mcor[upper.tri(mcor,diag = TRUE)] <- 100
mcor <- reshape2::melt(mcor)
mcor <- mcor %>% filter(value != 100)
mcor <- mcor %>% mutate(keep = ifelse(Var1=='fibrosis' | Var2=='fibrosis' | Var1=='NAS' | Var2=='NAS' | Var1=='sex' | Var2=='sex',TRUE,FALSE)) %>%
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(keep = ifelse((Var1=='fibrosis' & Var2=='NAS') | (Var1=='fibrosis' & Var2=='sex') | (Var1=='NAS' & Var2=='fibrosis') | (Var1=='NAS' & Var2=='sex') | (Var1=='sex' & Var2=='fibrosis') | (Var1=='sex' & Var2=='NAS'),
                       FALSE,TRUE))%>%
  filter(keep==TRUE) %>% select(-keep) 
mcor <- mcor %>% group_by(Var2) %>% mutate(avg_abs_cor = mean(abs(value))) %>% select(c('LV'='Var2'),avg_abs_cor) %>%
  unique()
lvs <- mcor$LV[order(-mcor$avg_abs_cor)]
lvs <- as.character(lvs[1:2])
lv1_var <- round(100 * var(Zh_plsr[,c(lvs[1])])/sum(apply(Xm,2,var)),2) 
lv2_var <- round(100 * var(Zh_plsr[,c(lvs[2])])/sum(apply(Xm,2,var)),2) 
invivo_plsr <- as.data.frame(Zh_plsr[,c(lvs[1],lvs[2])])
invivo_plsr <- cbind(invivo_plsr,as.data.frame(Yh))
invivo_plsr <- invivo_plsr %>% gather('phenotype','Score',-all_of(lvs))
colnames(invivo_plsr)[1:2] <- c('V1','V2')
ggplot(invivo_plsr,aes(x=V1,y=V2,color=Score)) +
  geom_point(size=5)+
  scale_color_viridis_c()+
  xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv1_var,'%)')) + ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv2_var,'%)')) +
  theme_pubr(base_size = 30,base_family = 'Arial')+
  theme(text = element_text(size=30,family = 'Arial'),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(linewidth = 4),
        legend.position = 'right')+
  facet_wrap(~phenotype)
## visualize plsr for backprojected data
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
ggplot(invivo_plsr_backprj,aes(x=V1,y=V2,color=Score)) +
  # geom_point(color='#FD0714',size=5)+
  geom_point(size=5)+
  scale_color_viridis_c()+
  xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv1_var,'%)')) + ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv2_var,'%)')) +
  theme_pubr(base_size = 30,base_family = 'Arial')+
  theme(text = element_text(size=30,family = 'Arial'),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(linewidth = 4),
        legend.position = 'right')+
  facet_wrap(~phenotype)

invivo_plsr_all <- rbind(invivo_plsr %>% mutate(type = 'Original human data'),
                         invivo_plsr_backprj %>% mutate(type = 'Truncated human data'))
invivo_plsr_all <- invivo_plsr_all %>% mutate(phenotype = ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
invivo_plsr_all <- invivo_plsr_all %>% mutate(Phenotype = Score) %>% select(-Score)
ggplot(invivo_plsr_all,aes(x=V1,y=V2,color=Phenotype)) + #%>% select(-phenotype,-Score) %>% unique()
  # geom_point(color='#FD0714',size=5)+
  geom_point(size=5)+
  scale_x_continuous(limits = c(-50,60))+
  scale_y_continuous(limits = c(-40,30))+
  scale_color_viridis_c()+
  xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv1_var,'%)')) + ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv2_var,'%)')) +
  theme_pubr(base_size = 30,base_family = 'Arial')+
  theme(text = element_text(size=30,family = 'Arial'),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(linewidth = 4),
        strip.background = element_blank(),
        legend.position = 'bottom')+
  facet_wrap(vars(type,phenotype),scales = 'free')

# Find extra basis
phi <- Wh %*% Bh
Wm_opt <- analytical_solution_opt(y=Yh,
                                  W_invitro = Wm,
                                  phi = phi)
Wm_tot <- cbind(Wm, Wm_opt)

Zh <- as.data.frame(Xh %*% Wm_opt)
Zh$NAS <- Yh[,1]
Zh$fibrosis <- Yh[,2]
Zh$sex <- Yh[,3]
Zh <- Zh %>% rownames_to_column('sample')
thetas <- seq(0,360,5)
phis <- seq(0,180,10)
df_proj <- data.frame()
df <- data.frame()
for (theta in thetas){
  for(p in phis){
    u <- c(cos(theta * pi / 180)*sin(p * pi / 180),sin(theta * pi / 180)*cos(p * pi / 180),cos(p * pi / 180))
    u <- as.matrix(u)
    Wproj <- Wm_opt %*% u
    Wtmp <- cbind(Wm,Wproj)
    Yhat <- cbind(1, Xh %*% Wtmp %*% t(Wtmp) %*% Wh) %*% rbind(apply(Yh,2,mean),Bh)
    corr_nas <- cor(Yhat[,1],Zh$NAS)
    corr_fib <- cor(Yhat[,2],Zh$fibrosis)
    confmat <- caret::confusionMatrix(data = as.factor(1*(Yhat[,3]>=0.5)),reference = as.factor(Yh[,3]))
    df <- rbind(df,
                data.frame(theta = theta,phi = p,phenotype = 'NAS',performance=corr_nas),
                data.frame(theta = theta,phi = p,phenotype = 'fibrosis',performance=corr_fib),
                data.frame(theta = theta,phi = p,phenotype = 'sex',performance=confmat$byClass['F1']))
  }
}

scatter_box_plot <- function(df,legend_title,
                             font_size =20,font_family = 'Arial',
                             point_shape = 21,point_size=2.8,point_stroke=1.2,
                             x_axis='LV extra 1',y_axis='LV extra 2',
                             box_y_width = 0.2,jitter_y_width=0.1,
                             jitter_x_height = 0.2,
                             theme_use = 'minimal',
                             dims = c('V1','V2'),
                             plotting=TRUE){
  if (theme_use=='bw'){
    th <- ggplot2::theme_bw(base_size = font_size, base_family = font_family)
  }else{
    th <- ggplot2::theme_minimal(base_size = font_size, base_family = font_family)
  }
  df <- df %>% select(all_of(dims),pheno)
  colnames(df)[1:2] <- c('V1','V2')
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
combined_plot_nas <- scatter_box_plot(Zh %>% mutate(pheno=NAS),'NAS',
                                      x_axis='LV extra 1',y_axis='LV extra 2',
                                      dims = c('V1','V2'))
combined_plot_fibrosis <- scatter_box_plot(Zh %>% mutate(pheno=fibrosis),'Fibrosis stage',
                                           x_axis='LV extra 1',y_axis='LV extra 2',
                                           dims = c('V1','V2'))
combined_plot_sex <- scatter_box_plot(Zh %>% mutate(pheno=sex),'sex',
                                      x_axis='LV extra 2',y_axis='LV extra 3',
                                      dims = c('V2','V3'))


df_mu_radial <- df%>% filter(phenotype!='sex') %>%  mutate(z=cos(phi*pi/180))  %>% filter(phi<=90) %>%
  group_by(theta,z) %>% mutate(mu = mean(performance)) %>% ungroup() %>% 
  select(-performance,-phenotype) %>% unique() %>%
  mutate(phenotype = 'average') %>% select(theta,phi,z,phenotype,c('performance'='mu')) %>%
  mutate(x = performance*cos(theta), y = performance*sin(theta))
df_radial <- df %>% filter(phenotype!='sex') %>%  mutate(z=cos(phi*pi/180))  %>% filter(phi<=90) %>% 
  mutate(x = performance*cos(theta), y = performance*sin(theta)) %>%
  mutate(phenotype = ifelse(phenotype == "fibrosis", "Fibrosis stage",phenotype)) #%>% mutate(theta = pi*theta/180)
plt_cor_radial <-   ggplot(df_radial %>% filter(theta<=90), 
                           aes(x = theta, y = performance, color = phenotype, group = phenotype)) + 
  geom_line(linewidth = 1.5) +
  geom_line(data = df_mu_radial %>% filter(theta<=90),
            aes(x = theta, y = performance),
            color='black',lwd=1.5,linetype='dashed',
            inherit.aes = FALSE)+
  scale_x_continuous(breaks = seq(0,90,15))+
  ylim(c(NA,1))+
  geom_hline(yintercept = 0.2,color='black',lwd=1)+
  geom_vline(xintercept = 90,color='black',lwd=1)+
  coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F, direction = 1) +
  labs(y = "Performance", x = 'LV extra 2') +
  ggtitle('LV extra 1')+
  theme_bw() +
  theme(text = element_text(size = 20),
        plot.title = element_text(vjust = -13,hjust = -0.1,size=18,face = 'bold'),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(hjust = 0.66,vjust = 5),
        axis.text.x = element_blank(),
        axis.title.x = element_text(vjust = 9,hjust = 1.05,size=18,face='bold'),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -1.5, unit = "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 1.2),
        panel.grid.minor = element_line(linewidth = 1.2)) +
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(~phi)
print(plt_cor_radial)

ggplot(df_radial %>% filter(theta<=90),aes(x=theta,y=performance,color=phenotype))+
  geom_line()+
  facet_wrap(~phi)

### Interohgate loadings--------------------
rownames(Wm_opt) <- rownames(Wm)
extra_basis_pathway_activity <- pathway_activity_interpretation(Wm_opt,
                                                                Wm)
W <- Wm_opt[,2:3]
colnames(W) <- c('V1','V2')
extra_basis_pathway_activity_v3 <- pathway_activity_interpretation(W,
                                                                Wm,
                                                                plotting = FALSE,
                                                                lim = 10)
extra_basis_pathway_activity_v3 <- extra_basis_pathway_activity_v3$figure[[2]] + ggtitle('LV extra 3')
print(extra_basis_pathway_activity_v3)

## hallmarks
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
p3 <- (ggplot(df_msig %>% filter(LV=='V3') %>% arrange(NES) %>%
                filter(padj<=0.2),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         # ggtitle('Hallmarks enriched in LV extra 1')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p3)

### interogate study stats per gender
ggboxplot(as.data.frame(Yh) %>% gather('phenotype','score',-sex),
          x='sex',y='score',add='jitter')+
  stat_compare_means(comparisons = list(c('1','0')))+
  facet_wrap(~phenotype)
