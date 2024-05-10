# load in packages
library(tidyverse)
library(DESeq2)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(ggradar)
library(pls)
library(pheatmap)
library(patchwork)
library(limma)
library(edgeR)
# library(dorothea)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(hgu133a.db)
library(rstatix)
library(fgsea)
library(topGO)
library(GO.db)
# library(progeny)
library(OmnipathR)
library(EGSEAdata)
# library(nichenetr)
source("../utils/plotting_functions.R")
source("functions_translation_jose.R")
source("CrossValidationUtilFunctions.R")
source('enrichment_calculations.R')

#### Load pre-processed data----------------------
# data <- readRDS("../data/preprocessed_NAFLD.rds")
# data_A <- data$data_A
# data_B <- data$data_B
# Y_A <- data$Y_A

# ## normalize
# ## I do not think I need any more common genes to interogate
# # Intersect datasets - keep common genes
# genes_common <- intersect(rownames(data_A), rownames(data_B))
# # Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
# data_A <- apply(data_A[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
# data_B <- apply(data_B[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
# # Remove features absent from at least 10% in each samples
# keep_gene <- (rowSums(data_A >= 1) >= 0.1*ncol(data_A)) &
#   (rowSums(data_B >= 1) >= 0.1*ncol(data_B))
# # Run PCA on each dataset. Not scaling PCA
# # Log and center dataset A and run PCA
# X_A <- log2(1 + data_A[keep_gene,])
# X_A <- t(X_A - rowMeans(X_A))
# # Log and center dataset B and run PCA
# Xm <- log2(1 + data_B[keep_gene,])
# Xm <- t(Xm - rowMeans(Xm))
# # Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
# # data_B <- apply(data_B, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
# # Remove features absent from at least 10% in each samples
# # keep_gene <- rowSums(data_B >= 1) >= 0.1*ncol(data_B)
# # Run PCA on each dataset. Not scaling PCA
# # Log and center dataset B and run PCA
# # keep_gene <- rownames(Woptim)[which(rownames(Woptim) %in% rownames(data_B))]
# # Xm <- log2(1 + data_B[keep_gene,])
# # Xm <- t(Xm - rowMeans(Xm))

dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver","Govaere")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
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

Woptim <- readRDS('../results/Wm_extra.rds')
Wtot <- readRDS('../results/Wm_total.rds')
Wm_combo<- readRDS('../results/Wm_combo.rds')
rownames(Woptim) <- rownames(Wm)
rownames(Wtot) <- rownames(Wm)
Zh <- as.data.frame(Xh %*% Woptim)
Zh$NAS <- Yh[,1]
Zh$fibrosis <- Yh[,2]
Zh <- Zh %>% rownames_to_column('sample')
w1 <- c(seq(-1,1,0.01))
w2 <- c(seq(-1,1,0.01))
r <- NULL
k <- 1
df <- data.frame()
# for (a in w1){
#   for (b in w2){
#     r[k] <- b/a
#     tmp <- Zh %>% mutate(u = a * V1 + b * V2)
#     corr_nas <- cor(tmp$u,tmp$NAS)
#     corr_fib <- cor(tmp$u,tmp$fibrosis)
#     df <- rbind(df,
#                 data.frame(w1=a,w2=b,phenotype = 'NAS',corr=corr_nas),
#                 data.frame(w1=a,w2=b,phenotype = 'fibrosis',corr=corr_fib))
#   }
# }
# thetas <- c(0,5,10,20,30,45,60,75,90,120,seq(140,240,20),270,seq(300,360,20))
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
# ggplot(df,aes(x=w1,y=w2,fill=corr)) +
#   scale_fill_gradient2(low = 'blue',high = 'red',mid='white',midpoint = 0, limits = c(-1, 1))+
#   geom_tile()+
#   facet_wrap(~phenotype) +
#   theme_pubr(base_size = 20,base_family = 'Arial')+
#   theme(text = element_text(size=20,family='Arial'),
#         legend.position = 'right')
# ggsave('../results/pc_loadings_scores_analysis/optimal_directions_combinations_corr.png',
#        height = 9,
#        width = 12,
#        units = 'in',
#        dpi=600)
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

# p2 <- ggplot(df %>% spread('phenotype','corr'),aes(x=NAS,y=fibrosis,fill=theta * pi / 180)) +
#   geom_point(shape=22)+
#   scale_fill_gradient(low = 'red',high = 'white',limits = c(0, 2*pi))+
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   theme_minimal(base_size = 20,base_family = 'Arial')+
#   theme(text = element_text(size=20,family='Arial'),
#         legend.position = 'right')

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
ggsave('../results/pc_loadings_scores_analysis/optimal_space_phenotypes_correlations_directions.png',
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)

(ggplot(Zh ,aes(x=V1,y=V2,colour=NAS))+
  geom_point()+
  scale_color_viridis_c()+
  theme_minimal(base_size=20,base_family = 'Arial')+
  theme(text= element_text(size=20,family = 'Arial'),
        legend.position = 'right')) +
  (ggplot(Zh ,aes(x=V1,y=V2,colour=fibrosis))+
                                        geom_point()+
                                        scale_color_viridis_c()+
                                        theme_minimal(base_size=20,base_family = 'Arial')+
                                        theme(text= element_text(size=20,family = 'Arial'),
                                              legend.position = 'right'))
ggsave('../results/projected_govaere_samples_on_extra_basis.png',
       height = 12,
       width = 16,
       units = 'in',
       dpi=600)

## Check out 
Woptim <- as.data.frame(Woptim) %>% rownames_to_column('gene')
Woptim$significant <- ''
Woptim$significant[order(-abs(Woptim$V1))[1:20]] <- Woptim$gene[order(-abs(Woptim$V1))[1:20]]
Woptim <- Woptim[order(Woptim$V1),]
Woptim$gene <- factor(Woptim$gene,levels = Woptim$gene)


### Perform PCA and get loadings----------------
# PCA_alldata <- prcomp(Xm, scale. = F, center = T)
loadings <- Wm
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,PC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$PC12))[1:10]] <- loadings$gene[order(-abs(loadings$PC12))[1:10]]
loadings <- loadings[order(loadings$PC12),]
loadings$gene <- factor(loadings$gene,levels = loadings$gene)
### Make barplot to look at top genes----------------
ggplot(Woptim,aes(x=gene,y=w_optimal,color = significant)) + geom_point() +
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
  xlab('genes') + ylab('optimal direction loadings')+
  scale_x_discrete(expand = c(0.1, 0.1))+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(text = element_text(family = 'Arial',size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none')
ggsave('../results/pc_loadings_scores_analysis/gene_woptimal_loadings.png',
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

### Infer statistical significance for gene loadings of PC12 and PC8----------------------------------
gene_loadings <- PCA_alldata$rotation
gene_loadings <- as.data.frame(gene_loadings) %>% rownames_to_column('gene')
null_loadings <- gene_loadings %>% gather('PC','null_act',-gene) %>% select(-PC)
gene_loadings$significant <- ''
gene_loadings_PC8 <- left_join(gene_loadings %>% select(gene,PC8,significant),
                               null_loadings) %>%
  group_by(gene) %>% mutate(p.value = ifelse(PC8>=0,
                                             sum(null_act>=PC8)/(ncol(gene_loadings)-2),
                                             sum(null_act<=PC8)/(ncol(gene_loadings)-2))) %>%
  ungroup()
adjustement_tmp <- gene_loadings_PC8 %>% select(gene,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'fdr')
adjustement_tmp <- adjustement_tmp %>% select(gene,p.adj)
gene_loadings_PC8 <- left_join(gene_loadings_PC8,adjustement_tmp)%>% select(-null_act) %>% unique()
gene_loadings_PC12 <- left_join(gene_loadings %>% select(gene,PC12,significant),
                                null_loadings) %>%
  group_by(gene) %>% mutate(p.value = ifelse(PC12>=0,
                                             sum(null_act>=PC12)/(ncol(gene_loadings)-2),
                                             sum(null_act<=PC12)/(ncol(gene_loadings)-2))) %>%
  ungroup()
adjustement_tmp <- gene_loadings_PC12 %>% select(gene,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'fdr')
adjustement_tmp <- adjustement_tmp %>% select(gene,p.adj)
gene_loadings_PC12 <- left_join(gene_loadings_PC12,adjustement_tmp) %>% select(-null_act) %>% unique()
gene_loadings_PC8$significant[order(-abs(gene_loadings_PC8$PC8))[1:10]] <- gene_loadings_PC8$gene[order(-abs(gene_loadings_PC8$PC8))[1:10]]
gene_loadings_PC8 <- gene_loadings_PC8[order(gene_loadings_PC8$PC8),]
gene_loadings_PC8$gene <- factor(gene_loadings_PC8$gene,levels = gene_loadings_PC8$gene)
gene_loadings_PC12$significant[order(-abs(gene_loadings_PC12$PC12))[1:10]] <- gene_loadings_PC12$gene[order(-abs(gene_loadings_PC12$PC12))[1:10]]
gene_loadings_PC12 <- gene_loadings_PC12[order(gene_loadings_PC12$PC12),]
gene_loadings_PC12$gene <- factor(gene_loadings_PC12$gene,levels = gene_loadings_PC12$gene)
# change significant
# gene_loadings_PC12 <- gene_loadings_PC12 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                        ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
# gene_loadings_PC8 <- gene_loadings_PC8 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                      ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
gene_loadings_PC12 <- gene_loadings_PC12 %>% mutate(statistical=ifelse(p.value<=0.01,'p-value<=0.01',
                                                                       ifelse(p.value<=0.05,'p-value<=0.05','p-value>0.05')))
gene_loadings_PC8 <- gene_loadings_PC8 %>% mutate(statistical=ifelse(p.value<=0.01,'p-value<=0.01',
                                                                     ifelse(p.value<=0.05,'p-value<=0.05','p-value>0.05')))

# Combine into one data frame
gene_loadings_plot_frame <- rbind(gene_loadings_PC8 %>% dplyr::rename(value = PC8) %>% mutate(PC='PC8'),
                                  gene_loadings_PC12 %>% dplyr::rename(value = PC12) %>% mutate(PC='PC12'))
gene_loadings_plot_frame <- gene_loadings_plot_frame %>% mutate(PC=paste0('used ',PC,' gene loadings'))

p <- ggplot(gene_loadings_plot_frame %>% dplyr::rename(`statistical threshold`=statistical),
            aes(x=value,y=-log10(p.value),color = `statistical threshold`)) + geom_point() +
  scale_color_manual(values = c('red','#fa8e8e','black'))+
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40,point.padding = 0.5)+
  xlab('gene loading') + ylab(expression(-log[10]('p-value'))) +
  geom_hline(yintercept=-log10(0.01), linetype="dashed",color = "#525252", size=1) +
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'top')+
  facet_wrap(~PC,scales = 'free_x')+
  guides(color = guide_legend(
    override.aes = list(
      linetype = NA,
      size = 3
    )
  ))
print(p)
ggsave('../results/pc_loadings_scores_analysis/gene_loadings_volcano.png',
       plot = p,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

p2 <- (ggplot(gene_loadings_PC8,aes(x=gene,y=PC8,color = statistical)) + geom_point() +
         # scale_color_gradient(high = 'red',low='white')+
         scale_color_manual(values = c('red','#fa8e8e','black'))+
         geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
         xlab('TFs') + ylab('PC8 loadings')+
         scale_x_discrete(expand = c(0.1, 0.1))+
         theme_pubr(base_family = 'Arial',base_size = 20)+
         theme(text = element_text(family = 'Arial',size=20),
               axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               legend.position = 'none')) +
  (ggplot(gene_loadings_PC12,aes(x=gene,y=PC12,color = statistical)) + geom_point() +
     scale_color_manual(values =c('red','#fa8e8e','black'))+
     geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
     xlab('TFs') + ylab('PC12 loadings')+
     scale_x_discrete(expand = c(0.1, 0.1))+
     theme_pubr(base_family = 'Arial',base_size = 20)+
     theme(text = element_text(family = 'Arial',size=20),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = 'none'))
print(p2)
ggsave('../results/pc_loadings_scores_analysis/gene_loadings_points_withstats.png',
       plot = p2,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)
### Infer Pathway activity with progeny----------------------------------------------------------------------------
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
Woptim <- readRDS('../results/Wm_extra.rds')
rownames(Woptim) <- rownames(Wm)
# library(doFuture)
# # parallel: set number of workers
# cores <- 16
# registerDoFuture()
# plan(multisession,workers = cores)
# iters <- 1000
# null_activity <- foreach(i = seq(iters)) %dopar% {
#   shuffled_genes <-  Wm_combo_selected[sample(nrow(Wm_combo_selected)),]
#   rownames(shuffled_genes) <- rownames(Wm_combo_selected)
#   tmp_paths <- decoupleR::run_viper(shuffled_genes, net_prog,minsize = 1,verbose = FALSE)
#   tmp_paths <- tmp_paths %>% select(source,condition,score) %>% spread('condition','score')
#   tmp_paths <- tmp_paths %>% select(source,NAS,fibrosis)
#   tmp_paths
# }
# null_activity <- do.call(rbind,null_activity)
# saveRDS(null_activity,'../results/pc_loadings_scores_analysis/null_activity_translatable_lvs.rds')
# translatable_progenies <- decoupleR::run_viper(Wm_combo_selected, net_prog,minsize = 1,verbose = FALSE)
# all_null_distribution_nas <- null_activity[,1]
# all_null_distribution_fibrosis <- null_activity[,2]
# null_activity <- null_activity %>% gather('condition','random_val',-source)
# translatable_progenies <- left_join(translatable_progenies,null_activity)
# translatable_progenies <- translatable_progenies %>% select(c('Pathway'='source'),condition,score,random_val)

# df_translatable_pathways_progeny <- translatable_progenies %>% group_by(condition,Pathway) %>%
#   mutate(bool=abs(random_val)>=abs(score)) %>% mutate(p_value_null = sum(bool)) %>%
#   ungroup() %>% select(-bool) %>% mutate(p_value_null=p_value_null/iters)
# df_pathways_progeny_weights_all  <- df_pathways_progeny_weights_all %>% 
#   mutate(p.adj = p_value_null*length(unique(df_pathways_progeny_weights_all$Pathway)))
# df_pathways_progeny_weights_all <- df_pathways_progeny_weights_all %>% group_by(Pathway) %>%
#   mutate(sig1 = sum(p.adj<=0.05)) %>% ungroup() %>%
#   group_by(direction) %>% mutate(sig2 = sum(p.adj<=0.05)) %>% ungroup()
extra_basis_paths <- decoupleR::run_viper(Woptim, net_prog,minsize = 1,verbose = TRUE) %>% select(-statistic)
PC_space_paths <- decoupleR::run_viper(Wm, net_prog,minsize = 1,verbose = TRUE)
PC_space_paths <- PC_space_paths %>% select(source,c('pc_activity'='score'))
PC_space_paths <- PC_space_paths$pc_activity
extra_basis_paths <- extra_basis_paths %>% select(-p_value)
# extra_basis_paths <- left_join(extra_basis_paths %>% select(-p_value),PC_space_paths)
colnames(extra_basis_paths)[1] <- 'Pathway'

# extra_basis_paths <- extra_basis_paths %>% group_by(condition,Pathway) %>%
#   mutate(bool=abs(PC_space_paths)>=abs(score)) %>% mutate(p_value_null = sum(bool)) %>%
#   ungroup() %>% select(-bool) %>% mutate(p_value_null=p_value_null/length(PC_space_paths))
extra_basis_paths <- extra_basis_paths %>% group_by(condition,Pathway) %>%
  mutate(p_value=sum(abs(PC_space_paths)>=abs(score))/length(PC_space_paths))
extra_basis_paths  <- extra_basis_paths %>%
  mutate(p.adj = p_value*length(unique(extra_basis_paths$Pathway))) %>%
  mutate(p.adj = ifelse(p.adj>1,1,p.adj))
# extra_basis_paths <- extra_basis_paths %>% group_by(Pathway) %>%
#   mutate(sig1 = sum(p.adj<0.1)) %>% ungroup() %>%
#   group_by(condition) %>% mutate(sig2 = sum(p.adj<0.1)) %>% ungroup()


# %>% 
#   mutate(Pathway=factor(Pathway,levels=extra_basis_paths$source[order(extra_basis_paths$score)]))

(ggplot(extra_basis_paths %>% select(c('activity'='score'),Pathway,p_value,condition) %>% 
         filter (condition=='V1'),
       aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity',color='black') +
  scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
  scale_x_continuous(n.breaks = 8,limits = c(-8,8))+
  ggtitle('Extra LV1')+
  ylab('Pathway')+
  geom_text(aes(label = ifelse(p_value <= 0.0001, "****",
                               ifelse(p_value <= 0.001,"***",
                                      ifelse(p_value<=0.01,"**",
                                             ifelse(p_value<=0.05,'*',
                                                    ifelse(p_value<=0.1,'\u2219',
                                                           'ns'))))),
                x = ifelse(activity < 0, activity - 0.2, activity + 0.2)),
            size = 6,
            color = 'black',
            angle=90) +                                   
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size = 20,family = 'Arial'),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))) +
  (ggplot(extra_basis_paths %>% select(c('activity'='score'),Pathway,p_value,condition) %>% 
            filter (condition=='V2'),
          aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity',color='black') +
     ggtitle('Extra LV2')+
     ylab('Pathway')+
     scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
     scale_x_continuous(n.breaks = 8,limits = c(-8,8))+
     geom_text(aes(label = ifelse(p_value <= 0.0001, "****",
                                  ifelse(p_value <= 0.001,"***",
                                         ifelse(p_value<=0.01,"**",
                                                ifelse(p_value<=0.05,'*',
                                                       ifelse(p_value<=0.1,'\u2219',
                                                              'ns'))))),
                   x = ifelse(activity < 0, activity - 0.2, activity + 0.2)),
               size = 6,
               color = 'black',
               angle=90) +                                   
     theme_pubr(base_size = 20,base_family = 'Arial')+
     theme(text = element_text(size = 20,family = 'Arial'),
           legend.position = 'right',
           plot.title = element_text(hjust = 0.5),
           axis.title.y = element_blank())) 
ggsave('../results/pc_loadings_scores_analysis/progenies_only_optimal_loadings_barplot.png',
       width = 12,
       height = 9,
       dpi = 600)


#### Run CARNIVAL--------------------------
library(CARNIVAL)
gene_loadings <- PCA_alldata$rotation[,paste0('PC',c(1,2,3,8,12))]
minNrOfGenes  <-  5
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
settings = list(verbose = TRUE, minsize = minNrOfGenes)
## read optimal W
Woptim <- readRDS('../results/Wm_opt_2.rds')
TF_activities = run_viper(cbind(Woptim,Woptim), dorotheaData, options =  settings)
TF_activities <- as.matrix(TF_activities[,1])
tfs <- rownames(TF_activities)
## load omnipath network
interactions <- import_omnipath_interactions()
interactions <- interactions  %>% mutate(kegg=grepl(pattern="KEGG",x=sources))%>% 
  filter(n_resources>=3 | kegg==T)  %>% 
  filter(n_references>=3 | kegg==T)  %>% 
  dplyr::select(-kegg) %>%
  dplyr::select(c('source'='source_genesymbol'),
                c('target'='target_genesymbol'),
                is_inhibition,is_stimulation) %>% unique()

interactions <- interactions %>% filter(!(is_inhibition==0 & is_stimulation==0)) %>% unique()
interactions <- interactions %>% mutate(interaction=ifelse(is_stimulation!=0,1,-1)) %>%
  dplyr::select(source,interaction,target) %>% unique()
### Run carnival
log_con <- file("../results/pc_loadings_scores_analysis/log_cavinval.txt", open="a")
for (j in 1:ncol(TF_activities)){
  
  cat(paste0('Iteration ',j,'/',ncol(TF_activities)), file = log_con, sep="\n")
  tf_activities <- TF_activities[which(rownames(TF_activities) %in% c(interactions$source,
                                                                      interactions$target)),
                                 1]
  tf_activities <- tf_activities[order(-abs(tf_activities))[1:20]]
  # Run carnival
  # YOU MUST FIND WHERE CPLEX OR CUROBI IS INSTALLED IN YOUR OWN COMPUTER
  # GurPath <- 'C:/gurobi1100/win64/bin/gurobi_cl.exe'
  # CplexPath <- 'C:/Program Files/IBM/ILOG/CPLEX_Studio_Community2211/cplex/bin/x64_win64/cplex.exe'
  CbcPath <- 'C:/Users/nmeim/Documents/CBC/bin/cbc.exe'
  carnivalOptions <- list("solver"="cbc",
                          "betaWeight"=0.2,
                          "solverPath"=CbcPath,
                          "timelimit"=1800,
                          "poolrelGap"=1e-4,
                          "lpFilename"="",
                          "outputFolder"="",
                          "cleanTmpFiles"=TRUE,
                          "keepLPFiles"=TRUE,
                          "poolCap"=100)
  # carnivalOptions <- defaultCplexCarnivalOptions()
  # carnivalOptions$timelimit <- 1200
  # carnivalOptions$solverPath <- CplexPath
  # Output dir
  Result_dir <- paste0("../results/pc_loadings_scores_analysis/CARNIVAL/Woptim")
  dir.create(Result_dir, showWarnings = FALSE)
  carnivalOptions$outputFolder <- Result_dir
  
  inverseCarnivalResults <- runInverseCarnival( measurements = tf_activities, 
                                                priorKnowledgeNetwork = interactions, 
                                                carnivalOptions = carnivalOptions)
  
  # Save interaction networks
  nets <- inverseCarnivalResults$sifAll
  nodes <- inverseCarnivalResults$attributesAll
  for (i in 1:length(nets)){
    t <- nets[[i]]
    t <- as.data.frame(t)
    t$Node1 <- as.character(t$Node1)
    t$Node2 <- as.character(t$Node2)
    t$Sign <- as.character(t$Sign)
    write_tsv(t,paste0(Result_dir,'/','interactions_1_model',i,'.tsv'))
    t <- nodes[[i]]
    t <- as.data.frame(t)
    t$Nodes <- as.character(t$Nodes)
    t$Activity <- as.numeric(t$Activity)
    write_delim(t,paste0(Result_dir,'/','nodesActivity_1_model',i,'.txt'),delim = '\t')
  }
  t <- as.data.frame(inverseCarnivalResults$weightedSIF)
  t$Node1 <- as.character(t$Node1)
  t$Node2 <- as.character(t$Node2)
  t$Sign <- as.character(t$Sign)
  t$Weight <- as.numeric(t$Weight)
  write_delim(t,paste0(Result_dir,'/','weightedModel_1.txt'),delim = '\t')
  t <- as.data.frame(inverseCarnivalResults[["nodesAttributes"]])
  write_delim(t,paste0(Result_dir,'/','nodesAttributes_1.txt'),delim = '\t')
  message(paste0('Finished ',colnames(TF_activities)[j]))
}
close(log_con)

### Perform Fisher test between KEGG pathways and nodes in the net
Result_dir <- paste0("../results/pc_loadings_scores_analysis/CARNIVAL/PC12")
nodes <- read.delim(paste0(Result_dir,'/','nodesAttributes_1.txt'))
nodes <- nodes %>% filter(AvgAct!=0)
net <- read.delim(paste0(Result_dir,'/','weightedModel_1.txt'))
net <- net %>% filter(Node1 %in% nodes$Node) %>% filter(Node2 %in% nodes$Node) %>% unique()
egsea.data(species = "human",returnInfo = TRUE)
pathways <- kegg.pathways$human
pathways <- pathways$kg.sets
entrez_ids <- mapIds(org.Hs.eg.db, keys = nodes$Node, column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
all_nodes <- unique(interactions$source,interactions$target)
all_entrez_ids <- mapIds(org.Hs.eg.db, keys = all_nodes, column = "ENTREZID", keytype = "SYMBOL")
all_entrez_ids <- unname(all_entrez_ids)
inds <- which(!is.na(all_entrez_ids))
all_entrez_ids <- all_entrez_ids[inds]

## Perform Fisher Exact Test for each KEGG pathway
fisherExactTest <- function(Set,selected_nodes,all_possible_nodes){
  a  <- sum(selected_nodes %in% Set)
  b <- length(selected_nodes) - a
  c <- length(Set) - a
  d <- length(all_possible_nodes) - a -b-c
  contigencyMat <- data.frame('selected'=c(a,b),
                              'not_selected'= c(c,d))
  rownames(contigencyMat) <- c('in_set','not_in_set')
  return(fisher.test(contigencyMat)$p.value)
}
# universe_nodes <- unique(c(unique(unname(unlist(pathways))),all_entrez_ids))
fisher_results <- lapply(pathways,fisherExactTest,entrez_ids,all_entrez_ids)
fisher_results <- as.data.frame(unlist(fisher_results))
colnames(fisher_results) <- 'p.value'
fisher_results <- fisher_results %>% rstatix::adjust_pvalue(method = 'bonferroni')
fisher_results <- fisher_results %>% rownames_to_column('pathway')
my_data_frame <- bind_rows(lapply(names(pathways), function(name) {
  data.frame(names = name, values = list(pathways[[name]]), stringsAsFactors = FALSE)
})) %>% gather('key','node',-names) %>% filter(!is.na(key))
my_data_frame <- distinct(my_data_frame %>% select(c('pathway'='names'),node))
fisher_results <- left_join(fisher_results,my_data_frame)
fisher_results <- fisher_results %>% filter(p.value.adj<0.05)
fisher_results <- fisher_results %>% filter(!is.na(node))
fishered_nodes <- mapIds(org.Hs.eg.db, keys = fisher_results$node, column = "SYMBOL", keytype = "ENTREZID")
fisher_results$node <- fishered_nodes
fisher_results <- fisher_results %>% filter(node %in% nodes$Node)
print(unique(fisher_results$pathway))
fisher_results <- fisher_results %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))
fisher_results <- fisher_results %>% group_by(pathway) %>% mutate(num_nodes = n_distinct(node)) %>% ungroup()
fisher_results <- fisher_results %>% filter(!grepl('cancer',pathway))
fisher_results <- fisher_results %>% filter(!grepl('infection',pathway))
fisher_results <- fisher_results %>% filter(grepl('signaling',pathway) | grepl('pathway',pathway))
fisher_results <- fisher_results %>% filter(num_nodes>3) %>% select(pathway,node) %>% unique()
colnames(fisher_results) <- c('pathway','Node')
nodes <- left_join(nodes,fisher_results)
write_delim(nodes,paste0(Result_dir,'/','edited_node_attributes.txt'),delim = '\t')


### Infer TF activity with Dorothea from loadings of PC12 and PC8--------------------------------------------------
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
Woptim <- readRDS('../results/Wm_extra.rds')
rownames(Woptim) <- rownames(Wm)
extra_basis_tfs <- decoupleR::run_viper(Woptim, dorotheaData,minsize = 5,verbose = TRUE) %>% select(-statistic)
PC_space_paths <- decoupleR::run_viper(Wm, dorotheaData,minsize = 5,verbose = TRUE)
PC_space_paths <- PC_space_paths %>% select(source,c('pc_activity'='score'))
PC_space_paths <- PC_space_paths$pc_activity
extra_basis_tfs <- extra_basis_tfs %>% select(-p_value)
# extra_basis_tfs <- left_join(extra_basis_tfs %>% select(-p_value),PC_space_paths)
colnames(extra_basis_tfs)[1] <- 'TF'

# extra_basis_tfs <- extra_basis_tfs %>% group_by(condition,TF) %>%
#   mutate(bool=abs(PC_space_paths)>=abs(score)) %>% mutate(p_value_null = sum(bool)) %>%
#   ungroup() %>% select(-bool) %>% mutate(p_value_null=p_value_null/length(PC_space_paths))
extra_basis_tfs <- extra_basis_tfs %>% group_by(condition,TF) %>%
  mutate(p_value=sum(abs(PC_space_paths)>=abs(score))/length(PC_space_paths))
extra_basis_tfs  <- extra_basis_tfs %>%
  mutate(p.adj = p_value*length(unique(extra_basis_tfs$TF))) %>%
  mutate(p.adj = ifelse(p.adj>1,1,p.adj))
# extra_basis_tfs <- extra_basis_tfs %>% group_by(TF) %>%
#   mutate(sig1 = sum(p.adj<0.1)) %>% ungroup() %>%
#   group_by(condition) %>% mutate(sig2 = sum(p.adj<0.1)) %>% ungroup()
extra_basis_tfs$significant <- NA  
extra_basis_tfs <- extra_basis_tfs %>% mutate(significant = ifelse(p_value<0.1,
                                                               TF,
                                                               significant))
(ggplot(extra_basis_tfs %>% select(c('activity'='score'),TF,p_value,condition,significant) %>% 
          filter (condition=='V1'),
        aes(x=reorder(TF,activity),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
    geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
    scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.01,1))+
    scale_y_continuous(n.breaks = 6,limits = c(-6,6))+
    ggtitle('Extra LV1')+
    xlab('TF')+
    theme_pubr(base_size = 20,base_family = 'Arial')+
    theme(text = element_text(size = 20,family = 'Arial'),
          legend.position = 'right',
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank())) +
  (ggplot(extra_basis_tfs %>% select(c('activity'='score'),TF,p_value,condition,significant) %>% 
            filter (condition=='V2'),
          aes(x=reorder(TF,activity),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
     geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
     scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.01,1))+
     scale_y_continuous(n.breaks = 6,limits = c(-6,6))+
     ggtitle('Extra LV2')+
     xlab('TF')+
     theme_pubr(base_size = 20,base_family = 'Arial')+
     theme(text = element_text(size = 20,family = 'Arial'),
           legend.position = 'right',
           plot.title = element_text(hjust = 0.5),
           axis.title.y = element_blank(),
           axis.text.x = element_blank())) 
ggsave('../results/pc_loadings_scores_analysis/TF_only_optimal_loadings_barplot.png',
       width = 12,
       height = 9,
       dpi = 600)

### Infer protein activity-----------------------------------------
interactions <- import_omnipath_interactions()
interactions <- interactions %>% mutate(kegg=grepl(pattern="KEGG",x=sources)) %>%
  mutate(signoir=grepl(pattern="SIGNOR",x=sources)) %>% filter(kegg==T | signoir==T)
interactions <- interactions  %>% filter(n_resources>1)  %>% #filter(n_references>=2)%>% 
  dplyr::select(c('source'='source_genesymbol'),
                c('target'='target_genesymbol'),
                is_inhibition,is_stimulation) %>% unique()
interactions <- interactions %>% mutate(interaction=ifelse(is_stimulation==1,
                                                           ifelse(is_inhibition!=1,1,0),
                                                           ifelse(is_inhibition!=0,-1,0))) %>%
  dplyr::select(source,interaction,target) %>% unique()
# interactions <- interactions %>% filter(interaction!=0)
minNrOfGenes  <-  5
settings = list(verbose = TRUE, minsize = minNrOfGenes)
# protein_activity = run_viper(PCA_alldata$rotation, 
#                              rbind(interactions %>% mutate(confidence = 'A') %>% 
#                                      select(c('tf'='source'),confidence,target,c('mor'='interaction')) %>%
#                                      filter(!(paste0(tf,' , ',target) %in% paste0(dorotheaData$tf,' , ',dorotheaData$target))),
#                                    dorotheaData) %>% unique(), 
#                              options =  settings)
protein_activity = run_viper(PCA_alldata$rotation, 
                             interactions %>% mutate(confidence = 'A') %>% 
                               select(c('tf'='source'),confidence,target,c('mor'='interaction')), 
                             options =  settings)
pheatmap(protein_activity[,paste0('PC',c(8,11,12,13,35))])
protein_activity <- as.data.frame(protein_activity) %>% rownames_to_column('prot')
# Build distribution of prot activities across all PCs and 
# check what prot activities you can anyway observe in the data
null_activity <- protein_activity %>% gather('PC','null_act',-prot) %>% select(-PC)
protein_activity$significant <- ''
protein_activity_PC8 <- left_join(protein_activity %>% select(prot,PC8,significant),
                                  null_activity) %>%
  group_by(prot) %>% mutate(p.value = ifelse(PC8>=0,
                                             sum(null_act>=PC8)/(ncol(protein_activity)-2),
                                             sum(null_act<=PC8)/(ncol(protein_activity)-2))) %>%
  ungroup()
adjustement_tmp <- protein_activity_PC8 %>% select(prot,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(prot,p.adj)
protein_activity_PC8 <- left_join(protein_activity_PC8,adjustement_tmp)%>% select(-null_act) %>% unique()
protein_activity_PC12 <- left_join(protein_activity %>% select(prot,PC12,significant),
                                   null_activity) %>%
  group_by(prot) %>% mutate(p.value = ifelse(PC12>=0,
                                             sum(null_act>=PC12)/(ncol(protein_activity)-2),
                                             sum(null_act<=PC12)/(ncol(protein_activity)-2))) %>%
  ungroup()
adjustement_tmp <- protein_activity_PC12 %>% select(prot,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(prot,p.adj)
protein_activity_PC12 <- left_join(protein_activity_PC12,adjustement_tmp) %>% select(-null_act) %>% unique()
protein_activity_PC8$significant[order(-abs(protein_activity_PC8$PC8))[1:20]] <- protein_activity_PC8$prot[order(-abs(protein_activity_PC8$PC8))[1:20]]
protein_activity_PC8 <- protein_activity_PC8[order(protein_activity_PC8$PC8),]
protein_activity_PC8$prot <- factor(protein_activity_PC8$prot,levels = protein_activity_PC8$prot)
protein_activity_PC12$significant[order(-abs(protein_activity_PC12$PC12))[1:20]] <- protein_activity_PC12$prot[order(-abs(protein_activity_PC12$PC12))[1:20]]
protein_activity_PC12 <- protein_activity_PC12[order(protein_activity_PC12$PC12),]
protein_activity_PC12$prot <- factor(protein_activity_PC12$prot,levels = protein_activity_PC12$prot)
# change significant
protein_activity_PC12 <- protein_activity_PC12 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                             ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
protein_activity_PC8 <- protein_activity_PC8 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                           ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))

# protein_activity_PC12 <- protein_activity_PC12 %>% mutate(statistical=ifelse(p.adj<=0.01,'p.adj<=0.01',
#                                                                              ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                                     ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1'))))
# protein_activity_PC8 <- protein_activity_PC8 %>% mutate(statistical=ifelse(p.adj<=0.01,'p.adj<=0.01',
#                                                                            ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                                   ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1'))))

# Combine into one data frame
protein_activity_plot_frame <- rbind(protein_activity_PC8 %>% dplyr::rename(value = PC8) %>% mutate(PC='PC8'),
                                     protein_activity_PC12 %>% dplyr::rename(value = PC12) %>% mutate(PC='PC12'))
protein_activity_plot_frame <- protein_activity_plot_frame %>% mutate(PC=paste0('used ',PC,' gene loadings'))
### Make barplot to look at top prots
# p <- (ggplot(protein_activity_PC8,aes(x=prot,y=PC8,color = statistical)) + geom_point() +
#   # scale_color_gradient(high = 'red',low='white')+
#     scale_color_manual(values = c('#fa8e8e','black'))+
#   geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
#   xlab('prots') + ylab('activity from PC8 loadings')+
#   scale_x_discrete(expand = c(0.1, 0.1))+
#   theme_pubr(base_family = 'Arial',base_size = 20)+
#   theme(text = element_text(family = 'Arial',size=20),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = 'none')) +
#   (ggplot(protein_activity_PC12,aes(x=prot,y=PC12,color = statistical)) + geom_point() +
#      scale_color_manual(values =c('red','#fa8e8e','black'))+
#      # scale_color_gradient(high = 'red',low='white')+
#      geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
#      xlab('prots') + ylab('activity from PC12 loadings')+
#      scale_x_discrete(expand = c(0.1, 0.1))+
#      theme_pubr(base_family = 'Arial',base_size = 20)+
#      theme(text = element_text(family = 'Arial',size=20),
#            axis.ticks.x = element_blank(),
#            axis.text.x = element_blank(),
#            legend.position = 'none'))
p <- ggplot(protein_activity_plot_frame %>% dplyr::rename(`statistical threshold`=statistical),
            aes(x=value,y=-log10(p.adj),color = `statistical threshold`)) + geom_point() +
  scale_color_manual(values = c('black'))+
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40,point.padding = 0.5)+
  xlab('activity') + ylab(expression(-log[10]('adjusted p-value'))) +
  geom_hline(yintercept=-log10(0.1), linetype="dashed",color = "#525252", size=1) +
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'top')+
  facet_wrap(~PC,scales = 'free_x')+
  guides(color = guide_legend(
    override.aes = list(
      linetype = NA,
      size = 3
    )
  ))
print(p)
ggsave('../results/pc_loadings_scores_analysis/protein_activity_from_loadings_volcano.png',
       plot = p,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

### Perform KEGG pathway analysis (GSEA)------------
# rownames(loadings) <- NULL
Woptim <- readRDS('../results/Wm_extra.rds')
rownames(Woptim) <- rownames(Wm)
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Woptim), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Woptim[inds,])
rownames(meas) <- entrez_ids
keggs <- fastenrichment(colnames(meas),
                        rownames(meas),
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
# df_keggs <- df_keggs %>% mutate(sig=ifelse(abs(NES)>1.5 & padj<0.1,'yes','no'))
# df_keggs <- df_keggs%>% mutate(label = ifelse(sig=='yes',pathway,NA))
p1 <- (ggplot(df_keggs %>% filter(abs(NES)>1) %>% filter(padj<0.5) %>% filter(LV=='V1') %>% 
                arrange(NES),aes(x=NES,y=reorder(pathway,-NES),fill=padj))+ 
         geom_bar(stat = 'identity',color='black',size=1.5) +
         scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         xlab('Normalized Enrichment Score') + ylab('KEGG Pathway')+
         ggtitle('KEEGs enriched in extra basis 1')+
         theme_pubr(base_family = 'Arial',base_size = 15)+
         theme(text = element_text(family = 'Arial',size=15),
               axis.text.y = element_text(size=13),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'right',
               legend.justification = "center"))+
  (ggplot(df_keggs %>% filter(abs(NES)>1)%>% filter(padj<0.25) %>% filter(LV=='V2') %>% 
            arrange(NES),aes(x=NES,y=reorder(pathway,-NES),fill=padj))+ 
     geom_bar(stat = 'identity',color='black',size=1.5) +
     scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
     xlab('Normalized Enrichment Score') +ylab('KEEG Pathway')+
     ggtitle('KEEGs enriched in extra basis 2')+
     theme_pubr(base_family = 'Arial',base_size = 15)+
     theme(text = element_text(family = 'Arial',size=15),
           axis.text.y = element_text(size=13),
           plot.title = element_text(hjust = 0.5),
           legend.key.size = unit(1.5, "lines"),
           legend.position = 'right',
           legend.justification = "center",
           axis.title.y = element_blank()))
print(p1)
ggsave('../results/pc_loadings_scores_analysis/keggs_on_optimal_loadings.png',
       plot=p1,
       width=12,
       height=9,
       units = 'in',
       dpi = 600)

### Perform GO Terms analysis (GSEA)------------
loadings <- PCA_alldata$rotation
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,PC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$PC12))[1:10]] <- loadings$gene[order(-abs(loadings$PC12))[1:10]]
loadings <- loadings[order(loadings$PC12),]
loadings$gene <- factor(loadings$gene,levels = loadings$gene)
meas <- loadings[,2]
# meas <- as.matrix(loadings[inds,1])
# rownames(meas) <- NULL
# rownames(meas) <- entrez_ids
names(meas) <- loadings$gene

## run gsea
genes <- factor(x = rep(1,length(meas)),levels = c(0,1))
names(genes) <- names(meas)

GOobject <- new("topGOdata",ontology = "BP", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                ID = "symbol", nodeSize = 10)
term.genes <- genesInTerm(GOobject, GOobject@graph@nodes)
gos <- fgseaSimple(pathways = term.genes,
                   stats=meas,
                   minSize=10,
                   maxSize=500,
                   nperm = 10000)
colnames(gos)[1] <- 'GO'
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
df_gos <- left_join(gos,go_annotations %>% filter(ontology=='BP') %>% dplyr::select(GO,`GO Terms`))
df_gos <- df_gos %>% mutate(sig=ifelse(abs(NES)>2 & padj<0.05,'yes','no'))
df_gos <- df_gos%>% mutate(label = ifelse(sig=='yes',`GO Terms`,NA))
p1 <- ggplot(df_gos %>% arrange(NES),aes(x=NES,y=-log10(padj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('GO Terms enriched in PC12 loadings')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

df_gos <- df_gos[order(df_gos$NES),]
df_gos$`GO Terms` <- factor(df_gos$`GO Terms`,levels = df_gos$`GO Terms`)
p2 <- ggplot(df_gos %>% arrange(NES) %>% filter(padj<0.05 & abs(NES)>2),aes(x=NES,y=`GO Terms`,fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(low = "red",high = "white",limits = c(min(df_gos$padj),0.05)) +
  xlab('Normalized Enrichment Score') +
  ggtitle('GO Terms enriched in PC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/pc_loadings_scores_analysis/go_on_pc12.png',
       plot=p2,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)

### Do GO Enrichment for Woptim
gos <- fastenrichment(colnames(Woptim),
                      rownames(Woptim),
                      Woptim,
                      enrichment_space = 'go_bp',
                      gene_id_type="symbol",
                      n_permutations = 10000,
                      order_columns=F)
go_nes <- as.data.frame(gos$NES$`NES GO BP`) %>% rownames_to_column('GO')
colnames(go_nes)[2] <- 'NES'
go_pval <- as.data.frame(gos$Pval$`Pval GO BP`) %>% rownames_to_column('GO')
colnames(go_pval)[2] <- 'padj'
df_gos <- left_join(go_nes,go_pval)
df_gos <- df_gos %>% mutate(GO=strsplit(GO,"_"))
df_gos <- df_gos %>% unnest(GO) %>% filter(!(GO %in% c("GO","FL1000","BP")))
df_gos <- df_gos %>% mutate(GO=as.character(GO))
df_gos <- left_join(df_gos,go_annotations)

### GSEA with Hallmarks MSIGDB------------------------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Woptim), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Woptim[inds,])
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

# df_msig <- df_msig %>% mutate(Hallmark=strsplit(Hallmark,"_"))
# df_msig <- df_msig %>% unnest(Hallmark) %>% filter(!(Hallmark %in% c("HALLMARK","FL1000","MSIG","H")))
# df_msig <- df_msig %>% mutate(Hallmark=as.character(Hallmark))
# df_msig <- df_msig %>% mutate(Hallmark=gsub('_'," ",Hallmark))
# df_msig_tmp <- df_msig %>% filter(PC=='PC12')
# df_msig$Hallmark <- factor(df_msig$Hallmark,levels = df_msig$Hallmark[order(df_msig$NES)])
p1 <- (ggplot(df_msig %>% filter(LV=='V1') %>% arrange(NES),aes(x=NES,y=reorder(Hallmark,-NES),fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
  xlab('Normalized Enrichment Score') + ylab('Hallmark')+
  ggtitle('Hallmarks enriched in extra basis 1')+
  theme_pubr(base_family = 'Arial',base_size = 15)+
  theme(text = element_text(family = 'Arial',size=15),
        axis.text.y = element_text(size=13),
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
           axis.text.y = element_text(size=13),
           plot.title = element_text(hjust = 0.5),
           legend.key.size = unit(1.5, "lines"),
           legend.position = 'right',
           legend.justification = "center",
           axis.title.y = element_blank()))
print(p1)
ggsave('../results/pc_loadings_scores_analysis/hallmark_on_optimal_loadings.png',
       plot=p1,
       width=12,
       height=9,
       units = 'in',
       dpi = 600)