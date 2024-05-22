# load in packages
library(tidyverse)
library(DESeq2)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(patchwork)
library(pls)
library(pheatmap)
library(patchwork)
library(limma)
library(edgeR)
library(dorothea)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(hgu133a.db)
library(factoextra)
source('vector_space_interpretation.R')
source('enrichment_calculations.R')

# Load data and do PCA----------------------------------------------------------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Hoang"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Hoang$metadata  %>% select(nafld_activity_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()
## get in vitro metadata and keep only control high NPC lean
invitro_metadata <- as.data.frame(data_list[[target_dataset]]$metadata)
invitro_metadata <- invitro_metadata[,1:13] %>% filter(Number_of_cues==0) %>% filter(NPC=='high') %>% filter(Background=='lean')
Xm <- Xm[invitro_metadata$sampleName,]

## Load perturbation and perturb
dx_lean <- data.table::fread('../results/optimized_mps/dx_lean_govaere_kostrzewski.csv') %>% select(-V1)
Xper <- as.matrix(Xm) + do.call(rbind, replicate(nrow(Xm), as.matrix(dx_lean), simplify=FALSE))
X_all <- rbind(Xm,Xper)
conditions <- c(rep('control',nrow(Xm)),rep('perturbed',nrow(Xm)))
rownames(X_all) <- conditions
pca <- prcomp(X_all,scale. = FALSE, center = TRUE)
fviz_screeplot(pca)
df_pcs <- as.data.frame(pca$x)
df_pcs$condition <- conditions
ggplot(df_pcs,aes(x=PC1,y=PC2,color=condition)) +
  geom_point(size=1.5)+
  theme_pubr(base_size=20,base_family = 'Arial')
w <- as.matrix(pca$rotation[,1])
colnames(w) <- 'P1'
plot_extra_gene_loadings_pc1 <- plot_gene_loadings(loadings = w,
                                                   selection='P1',
                                                   y_lab = 'gene weight')
optimal_var_pathway_activity <- pathway_activity_interpretation(w,
                                                                Wm,
                                                                plotting = FALSE)
optimal_var_pathway_activity <- optimal_var_pathway_activity[[2]]
ggplot(optimal_var_pathway_activity %>% select(c('activity'='score'),Pathway,p_value,condition) %>% 
         filter (condition=='P1'),
       aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity',color='black') +
  scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
  # scale_x_continuous(n.breaks = 8,limits = c(-8,8))+
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
  theme_pubr(base_size = 18,base_family = 'Arial')+
  theme(text = element_text(size = 18,family = 'Arial'),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))
ggsave(paste0('../results/optimized_mps/progenies_optimal_variance_paths_',
              tolower(target_dataset),
              '_barplot.png'),
       width = 12,
       height = 9,
       dpi = 600)

### Perform dorothea analysis using directly dXs ---------------
# dex_data <- top.table_lean %>% select(gene,logFC) %>% column_to_rownames('gene')
dex_data <- t(dx_lean)
### get high-confidence dorothea regulon
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source'
TF_activities = decoupleR::run_viper(dex_data, dorotheaData,minsize = 5)
TF_activities <- TF_activities%>%
  mutate(significance = ifelse(p_value<0.01,ifelse(score>1.5,'upregulated',ifelse(score<(-1.5),'downregulated',NA)),NA))
TF_activities <- TF_activities %>% mutate(label = ifelse(!is.na(significance),source,NA))
ggplot(TF_activities,aes(x=score,y=-log10(p_value),color=significance)) + geom_point() +
  ylab(expression(-log[10]('p-value'))) +
  scale_color_manual(values=c("#EC2326","#008c41"))+
  geom_hline(yintercept=-log10(0.01), linetype="dashed",color = "#525252", size=0.5) +
  geom_vline(xintercept=-1.5, linetype="dashed",color = "#EC2326", size=1) +
  geom_vline(xintercept=1.5, linetype="dashed",color = "#008c41", size=1)+ 
  geom_label_repel(aes(label=label),size=6)+
  theme(text = element_text(size=18,family = 'Arial'),
        plot.title = element_text(hjust=0.5))
ggsave(paste0('../results/optimized_mps/perturbed_tf_activity_',
              tolower(target_dataset),
              '_barplot.png'),
       width = 12,
       height = 9,
       dpi = 600)

### Perform Progeny analysis using t-statistic OR directly dXs ---------------
# dex_data <- top.table_lean %>% select(gene,logFC) %>% column_to_rownames('gene')
dex_data <- t(dx_lean)
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
path_activities = decoupleR::run_viper(dex_data, net_prog,minsize = 1)
# path_activities <- pathway_activity_interpretation(dex_data,
#                                                    Wm,
#                                                    plotting = FALSE)
# path_activities <- path_activities[[2]]
ggplot(path_activities %>% select(c('activity'='score'),source,p_value) ,
       aes(x=activity,y=reorder(source,activity),fill=activity)) + geom_bar(stat = 'identity',color='black') +
  scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
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
  theme_pubr(base_size = 18,base_family = 'Arial')+
  theme(text = element_text(size = 18,family = 'Arial'),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))
ggsave(paste0('../results/optimized_mps/perturbed_pathway_activity_',
              tolower(target_dataset),
              '_barplot.png'),
       width = 12,
       height = 9,
       dpi = 600)

# use limma with lean transformed data-------------------------------------
genes_var <- apply(X_all,2,var)
inds <- which(genes_var>1e-6)
design <- model.matrix(~0+conditions)
colnames(design) <- c('control','perturbed')
X_all <- X_all[,inds]
fit_lean <- lmFit(t(X_all),design,method = 'ls')
# plotSA(fit_lean, main="Final model: Mean-variance trend")
contr_lean <- makeContrasts(contrasts = 'perturbed-control',levels = colnames(coef(fit_lean)))
contr_fit_lean <- contrasts.fit(fit_lean, contr_lean)
fit_ebayes_lean <- eBayes(contr_fit_lean,robust=TRUE,trend=TRUE)
# plotSA(fit_ebayes_lean, main="Final model: Mean-variance trend")
top.table_lean <- topTable(fit_ebayes_lean,number = ncol(X_all))
# Make adjustements to get the desired plots
top.table_lean <- top.table_lean %>% mutate(significance = ifelse(adj.P.Val<0.05,ifelse(logFC>1.5,'upregulated',ifelse(logFC<(-1.5),'downregulated',NA)),NA))
top.table_lean <- top.table_lean %>% rownames_to_column('gene')
top.table_lean <- top.table_lean %>% mutate(label = ifelse(!is.na(significance),gene,NA))
ggplot(top.table_lean,aes(x=logFC,y=-log10(adj.P.Val),color=significance)) + geom_point() +
  xlab(expression(log[2]('Fold Change'))) + ylab(expression(-log[10]('adjusted p-value'))) +
  scale_color_manual(values=c("#EC2326","#008c41"))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "#525252", size=0.5) +
  geom_vline(xintercept=-1.5, linetype="dashed",color = "#EC2326", size=1) +
  geom_vline(xintercept=1.5, linetype="dashed",color = "#008c41", size=1)+ 
  # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
  geom_label_repel(aes(label=label),size=6)+
  ggtitle('Proposed perturbations for lean-derived samples')+
  theme(text = element_text(size=18,family = 'Arial'),
        plot.title = element_text(hjust=0.5))
ggsave('../results/optimized_mps/lean_volcano.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
