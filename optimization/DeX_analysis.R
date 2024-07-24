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
library(msigdbr)
source('vector_space_interpretation.R')
source('functions_translation.R')
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
tf_activities <- TF_activity_interpretation(dex_data,Wm,dorotheaData,plotting=TRUE,lim=10)
p_tf <- tf_activities$figure[[1]] + ggtitle(expression('Transcription factors driving gene '*Delta*'X'))
print(p_tf)
ggsave(paste0('../figures/perturbed_tf_activity_',
              tolower(target_dataset),
              '_dorothea.png'),
       plot = p_tf,
       width = 12,
       height = 6,
       dpi = 600)
ggsave(paste0('../figures/perturbed_tf_activity_',
              tolower(target_dataset),
              '_dorothea.eps'),
       plot = p_tf,
       device = cairo_ps,
       width = 12,
       height = 6,
       dpi = 600)

### Perform Progeny analysis using t-statistic OR directly dXs ---------------
# dex_data <- top.table_lean %>% select(gene,logFC) %>% column_to_rownames('gene')
dex_data <- t(dx_lean)
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
# path_activities = decoupleR::run_viper(dex_data, net_prog,minsize = 1)
path_activities <- pathway_activity_interpretation(dex_data,
                                                   Wm,
                                                   plotting = TRUE)
path_activities <- path_activities$figure[[1]]
path_activities <- path_activities + ggtitle('Patwhay activity') +
  theme(text = element_text(size = 18,family = 'Arial'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=19),
        # plot.title = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none')
print(path_activities)

### Hallmarks GSEA using dX--------------------------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(dex_data), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(dex_data[inds,])
rownames(meas) <- entrez_ids
msig <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes <- as.data.frame(msig$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark')  #%>% gather('PC','NES',-Hallmark)
colnames(msig_nes) <- c('Hallmark','NES')
msig_pval <- as.data.frame(msig$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')#%>% gather('PC','padj',-Hallmark)
colnames(msig_pval) <- c('Hallmark','padj')
df_msig <- left_join(msig_nes,msig_pval)
df_msig <- df_msig %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))
df_msig <- df_msig %>% mutate(Hallmark = str_replace_all(Hallmark,'_'," "))
p_msig <- ggplot(df_msig %>% #mutate(Hallmark=tolower(Hallmark)) %>% 
                   arrange(NES) %>% filter(padj<=0.1),
                 aes(x=NES,y=reorder(Hallmark,NES),fill=NES))+ 
             geom_bar(stat = 'identity') +
             # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
             scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
             xlab('Normalized Enrichment Score') + ylab('Hallmark')+
             theme_minimal(base_family = 'Arial',base_size = 18)+
             theme(text = element_text(family = 'Arial',size=18),
                   axis.text.y = element_text(size=15),
                   axis.title.y = element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   # legend.key.size = unit(1.5, "lines"),
                   legend.position = 'none')
print(p_msig)
ggsave(paste0('../figures/hallmarks_',
              tolower(ref_dataset),'_',
              tolower(target_dataset),
              '_max_var.png'),
       plot=p_msig,
       width=9,
       height=9,
       units = 'in',
       dpi = 600)

combined_plot <- path_activities + p_msig#+
  # plot_annotation(title = "Optimized MPS to maximize the total human variace captured",
  #                 theme = theme(plot.title = element_text(size = 18, hjust = 0.5, family = 'Arial')))
print(combined_plot)

ggsave(paste0('../figures/perturbed_pathway_and_msig_',
              tolower(target_dataset),
              '_barplot.png'),
       plot=combined_plot,
       width = 18,
       height = 9,
       dpi = 600)
setEPS()
postscript(paste0('../figures/perturbed_pathway_and_msig_',
                  tolower(target_dataset),
                  '_barplot.eps'),
           width = 18,
           height = 9)
print(combined_plot)
dev.off()

# Interrogate dX for sex-related genesets enrichment-------------------------------------
# Retrieve all gene sets
# List available libraries
msigdb_gene_sets <- msigdbr(species = "Homo sapiens")
# Filter for gene sets related to sex differences
sex_diff_gene_sets <- msigdb_gene_sets %>%
  filter(grepl("sex|gender", gs_name, ignore.case = TRUE))
all_sex_genesets <- sex_diff_gene_sets %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()
all_genesets <- msigdb_gene_sets %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()
### Gene-set enrichment witll MSIG data
msig_sex_diff <- apply(dex_data,
                       MARGIN = 2,
                       fgsea,pathways = all_sex_genesets,
                       minSize=10,
                       maxSize=500,
                       nperm = 10000)
df_msig_sex <- msig_sex_diff[[1]]
df_msig_sex <- df_msig_sex %>% select(c('geneset'='pathway'),ES,NES,pval,padj) %>% unique()
df_msig_sex <- left_join(df_msig_sex,distinct(msigdb_gene_sets %>% select(c('geneset'='gs_name'),gs_cat,gs_subcat)))
df_msig_sex <- df_msig_sex %>% mutate(type = ifelse(geneset %in% sex_diff_gene_sets$gs_name,'sex diff','other')) 
df_msig_sex <- df_msig_sex %>% 
  mutate(geneset = ifelse(gs_subcat=='GO:BP',substr(geneset, nchar('GOBP_')+1, nchar(geneset)),
ifelse(gs_subcat=='GO:CC',substr(geneset, nchar('GOCC_')+1, nchar(geneset)),
       ifelse(gs_subcat=='HPO',substr(geneset, nchar('HP_')+1, nchar(geneset)),
              ifelse(gs_subcat=='CP:WIKIPATHWAYS',substr(geneset, nchar('WP_')+1, nchar(geneset)),
                     ifelse(gs_subcat=='IMMUNESIGDB',substr(geneset, nchar('GSE18281_')+1, nchar(geneset)),
                            ifelse(gs_subcat=='GO:MF',substr(geneset, nchar('GOMF_')+1, nchar(geneset)),
                                   geneset)))))))%>% 
  mutate(geneset = str_replace_all(geneset,'_'," "))
# df_msig_sex <- df_msig_sex %>% mutate(geneset=tolower(geneset))
# df_msig_sex <- df_msig_sex %>% filter(padj<0.25) %>% filter(abs(NES)>=0.5)
p_msig_sex_diff <- ggplot(df_msig_sex %>% 
                   arrange(NES), #%>% filter(padj<=0.1)
                 aes(x=NES,y=reorder(geneset,NES),fill=padj))+ 
  geom_bar(stat = 'identity') +
  scale_fill_gradient(trans='log10',low='indianred',high = 'whitesmoke',limits = c(0.01,1)) +
  # scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
  geom_text(aes(label = ifelse(geneset %in% c('SOMATIC SEX DETERMINATION','MALE SEX DIFFERENTIATION'),
                               paste0('p-adj.=',round(padj,3)),NA),
                x = ifelse(NES < 0, 0.3, -0.3)),
            size = 6,
            color = 'black',
            angle=0) + 
  xlab('Normalized Enrichment Score') + ylab('geneset')+
  ggtitle('Enrichment of sex-related genesets from the MSIG database')+
  labs(fill = "p-adj.")+
  theme_minimal(base_family = 'Arial',base_size = 18)+
  theme(text = element_text(family = 'Arial',size=18),
        axis.text.y = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        # legend.key.size = unit(1.5, "lines"),
        legend.position = 'right')
print(p_msig_sex_diff)
ggsave(paste0('../figures/perturbed_msig_sex_',
              tolower(target_dataset),
              '_barplot.png'),
       plot=p_msig_sex_diff,
       width = 16,
       height = 8,
       dpi = 600)
ggsave(paste0('../figures/perturbed_msig_sex_',
              tolower(target_dataset),
              '_barplot.eps'),
       plot=p_msig_sex_diff,
       device = cairo_ps,
       width = 16,
       height = 8,
       dpi = 600)
