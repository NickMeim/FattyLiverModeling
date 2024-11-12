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
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
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
## get in vitro metadata and keep only control high NPC lean
invitro_metadata <- as.data.frame(data_list[[target_dataset]]$metadata)
invitro_metadata <- invitro_metadata[,1:13] %>% filter(Number_of_cues==0) %>% filter(NPC=='high') %>% filter(Background=='lean')
# Xm <- Xm[invitro_metadata$sampleName,]

## Load perturbation
dx_lean <- data.table::fread('../results/optimized_mps/dx_lean_govaere_kostrzewski_corrected.csv') %>% select(-V1)

### Visualize genes in dx ---------------
dex_data <- t(dx_lean)
plot_genes_dx <- plot_gene_loadings(loadings = dex_data,
                                    selection='V1',
                                    y_lab = expression('proposed '*Delta*'X on lean samples'),
                                    top=20)
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
# Save
saveRDS(df_msig, file = paste0("../results/hallmark_enrichment_", tolower(target_dataset),"_dx_lean.rds"))
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