library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(pROC)
source('functions_translation.R')
source('CrossValidationUtilFunctions.R')
source("../utils/plotting_functions.R")
source("vector_space_interpretation.R")
source('enrichment_calculations.R')

# function to calculate quickly the Jaccard similarity
jaccard <- function(M=NULL,x=NULL,y=NULL){
  # First calculate for positive sign
  if (is.null(M)){
    xp <- 1* (x==1)
    yp <- 1* (y==1)
    intersection_p  <-  sum(xp*yp)
    union_p  <-  sum(xp) + sum(yp) - intersection_p
  }else{
    Mp <- 1*(M==1)
    intersection_p <- Mp %*% t(Mp)
    union_p <- apply(as.matrix(rowSums(Mp)),1,'+',as.matrix(rowSums(Mp))) - intersection_p
  }
  Jp <- intersection_p/union_p
  
  # Then calculate negatives
  if (is.null(M)){
    xn <- 1* (x==(-1))
    yn <- 1* (y==(-1))
    intersection_n  <-  sum(xn*yn)
    union_n  <-  sum(xn) + sum(yn) - intersection_n
  }else{
    Mn <- 1*(M==(-1))
    intersection_n <- Mn %*% t(Mn)
    union_n <- apply(as.matrix(rowSums(Mn)),1,'+',as.matrix(rowSums(Mn))) - intersection_n
  }
  Jn <- intersection_n/union_n
  J <- 0.5 * (Jp + Jn)
  return(J)
}

jaccard_2mats <- function(M1,M2){
  # First calculate for positive sign
  M1p <- 1*(M1==1)
  M2p <- 1*(M2==1)
  intersection_p <- M1p %*% t(M2p)
  union_p <- apply(as.matrix(rowSums(M1p)),1,'+',as.matrix(rowSums(M2p))) - intersection_p
  Jp <- intersection_p/union_p
  
  # Then calculate negatives
  M1n <- 1*(M1==(-1))
  M2n <- 1*(M2==(-1))
  intersection_n <- M1n %*% t(M2n)
  union_n <- apply(as.matrix(rowSums(M1n)),1,'+',as.matrix(rowSums(M2n))) - intersection_n
  Jn <- intersection_n/union_n
  J <- 0.5 * (Jp + Jn)
  return(J)
}

### Load extra basis and dX from variance optimization and infer TF activity, Hallmark enrichment, KEGG enrichment,GO Terms enrichment------------------------
Wm_opt <- readRDS('../results/Wm_kostrzewski_extra.rds')
dx_lean <- data.table::fread('../results/optimized_mps/dx_lean_govaere_kostrzewski.csv') %>% select(-V1)
W_all <- cbind(Wm_opt,t(dx_lean))

### TF enrichment
dorotheaData <- read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter <- is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData <- dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
minNrOfGenes  <-  5
TF_activities_extra <- decoupleR::run_viper(W_all, dorotheaData,minsize = minNrOfGenes,verbose = FALSE)
TF_activities_extra_filt <- TF_activities_extra %>% filter(p_value<0.01) %>% filter(abs(score)>0.5)
TF_activities_extra_filt <- TF_activities_extra %>% mutate(score = ifelse(p_value>=0.01,0,score)) %>% 
  mutate(score = ifelse(abs(score)<=0.5,0,score)) %>% 
  mutate(sig = ifelse(score<0,-1,
                      ifelse(score>0,1,0))) %>% 
  select(c('sig_id'='condition'),c('TF'='source'),sig)
TF_activities_extra_filt <- as.matrix(TF_activities_extra_filt %>% spread('TF','sig') %>% column_to_rownames('sig_id'))
TF_activities_extra <- as.matrix(TF_activities_extra %>% select(condition,source,score) %>% spread('source','score')%>% 
                                   column_to_rownames('condition'))

### Hallmark enrichment
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(W_all), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(W_all[inds,])
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
df_msig_extra_filt <- df_msig %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(c('sig_id'='LV'),Hallmark,sig)
df_msig_extra_filt <- as.matrix(df_msig_extra_filt %>% spread('Hallmark','sig') %>% column_to_rownames('sig_id'))
df_msig_extra <- as.matrix(df_msig %>% select(LV,Hallmark,NES) %>% spread('Hallmark','NES') %>% column_to_rownames('LV'))

### KEGG enrichment
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(W_all), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(W_all[inds,])
rownames(meas) <- entrez_ids
kegg <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'kegg',
                       n_permutations = 10000,
                       order_columns=F)
kegg_nes <- as.data.frame(kegg$NES$`NES KEGG`) %>% rownames_to_column('KEGG') 
kegg_nes <- kegg_nes %>% gather('LV','NES',-KEGG)
kegg_pval <- as.data.frame(kegg$Pval$`Pval KEGG`) %>% rownames_to_column('KEGG')
kegg_pval <- kegg_pval %>% gather('LV','padj',-KEGG)
df_kegg <- left_join(kegg_nes,kegg_pval)
df_kegg <- df_kegg %>% mutate(KEGG=strsplit(KEGG,"_"))
df_kegg <- df_kegg %>% unnest(KEGG) %>% filter(!(KEGG %in% c("KEGG","FL1000")))
df_kegg <- df_kegg %>% mutate(KEGG=as.character(KEGG))
df_kegg <- df_kegg %>% mutate(KEGG=substr(KEGG, 9, nchar(KEGG)))
df_kegg_extra_filt <- df_kegg %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(c('sig_id'='LV'),KEGG,sig)
df_kegg_extra_filt <- as.matrix(df_kegg_extra_filt %>% spread('KEGG','sig') %>% column_to_rownames('sig_id'))
df_kegg_extra <- as.matrix(df_kegg %>% select(LV,KEGG,NES) %>% spread('KEGG','NES') %>% column_to_rownames('LV'))

## GO Terms BP
gos <- fastenrichment(colnames(W_all),
                      rownames(W_all),
                      W_all,
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
df_gos_bp_extra_filt <- df_gos %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(c('sig_id'='LV'),GO,sig)
df_gos_bp_extra_filt <- as.matrix(df_gos_bp_extra_filt %>% spread('GO','sig') %>% column_to_rownames('sig_id'))
df_gos_bp_extra <- as.matrix(df_gos %>% select(LV,GO,NES) %>% spread('GO','NES') %>% column_to_rownames('LV'))

## GO Terms CC
gos <- fastenrichment(colnames(W_all),
                      rownames(W_all),
                      W_all,
                      enrichment_space = 'go_cc',
                      gene_id_type = 'symbol',
                      n_permutations = 10000,
                      order_columns=F)
go_nes <- as.data.frame(gos$NES$`NES GO CC`) %>% rownames_to_column('GO')
go_nes <- go_nes %>% gather('LV','NES',-GO)
go_pval <- as.data.frame(gos$Pval$`Pval GO CC`) %>% rownames_to_column('GO')
go_pval <- go_pval %>% gather('LV','padj',-GO)
df_gos <- left_join(go_nes,go_pval)
df_gos <- df_gos %>% mutate(GO=strsplit(GO,"_"))
df_gos <- df_gos %>% unnest(GO) %>% filter(!(GO %in% c("GO","FL1000","CC")))
df_gos <- df_gos %>% mutate(GO=as.character(GO))
df_gos_cc_extra_filt <- df_gos %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(c('sig_id'='LV'),GO,sig)
df_gos_cc_extra_filt <- as.matrix(df_gos_cc_extra_filt %>% spread('GO','sig') %>% column_to_rownames('sig_id'))
df_gos_cc_extra <- as.matrix(df_gos %>% select(LV,GO,NES) %>% spread('GO','NES') %>% column_to_rownames('LV'))

### Load L1000 gex and enrichment data to calculate Jaccard and cosine similarities---------------------
# Load the L1000 dataset
drug_ligand_ex <- readRDS("../data/l1000_drugs_ligands_expression.rds")
gene_info <- readRDS("../data/l1000_geneInfo.rds")
print(all(gene_info$gene_id==rownames(drug_ligand_ex)))
rownames(drug_ligand_ex) <- gene_info$gene_symbol
cell_info <- readRDS("../data/l1000_cell_line_info_ligands_drugs.rds")
# cell_info <- cell_info %>% filter(cell_lineage=='liver')
sigInfo <- readRDS("../data/l1000_meta_data.rds")
# sigInfo <- sigInfo %>% filter(cell_iname %in% cell_info$cell_iname)
drug_ligand_ex <- drug_ligand_ex[,sigInfo$sig_id]
sigInfo <- sigInfo %>% mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

### Run progeny to infer pathway activity and keep only drugs that can affect these pathways in some cell line
### The pathways of interest are JAK-STAT (or NFKb) at opposite sign as p53
paths <- c('JAK-STAT','NFkB','p53')
combined_sign_jakstat <- -1
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
L1000_progenies <- decoupleR::run_viper(drug_ligand_ex, net_prog,minsize = 1,verbose = TRUE) %>% select(-statistic)
L1000_progenies_filt <- L1000_progenies %>% filter(source %in% paths) %>% group_by(condition) %>% mutate(min_p = min(p_value)) %>%
  ungroup() %>% filter(min_p<0.01) %>% select(-min_p)
L1000_progenies_filt_jax_nfkb <- L1000_progenies_filt %>% filter(source %in% paths[1:2]) %>% 
  select(-p_value) %>% spread('source','score') %>%
  group_by(condition) %>% mutate(s = sign(`JAK-STAT` * NFkB)) %>% ungroup() %>%
  gather('source','score',-condition,-s) %>% filter(s>0)
L1000_progenies_filt <- L1000_progenies_filt %>% filter(condition %in% L1000_progenies_filt_jax_nfkb$condition) %>%
  filter(source %in% c('JAK-STAT','p53')) %>%# we care more about JAK-STAT and p53 and since we made sure NFkb has the same sign as JAK-STAT we are going to ingore it
  select(-p_value) %>% spread('source','score') %>%
  group_by(condition) %>% mutate(s = sign(`JAK-STAT` * p53)) %>% ungroup() %>%
  gather('source','score',-condition,-s) %>% filter(s<0)
L1000_progenies_filt <- L1000_progenies %>% filter(condition %in% L1000_progenies_filt$condition)
sigs_selected <- unique(L1000_progenies_filt$condition)

sigInfo <- sigInfo %>% filter(sig_id %in% sigs_selected)
drug_ligand_ex <- drug_ligand_ex[,sigInfo$sig_id]

### Run viper to infer TF activity
dorotheaData <- read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter <- is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData <- dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
minNrOfGenes  <-  5
TF_activities <- decoupleR::run_viper(drug_ligand_ex, dorotheaData,minsize = minNrOfGenes,verbose = FALSE)
mat <- t(as.matrix(TF_activities %>% select(condition,source,score) %>% 
                     spread('source','score') %>% column_to_rownames('condition')))
df_tf_cosine <- lsa::cosine(mat)
df_tf_cosine[lower.tri(df_tf_cosine,diag = T)]<- -100
df_tf_cosine <- reshape2::melt(df_tf_cosine)
df_tf_cosine <- df_tf_cosine %>% filter(value != -100)
gc()
df_tf_cosine <- left_join(df_tf_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                          by = c("Var1"="sig_id"))
df_tf_cosine <- left_join(df_tf_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                          by = c("Var2"="sig_id"))
df_tf_cosine <- df_tf_cosine %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
df_tf_cosine <- df_tf_cosine %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
df_tf_cosine <- df_tf_cosine %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
gc()
# ROC analysis
roc_object <- roc(df_tf_cosine$is_same_drug_cell, df_tf_cosine$value)
plot(roc_object)
optimal_threshold_cosine_tfs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
# Print the optimal threshold
print(optimal_threshold_cosine_tfs)
y_pred <- 1*(df_tf_cosine$value>optimal_threshold_cosine_tfs$threshold)
y <- df_tf_cosine$is_same_drug_cell
confusion <- confusionMatrix(data=factor(y_pred,levels = c(1,0)), 
                             reference = factor(y,levels = c(1,0)))
print(confusion)
ggviolin(df_tf_cosine,x='is_duplicate',y='value',fill = 'is_duplicate',draw_quantiles = TRUE)
ggviolin(df_tf_cosine,x='is_same_drug_cell',y='value',fill = 'is_same_drug_cell',draw_quantiles = TRUE)
# repeat for Jaccard
TF_activities_sig <- TF_activities %>% mutate(score = ifelse(p_value>=0.01,0,score)) %>% 
  mutate(score = ifelse(abs(score)<=0.5,0,score)) %>% 
  mutate(sig = ifelse(score<0,-1,
                      ifelse(score>0,1,0))) %>% 
  select(c('sig_id'='condition'),c('TF'='source'),sig)
TF_activities_sig <- as.matrix(TF_activities_sig %>% spread('TF','sig') %>% column_to_rownames('sig_id'))
jaccard_similarity <-  jaccard(TF_activities_sig)
jaccard_similarity[lower.tri(jaccard_similarity,diag = T)]<- -100
df_jaccard_tfs <- reshape2::melt(jaccard_similarity)
df_jaccard_tfs <- df_jaccard_tfs %>% filter(value != -100)
gc()
df_jaccard_tfs <- left_join(df_jaccard_tfs,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                            by = c("Var1"="sig_id"))
df_jaccard_tfs <- left_join(df_jaccard_tfs,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                            by = c("Var2"="sig_id"))
df_jaccard_tfs <- df_jaccard_tfs %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
df_jaccard_tfs <- df_jaccard_tfs %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
df_jaccard_tfs <- df_jaccard_tfs %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
gc()
# ROC analysis
roc_object <- roc(df_jaccard_tfs$is_same_drug_cell, df_jaccard_tfs$value)
plot(roc_object)
optimal_threshold_jaccard_tfs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
# Print the optimal threshold
print(optimal_threshold_jaccard_tfs)
y_pred <- 1*(df_jaccard_tfs$value>optimal_threshold_jaccard_tfs$threshold)
y <- df_jaccard_tfs$is_same_drug_cell
confusion <- confusionMatrix(data=factor(y_pred,levels = c(1,0)), 
                             reference = factor(y,levels = c(1,0)))
print(confusion)
ggviolin(df_jaccard_tfs,x='is_duplicate',y='value',fill = 'is_duplicate',add='jitter')
ggviolin(df_jaccard_tfs,x='is_same_drug_cell',y='value',fill = 'is_same_drug_cell',add='jitter')

### Load Hallmarks results


### Load the all the clinilca datasets and generate a distribution similarity scores-----------------------
### in terms of affected TFs, pathway activities, hallmark genesets, KEGGs, GO Terms
dataset_names <- c("Govaere",'Hoang' ,"Kostrzewski", "Wang", "Feaver")
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

### Identify external perturbations using L1000 data--------------
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))

### Infer pathway activity and keep drugs that affect desired pathways
paths <- c('JAK-STAT','p53')
combined_sign <- -1
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
L1000_progenies <- decoupleR::run_viper(drug_ligand_ex, net_prog,minsize = 1,verbose = TRUE) %>% select(-statistic)
L1000_progenies_filt <- L1000_progenies %>% filter(source %in% paths) %>% filter(p_value<=0.05)
L1000_progenies_filt <- L1000_progenies_filt %>% spread('source','score') %>% mutate(s = sign(`JAK-STAT` * p53)) %>%
  filter(s ==(-1) | is.na(s)) %>% select(-s,-p_value) %>% gather('pathway','score',-condition)
L1000_progenies_filt <- L1000_progenies_filt %>% filter(!is.na(score))
colnames(L1000_progenies_filt)[1] <- 'sig_id'
L1000_progenies_filt <- left_join(L1000_progenies_filt,sigInfo) %>% 
  select(sig_id,cmap_name,pert_idose,pert_itime,pert_type,cell_iname,pathway,score) %>% unique()
L1000_progenies_filt <- L1000_progenies_filt %>% group_by(cmap_name,pert_idose,pert_itime,pathway) %>% 
  mutate(mu_score = mean(score)) %>% ungroup() %>% spread('pathway','mu_score') %>% mutate(s = sign(`JAK-STAT` * p53)) %>%
  filter(s ==(-1) | is.na(s)) %>% select(-s) %>% 
  gather('pathway','mu_score',-sig_id,-cmap_name,-pert_idose,-pert_itime,-pert_type,-cell_iname,-score)
L1000_progenies_filt <- L1000_progenies_filt %>% filter(abs(mu_score)>1)
L1000_progenies_filt_liver <- L1000_progenies_filt %>% filter(cell_iname %in% c('HEPG2','HUH7')) %>%
  select(-cell_iname,-mu_score) %>% unique() %>% group_by(cmap_name,pert_idose,pert_itime,pathway) %>% 
  mutate(mu_score = mean(score)) %>% ungroup() %>% spread('pathway','mu_score') %>% mutate(s = sign(`JAK-STAT` * p53)) %>%
  filter(s ==(-1) | is.na(s)) %>% select(-s) %>% 
  gather('pathway','mu_score',-sig_id,-cmap_name,-pert_idose,-pert_itime,-pert_type,-score)
L1000_progenies_filt_liver <- L1000_progenies_filt_liver %>% filter(abs(mu_score)>1)
L1000_progenies_filt <- L1000_progenies_filt %>% filter(cmap_name %in% L1000_progenies_filt_liver$cmap_name)  
L1000_selected <- L1000_progenies_filt %>% group_by(pathway) %>%
  slice_max(order_by = abs(mu_score), n = 10) %>% ungroup() %>% 
  select(cmap_name,pert_idose,pert_itime,pert_type,pathway,mu_score) %>% unique()
L1000_selected <- L1000_selected %>% group_by(cmap_name) %>% mutate(dose_counts = n_distinct(pert_idose)) %>% 
  mutate(time_counts = n_distinct(pert_itime)) %>% ungroup()
# L1000_selected_liver <- L1000_progenies_filt_liver %>% group_by(pathway) %>%
#   slice_max(order_by = abs(mu_score), n = 10) %>% ungroup() %>% unique()
## Project L1000 data into extra basis
# drug_ligand_ex_subset <- drug_ligand_ex[which(rownames(drug_ligand_ex) %in% rownames(Wm_opt)),]
# Wm_opt_subset <- Wm_opt[which(rownames(Wm_opt) %in% rownames(drug_ligand_ex_subset)),] 
# drug_ligand_ex_subset <- drug_ligand_ex_subset[rownames(Wm_opt_subset),]
# Z_l1000 <- as.data.frame(t(drug_ligand_ex_subset) %*% Wm_opt_subset)
# Z_l1000 <- Z_l1000 %>% rownames_to_column('sample')
# Zh <- as.data.frame(Xh %*% Wm_opt)
# Zh$NAS <- Yh[,1]
# Zh$fibrosis <- Yh[,2]
# Zh <- Zh %>% rownames_to_column('sample')
# Z_l1000$NAS <- NA
# Z_l1000$fibrosis <- NA
# Z_all <- rbind(Zh,Z_l1000)
# Z_all <- Z_all %>% mutate(NAS=ifelse(is.na(NAS),'L1000',
#                                            ifelse(NAS<3,'Not NASH',
#                                                   ifelse(NAS<=4,'Borderline NASH',
#                                                          'Definate NASH'))))
# Z_all <- Z_all %>% mutate(fibrosis=ifelse(is.na(fibrosis),'L1000',fibrosis))
# (ggplot(Z_all ,aes(x=V1,y=V2,colour=NAS))+
#     geom_point()+
#     # scale_color_viridis_d()+
#     theme_minimal(base_size=20,base_family = 'Arial')+
#     theme(text= element_text(size=20,family = 'Arial'),
#           legend.position = 'right')) +
#   (ggplot(Z_all ,aes(x=V1,y=V2,colour=fibrosis))+
#      geom_point()+
#      # scale_color_viridis_d()+
#      theme_minimal(base_size=20,base_family = 'Arial')+
#      theme(text= element_text(size=20,family = 'Arial'),
#            legend.position = 'right'))
# ggsave(paste0('results/projected_L1000_',
#               tolower(ref_dataset),
#               '_samples_on_extra_basis',
#               tolower(target_dataset),
#               '.png'),
#        height = 12,
#        width = 16,
#        units = 'in',
#        dpi=600)
# 
# ## Calculate distances between L1000 and projected human samples
# mat <- as.matrix(Z_all %>% column_to_rownames('sample') %>% select(V1,V2))
# distances <- as.matrix(dist(mat, method = "euclidean", diag = TRUE, upper = TRUE))
# distances <- distances[Z_l1000$sample,Zh$sample]
# distances <- as.data.frame(distances) %>% rownames_to_column('sig_id') %>% gather('sample','dist',-sig_id)
# distances <- left_join(distances,Zh %>% select(-V1,-V2) %>% unique())
# # distances <- distances %>%mutate(NAS=ifelse(NAS<3,'Not NASH',
# #                                             ifelse(NAS<=4,'Borderline NASH',
# #                                                    'Definate NASH')))
# distances_filt_nas <- left_join(distances,sigInfo) %>% select(cmap_name,pert_type,sample,dist,NAS,fibrosis) %>% unique()
# distances_filt_nas <- distances_filt_nas %>% group_by(cmap_name,NAS) %>% mutate(per_sig_mean_nas_dist = mean(dist)) %>%
#   ungroup() %>% group_by(cmap_name) %>% mutate(min_nas_dist = min(per_sig_mean_nas_dist)) %>% ungroup() %>%
#   filter(per_sig_mean_nas_dist==min_nas_dist) %>% select(cmap_name,pert_type,NAS,per_sig_mean_nas_dist) %>% unique()
# # distances_filt_nas <- distances_filt_nas %>% group_by(NAS) %>% mutate(counts = n_distinct(sig_id)) %>% ungroup()
# distances_filt_nas <- distances_filt_nas %>% group_by(NAS) %>%
#   slice_max(order_by = per_sig_mean_nas_dist, n = 25)