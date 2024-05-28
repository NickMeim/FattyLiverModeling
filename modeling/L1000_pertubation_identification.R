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
  sum1 <- rowSums(M1p)
  sum2 <- rowSums(M2p)
  # Create matrices of row sums
  sum1_matrix <- matrix(sum1, nrow = length(sum1), ncol = length(sum2), byrow = FALSE)
  sum2_matrix <- matrix(sum2, nrow = length(sum1), ncol = length(sum2), byrow = TRUE)
  union_p <- sum1_matrix + sum2_matrix - intersection_p
  Jp <- intersection_p/union_p
  
  # Then calculate negatives
  M1n <- 1*(M1==(-1))
  M2n <- 1*(M2==(-1))
  intersection_n <- M1n %*% t(M2n)
  sum1 <- rowSums(M1n)
  sum2 <- rowSums(M2n)
  # Create matrices of row sums
  sum1_matrix <- matrix(sum1, nrow = length(sum1), ncol = length(sum2), byrow = FALSE)
  sum2_matrix <- matrix(sum2, nrow = length(sum1), ncol = length(sum2), byrow = TRUE)
  union_n <- sum1_matrix + sum2_matrix - intersection_n
  Jn <- intersection_n/union_n
  J <- 0.5 * (Jp + Jn)
  return(J)
}

cosine_2mats <- function(A, B) {
  # Normalize rows of A
  norm_A <- sqrt(rowSums(A^2))
  A_normalized <- A / norm_A
  
  # Normalize rows of B
  norm_B <- sqrt(rowSums(B^2))
  B_normalized <- B / norm_B
  
  # Compute cosine similarity
  similarity <- A_normalized %*% t(B_normalized)
  
  return(similarity)
}


### Load extra basis and dX from variance optimization and infer:
### TF activity, Hallmark enrichment, KEGG enrichment,GO Terms enrichment------------------------
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

### Load the L1000 dataset---------------------------------------------------
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
duplicates <- sigInfo %>% group_by(duplIdentifier) %>% mutate(counts = n()) %>% ungroup() %>% filter(counts>1)
dupl_sigs <- unique(duplicates$sig_id)
no_dupl_sigs <- unique(sigInfo$sig_id[!(sigInfo$sig_id %in% dupl_sigs)])
sameDrugCell <- sigInfo %>% group_by(cmap_name,cell_iname) %>% mutate(counts = n()) %>% ungroup() %>% filter(counts>1)
sameDrugCell_sigs <- unique(sameDrugCell$sig_id)
no_sameDrugCell_sigs <- unique(sigInfo$sig_id[!(sigInfo$sig_id %in% dupl_sigs)])
sameDrug <- sigInfo %>% group_by(cmap_name) %>% mutate(counts = n()) %>% ungroup() %>% filter(counts>1)
sameDrug_sigs <- unique(sameDrug$sig_id)
no_sameDrug_sigs <- unique(sigInfo$sig_id[!(sigInfo$sig_id %in% sameDrug_sigs)])
gc()
### Run progeny to infer pathway activity and keep only drugs that can affect these pathways in some cell line------------
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

# sigInfo <- sigInfo %>% filter(sig_id %in% sigs_selected)
# drug_ligand_ex <- drug_ligand_ex[,sigInfo$sig_id]

### Repeat with IFNA,G VS E2F targets,G2M checkpoint hallmarks---------------------------------------
halls <- c('G2M_CHECKPOINT','E2F_TARGETS','INTERFERON_ALPHA_RESPONSE','INTERFERON_GAMMA_RESPONSE')
l1000_hallmarks <- readRDS('../../../L1000_2021_11_23/GSEA/df_hallmarks_ligand_drugs_10k_filtered.rds')
l1000_hallmarks <- l1000_hallmarks %>% filter(sigID %in% sigs_selected) %>% filter(pathway %in% halls)
l1000_hallmarks_filt <- l1000_hallmarks %>% group_by(sigID) %>% mutate(min_p = min(padj)) %>%
  ungroup() %>% filter(min_p<=0.05) %>% select(-min_p)
l1000_hallmarks_filt_IFN <- l1000_hallmarks_filt %>% filter(pathway %in% halls[3:4]) %>% 
  select(-padj) %>% spread('pathway','NES') %>%
  group_by(sigID) %>% mutate(s = sign(INTERFERON_ALPHA_RESPONSE * INTERFERON_GAMMA_RESPONSE)) %>% ungroup() %>%
  gather('pathway','NES',-sigID,-s) %>% filter(s>0)
l1000_hallmarks_filt_E2FG2M <- l1000_hallmarks_filt %>% filter(pathway %in% halls[1:2]) %>% 
  select(-padj) %>% spread('pathway','NES') %>%
  group_by(sigID) %>% mutate(s = sign(E2F_TARGETS * G2M_CHECKPOINT)) %>% ungroup() %>%
  gather('pathway','NES',-sigID,-s) %>% filter(s>0)
l1000_hallmarks_filt <- l1000_hallmarks_filt %>%
  filter((sigID %in% l1000_hallmarks_filt_IFN$sigID) & (sigID %in% l1000_hallmarks_filt_E2FG2M$sigID)) %>%
  filter(pathway %in% c('E2F_TARGETS','INTERFERON_ALPHA_RESPONSE')) %>%
  select(-padj) %>% spread('pathway','NES') %>%
  group_by(sigID) %>% mutate(s = sign(E2F_TARGETS * INTERFERON_ALPHA_RESPONSE)) %>% ungroup() %>%
  gather('pathway','NES',-sigID,-s) %>% filter(s<0)
l1000_hallmarks_filt <- l1000_hallmarks %>% filter(sigID %in% l1000_hallmarks_filt$sigID)
sigs_selected <- unique(l1000_hallmarks_filt$sigID)

### save selected L1000 enrichments combined with extra basis
l1000_hallmarks <- readRDS('../../../L1000_2021_11_23/GSEA/df_hallmarks_ligand_drugs_10k_filtered.rds')
l1000_gos <- readRDS('../../../L1000_2021_11_23/GSEA/df_gobp_ligand_drugs_10k_filtered.rds')
l1000_keggs <- readRDS('../../../L1000_2021_11_23/GSEA/df_gsea_keggs_ligand_drugs_10k_filtered.rds')
dorotheaData <- read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter <- is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData <- dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
minNrOfGenes  <-  5
TF_activities <- decoupleR::run_viper(drug_ligand_ex[,sigs_selected], dorotheaData,minsize = minNrOfGenes,verbose = FALSE)
TF_activities <- TF_activities %>% select(source,condition,score)%>% spread('condition','score') %>% 
  filter(source %in% colnames(TF_activities_extra)) %>% column_to_rownames('source')
l1000_hallmarks <- l1000_hallmarks %>% filter(sigID %in% sigs_selected) %>% select(-padj) %>% spread('sigID','NES') %>% 
  filter(pathway %in% colnames(df_msig_extra)) %>% column_to_rownames('pathway')
l1000_gos <- l1000_gos %>% filter(sigID %in% sigs_selected)%>% select(-padj) %>% spread('sigID','NES') %>% 
  filter(pathway %in% colnames(df_gos_bp_extra)) %>% column_to_rownames('pathway')
l1000_keggs <- l1000_keggs %>% filter(sigID %in% sigs_selected)%>% select(-padj) %>% spread('sigID','NES') %>% 
  mutate(pathway=substr(pathway, 9, nchar(pathway)))%>% 
  filter(pathway %in% colnames(df_kegg_extra)) %>% column_to_rownames('pathway')

all(rownames(t(df_msig_extra)[rownames(l1000_hallmarks),])==rownames(l1000_hallmarks))
all(rownames(t(df_gos_bp_extra)[rownames(l1000_gos),])==rownames(l1000_gos))
all(rownames(t(df_kegg_extra)[rownames(l1000_keggs),])==rownames(l1000_keggs))
all(rownames(t(TF_activities_extra)[rownames(TF_activities),])==rownames(TF_activities))

combo_hallmarks <- as.matrix(cbind(t(df_msig_extra)[rownames(l1000_hallmarks),],l1000_hallmarks))
combo_gos <- as.matrix(cbind(t(df_gos_bp_extra)[rownames(l1000_gos),],l1000_gos))
combo_keggs <- as.matrix(cbind(t(df_kegg_extra)[rownames(l1000_keggs),],l1000_keggs))
combo_tfs <- as.matrix(cbind(t(TF_activities_extra)[rownames(TF_activities),],TF_activities))
# saveRDS(combo_hallmarks,'../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_hallmarks.rds')
# saveRDS(combo_gos,'../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_gos.rds')
# saveRDS(combo_keggs,'../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_keggs.rds')
# saveRDS(combo_tfs,'../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_tfs.rds')

### Run viper to infer TF activity from L1000 and find threshold to consider neighbors---------------------------------------
dorotheaData <- read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter <- is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData <- dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
minNrOfGenes  <-  5
TF_activities <- decoupleR::run_viper(drug_ligand_ex, dorotheaData,minsize = minNrOfGenes,verbose = FALSE)
optimals <- NULL
for (i in 1:30){
  # mat <- t(as.matrix(TF_activities %>% select(condition,source,score) %>%
  #                      filter(condition %in% c(dupl_sigs,sample(no_dupl_sigs,length(dupl_sigs)))) %>%
  #                      spread('source','score') %>% column_to_rownames('condition')))
  mat <- t(as.matrix(TF_activities %>% select(condition,source,score) %>%
                       filter(condition %in% unique(c(sameDrugCell_sigs,
                                                      sample(unique(no_sameDrugCell_sigs),
                                                             length(sameDrugCell_sigs),
                                                             replace = FALSE)))) %>%
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
  optimal_threshold_cosine_tfs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_cosine_tfs$threshold
  print(paste0('Finished iteration ',i))
  # y_pred <- 1*(df_tf_cosine$value>optimal_threshold_cosine_tfs$threshold)
  # y <- df_tf_cosine$is_same_drug_cell
  # confusion <- confusionMatrix(data=factor(y_pred,levels = c(1,0)),
  #                              reference = factor(y,levels = c(1,0)))
  # print(confusion)
  # ggviolin(df_tf_cosine,x='is_duplicate',y='value',fill = 'is_duplicate',draw_quantiles = TRUE)
  # ggviolin(df_tf_cosine,x='is_same_drug_cell',y='value',fill = 'is_same_drug_cell',draw_quantiles = TRUE)
}
# Print the optimal threshold
hist(optimals)
optimal_threshold_cosine_tfs <- mean(optimals)
print(optimal_threshold_cosine_tfs)

# repeat for Jaccard
TF_activities_sig <- TF_activities %>% mutate(score = ifelse(p_value>=0.01,0,score)) %>% 
  mutate(score = ifelse(abs(score)<=0.5,0,score)) %>% 
  mutate(sig = ifelse(score<0,-1,
                      ifelse(score>0,1,0))) %>% 
  select(c('sig_id'='condition'),c('TF'='source'),sig)
optimals <- NULL
for (i in 1:30){
  mat <- as.matrix(TF_activities_sig %>%
                     filter(sig_id %in% unique(c(sameDrugCell_sigs,
                                                    sample(unique(no_sameDrugCell_sigs),
                                                           length(sameDrugCell_sigs),
                                                           replace = FALSE)))) %>%
                       spread('TF','sig') %>% column_to_rownames('sig_id'))
  jaccard_similarity <-  jaccard(mat)
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
  # plot(roc_object)
  optimal_threshold_jaccard_tfs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_jaccard_tfs$threshold
  print(paste0('Finished iteration ',i))
}
hist(optimals)
optimal_threshold_jaccard_tfs <- mean(optimals)
print(optimal_threshold_jaccard_tfs)

### Load Hallmarks results from L1000 and find threshold to consider neighbors---------------------------------------
l1000_hallmarks <- readRDS('../../../L1000_2021_11_23/GSEA/df_hallmarks_ligand_drugs_10k_filtered.rds')
optimals <- NULL
for (i in 1:30){
  mat <- t(as.matrix(l1000_hallmarks %>% select(sigID,pathway,NES) %>% 
                       filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                   sample(unique(no_sameDrugCell_sigs),
                                                          length(sameDrugCell_sigs),
                                                          replace = FALSE)))) %>%
                       spread('pathway','NES') %>% column_to_rownames('sigID')))
  df_hallmark_cosine <- lsa::cosine(mat)
  df_hallmark_cosine[lower.tri(df_hallmark_cosine,diag = T)]<- -100
  df_hallmark_cosine <- reshape2::melt(df_hallmark_cosine)
  df_hallmark_cosine <- df_hallmark_cosine %>% filter(value != -100)
  gc()
  df_hallmark_cosine <- left_join(df_hallmark_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                            by = c("Var1"="sig_id"))
  df_hallmark_cosine <- left_join(df_hallmark_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                            by = c("Var2"="sig_id"))
  df_hallmark_cosine <- df_hallmark_cosine %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_hallmark_cosine <- df_hallmark_cosine %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_hallmark_cosine <- df_hallmark_cosine %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_hallmark_cosine$is_same_drug_cell, df_hallmark_cosine$value)
  optimal_threshold_cosine_hallmakrs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_cosine_hallmakrs$threshold
  print(paste0('Finished iteration ',i))
}
# Print the optimal threshold
hist(optimals)
optimal_threshold_cosine_hallmakrs <- mean(optimals)
print(optimal_threshold_cosine_hallmakrs)

# repeat for Jaccard
hallmarks_sig <- l1000_hallmarks %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(sigID,pathway,sig)
optimals <- NULL
for (i in 1:30){
  mat <- as.matrix(hallmarks_sig %>% 
                     filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                sample(unique(no_sameDrugCell_sigs),
                                                       length(sameDrugCell_sigs),
                                                       replace = FALSE)))) %>%
                     spread('pathway','sig') %>% column_to_rownames('sigID'))
  jaccard_similarity <-  jaccard(mat)
  jaccard_similarity[lower.tri(jaccard_similarity,diag = T)]<- -100
  df_jaccard_hallmarks <- reshape2::melt(jaccard_similarity)
  df_jaccard_hallmarks <- df_jaccard_hallmarks %>% filter(value != -100)
  gc()
  df_jaccard_hallmarks <- left_join(df_jaccard_hallmarks,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                              by = c("Var1"="sig_id"))
  df_jaccard_hallmarks <- left_join(df_jaccard_hallmarks,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                              by = c("Var2"="sig_id"))
  df_jaccard_hallmarks <- df_jaccard_hallmarks %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_jaccard_hallmarks <- df_jaccard_hallmarks %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_jaccard_hallmarks <- df_jaccard_hallmarks %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_jaccard_hallmarks$is_same_drug_cell, df_jaccard_hallmarks$value)
  # plot(roc_object)
  optimal_threshold_jaccard_hallmarks <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_jaccard_hallmarks$threshold
  print(paste0('Finished iteration ',i))
}
hist(optimals)
optimal_threshold_jaccard_hallmarks <- mean(optimals)
print(optimal_threshold_jaccard_hallmarks)


### Load KEGGS results from L1000 and find threshold to consider neighbors---------------------------------------
l1000_keggs <- readRDS('../../../L1000_2021_11_23/GSEA/df_gsea_keggs_ligand_drugs_10k_filtered.rds')
l1000_keggs <- l1000_keggs %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))
optimals <- NULL
for (i in 1:30){
  mat <- t(as.matrix(l1000_keggs %>% select(sigID,pathway,NES) %>% 
                       filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                  sample(unique(no_sameDrugCell_sigs),
                                                         length(sameDrugCell_sigs),
                                                         replace = FALSE)))) %>%
                       spread('pathway','NES') %>% column_to_rownames('sigID')))
  df_kegg_cosine <- lsa::cosine(mat)
  df_kegg_cosine[lower.tri(df_kegg_cosine,diag = T)]<- -100
  df_kegg_cosine <- reshape2::melt(df_kegg_cosine)
  df_kegg_cosine <- df_kegg_cosine %>% filter(value != -100)
  gc()
  df_kegg_cosine <- left_join(df_kegg_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                  by = c("Var1"="sig_id"))
  df_kegg_cosine <- left_join(df_kegg_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                  by = c("Var2"="sig_id"))
  df_kegg_cosine <- df_kegg_cosine %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_kegg_cosine <- df_kegg_cosine %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_kegg_cosine <- df_kegg_cosine %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_kegg_cosine$is_same_drug_cell, df_kegg_cosine$value)
  optimal_threshold_cosine_keggs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_cosine_keggs$threshold
  print(paste0('Finished iteration ',i))
}
# Print the optimal threshold
hist(optimals)
optimal_threshold_cosine_keggs <- mean(optimals)
print(optimal_threshold_cosine_keggs)

# repeat for Jaccard
kegg_sig <- l1000_keggs %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(sigID,pathway,sig)
optimals <- NULL
for (i in 1:30){
  mat <- as.matrix(kegg_sig %>% 
                     filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                sample(unique(no_sameDrugCell_sigs),
                                                       length(sameDrugCell_sigs),
                                                       replace = FALSE)))) %>%
                     spread('pathway','sig') %>% column_to_rownames('sigID'))
  df_jaccard_keggs <-  jaccard(mat)
  df_jaccard_keggs[lower.tri(df_jaccard_keggs,diag = T)]<- -100
  df_jaccard_keggs <- reshape2::melt(df_jaccard_keggs)
  df_jaccard_keggs <- df_jaccard_keggs %>% filter(value != -100)
  gc()
  df_jaccard_keggs <- left_join(df_jaccard_keggs,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                    by = c("Var1"="sig_id"))
  df_jaccard_keggs <- left_join(df_jaccard_keggs,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                    by = c("Var2"="sig_id"))
  df_jaccard_keggs <- df_jaccard_keggs %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_jaccard_keggs <- df_jaccard_keggs %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_jaccard_keggs <- df_jaccard_keggs %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_jaccard_keggs$is_same_drug_cell, df_jaccard_keggs$value)
  # plot(roc_object)
  optimal_threshold_jaccard_keggs <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_jaccard_keggs$threshold
  print(paste0('Finished iteration ',i))
}
hist(optimals)
optimal_threshold_jaccard_keggs <- mean(optimals)
print(optimal_threshold_jaccard_keggs)

### Load GO BP results from L1000 and find threshold to consider neighbors---------------------------------------
l1000_gos_bp <- readRDS('../../../L1000_2021_11_23/GSEA/df_gobp_ligand_drugs_10k_filtered.rds')
optimals <- NULL
for (i in 1:30){
  mat <- t(as.matrix(l1000_gos_bp %>% select(sigID,pathway,NES) %>% 
                       filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                  sample(unique(no_sameDrugCell_sigs),
                                                         length(sameDrugCell_sigs),
                                                         replace = FALSE)))) %>%
                       spread('pathway','NES') %>% column_to_rownames('sigID')))
  df_gos_bp_cosine <- lsa::cosine(mat)
  df_gos_bp_cosine[lower.tri(df_gos_bp_cosine,diag = T)]<- -100
  df_gos_bp_cosine <- reshape2::melt(df_gos_bp_cosine)
  df_gos_bp_cosine <- df_gos_bp_cosine %>% filter(value != -100)
  gc()
  df_gos_bp_cosine <- left_join(df_gos_bp_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                              by = c("Var1"="sig_id"))
  df_gos_bp_cosine <- left_join(df_gos_bp_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                              by = c("Var2"="sig_id"))
  df_gos_bp_cosine <- df_gos_bp_cosine %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_gos_bp_cosine <- df_gos_bp_cosine %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_gos_bp_cosine <- df_gos_bp_cosine %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_gos_bp_cosine$is_same_drug_cell, df_gos_bp_cosine$value)
  optimal_threshold_cosine_gos_bp <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_cosine_gos_bp$threshold
  print(paste0('Finished iteration ',i))
}
# Print the optimal threshold
hist(optimals)
optimal_threshold_cosine_gos_bp <- mean(optimals)
print(optimal_threshold_cosine_gos_bp)

# repeat for Jaccard
gos_bp_sig <- l1000_gos_bp %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(sigID,pathway,sig)
optimals <- NULL
for (i in 1:30){
  mat <- as.matrix(gos_bp_sig %>% 
                     filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                sample(unique(no_sameDrugCell_sigs),
                                                       length(sameDrugCell_sigs),
                                                       replace = FALSE)))) %>%
                     spread('pathway','sig') %>% column_to_rownames('sigID'))
  df_jaccard_gos_bp <-  jaccard(mat)
  df_jaccard_gos_bp[lower.tri(df_jaccard_gos_bp,diag = T)]<- -100
  df_jaccard_gos_bp <- reshape2::melt(df_jaccard_gos_bp)
  df_jaccard_gos_bp <- df_jaccard_gos_bp %>% filter(value != -100)
  gc()
  df_jaccard_gos_bp <- left_join(df_jaccard_gos_bp,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                by = c("Var1"="sig_id"))
  df_jaccard_gos_bp <- left_join(df_jaccard_gos_bp,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                by = c("Var2"="sig_id"))
  df_jaccard_gos_bp <- df_jaccard_gos_bp %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_jaccard_gos_bp <- df_jaccard_gos_bp %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_jaccard_gos_bp <- df_jaccard_gos_bp %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_jaccard_gos_bp$is_same_drug_cell, df_jaccard_gos_bp$value)
  # plot(roc_object)
  optimal_threshold_jaccard_gos_bp <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_jaccard_gos_bp$threshold
  print(paste0('Finished iteration ',i))
}
hist(optimals)
optimal_threshold_jaccard_gos_bp <- mean(optimals)
print(optimal_threshold_jaccard_gos_bp)
gc()

### Load GO CC results from L1000 and find threshold to consider neighbors---------------------------------------
l1000_gos_cc <- readRDS('../../../L1000_2021_11_23/GSEA/df_gocc_ligand_drugs_10k_filtered.rds')
optimals <- NULL
for (i in 1:30){
  mat <- t(as.matrix(l1000_gos_cc %>% select(sigID,pathway,NES) %>% 
                       filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                  sample(unique(no_sameDrugCell_sigs),
                                                         length(sameDrugCell_sigs),
                                                         replace = FALSE)))) %>%
                       spread('pathway','NES') %>% column_to_rownames('sigID')))
  df_gos_bp_cosine <- lsa::cosine(mat)
  df_gos_bp_cosine[lower.tri(df_gos_bp_cosine,diag = T)]<- -100
  df_gos_bp_cosine <- reshape2::melt(df_gos_bp_cosine)
  df_gos_bp_cosine <- df_gos_bp_cosine %>% filter(value != -100)
  gc()
  df_gos_bp_cosine <- left_join(df_gos_bp_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                              by = c("Var1"="sig_id"))
  df_gos_bp_cosine <- left_join(df_gos_bp_cosine,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                              by = c("Var2"="sig_id"))
  df_gos_bp_cosine <- df_gos_bp_cosine %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_gos_bp_cosine <- df_gos_bp_cosine %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_gos_bp_cosine <- df_gos_bp_cosine %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_gos_bp_cosine$is_same_drug_cell, df_gos_bp_cosine$value)
  optimal_threshold_cosine_gos_cc <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_cosine_gos_cc$threshold
  print(paste0('Finished iteration ',i))
}
# Print the optimal threshold
hist(optimals)
optimal_threshold_cosine_gos_cc <- mean(optimals)
print(optimal_threshold_cosine_gos_cc)

# repeat for Jaccard
gos_cc_sig <- l1000_gos_cc %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(sigID,pathway,sig)
optimals <- NULL
for (i in 1:30){
  mat <- as.matrix(gos_cc_sig %>% 
                     filter(sigID %in% unique(c(sameDrugCell_sigs,
                                                sample(unique(no_sameDrugCell_sigs),
                                                       length(sameDrugCell_sigs),
                                                       replace = FALSE)))) %>%
                     spread('pathway','sig') %>% column_to_rownames('sigID'))
  df_jaccard_gos_cc <-  jaccard(mat)
  df_jaccard_gos_cc[lower.tri(df_jaccard_gos_cc,diag = T)]<- -100
  df_jaccard_gos_cc <- reshape2::melt(df_jaccard_gos_cc)
  df_jaccard_gos_cc <- df_jaccard_gos_cc %>% filter(value != -100)
  gc()
  df_jaccard_gos_cc <- left_join(df_jaccard_gos_cc,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                by = c("Var1"="sig_id"))
  df_jaccard_gos_cc <- left_join(df_jaccard_gos_cc,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier),
                                by = c("Var2"="sig_id"))
  df_jaccard_gos_cc <- df_jaccard_gos_cc %>% mutate(is_duplicate = 1*(duplIdentifier.x==duplIdentifier.y))
  df_jaccard_gos_cc <- df_jaccard_gos_cc %>% mutate(is_same_drug = 1*(cmap_name.x==cmap_name.y))
  df_jaccard_gos_cc <- df_jaccard_gos_cc %>% mutate(is_same_drug_cell = 1*(cmap_name.x==cmap_name.y)*(cell_iname.x==cell_iname.y))
  gc()
  # ROC analysis
  roc_object <- roc(df_jaccard_gos_cc$is_same_drug_cell, df_jaccard_gos_cc$value)
  # plot(roc_object)
  optimal_threshold_jaccard_gos_cc <- coords(roc_object, "best", ret = "threshold", best.method = "youden")
  optimals[i] <- optimal_threshold_jaccard_gos_cc$threshold
  print(paste0('Finished iteration ',i))
}
hist(optimals)
optimal_threshold_jaccard_gos_cc <- mean(optimals)
print(optimal_threshold_jaccard_gos_cc)

### Save all L1000 thresholds--------------
# saveRDS(optimal_threshold_cosine_tfs,'../results/optimal_threshold_cosine_tfs_l1000.rds')
# saveRDS(optimal_threshold_cosine_keggs,'../results/optimal_threshold_cosine_keggs_l1000.rds')
# saveRDS(optimal_threshold_cosine_hallmakrs,'../results/optimal_threshold_cosine_hallmakrs_l1000.rds')
# saveRDS(optimal_threshold_cosine_gos_bp,'../results/optimal_threshold_cosine_gos_bp_l1000.rds')
# saveRDS(optimal_threshold_cosine_gos_cc,'../results/optimal_threshold_cosine_gos_cc_l1000.rds')
# 
# saveRDS(optimal_threshold_jaccard_tfs,'../results/optimal_threshold_jaccard_tfs_l1000.rds')
# saveRDS(optimal_threshold_jaccard_keggs,'../results/optimal_threshold_jaccard_keggs_l1000.rds')
# saveRDS(optimal_threshold_jaccard_hallmarks,'../results/optimal_threshold_jaccard_hallmakrs_l1000.rds')
# saveRDS(optimal_threshold_jaccard_gos_bp,'../results/optimal_threshold_jaccard_gos_bp_l1000.rds')
# saveRDS(optimal_threshold_jaccard_gos_cc,'../results/optimal_threshold_jaccard_gos_cc_l1000.rds')

### Get L1000 hits-----------------
sigInfo <- sigInfo %>% filter(sig_id %in% sigs_selected)
drug_ligand_ex <- drug_ligand_ex[,sigInfo$sig_id]

## keep only data from L1000 perturbations that affect these pathways
l1000_keggs <- l1000_keggs %>% filter(sigID %in% sigs_selected)
l1000_hallmarks <- l1000_hallmarks %>% filter(sigID %in% sigs_selected)
l1000_gos_bp <- l1000_gos_bp %>% filter(sigID %in% sigs_selected)
l1000_gos_cc <- l1000_gos_cc %>% filter(sigID %in% sigs_selected)
TF_activities <- TF_activities %>% filter(condition %in% sigs_selected)
gc()

# saveRDS(l1000_keggs,'../../archived github/l1000_keggs_filtered.rds')
# saveRDS(l1000_hallmarks,'../../archived github/l1000_hallmarks_filtered.rds')
# saveRDS(l1000_gos_bp,'../../archived github/l1000_gos_pb_filtered.rds')
# saveRDS(l1000_gos_cc,'../../archived github/l1000_gos_cc_filtered.rds')
# saveRDS(TF_activities,'../../archived github/TF_activities_filtered.rds')

# Calculate similarity with extra basis and optimized variance
# Remember direction of similarity does not matter !
# Even reverse down,up elements or cosine -1 are fine for the purpose of this exercise !
# So the found threshold is for the same exact direction but also:
# cosine less than -threshold is ok and jaccard similarity can be calculated with the reversed bottom-up elements
optimal_threshold_cosine_tfs <- readRDS('../results/optimal_threshold_cosine_tfs_l1000.rds')
optimal_threshold_cosine_keggs <- readRDS('../results/optimal_threshold_cosine_keggs_l1000.rds')
optimal_threshold_cosine_hallmakrs <- readRDS('../results/optimal_threshold_cosine_hallmakrs_l1000.rds')
optimal_threshold_cosine_gos_bp <- readRDS('../results/optimal_threshold_cosine_gos_bp_l1000.rds')
optimal_threshold_cosine_gos_cc <- readRDS('../results/optimal_threshold_cosine_gos_cc_l1000.rds')
optimal_threshold_jaccard_tfs <- readRDS('../results/optimal_threshold_jaccard_tfs_l1000.rds')
optimal_threshold_jaccard_keggs <- readRDS('../results/optimal_threshold_jaccard_keggs_l1000.rds')
optimal_threshold_jaccard_hallmarks <- readRDS('../results/optimal_threshold_jaccard_hallmakrs_l1000.rds')
optimal_threshold_jaccard_gos_bp <- readRDS('../results/optimal_threshold_jaccard_gos_bp_l1000.rds')
optimal_threshold_jaccard_gos_cc <- readRDS('../results/optimal_threshold_jaccard_gos_cc_l1000.rds')
# l1000_keggs <- readRDS('../../archived github/l1000_keggs_filtered.rds')
# l1000_hallmarks <- readRDS('../../archived github/l1000_hallmarks_filtered.rds')
# l1000_gos_bp <- readRDS('../../archived github/l1000_gos_pb_filtered.rds')
# l1000_gos_cc <- readRDS('../../archived github/l1000_gos_cc_filtered.rds')
# TF_activities <- readRDS('../../archived github/TF_activities_filtered.rds')

### First TF activity similarities
### Calculate Jaccard similarity
TF_activities_sig <- TF_activities %>% mutate(score = ifelse(p_value>=0.01,0,score)) %>% 
  mutate(score = ifelse(abs(score)<=0.5,0,score)) %>% 
  mutate(sig = ifelse(score<0,-1,
                      ifelse(score>0,1,0))) %>% 
  select(c('sig_id'='condition'),c('TF'='source'),sig)
M_L1000 <- as.matrix(TF_activities_sig %>%  spread('TF','sig') %>% column_to_rownames('sig_id'))
M2 <- TF_activities_extra_filt[,colnames(M_L1000)]
print(all(colnames(M_L1000)==colnames(M2)))
jaccard_sim <- jaccard_2mats(M_L1000,M2)
jaccard_sim_reversed <- jaccard_2mats(-1*M_L1000,M2)
jaccard_sim_reversed <- as.data.frame(jaccard_sim_reversed) %>% rownames_to_column('sig_id') %>% gather('LV','jaccard',-sig_id) 
jaccard_sim <- as.data.frame(jaccard_sim) %>% rownames_to_column('sig_id') %>% gather('LV','jaccard',-sig_id) 
jaccard_sim <- left_join(jaccard_sim,jaccard_sim_reversed,by=c('sig_id','LV'))
jaccard_sim <- jaccard_sim %>% group_by(sig_id,LV) %>% mutate(jaccard= max(c(jaccard.x,jaccard.y))) %>% ungroup()
jaccard_sim <- jaccard_sim %>% select(sig_id,LV,jaccard)
### Calculate cosine similarity
M_L1000 <- as.matrix(TF_activities %>% select(condition,source,score) %>% spread('source','score') %>% 
                       column_to_rownames('condition'))
M2 <- TF_activities_extra[,colnames(M_L1000)]
print(all(colnames(M_L1000)==colnames(M2)))
cos_sim <- cosine_2mats(M_L1000,M2)
cos_sim <- as.data.frame(cos_sim) %>% rownames_to_column('sig_id') %>% gather('LV','cosine',-sig_id)
cos_sim <- left_join(cos_sim,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,pert_type))

TF_similarities <- left_join(cos_sim,jaccard_sim,by=c('sig_id','LV'))

### Hallmarks similarities
### Calculate Jaccard similarity
hallmarks_sig <- l1000_hallmarks %>% mutate(NES = ifelse(padj>=0.05,0,NES)) %>% 
  mutate(NES = ifelse(abs(NES)<=0.5,0,NES)) %>% 
  mutate(sig = ifelse(NES<0,-1,
                      ifelse(NES>0,1,0))) %>% 
  select(c('sig_id'='sigID'),c('geneset'='pathway'),sig)
M_L1000 <- as.matrix(hallmarks_sig %>%  spread('geneset','sig') %>% column_to_rownames('sig_id'))
M2 <- df_msig_extra_filt[,colnames(M_L1000)]
print(all(colnames(M_L1000)==colnames(M2)))
jaccard_sim <- jaccard_2mats(M_L1000,M2)
jaccard_sim_reversed <- jaccard_2mats(-1*M_L1000,M2)
jaccard_sim_reversed <- as.data.frame(jaccard_sim_reversed) %>% rownames_to_column('sig_id') %>% gather('LV','jaccard',-sig_id) 
jaccard_sim <- as.data.frame(jaccard_sim) %>% rownames_to_column('sig_id') %>% gather('LV','jaccard',-sig_id) 
jaccard_sim <- left_join(jaccard_sim,jaccard_sim_reversed,by=c('sig_id','LV'))
jaccard_sim <- jaccard_sim %>% group_by(sig_id,LV) %>% mutate(jaccard= max(c(jaccard.x,jaccard.y))) %>% ungroup()
jaccard_sim <- jaccard_sim %>% select(sig_id,LV,jaccard)
### Calculate cosine similarity
M_L1000 <- as.matrix(l1000_hallmarks %>% select(sigID,pathway,NES) %>% spread('pathway','NES') %>% 
                       column_to_rownames('sigID'))
M2 <- df_msig_extra[,colnames(M_L1000)]
print(all(colnames(M_L1000)==colnames(M2)))
cos_sim <- cosine_2mats(M_L1000,M2)
cos_sim <- as.data.frame(cos_sim) %>% rownames_to_column('sig_id') %>% gather('LV','cosine',-sig_id)
cos_sim <- left_join(cos_sim,sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,pert_type))
