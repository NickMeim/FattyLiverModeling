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

### Get L1000 hits-----------------------------------------------------------------------------------------
# combo_hallmarks <- readRDS('../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_hallmarks.rds')
# combo_gos<- readRDS('../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_gos.rds')
# combo_keggs<- readRDS('../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_keggs.rds')
# combo_tfs<- readRDS('../../../L1000_2021_11_23/combo_L1000_liverExtraBasis_tfs.rds')
sigInfo <- sigInfo %>% filter(sig_id %in% sigs_selected)
drug_ligand_ex <- drug_ligand_ex[,sigInfo$sig_id]

### Load pre-calculated GSEA-distances (which used the above combined data)
# mean_dist_hallmarks <- readRDS('../../../L1000_2021_11_23/GSEA/hallmarks_gsea_dist_l1000_with_extra_basis.rds')
# mean_dist_hallmarks <- mean_dist_hallmarks[-which(rownames(mean_dist_hallmarks) %in% c('V1','V2','V3')),c('V1','V2','V3')]
mean_dist_keggs <- readRDS('../../../L1000_2021_11_23/GSEA/keggs_gsea_dist_l1000_with_extra_basis.rds')
### Convert matrix into data frame
# Keep only unique (non-self) pairs
mean_dist_keggs[lower.tri(mean_dist_keggs,diag = T)] <- -100
dist_keggs <- reshape2::melt(mean_dist_keggs)
dist_keggs <- dist_keggs %>% filter(value != -100)
colnames(dist_keggs)[3] <- 'keggs_dist'
l1000_dist_keggs <- left_join(dist_keggs,
                        sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,duplIdentifier,dupl_counts),
                        by = c("Var1"="sig_id"))
l1000_dist_keggs <- left_join(l1000_dist_keggs,
                        sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,duplIdentifier,dupl_counts),
                        by = c("Var2"="sig_id"))
l1000_dist_keggs <- l1000_dist_keggs %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y)) %>% filter(!is.na(is_duplicate))
dist_keggs <- dist_keggs %>% filter(Var1 %in% c('V1','V2','V3')) %>% filter(Var2 %in% sigInfo$sig_id)

mean_dist_gos <- readRDS('../../../L1000_2021_11_23/GSEA/gobp_gsea_dist_l1000_with_extra_basis.rds')
### Convert matrix into data frame
# Keep only unique (non-self) pairs
mean_dist_gos[lower.tri(mean_dist_gos,diag = T)] <- -100
dist_gos <- reshape2::melt(mean_dist_gos)
dist_gos <- dist_gos %>% filter(value != -100)
colnames(dist_gos)[3] <- 'gos_dist'
l1000_dist_gos <- left_join(dist_gos,
                              sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,duplIdentifier,dupl_counts),
                              by = c("Var1"="sig_id"))
l1000_dist_gos <- left_join(l1000_dist_gos,
                              sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,duplIdentifier,dupl_counts),
                              by = c("Var2"="sig_id"))
l1000_dist_gos <- l1000_dist_gos %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y)) %>% filter(!is.na(is_duplicate))
dist_gos <- dist_gos %>% filter(Var1 %in% c('V1','V2','V3')) %>% filter(Var2 %in% sigInfo$sig_id)

mean_dist_tfs <- readRDS('../../../L1000_2021_11_23/GSEA/tfs_gsea_dist_l1000_with_extra_basis.rds')
### Convert matrix into data frame
# Keep only unique (non-self) pairs
mean_dist_tfs[lower.tri(mean_dist_tfs,diag = T)] <- -100
dist_tfs <- reshape2::melt(mean_dist_tfs)
dist_tfs <- dist_tfs %>% filter(value != -100)
colnames(dist_tfs)[3] <- 'tfs_dist'
l1000_dist_tfs <- left_join(dist_tfs,
                            sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,duplIdentifier,dupl_counts),
                            by = c("Var1"="sig_id"))
l1000_dist_tfs <- left_join(l1000_dist_tfs,
                            sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,duplIdentifier,dupl_counts),
                            by = c("Var2"="sig_id"))
l1000_dist_tfs <- l1000_dist_tfs %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y)) %>% filter(!is.na(is_duplicate))
dist_tfs <- dist_tfs %>% filter(Var1 %in% c('V1','V2','V3')) %>% filter(Var2 %in% sigInfo$sig_id)

## all together
all_dists <- left_join(dist_tfs,dist_gos)
all_dists <- left_join(all_dists,dist_keggs)
all_dists <- left_join(all_dists,sigInfo %>% select(sig_id,cmap_name,pert_idose,pert_itime,cell_iname,pert_type),by=c('Var2'='sig_id'))
all_dists <- all_dists %>% select(Var1,Var2,cmap_name,pert_idose,pert_itime,pert_type,tfs_dist,gos_dist,keggs_dist)

### Calculate dist threshold based on L1000 duplicates
l1000_dist_tfs <- l1000_dist_tfs %>%
  mutate(`same drug and cell` = (paste0(cmap_name.x,'_',cell_iname.x)==paste0(cmap_name.y,'_',cell_iname.y)))
l1000_dist_tfs <- l1000_dist_tfs %>% group_by(`same drug and cell`) %>% mutate(mu = mean(tfs_dist)) %>% mutate(std = sd(tfs_dist)) %>% ungroup()
thresh_tfs <- l1000_dist_tfs$mu[l1000_dist_tfs$`same drug and cell`][1] + 3*l1000_dist_tfs$std[l1000_dist_tfs$`same drug and cell`][1]
ggviolin(l1000_dist_tfs,
         x='same drug and cell',y='tfs_dist',fill = 'same drug and cell') +
  geom_boxplot(aes(fill = `same drug and cell`),width = 0.15,outliers = FALSE)+
  geom_hline(yintercept = thresh_tfs,
             linetype = 'dashed',color='black',lwd=1)+
  ylab('GSEA-based distance') + xlab('')+
  ggtitle('Distance distributions of same drug and cell line signatures')+
  theme(text = element_text(family='Arial',size=20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        legend.position = 'bottom')
ggsave('../results/l1000_same_drug_cell_separation_gsea_dist_tfs.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)

l1000_dist_keggs <- l1000_dist_keggs %>%
  mutate(`same drug and cell` = (paste0(cmap_name.x,'_',cell_iname.x)==paste0(cmap_name.y,'_',cell_iname.y)))
l1000_dist_keggs <- l1000_dist_keggs %>% group_by(`same drug and cell`) %>% mutate(mu = mean(keggs_dist)) %>% mutate(std = sd(keggs_dist)) %>% ungroup()
thresh_keggs <- l1000_dist_keggs$mu[l1000_dist_keggs$`same drug and cell`][1] + 3*l1000_dist_keggs$std[l1000_dist_keggs$`same drug and cell`][1]
ggviolin(l1000_dist_keggs,
         x='same drug and cell',y='keggs_dist',fill = 'same drug and cell') +
  geom_boxplot(aes(fill = `same drug and cell`),width = 0.15,outliers = FALSE)+
  geom_hline(yintercept = thresh_keggs,
             linetype = 'dashed',color='black',lwd=1)+
  ylab('GSEA-based distance') + xlab('')+
  ggtitle('Distance distributions of same drug and cell line signatures')+
  theme(text = element_text(family='Arial',size=20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        legend.position = 'bottom')
ggsave('../results/l1000_same_drug_cell_separation_gsea_dist_keggs.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)

l1000_dist_gos <- l1000_dist_gos %>%
  mutate(`same drug and cell` = (paste0(cmap_name.x,'_',cell_iname.x)==paste0(cmap_name.y,'_',cell_iname.y)))
l1000_dist_gos <- l1000_dist_gos %>% group_by(`same drug and cell`) %>% mutate(mu = mean(gos_dist)) %>% mutate(std = sd(gos_dist)) %>% ungroup()
thresh_gos <- l1000_dist_gos$mu[l1000_dist_gos$`same drug and cell`][1] + 3*l1000_dist_gos$std[l1000_dist_gos$`same drug and cell`][1]
ggviolin(l1000_dist_gos,
         x='same drug and cell',y='gos_dist',fill = 'same drug and cell') +
  geom_boxplot(aes(fill = `same drug and cell`),width = 0.15,outliers = FALSE)+
  geom_hline(yintercept = thresh_gos,
             linetype = 'dashed',color='black',lwd=1)+
  ylab('GSEA-based distance') + xlab('')+
  ggtitle('Distance distributions of same drug and cell line signatures')+
  theme(text = element_text(family='Arial',size=20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        legend.position = 'bottom')
ggsave('../results/l1000_same_drug_cell_separation_gsea_dist_gos.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)

### Keep L1000 - extra basis hits
all_dists_filt <- all_dists %>% filter(tfs_dist<thresh_tfs | tfs_dist>(2-thresh_tfs)) %>% 
  filter(keggs_dist<thresh_keggs| keggs_dist>(2-thresh_keggs)) %>%
  filter(gos_dist<thresh_gos| gos_dist>(2-thresh_gos))
### Drug-target information
egsea.data()
#jak_stat_signaling
jak_stat_signaling_targets <- kegg.pathways$human$kg.sets[grepl('Jak',names(kegg.pathways$human$kg.sets))]
jak_stat_signaling_targets <- do.call('c',jak_stat_signaling_targets)
jak_stat_signaling_targets <- unique(jak_stat_signaling_targets)
jak_stat_signaling_targets <- mapIds(org.Hs.eg.db, keys = jak_stat_signaling_targets, column = "SYMBOL", keytype = "ENTREZID")
jak_stat_signaling_targets <- unname(jak_stat_signaling_targets)
inds <- which(!is.na(jak_stat_signaling_targets))
jak_stat_signaling_targets <- jak_stat_signaling_targets[inds]
#p53_signaling
p53_signaling_targets <- kegg.pathways$human$kg.sets[grepl('p53',names(kegg.pathways$human$kg.sets))]
p53_signaling_targets <- do.call('c',p53_signaling_targets)
p53_signaling_targets <- unique(p53_signaling_targets)
p53_signaling_targets <- mapIds(org.Hs.eg.db, keys = p53_signaling_targets, column = "SYMBOL", keytype = "ENTREZID")
p53_signaling_targets <- unname(p53_signaling_targets)
inds <- which(!is.na(p53_signaling_targets))
p53_signaling_targets <- p53_signaling_targets[inds]
#nfkb_signaling
nfkb_signaling_targets <- kegg.pathways$human$kg.sets[grepl('NF-kappa B',names(kegg.pathways$human$kg.sets))]
nfkb_signaling_targets <- do.call('c',nfkb_signaling_targets)
nfkb_signaling_targets <- unique(nfkb_signaling_targets)
nfkb_signaling_targets <- mapIds(org.Hs.eg.db, keys = nfkb_signaling_targets, column = "SYMBOL", keytype = "ENTREZID")
nfkb_signaling_targets <- unname(nfkb_signaling_targets)
inds <- which(!is.na(nfkb_signaling_targets))
nfkb_signaling_targets <- nfkb_signaling_targets[inds]
## HIF-1 signaling
hif1_signaling_targets <-  kegg.pathways$human$kg.sets[grepl('HIF',names(kegg.pathways$human$kg.sets))]
hif1_signaling_targets <- do.call('c',hif1_signaling_targets)
hif1_signaling_targets <- unique(hif1_signaling_targets)
hif1_signaling_targets <- mapIds(org.Hs.eg.db, keys = hif1_signaling_targets, column = "SYMBOL", keytype = "ENTREZID")
hif1_signaling_targets <- unname(hif1_signaling_targets)
inds <- which(!is.na(hif1_signaling_targets))
hif1_signaling_targets <- hif1_signaling_targets[inds]
## Insulin signaling
insulin_signaling_targets <-  kegg.pathways$human$kg.sets[grepl('Insulin signaling',names(kegg.pathways$human$kg.sets))]
insulin_signaling_targets <- do.call('c',insulin_signaling_targets)
insulin_signaling_targets <- unique(insulin_signaling_targets)
insulin_signaling_targets <- mapIds(org.Hs.eg.db, keys = insulin_signaling_targets, column = "SYMBOL", keytype = "ENTREZID")
insulin_signaling_targets <- unname(insulin_signaling_targets)
inds <- which(!is.na(insulin_signaling_targets))
insulin_signaling_targets <- insulin_signaling_targets[inds]

pathways_info <- rbind(data.frame(pathway = rep('Insulin signaling',length(insulin_signaling_targets)),
                            target = insulin_signaling_targets),
                       data.frame(pathway = rep('HIF-1 signaling',length(hif1_signaling_targets)),
                                  target = hif1_signaling_targets),
                       data.frame(pathway = rep('NFkB signaling',length(nfkb_signaling_targets)),
                                  target = nfkb_signaling_targets),
                       data.frame(pathway = rep('P53 signaling',length(p53_signaling_targets)),
                                  target = p53_signaling_targets),
                       data.frame(pathway = rep('JAK-STAT signaling',length(jak_stat_signaling_targets)),
                                  target = jak_stat_signaling_targets))

pert_info <- read.delim('../../../L1000_2021_11_23/compoundinfo_beta.txt')
pert_info <- left_join(pert_info,pathways_info)
# pert_info <- pert_info %>% 
#   mutate(target = ifelse(target %in% c('INSR','TYK2',
#                                        jak_stat_signaling_targets[grepl('IL',jak_stat_signaling_targets)],
#                                        jak_stat_signaling_targets[grepl('STAT',jak_stat_signaling_targets)],
#                                        jak_stat_signaling_targets[grepl('SOC',jak_stat_signaling_targets)],
#                                        jak_stat_signaling_targets[grepl('PIAS',jak_stat_signaling_targets)],
#                                        "IRF9") | grepl('IFN',target) | grepl('JAK',target),target,
#                          ifelse(target %in% jak_stat_signaling_targets,'JAK-STAT signaling',
#                                 ifelse(target %in% p53_signaling_targets,'P53 sigaling',
#                                        ifelse(target %in% nfkb_signaling_targets,'NFkB sigaling','other')))))
pert_info <- pert_info %>% select(pert_id,cmap_name,target,pathway,target) %>% unique()
pert_info <- pert_info %>% mutate(pathway = ifelse(target=='AR','AR',
                                                   ifelse(target=='EGFR','EGFR',pathway)))
pert_info <- pert_info %>% filter(target!='other')
pert_info <- pert_info %>% filter(pathway!='')
pert_info <- pert_info %>% filter(!is.na(pathway))
target_info <- left_join(sigInfo %>% select(sig_id,pert_id,cmap_name,pert_type),
                         pert_info)
# filter to keep relevent experiments---
target_info <- target_info %>% filter(sig_id %in% c(all_dists_filt$Var1,all_dists_filt$Var2))
# filter to keep relevent experiments---
target_info <- target_info %>% mutate(pathway = ifelse(is.na(pathway),
                                                      ifelse(pert_type=='trt_lig',
                                                             ifelse(grepl('IL',cmap_name) | grepl('IFN',cmap_name),'JAK-STAT signaling',
                                                                    'other'),'other'),
                                                      pathway))
target_info <- target_info %>% mutate(pathway = ifelse(cmap_name=='INS','Insulin signaling',pathway))

all_dists_filt2 <-left_join(all_dists_filt,
                           target_info %>% select(-pert_id) %>% unique(),
                           by=c('Var2'='sig_id','cmap_name','pert_type'))
all_dists_filt2 <- all_dists_filt2 %>% mutate(dist = (tfs_dist+gos_dist+keggs_dist)/3) %>%
  group_by(Var1,cmap_name) %>% mutate(avg_dist = mean(dist)) %>% mutate(score=avg_dist-1) %>%
  ungroup()
all_dists_filt2 <- all_dists_filt2 %>% mutate(type=pert_type)
top_show <- all_dists_filt2 %>% group_by(Var1) %>% mutate(ranking = rank(avg_dist)) %>% mutate(max_rank = max(ranking)) %>% ungroup()
top_show <- rbind(top_show %>% filter(Var1=='V1') %>% filter(ranking>=(max_rank-5) | ranking<=5),
                  top_show %>% filter(Var1=='V2') %>% filter(ranking>=(max_rank-5) | ranking<=5),
                  top_show %>% filter(Var1=='V3') %>% filter(ranking>=(max_rank-5) | ranking<=5))
top_show <- distinct(top_show %>% select(Var1,cmap_name)) %>% mutate(significant=cmap_name)
all_dists_filt2 <- left_join(all_dists_filt2,top_show)
all_dists_filt2 <- all_dists_filt2 %>% 
  mutate(significant = ifelse(grepl('IFN',cmap_name) | grepl('IL',cmap_name) | grepl('INS',cmap_name),cmap_name,significant))

### Visualize
num_colors <- length(unique(all_dists_filt2$pathway))
palette_colors <- scales::hue_pal()(num_colors)
# tmp <- all_dists_filt2 %>%mutate(target = ifelse(target=='other',NA,target))
color_mapping <- setNames(palette_colors, unique(all_dists_filt2$pathway[all_dists_filt2$pathway!='other']))
p1 <- ggplot(all_dists_filt2 %>% filter(Var1=='V1') %>% select(cmap_name,pathway,avg_dist,significant) %>% unique() %>% 
               mutate(pathway = ifelse(pathway=='other',NA,pathway)) %>% mutate(s = ifelse(is.na(pathway),'big','small')),
             aes(x=as.numeric(reorder(cmap_name,avg_dist)),y=avg_dist)) +
  geom_point(aes(size=s),alpha=0.7,color = '#CD5C5C') +
  geom_text_repel(aes(label=significant,color = pathway),show.legend = TRUE,size=8,max.overlaps=60,box.padding = 1)+
  scale_size_manual(values = c(0.7,2))+
  xlab('Rank') + ylab('score')+
  ggtitle('L1000 perturbations for LV extra 1')+
  guides(size='none')+
  labs(color='target')+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        panel.grid.major = element_line(),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = 'right')
print(p1)
ggsave('../figures/l1000_for_extra_lv1.png',
       plot=p1,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
ggsave('../figures/l1000_for_extra_lv1.eps',
       plot=p1,
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
p2 <- ggplot(all_dists_filt2 %>% filter(Var1=='V2') %>% select(cmap_name,pathway,avg_dist,significant) %>% unique() %>% 
               mutate(pathway = ifelse(pathway=='other',NA,pathway))%>% mutate(s = ifelse(is.na(pathway),'big','small')),
             aes(x=as.numeric(reorder(cmap_name,avg_dist)),y=avg_dist)) +
        geom_point(aes(size=s),alpha=0.7,color = '#CD5C5C') +
        geom_text_repel(aes(label=significant,color = pathway),show.legend = TRUE,size=8,max.overlaps=60,box.padding = 1)+
        scale_size_manual(values = c(0.7,2))+
        xlab('Rank') + ylab('score')+
        ggtitle('L1000 perturbations for LV extra 2')+
        guides(size='none')+
        labs(color='target')+
        theme_pubr(base_family = 'Arial',base_size = 24)+
        theme(text = element_text(family = 'Arial',size=24),
              panel.grid.major = element_line(),
              # axis.ticks.x = element_blank(),
              # axis.text.x = element_blank(),
              plot.title = element_blank(),
              legend.position = 'right')
print(p2)
ggsave('../figures/l1000_for_extra_lv2.png',
       plot=p2,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
ggsave('../figures/l1000_for_extra_lv2.eps',
       plot=p2,
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
p3 <- ggplot(all_dists_filt2 %>% filter(Var1=='V3') %>% select(cmap_name,pathway,avg_dist,significant) %>% unique() %>% 
            mutate(pathway = ifelse(pathway=='other',NA,pathway))%>% mutate(s = ifelse(is.na(pathway),'big','small')),
          aes(x=as.numeric(reorder(cmap_name,avg_dist)),y=avg_dist)) +
     geom_point(aes(size=s),alpha=0.7,color = '#CD5C5C') +
     geom_text_repel(aes(label=significant,color = pathway),show.legend = TRUE,size=8,max.overlaps=60,box.padding = 1)+
     xlab('Rank') + ylab('score')+
     ggtitle('L1000 perturbations for direction of maximized variance')+
     guides(size='none')+
     labs(color='target')+
     scale_size_manual(values = c(0.7,2))+
     theme_pubr(base_family = 'Arial',base_size = 24)+
     theme(text = element_text(family = 'Arial',size=24),
           panel.grid.major = element_line(),
           # axis.ticks.x = element_blank(),
           # axis.text.x = element_blank(),
           plot.title = element_blank(),
           legend.position = 'right')
print(p3)
ggsave('../figures/l1000_for_extra_lv3.png',
       plot=p3,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
ggsave('../figures/l1000_for_extra_lv3.eps',
       plot=p3,
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
