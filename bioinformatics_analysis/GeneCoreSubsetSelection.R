library(AnnotationDbi)
library(tidyverse)
library(fgsea)
library(gage)
library(EGSEAdata)
library(org.Hs.eg.db)
library(cmapR)
library(rhdf5)
library(topGO)
library(GO.db)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)
source('enrichment_calculations.R')
egsea.data(species = "human",returnInfo = TRUE)
library(igraph)

geos <- c("GSE126848","GSE130970","GSE134422","GSE135251","GSE162694")
all_significant_keggs <- readRDS('../results/bioinformatics_results/all_significant_keggs.rds')

### Initialize KEGGS, TFs, msigdb genesets annotation and GO Terms annotation-----------------------
all_msig_genesets <- NULL
msig_overview <- t(as.data.frame(msigdb[[1]]))
for (i in 2:(length(msigdb)-1)) {
  msig_overview <- rbind(msig_overview,t(as.data.frame(msigdb[[i]])))
}
msig_overview <- as.data.frame(msig_overview)
rownames(msig_overview) <- msig_overview$STANDARD_NAME
hallmark <- list()
c1 <- list()
c2 <- list()
c3 <- list()
c4 <- list()
c5 <- list()
c6 <- list()
c7 <- list()
h <- 1
c1_ind <- 1
c2_ind <- 1
c3_ind <- 1
c4_ind <- 1
c5_ind <- 1
c6_ind <- 1
c7_ind <- 1
for (i in 1:(length(msigdb)-1)) {
  a <- msigdb[[i]]
  
  if (a["CATEGORY_CODE"] == "H") {
    nm <- a["STANDARD_NAME"]
    hallmark <- c(hallmark,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(hallmark)[h] <- substr(nm,nchar('HALLMARK_')+1,nchar(nm))
    h <- h+1
  }
  if (a["CATEGORY_CODE"] == "C1") {
    nm <- a["STANDARD_NAME"]
    c1 <- c(c1,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c1)[c1_ind] <- nm
    c1_ind <- c1_ind+1
  }
  if (a["CATEGORY_CODE"] == "C2") {
    nm <- a["STANDARD_NAME"]
    c2 <- c(c2,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c2)[c2_ind] <- nm
    c2_ind <- c2_ind+1
  }
  if (a["CATEGORY_CODE"] == "C3") {
    nm <- a["STANDARD_NAME"]
    c3 <- c(c3,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c3)[c3_ind] <- nm
    c3_ind <- c3_ind+1
  }
  if (a["CATEGORY_CODE"] == "C4") {
    nm <- a["STANDARD_NAME"]
    c4 <- c(c4,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c4)[c4_ind] <- nm
    c4_ind <- c4_ind+1
  }
  if (a["CATEGORY_CODE"] == "C5") {
    nm <- a["STANDARD_NAME"]
    c5 <- c(c5,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c5)[c5_ind] <- nm
    c5_ind <- c5_ind+1
  }
  if (a["CATEGORY_CODE"] == "C6") {
    nm <- a["STANDARD_NAME"]
    c6 <- c(c6,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c6)[c6_ind] <- nm
    c6_ind <- c6_ind+1
  }
  if (a["CATEGORY_CODE"] == "C7") {
    nm <- a["STANDARD_NAME"]
    c7 <- c(c7,str_split(a["MEMBERS_EZID"],pattern = ","))
    names(c7)[c7_ind] <- nm
    c7_ind <- c7_ind+1
  }
}
all_msig_genesets <- c(all_msig_genesets,hallmark)
all_msig_genesets <- c(all_msig_genesets,c1)
all_msig_genesets <- c(all_msig_genesets,c2)
all_msig_genesets <- c(all_msig_genesets,c3)
all_msig_genesets <- c(all_msig_genesets,c4)
all_msig_genesets <- c(all_msig_genesets,c5)
all_msig_genesets <- c(all_msig_genesets,c6)
all_msig_genesets <- c(all_msig_genesets,c7)

### Load appropriate regulon
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
minNrOfGenes  <-  5

### Load KEGGs annotation
pathways <- kegg.pathways$human
pathways <- pathways$kg.sets
new_names <- NULL
for (nm in names(pathways)){
  new_names <- c(new_names,substr(nm,10,nchar(nm)))
}
names(pathways) <- new_names

## Load all DESeq results-------------------------------------
# df_all <- data.frame()
df_msigdbs_all <- data.frame()
df_gos_all <- data.frame()
df_tfs <- data.frame()
counter <- 1
df_sig_all <- data.frame()
all_genes <- NULL
for (geo in geos){
  list_deseq<- readRDS(paste0('../results/bioinformatics_results/',geo,'DESeqDF.rds'))
  for (i in 1:length(list_deseq)){
    if (is.null(list_deseq[[i]]$gene_names)){
      df_sig <- as.data.frame(list_deseq[[i]]) %>% rownames_to_column('gene') %>% filter(abs(log2FoldChange)>1.5) %>% filter(padj<0.05) %>% 
        dplyr::select(gene) %>% unique()
      genes <- as.data.frame(list_deseq[[i]]) %>% rownames_to_column('gene')
      genes <- unique(genes$gene)
    } else{
      df_sig <- as.data.frame(list_deseq[[i]]) %>% filter(abs(log2FoldChange)>1.5) %>% filter(padj<0.05) %>% 
        dplyr::select(c('gene'='gene_names')) %>% unique()
      genes <- as.data.frame(list_deseq[[i]])$gene_names
      genes <- unique(genes)
    }
    df_sig <- df_sig %>% mutate(dataset = geo)
    df_sig_all <- rbind(df_sig_all,df_sig) %>% unique()
    all_genes <- c(all_genes,genes)
    
    df <- list_deseq[[i]] %>% select(stat)
    colnames(df) <- names(list_deseq)[i]
    if (i==1){
      df_deseq <- df
    }else{
      df_deseq <- cbind(df_deseq,df)
    }
  }
  if (is.null(list_deseq[[1]]$gene_names)){
    rownames(df_deseq) <- rownames(list_deseq[[1]])
  }else{
    rownames(df_deseq) <-  list_deseq[[1]]$gene_names
  }
  ### Infer TF activities
  TF_activities = decoupleR::run_viper(df_deseq, dorotheaData,minsize = minNrOfGenes,verbose = FALSE)
  TF_activities <- left_join(TF_activities %>% filter(p_value<0.01) %>% filter(abs(score)>0.5) %>% 
    select(source) %>% unique(),dorotheaData) %>% mutate(dataset = geo)
  df_tfs <- rbind(df_tfs,TF_activities)
  
  ## do go enrichment
  go_enrichment <- fastenrichment(colnames(df_deseq),
                                  rownames(df_deseq),
                                  as.matrix(df_deseq),
                                  enrichment_space = c("go_bp","go_cc","go_mf"),
                                  gene_id_type="symbol",
                                  n_permutations=5000)
  
  if (ncol(df_deseq)<2) {
    go_nes <- rbind(as.matrix(go_enrichment$NES$`NES GO BP`),
                    as.matrix(go_enrichment$NES$`NES GO MF`),
                    as.matrix(go_enrichment$NES$`NES GO CC`))
    colnames(go_nes) <- 'NES'
    go_nes <- as.data.frame(go_nes) %>% rownames_to_column('GO') %>% 
      mutate(category = ifelse(grepl('BP',GO),'BP',
                               ifelse(grepl('CC',GO),'CC','MF'))) %>% 
      mutate(GO=substr(GO, nchar('FL1000_GO_BP_')+1, nchar(GO)))
    go_pval <- rbind(as.matrix(go_enrichment$Pval$`Pval GO BP`),
                     as.matrix(go_enrichment$Pval$`Pval GO MF`),
                     as.matrix(go_enrichment$Pval$`Pval GO CC`))
    colnames(go_pval) <- 'padj'
    go_pval <- as.data.frame(go_pval) %>% rownames_to_column('GO') %>%
      mutate(category = ifelse(grepl('BP',GO),'BP',
                               ifelse(grepl('CC',GO),'CC','MF'))) %>% 
      mutate(GO=substr(GO, nchar('FL1000_GO_BP_')+1, nchar(GO)))
  } else{
    go_nes <- do.call(rbind,go_enrichment$NES)
    go_nes <- as.data.frame(go_nes) %>% rownames_to_column('GO') %>% gather('condition','NES',-GO) %>%
      mutate(category = ifelse(grepl('BP',GO),'BP',
                               ifelse(grepl('CC',GO),'CC','MF'))) %>% 
      mutate(GO=substr(GO, nchar('FL1000_GO_BP_')+1, nchar(GO)))
    go_pval <- do.call(rbind,go_enrichment$Pval)
    go_pval <- as.data.frame(go_pval) %>% rownames_to_column('GO') %>% gather('condition','padj',-GO) %>%
      mutate(category = ifelse(grepl('BP',GO),'BP',
                               ifelse(grepl('CC',GO),'CC','MF'))) %>% 
      mutate(GO=substr(GO, nchar('FL1000_GO_BP_')+1, nchar(GO)))
  }
  df_gos <- left_join(go_nes,go_pval)
  df_gos <- df_gos %>% filter(padj<0.05) %>% filter(abs(NES)>0.5)
  df_gos <- df_gos %>% select(c('gene_set'='GO')) %>% mutate(dataset=geo) %>% unique()
  ###Repeat for msigdb genesets
  entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(df_deseq), column = "ENTREZID", keytype = "SYMBOL")
  entrez_ids <- unname(entrez_ids)
  inds <- which(!is.na(entrez_ids))
  entrez_ids <- entrez_ids[inds]
  meas <- as.matrix(df_deseq[inds,])
  colnames(meas) <- colnames(df_deseq)
  rownames(meas) <- entrez_ids
  msigdb_enrichment <- fastenrichment(colnames(df_deseq),
                                      entrez_ids,
                                      meas,
                                      enrichment_space = c("msig_db_h","msig_db_c1",
                                                           "msig_db_c2","msig_db_c3",
                                                           "msig_db_c4","msig_db_c5",
                                                           "msig_db_c6","msig_db_c7"),
                                      gene_id_type="entrez",
                                      n_permutations=5000)
  if (ncol(df_deseq)<2) {
    msigdb_nes <- rbind(as.matrix(msigdb_enrichment$NES$`NES MSIG Hallmark`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C1`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C2`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C3`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C4`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C5`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C6`),
                    as.matrix(msigdb_enrichment$NES$`NES MSIG C7`))
    colnames(msigdb_nes) <- 'NES'
    msigdb_nes <- as.data.frame(msigdb_nes) %>% rownames_to_column('msigdb') %>%
      mutate(category = ifelse(grepl('FL1000_MSIG_C',msigdb),substr(msigdb,nchar('FL1000_MSIG_C'),nchar('FL1000_MSIG_C')+1),
                               'HALLMARKS'))  %>% 
      mutate(msigdb=ifelse(grepl('FL1000_MSIG_C',msigdb),
                           substr(msigdb, nchar('FL1000_MSIG_C1_')+1, nchar(msigdb)),
                           substr(msigdb, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(msigdb))))
    msigdb_pval <- rbind(as.matrix(msigdb_enrichment$Pval$`Pval MSIG Hallmark`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C1`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C2`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C3`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C4`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C5`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C6`),
                        as.matrix(msigdb_enrichment$Pval$`Pval MSIG C7`))
    colnames(msigdb_pval) <- 'padj'
    msigdb_pval <- as.data.frame(msigdb_pval) %>% rownames_to_column('msigdb') %>%
      mutate(category = ifelse(grepl('FL1000_MSIG_C',msigdb),substr(msigdb,nchar('FL1000_MSIG_C'),nchar('FL1000_MSIG_C')+1),
                               'HALLMARKS'))  %>% 
      mutate(msigdb=ifelse(grepl('FL1000_MSIG_C',msigdb),
                           substr(msigdb, nchar('FL1000_MSIG_C1_')+1, nchar(msigdb)),
                           substr(msigdb, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(msigdb))))
  } else{
    msigdb_nes <- do.call(rbind,msigdb_enrichment$NES)
    msigdb_nes <- as.data.frame(msigdb_nes) %>% rownames_to_column('msigdb') %>% gather('condition','NES',-msigdb) %>%
      mutate(category = ifelse(grepl('FL1000_MSIG_C',msigdb),substr(msigdb,nchar('FL1000_MSIG_C'),nchar('FL1000_MSIG_C')+1),
                               'HALLMARKS'))  %>% 
      mutate(msigdb=ifelse(grepl('FL1000_MSIG_C',msigdb),
                           substr(msigdb, nchar('FL1000_MSIG_C1_')+1, nchar(msigdb)),
                           substr(msigdb, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(msigdb))))
    msigdb_pval <- do.call(rbind,msigdb_enrichment$Pval)
    msigdb_pval <- as.data.frame(msigdb_pval) %>% rownames_to_column('msigdb') %>% gather('condition','padj',-msigdb) %>%
      mutate(category = ifelse(grepl('FL1000_MSIG_C',msigdb),substr(msigdb,nchar('FL1000_MSIG_C'),nchar('FL1000_MSIG_C')+1),
                               'HALLMARKS'))  %>% 
      mutate(msigdb=ifelse(grepl('FL1000_MSIG_C',msigdb),
                           substr(msigdb, nchar('FL1000_MSIG_C1_')+1, nchar(msigdb)),
                           substr(msigdb, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(msigdb))))
  }
  df_msigdbs <- left_join(msigdb_nes,msigdb_pval)
  df_msigdbs <- df_msigdbs %>% filter(padj<0.05) %>% filter(abs(NES)>0.5)
  df_msigdbs <- df_msigdbs %>% select(c('gene_set'='msigdb')) %>% mutate(dataset=geo) %>% unique()

  #combine
  df_msigdbs_all <- rbind(df_msigdbs_all,df_msigdbs)
  df_gos_all <- rbind(df_gos_all,df_gos)
  
  # saveRDS(df_msigdbs_all,'../results/bioinformatics_results/all_significant_msigs.rds')
  # saveRDS(df_gos_all,'../results/bioinformatics_results/all_significant_gos.rds')
  # saveRDS(df_tfs,'../results/bioinformatics_results/all_significant_tfs.rds')
  # saveRDS(df_sig_all,'../results/bioinformatics_results/all_significant_genes.rds')
  
  print(paste0('Finised dataset ',counter,' out of ',length(geos)))
  counter <- counter+1
}
all_genes <- unique(all_genes)

## Filter genesets and genes that appear in multiple datasets-------------
df_msigdbs_all <- readRDS('../results/bioinformatics_results/all_significant_msigs.rds')
df_gos_all <- readRDS('../results/bioinformatics_results/all_significant_gos.rds')
df_tfs <- readRDS('../results/bioinformatics_results/all_significant_tfs.rds')
df_sig_all <- readRDS('../results/bioinformatics_results/all_significant_genes.rds')
df_kegg_filt <- all_significant_keggs %>% group_by(pathway) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=3) %>% 
  select(-dataset,-counts) %>% unique()
df_gos_filt <- df_gos_all %>% group_by(gene_set) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=3) %>% 
  select(-dataset,-counts) %>% unique()
df_msigdbs_filt <- df_msigdbs_all %>% group_by(gene_set) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=3) %>% 
  select(-dataset,-counts) %>% unique()
df_tfs_filt <- df_tfs %>% group_by(source) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=3) %>% 
  select(-dataset,-counts) %>% unique()
df_genes_filt <- df_sig_all %>% group_by(gene) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=3) %>% 
  dplyr::select(-dataset,-counts) %>% unique()
  
## annotate gene sets and get important genes-------------------
## Annotate MSIG
df_msig_annotation <- data.frame()
all_msig_genesets_filt <- all_msig_genesets[names(all_msig_genesets) %in% df_msigdbs_filt$gene_set]
for (i in 1:length(all_msig_genesets)){
  geneset <- all_msig_genesets[[i]]
  df_msig_annotation <- rbind(df_msig_annotation,
                              data.frame(gene_set = names(all_msig_genesets)[[i]],
                                         genes = geneset))
}
df_msig_genes <- distinct(df_msig_annotation %>% filter(gene_set %in% df_msigdbs_filt$gene_set) %>% select(genes))
gene_ids <- mapIds(org.Hs.eg.db, keys = df_msig_genes$genes, column = "SYMBOL", keytype = "ENTREZID")
gene_ids <- unname(gene_ids)
inds <- which(!is.na(gene_ids))
gene_ids <- gene_ids[inds]
df_msig_genes <- gene_ids
df_msig_genes <- unique(df_msig_genes)
saveRDS(df_msig_genes,'../results/clinical_msig_genes.rds')

## Annotate KEGGS
pathways <- kegg.pathways$human
pathways <- pathways$kg.sets
pathways <- pathways[df_kegg_filt$pathway]
kegg_genes <- c()
for (set in pathways){
  kegg_genes <- c(kegg_genes,set)
}
kegg_genes <- unique(kegg_genes)
kegg_gene_symbols <- mapIds(org.Hs.eg.db, keys = kegg_genes, column = "SYMBOL", keytype = "ENTREZID")
kegg_gene_symbols <- unname(kegg_gene_symbols)
inds <- which(!is.na(kegg_gene_symbols))
kegg_gene_symbols <- kegg_gene_symbols[inds]
saveRDS(kegg_gene_symbols,'../results/cliniacal_kegg_gene.rds')

### Save genes
saveRDS(df_genes_filt$gene,'../results/clinical_significant_genes.rds')

## Annotate GOs
genes_annot <- factor(x = rep(1,length(all_genes)),levels = c(0,1))
names(genes_annot) <- all_genes
GOobject <- new("topGOdata",ontology = 'BP', allGenes = genes_annot, annot=annFUN.org, mapping="org.Hs.eg.db", 
                ID = 'symbol', nodeSize = 10)
term.genes_bp <- genesInTerm(GOobject, GOobject@graph@nodes)

GOobject <- new("topGOdata",ontology = 'MF', allGenes = genes_annot, annot=annFUN.org, mapping="org.Hs.eg.db", 
                ID = 'symbol', nodeSize = 10)
term.genes_mf <- genesInTerm(GOobject, GOobject@graph@nodes)

GOobject <- new("topGOdata",ontology = 'CC', allGenes = genes_annot, annot=annFUN.org, mapping="org.Hs.eg.db", 
                ID = 'symbol', nodeSize = 10)
term.genes_cc <- genesInTerm(GOobject, GOobject@graph@nodes)

all_go_terms_annotations <- c(term.genes_bp,term.genes_mf,term.genes_cc)

all_go_terms_annotations <- all_go_terms_annotations[df_gos_filt$gene_set]
go_genes <- c()
for (set in all_go_terms_annotations){
  go_genes <- c(go_genes,set)
}
go_genes <- unique(go_genes)
saveRDS(go_genes,'../results/cliniacal_go_gene.rds')

### Combine all to get a core gene signature
regulatory_genes <- unique(c(df_tfs_filt$source,df_tfs_filt$target))
core_genes <- c(go_genes,kegg_gene_symbols,df_msig_genes,df_genes_filt$gene,regulatory_genes)
core_genes <- unique(core_genes)
saveRDS(core_genes,'../results/core_genes_clinical.rds')

### Get signaling net and augment core genes to keep genes whose products are in the path-------------
### from potential ligands to all the other core genes
core_genes <- readRDS('../results/core_genes_clinical.rds')
library(OmnipathR)
interactions <- import_omnipath_interactions()
all_ligands <- curated_ligand_receptor_interactions()
# all_ligands <- all_ligands %>% filter(n_resources>=3)  %>% filter(n_references>=3)
interactions <- interactions  %>% #mutate(kegg=grepl(pattern="KEGG",x=sources))%>%
  # filter(n_resources>=2 | kegg==T | source %in% all_ligands$source | target %in% all_ligands$target)  %>%
  # filter(n_references>=2 | kegg==T | source %in% all_ligands$source | target %in% all_ligands$target)  %>%
  # dplyr::select(-kegg) %>%
  dplyr::select(c('source'='source_genesymbol'),
                c('target'='target_genesymbol'),
                is_inhibition,is_stimulation) %>% unique()

### Get potential ligands and get all TFs
df_tfs <- readRDS('../results/bioinformatics_results/all_significant_tfs.rds')
df_tfs_filt <- df_tfs %>% group_by(source) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=3) %>% 
  select(-dataset,-counts) %>% unique()
all_tfs <- unique(df_tfs_filt$source)
all_receptors <- all_ligands$target_genesymbol
all_ligands <- all_ligands$source_genesymbol

### keep selected stuff
all_ligands <- all_ligands[which(all_ligands %in% interactions$source | all_ligands %in% interactions$target)]
all_tfs <-  all_tfs[which(all_tfs %in% c(interactions$source,interactions$target))]
all_receptors <-  all_receptors[which(all_receptors %in% c(interactions$source,interactions$target))]
core_genes_kept <- core_genes[which(core_genes %in% c(interactions$source,interactions$target))]

### Find the subgraph
g = graph.data.frame(interactions[,c('source','target')])
needNodes <- character()
## loop through all nodes and calculate the path to each other node
pb = txtProgressBar(min = 0, max = length(all_ligands), initial = 0) 
for(i in 1:length(all_ligands)){
  paths <- shortest_paths(g, from=all_ligands[i], 
                          to=unique(c(core_genes_kept,all_tfs)),
                          mode="all")
  needNodes <- unique(c(needNodes, unlist(lapply(paths$vpath, names))))
  saveRDS(needNodes,'needNodes_filt3.rds')
  saveRDS(i,'i.rds')
  setTxtProgressBar(pb,i)
}
# saveRDS(needNodes,'needNodes.rds')
## subset the graph
subGr <- induced_subgraph(g, vids=needNodes)
interactions2keep <- igraph::as_data_frame(subGr)
interactions2keep <- distinct(interactions2keep)

## Select final core genes and make sure they are in all human samples--------------
final_core <- unique(c(core_genes,interactions2keep$from,interactions2keep$to,all_tfs))
saveRDS(final_core,'final_core.rds')
final_core <- readRDS('final_core.rds')
gc()

### Load all human datasets and keep only core genes appearing to all of them
human_genes_list <- NULL
counter <- 1
for (geo in geos){
  list_deseq<- readRDS(paste0('../results/bioinformatics_results/',geo,'DESeqDF.rds'))
  for (i in 1:length(list_deseq)){
    if (is.null(list_deseq[[i]]$gene_names)){
      df_sig <- as.data.frame(list_deseq[[i]]) %>% rownames_to_column('gene') %>% filter(abs(log2FoldChange)>1.5) %>% filter(padj<0.05) %>% 
        dplyr::select(gene) %>% unique()
      genes <- as.data.frame(list_deseq[[i]]) %>% rownames_to_column('gene')
      genes <- unique(genes$gene)
    } else{
      df_sig <- as.data.frame(list_deseq[[i]]) %>% filter(abs(log2FoldChange)>1.5) %>% filter(padj<0.05) %>% 
        dplyr::select(c('gene'='gene_names')) %>% unique()
      genes <- as.data.frame(list_deseq[[i]])$gene_names
      genes <- unique(genes)
    }
  }
  human_genes_list[[counter]] <- genes
  counter = counter+1
}
common_across_human_datasets <- Reduce(intersect,human_genes_list)
final_core_common <- final_core[which(final_core %in% common_across_human_datasets)]
saveRDS(final_core_common,'../preprocessing/common_core_genes.rds')

### Explore how much of the human variance and in vitro variance are these genes capturing-----------
final_core_common <- readRDS('../preprocessing/common_core_genes.rds')
expression_matrix <- readRDS('../data/ARCHS4/FattyLiver_expression_matrix.rds')
meta_data <- read.delim('../data/ARCHS4/FattyLiver_meta_data.tsv',row.names = 1)
final_core_common <- final_core_common[which(final_core_common %in% rownames(expression_matrix))]

explained_var <- NULL
explained_var_rand <- NULL
k <- 1
for (i in 1:length(geos)){
  geo <- geos[i]
  tmp <- meta_data %>% filter(series_id==geo)
  X <- expression_matrix[,tmp$sample]
  genes <- rownames(X)
  X <- apply(X, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
  X <- log10(X+1)
  rownames(X) <- genes
  var_X <- sum(apply(X,1,var))
  
  X_subset <- X[final_core_common,]
  var_subset <- sum(apply(X_subset,1,var))
  
  for (j in 1:100){
    X_subset_rand <- X[sample.int(n = nrow(X),size=length(final_core_common)),]
    var_subset_rand <- sum(apply(X_subset_rand,1,var))
    explained_var_rand[k] <-  var_subset_rand/var_X
    k <- k+1
  }
  
  explained_var[i] <-  var_subset/var_X
  
  print(paste0('Finished dataset ',i,' out of ',length(geos)))
}

hist(explained_var)
hist(explained_var_rand)
wilcox.test(explained_var,explained_var_rand)

df_res <- rbind(data.frame(variance = explained_var_rand , geneset=rep('random',length(explained_var_rand))),
                data.frame(variance = explained_var , geneset=rep('selected',length(explained_var))))
ggplot(df_res,aes(x=geneset,y=100*variance,color = geneset)) +
  geom_jitter(size=2)+
  ylab('explained variance (%)')+
  stat_compare_means(method = 'wilcox.test',
                     tip.length=0.05,
                     label.x = 1.35,
                     size=5)+
  theme_pubr(base_family = 'Arial',base_size=24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'none')
ggsave('../results/bioinformatics_results/explained_variance_using_core_set.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)

### Repeat for in-vitro models
# Load gene expression
data_wang <- read.delim2('../data/GSE166256_kallisto_counts.tsv') %>% column_to_rownames('Gene')
data_wang <- as.matrix(data_wang)
data_wang <- apply(data_wang,MARGIN = c(1,2),as.numeric)
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_mps <- data$data_B
data_feaver <- read.delim2('../data/GSE89063_expression_matrix.tsv') # %>% column_to_rownames('X')
data_feaver <- aggregate(data_feaver[,2:ncol(data_feaver)],by = list(data_feaver$X),median)
gc()
data_feaver <- data_feaver %>% column_to_rownames('Group.1')

### Pre-process MPS
data_mps  <- apply(data_mps, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_mps <- log10(1 + data_mps)
### Pre-process Wang
data_wang  <- apply(data_wang, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_wang <- log10(1 + data_wang)
## Pre-process feaver
data_feaver  <- apply(data_feaver, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_feaver <- log10(1 + data_feaver)

invitro_datasets <- list(data_feaver,data_mps,data_wang)
explained_var <- NULL
explained_var_rand <- NULL
k <- 1
for (i in 1:length(invitro_datasets)){
  X <- invitro_datasets[[i]]
  var_X <- sum(apply(X,1,var))
  
  if (all(final_core_common %in% rownames(X))){
    keep_genes <- final_core_common
  }else{
    keep_genes <- final_core_common[which(final_core_common %in% rownames(X))]
  }
  X_subset <- X[keep_genes,]
  var_subset <- sum(apply(X_subset,1,var))
  for (j in 1:100){
    X_subset_rand <- X[sample.int(n = nrow(X),size=length(keep_genes)),]
    var_subset_rand <- sum(apply(X_subset_rand,1,var))
    explained_var_rand[k] <-  var_subset_rand/var_X
    k <- k+1
  }
  
  explained_var[i] <-  var_subset/var_X
  
  print(paste0('Finished dataset ',i,' out of ',length(invitro_datasets)))
}

hist(explained_var)
hist(explained_var_rand)
wilcox.test(explained_var,explained_var_rand)

df_res_invitro <- rbind(data.frame(variance = explained_var_rand , geneset=rep('random',length(explained_var_rand))),
                data.frame(variance = explained_var , geneset=rep('selected',length(explained_var))))
ggplot(df_res_invitro,aes(x=geneset,y=100*variance,color = geneset)) +
  geom_jitter(size=2)+
  ylab('explained variance (%)')+
  stat_compare_means(method = 'wilcox.test',
                     tip.length=0.05,
                     label.x = 1.35,
                     size=5)+
  theme_pubr(base_family = 'Arial',base_size=24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'none')
ggsave('../results/bioinformatics_results/explained_variance_invitro_using_core_set.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
