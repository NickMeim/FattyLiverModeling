library(tidyverse)
source('distance_scores.R')
library(doFuture)
# parallel: set number of workers
cores <- 8
registerDoFuture()
plan(multisession,workers = cores)


# Access command-line arguments
# args <- commandArgs(trailingOnly = TRUE)
# args <- as.numeric(args)
args <- 4

list_types <- c('hallmarks','keggs','gobp','tfs')
print(paste0('Feature selected : ',list_types[args[1]]))
type_feats <- list_types[args[1]]


if (type_feats=='hallmarks'){
  l1000_feature_matrix <- readRDS('../../../L1000_2021_11_23/GSEA/df_hallmarks_ligand_drugs_10k_filtered.rds')
  thresholds <- c(6,10,15,20)
}else if (type_feats=='gobp'){
  l1000_feature_matrix <- readRDS('../../../L1000_2021_11_23/GSEA/df_gobp_ligand_drugs_10k_filtered.rds')
  thresholds <- c(10,20,30,50,50)
}else if (type_feats=='keggs'){
  l1000_feature_matrix <- readRDS('../../../L1000_2021_11_23/GSEA/df_gsea_keggs_ligand_drugs_10k_filtered.rds')
  thresholds <- c(5,10,20,30,50,50)
}else if (type_feats=='tfs'){
  thresholds <- c(6,10,15,20,30)
  drug_ligand_ex <- readRDS("../data/l1000_drugs_ligands_expression.rds")
  gene_info <- readRDS("../data/l1000_geneInfo.rds")
  rownames(drug_ligand_ex) <- gene_info$gene_symbol
  cell_info <- readRDS("../data/l1000_cell_line_info_ligands_drugs.rds")
  sigInfo <- readRDS("../data/l1000_meta_data.rds")
  drug_ligand_ex <- drug_ligand_ex[,sigInfo$sig_id]
  dorotheaData <- read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
  confidenceFilter <- is.element(dorotheaData$confidence, c('A', 'B'))
  dorotheaData <- dorotheaData[confidenceFilter,]
  colnames(dorotheaData)[1] <- 'source' 
  minNrOfGenes  <-  5
  l1000_feature_matrix <- decoupleR::run_viper(drug_ligand_ex, dorotheaData,minsize = minNrOfGenes,verbose = FALSE)
  l1000_feature_matrix <- l1000_feature_matrix %>% select(c('sigID' = 'condition'),c('pathway'='source'),c('NES'='score'))
}

l1000_feature_matrix <- as.matrix(l1000_feature_matrix %>% select(sigID,pathway,NES) %>% spread('sigID','NES') %>% 
                                    column_to_rownames('pathway'))

dist_all_dupls <- NULL
print('Begin calculating GSEA distance...')
### calculate distances
dist_all_dupls  <- foreach(thres = thresholds) %dopar% {
  distance_scores(num_table = l1000_feature_matrix ,threshold_count = thres,names = colnames(l1000_feature_matrix))
}
print('Finished calculating GSEA distance')

distance <- do.call(cbind,dist_all_dupls)
distance <- array(distance,c(dim=dim(dist_all_dupls[[1]]),length(dist_all_dupls)))
mean_dist <- apply(distance, c(1,2), mean, na.rm = TRUE)
colnames(mean_dist) <- colnames(l1000_feature_matrix)
rownames(mean_dist) <- colnames(l1000_feature_matrix)
print('Begin saving GSEA distance...')
saveRDS(mean_dist,paste0(type_feats,'_gsea_dist_drugs_ligands_filtered_10k.rds'))