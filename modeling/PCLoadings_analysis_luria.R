# load in packages
library(tidyverse)
library(dorothea)
library(AnnotationDbi)
library(OmnipathR)
library(CARNIVAL)

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

#### Load pre-processed data----------------------
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_A <- data$data_A
data_B <- data$data_B
Y_A <- data$Y_A

## normalize
# Intersect datasets - keep common genes
genes_common <- intersect(rownames(data_A), rownames(data_B))
# Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
data_A <- apply(data_A[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_B <- apply(data_B[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
# Remove features absent from at least 10% in each samples
keep_gene <- (rowSums(data_A >= 1) >= 0.1*ncol(data_A)) &
  (rowSums(data_B >= 1) >= 0.1*ncol(data_B))
# Run PCA on each dataset. Not scaling PCA
# Log and center dataset A and run PCA
X_A <- log2(1 + data_A[keep_gene,])
X_A <- t(X_A - rowMeans(X_A))
# Log and center dataset B and run PCA
X_B <- log2(1 + data_B[keep_gene,])
X_B <- t(X_B - rowMeans(X_B))

### Perform PCA and get loadings----------------
PCA_alldata <- prcomp(X_B, scale. = F, center = T)

#### Run CARNIVAL--------------------------
gene_loadings <- PCA_alldata$rotation
minNrOfGenes  <-  5
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities = run_viper(gene_loadings, dorotheaData, options =  settings)
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
j <- args
print(paste0('Iteration ',j,'/',ncol(TF_activities)))
tf_activities <- TF_activities[which(rownames(TF_activities) %in% c(interactions$source,
                                                                    interactions$target)),
                               colnames(TF_activities)[j]]
# Run carnival
CbcPath <- '../../CBC/bin/cbc'
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

# Output dir
Result_dir <- paste0("../results/pc_loadings_scores_analysis/",colnames(TF_activities)[j])
dir.create(Result_dir, showWarnings = FALSE)
carnivalOptions$outputFolder <- Result_dir
carnivalOptions$workdir <- Result_dir

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
