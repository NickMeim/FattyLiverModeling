# load in packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(pls)
library(pheatmap)
library(patchwork)
library(caret)
library(factoextra)

# Load data of each dataset----------------------------------------------------------
mps_metadata <- read.table('../data/GSE168285/GSE168285_sample_meta_data.txt',header = TRUE, sep = "\t") %>% select(-X)
wang_metadata <- data.table::fread('../data/METADATA_Wangetal2020.txt',header = TRUE, sep = "\t")
feaver_metadata <- read.table('../data/GSE89063_metadata.txt',header = TRUE, sep = ",")# %>% select(-X)

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

#Load human data
data_human <- data$data_A

# Pre-process
genes_common_mps <- reduce(list(rownames(data_mps),rownames(data_human)),intersect)
# data_mps <- data_mps[rowSums(data_mps >= 1) >= 0.1*ncol(data_mps),]
data_human_mps <-  data_human[genes_common_mps,]
data_mps <-  data_mps[genes_common_mps,]
keep_gene_mps <- (rowSums(data_mps >= 1) >= 0.1*ncol(data_mps)) &
  (rowSums(data_human_mps >= 1) >= 0.1*ncol(data_human_mps))
data_human_mps <-  data_human_mps[keep_gene_mps,]
data_mps <-  data_mps[keep_gene_mps,]
data_mps  <- apply(data_mps, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_mps <- log10(1 + data_mps)
data_mps <- data_mps[apply(data_mps,1,sd)>1e-2,]
data_mps <- scale(t(data_mps),scale = FALSE)
hist(data_mps)
data_human_mps  <- apply(data_human_mps, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_human_mps <- log10(1 + data_human_mps)
data_human_mps <- data_human_mps[apply(data_human_mps,1,sd)>1e-2,]
data_human_mps <- scale(t(data_human_mps),scale = FALSE)
hist(data_human_mps)
genes_common_mps <- reduce(list(colnames(data_mps),colnames(data_human_mps)),intersect)
data_human_mps <-  data_human_mps[,genes_common_mps]
data_mps <-  data_mps[,genes_common_mps]
# Dimensionality reduction using human PCs
pca_human_mps <- prcomp(data_human_mps,scale. = FALSE)
fviz_screeplot(pca_human_mps, addlabels = TRUE,n = 20)
perc_var <- cumsum(pca_human_mps$sdev^2)/sum(pca_human_mps$sdev^2)
pc2keep <- which(perc_var>=0.9)[1] 
W1 <- pca_human_mps$rotation[,1:pc2keep]
data_mps_projected <- data_mps %*% W1
hist(data_mps_projected)

### Pre-process Wang
genes_common_wang <- reduce(list(rownames(data_wang),rownames(data_human)),intersect)
data_human_wang <-  data_human[genes_common_wang,]
data_wang <-  data_wang[genes_common_wang,]
keep_gene_wang <- (rowSums(data_wang >= 1) >= 0.1*ncol(data_wang)) &
  (rowSums(data_human_wang >= 1) >= 0.1*ncol(data_human_wang))
data_human_wang <-  data_human_wang[keep_gene_wang,]
data_wang <-  data_wang[keep_gene_wang,]
data_wang  <- apply(data_wang, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_wang <- log10(1 + data_wang)
data_wang <- data_wang[apply(data_wang,1,sd)>1e-2,]
data_wang <- scale(t(data_wang),scale = FALSE)
hist(data_wang)
data_human_wang  <- apply(data_human_wang, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_human_wang <- log10(1 + data_human_wang)
data_human_wang <- data_human_wang[apply(data_human_wang,1,sd)>1e-2,]
data_human_wang <- scale(t(data_human_wang),scale = FALSE)
hist(data_human_wang)
genes_common_wang <- reduce(list(colnames(data_wang),colnames(data_human_wang)),intersect)
data_human_wang <-  data_human_wang[,genes_common_wang]
data_wang <-  data_wang[,genes_common_wang]
# Dimensionality reduction using human PCs
pca_human_wang <- prcomp(data_human_wang,scale. = FALSE)
fviz_screeplot(pca_human_wang, addlabels = TRUE,n = 20)
perc_var <- cumsum(pca_human_wang$sdev^2)/sum(pca_human_wang$sdev^2)
pc2keep <- which(perc_var>=0.9)[1] 
W2 <- pca_human_wang$rotation[,1:pc2keep]
data_wang_projected <- data_wang %*% W2
hist(data_wang_projected)


## Preprocess feaver
genes_common_feaver <- reduce(list(rownames(data_feaver),rownames(data_human)),intersect)
data_human_feaver <-  data_human[genes_common_feaver,]
data_feaver <-  data_feaver[genes_common_feaver,]
keep_gene_feaver <- (rowSums(data_feaver >= 1) >= 0.1*ncol(data_feaver)) &
  (rowSums(data_human_feaver >= 1) >= 0.1*ncol(data_human_feaver))
data_human_feaver <-  data_human_feaver[keep_gene_feaver,]
data_feaver <-  data_feaver[keep_gene_feaver,]
data_feaver  <- apply(data_feaver, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_feaver <- log10(1 + data_feaver)
data_feaver <- data_feaver[apply(data_feaver,1,sd)>1e-2,]
data_feaver <- scale(t(data_feaver),scale = FALSE)
hist(data_feaver)
data_human_feaver  <- apply(data_human_feaver, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_human_feaver <- log10(1 + data_human_feaver)
data_human_feaver <- data_human_feaver[apply(data_human_feaver,1,sd)>1e-2,]
data_human_feaver <- scale(t(data_human_feaver),scale = FALSE)
hist(data_human_feaver)
genes_common_feaver <- reduce(list(colnames(data_feaver),colnames(data_human_feaver)),intersect)
data_human_feaver <-  data_human_feaver[,genes_common_feaver]
data_feaver <-  data_feaver[,genes_common_feaver]
# Dimensionality reduction using human PCs
pca_human_feaver <- prcomp(data_human_feaver,scale. = FALSE)
fviz_screeplot(pca_human_feaver, addlabels = TRUE,n = 20)
perc_var <- cumsum(pca_human_feaver$sdev^2)/sum(pca_human_feaver$sdev^2)
pc2keep <- which(perc_var>=0.9)[1] 
W3 <- pca_human_feaver$rotation[,1:pc2keep]
data_feaver_projected <- data_feaver %*% W3
hist(data_feaver_projected)

all_data <- list(data_mps,data_mps_projected,data_human_mps,mps_metadata,
                 data_feaver,data_feaver_projected,data_human_feaver,feaver_metadata,
                 data_wang,data_wang_projected,data_human_wang,wang_metadata,
                 data$Y_A)
names(all_data) <- c('data_mps','data_mps_projected','data_human_mps','mps_metadata',
                     'data_feaver','data_feaver_projected','data_human_feaver','feaver_metadata',
                     'data_wang','data_wang_projected','data_human_wang','wang_metadata',
                     'NAS')
saveRDS(all_data,'../preprocessing/all_invitro_data.rds')

### Use SVMs to predict conditions-------------
all_data <- readRDS('../preprocessing/all_invitro_data.rds')
# data_mps <- all_data$data_mps
data_mps_projected <- all_data$data_mps_projected
# data_human_mps <- all_data$data_human_mps
mps_metadata <- all_data$mps_metadata
# data_feaver <- all_data$data_feaver
data_feaver_projected <- all_data$data_feaver_projected
# data_human_feaver <- all_data$data_human_feaver
feaver_metadata <- all_data$feaver_metadata
# data_wang<- all_data$data_wang
data_wang_projected<- all_data$data_wang_projected
# data_human_wang<- all_data$data_human_wang
wang_metadata<- all_data$wang_metadata

### First do that for predicting wang media condition per scaffold
scaffolds <- unique(wang_metadata$scaffold)
for (scaff in scaffolds){
  metadata <- wang_metadata %>% filter(scaffold==scaff)
  X <- data_wang_projected[paste0('X',metadata$filename),]
  Y <- metadata$media
  # Y <- ifelse(Y=="enhanced liver culture media",1,0)
  # Y <- factor(Y,levels = c(0,1))
  Y <- as.matrix(Y)
  rownames(Y) <- rownames(X)
  colnames(Y) <- 'media'
  data_scaff <- cbind(as.data.frame(X),as.data.frame(Y))
  preds <- NULL
  for (i in 1:nrow(data_scaff)){
    inds <- seq(1,nrow(data_scaff))
    inds <- inds[which(inds!=i)]
    mdl <- train(media ~ ., data = data_scaff[inds,], method = 'rf')
    yhat <- predict(mdl,newdata = data_scaff[i,])
    preds <- c(preds,yhat)
  }
  conf <- confusionMatrix()
  for (iter in 1:10){
    preds_shuffled <- NULL
    for (i in 1:nrow(data_scaff)){
      inds <- seq(1,nrow(data_scaff))
      inds <- inds[which(inds!=i)]
      mdl <- train(media ~ ., data = data_scaff[inds,], method = 'svmLinear')
      yhat <- predict(mdl,newdata = data_scaff[i,])
      preds <- c(preds,yhat)
    }
  }
}






