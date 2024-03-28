library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(pheatmap)
library(corrplot)
library(caret)
library(infotheo)
library(glmnet)
library(matlib)
library(factoextra)
source("aux_functions.R")
source("functions_translation_nikos.R")

### Load all raw data and re-do PCA---------------------------
# Load all data
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_mps <- data$data_B
data_human <- data$data_A
NAFLD_score <- data$Y_A
data_wang <- read.delim2('../data/GSE166256_kallisto_counts.tsv') %>% column_to_rownames('Gene')
data_wang <- as.matrix(data_wang)
data_wang <- apply(data_wang,MARGIN = c(1,2),as.numeric)
genes_common <- reduce(list(rownames(data_mps), rownames(data_human),rownames(data_wang)),intersect)

# Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
data_human <- apply(data_human[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_mps <- apply(data_mps[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_wang <- apply(data_wang[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
# Remove features absent from at least 10% in each samples
# keep_gene <- (rowSums(data_human >= 1) >= 0.1*ncol(data_human)) &
#   (rowSums(data_mps >= 1) >= 0.1*ncol(data_mps)) &
#   (rowSums(data_wang >= 1) >= 0.1*ncol(data_wang))
# data_human <- log2(1 + data_human[keep_gene,])
data_human <- log10(1 + data_human)
data_human <- t(data_human - rowMeans(data_human))
# data_mps <- log2(1 + data_mps[keep_gene,])
data_mps <- log10(1 + data_mps)
data_mps <- t(data_mps - rowMeans(data_mps))
# data_wang <- log2(1 + data_wang[keep_gene,])
data_wang <- log10(1 + data_wang)
data_wang <- t(data_wang - rowMeans(data_wang))

varA <-  sum(apply(data_human,2,sd)^2)
# perform PCA in each dataset
# first for MPS
pca_mps <- prcomp(data_mps,center = TRUE,scale. = FALSE)
# fviz_eig(pca_mps)
Z1 <- data_human %*% pca_mps$rotation
varAB <- apply(Z1,2,sd)^2
perc <- 100*varAB/varA
names(perc) <- colnames(Z1)
perc <- perc[order(-perc)][1:15]
perc <- as.data.frame(perc)
perc <- perc %>% rownames_to_column('PC')
p1 <- ggbarplot(perc,x='PC',y='perc',fill='purple') + ylab('% variance explained') + ggtitle('Percentage of variance explained of the human data')+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
Zb <- pca_mps$x
varB <- sum(apply(data_mps,2,sd)^2)
varPCB <- apply(Zb,2,sd)^2
perc2 <- 100*varPCB/varB
perc2 <- perc2[perc$PC]
perc2 <- as.data.frame(perc2)
perc2 <- perc2 %>% rownames_to_column('PC')
p2 <- ggbarplot(perc2,x='PC',y='perc2',fill='blue') + ylab('% variance explained') + ggtitle('Percentage of variance explained of the MPS data')+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
p <- p1/p2
# print(p)

# secondly for Wang
pca_wang <- prcomp(data_wang,center = TRUE,scale. = FALSE)
# fviz_eig(pca_wang)
Z2 <- data_human %*% pca_wang$rotation
varAC <- apply(Z2,2,sd)^2
perc <- 100*varAC/varA
names(perc) <- colnames(Z2)
perc <- perc[order(-perc)][1:15]
perc <- as.data.frame(perc)
perc <- perc %>% rownames_to_column('PC')
p1 <- ggbarplot(perc,x='PC',y='perc',fill='purple') + ylab('% variance explained') + ggtitle('Percentage of variance explained of the human data')+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
perc2 <- 100*varPCB/varB
perc2 <- perc2[perc$PC]
perc2 <- as.data.frame(perc2)
perc2 <- perc2 %>% rownames_to_column('PC')
p2 <- ggbarplot(perc2,x='PC',y='perc2',fill='blue') + ylab('% variance explained') + ggtitle('Percentage of variance explained of the MPS data')+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
p <- p1/p2
# print(p)

### I want to project from relevent Wang to non-relevent MPS
load('../results/two_MPS_model.RData')
non_interesting_mps <- colnames(features_joint)[grep('MPS',colnames(features_joint))]
selected_wang <- colnames(features_joint)[grep('Wang',colnames(features_joint))]
# 18 PCs for each dataset should be enough
# Z1 <- Z1[,1:ncol(Z1)]
colnames(Z1) <- paste0(colnames(Z1),'_MPS')
# Z1 <- Z1[,!(colnames(Z1) %in% non_interesting_mps)]
Z1 <- Z1[,1:ncol(Z2)]
# Z2 <- Z2[,1:ncol(Z2)]
colnames(Z2) <- paste0(colnames(Z2),'_Wang')
# Z2 <- Z2[,selected_wang] 
  
data_combined <- cbind(Z1,Z2)
data_combined <- as.data.frame(data_combined)
data_combined$NAS <-  data$Y_A 
L1 <- pca_mps$rotation
colnames(L1) <- paste0(colnames(L1),'_MPS')
L1 <- L1[,colnames(Z1)]
L2 <- pca_wang$rotation
colnames(L2) <- paste0(colnames(L2),'_Wang')
L2 <- L2[,colnames(Z2)]

# ### Load selected MPS data---------------
# load('../results/two_MPS_model.RData')
# data_combined <- as.data.frame(features_joint)
# data_combined$NAS <- NAFLD_score
# Z1 <- features_joint[,grep('MPS',colnames(features_joint))]
# Z2 <- features_joint[,!grepl('MPS',colnames(features_joint))]


### Fit linear regression model and get coefficients-------------------
fit <- lm(NAS~.,data_combined)
print(paste0('MAE=',mean(abs(NAFLD_score - predict(fit)))))
print(paste0('Kendall tau=',cor(NAFLD_score,predict(fit),method = 'kendall')))
print(paste0('r=',cor(NAFLD_score,predict(fit))))
print(paste0('R2=',cor(NAFLD_score,predict(fit))^2))
summary(fit)

A1 <- fit$coefficients[grep('MPS',names(fit$coefficients))]
A2 <- fit$coefficients[!grepl('MPS',names(fit$coefficients)) & !grepl('(Intercept)',names(fit$coefficients))]
b <- fit$coefficients['(Intercept)']
Y <- as.matrix(NAFLD_score)
# L1 = loadings_joint[,grepl('MPS',colnames(features_joint))]
# L2 = loadings_joint[,!grepl('MPS',colnames(features_joint))]

### Use PRECISE approach: Optimal P = L2 * L1^T (if dim(L) = factors x genes-----------------------------
print(dim(L1))
print(dim(L2))
P  <-  t(L1)  %*% L2
max(P)
colnames(P) <- colnames(Z2)
rownames(P) <- colnames(Z1)
dim(P)

max_sim <- apply(P[,selected_wang],1,max)
max_sim <- max_sim[order(-max_sim)]
print(max_sim)
M <- t(P[,selected_wang])

### Find most similar loadings between two systems
df_similar <- as.data.frame(P) %>% rownames_to_column('MPS') %>%
  gather('Wang','similarity',-'MPS') %>% mutate(MPS=strsplit(MPS,'_')) %>% unnest(MPS) %>% filter(MPS!='MPS') %>%
  mutate(Wang=strsplit(Wang,'_')) %>% unnest(Wang) %>% filter(Wang!='Wang') %>%
  group_by(Wang) %>% mutate(max_sim = max(similarity)) %>% ungroup() %>%
  filter(similarity==max_sim) %>% select(-max_sim) %>% unique() %>%
  filter(Wang %in% str_split_fixed(selected_wang,'_',2)[,1])


L2_hat <- L1 %*% t(P)
L1_hat <- L2 %*% P
hist(diag(cor(L2_hat,L2)),main= 'Pearson correlation of between true and reconstructed loadings')
print(mean(diag(cor(L2_hat,L2))))
hist(diag(cor(L2_hat,L2,method = 'spearman')),main='Spearman correlation of between true and reconstructed loadings')
print(mean(diag(cor(L2_hat,L2,method = 'spearman'))))

print(diag(cor(L1_hat,L1)))
print(mean(diag(cor(L1_hat,L1))))

# evaluate Z1, Z2 reconstruction
Z1_hat <- data_human %*% L1_hat
print(mean(diag(cor(Z1_hat,Z1))))
hist(diag(cor(Z1_hat,Z1)))

Z2_hat <- data_human %*% L2_hat
print(mean(diag(cor(Z2_hat,Z2))))
hist(diag(cor(Z2_hat,Z2)))


#### SVD Decomposition of similarity matrix--------------------------------------
P  <-  t(L2) %*% L1[,1:10]
decomposed <- svd(P)
S <- L2 %*% decomposed$u
colnames(S) <- paste0('S',seq(1,ncol(S)))
Tau <- L1[,1:10] %*% decomposed$v
colnames(Tau) <- paste0('T',seq(1,ncol(S)))
D <- decomposed$d
mm <- t(S) %*% Tau
genes_2 <- names(S[order(S[,1]),1])
genes_1 <- names(Tau[order(Tau[,1]),1])
no_top_bot <- 2000
top_1 <- genes_1[c(seq(1,no_top_bot),seq((length(genes_1)-no_top_bot),length(genes_1)))]
top_2 <- genes_2[c(seq(1,no_top_bot),seq((length(genes_2)-no_top_bot),length(genes_2)))]
sim <- length(intersect(top_1,top_2))/length(union(top_1,top_2))

###evaluate precise approach------------------------------------------------
Yhat <- Z1 %*% A1 + Z1 %*% t(P) %*% A2 + b
# if (max(Yhat)>max(NAFLD_score)){
#   Yhat <- Yhat - (max(Yhat)-max(NAFLD_score))
# }
# if (min(Yhat)<min(NAFLD_score)){
#   Yhat <- Yhat + abs((min(Yhat)-min(NAFLD_score)))
# }
print(paste0('MAE=',mean(abs(NAFLD_score - Yhat[,1]))))
print(paste0('Kendall tau=',cor(NAFLD_score,Yhat[,1],method = 'kendall')))
print(paste0('r=',cor(NAFLD_score,Yhat[,1])))
print(paste0('R2=',cor(NAFLD_score,Yhat[,1])^2))

print(paste0('MAE=',mean(abs(predict(fit) - Yhat[,1]))))
print(paste0('Kendall tau=',cor(predict(fit),Yhat[,1],method = 'kendall')))
print(paste0('r=',cor(predict(fit),Yhat[,1])))
print(paste0('R2=',cor(predict(fit),Yhat[,1])^2))

# Repeat for predicting L1
Yhat <- Z2 %*% P %*% A1 + Z2 %*% A2 + b
print(paste0('MAE=',mean(abs(NAFLD_score - Yhat[,1]))))
print(paste0('Kendall tau=',cor(NAFLD_score,Yhat[,1],method = 'kendall')))
print(paste0('r=',cor(NAFLD_score,Yhat[,1])))
print(paste0('R2=',cor(NAFLD_score,Yhat[,1])^2))
print(paste0('MAE=',mean(abs(predict(fit) - Yhat[,1]))))
print(paste0('Kendall tau=',cor(predict(fit),Yhat[,1],method = 'kendall')))
print(paste0('r=',cor(predict(fit),Yhat[,1])))
print(paste0('R2=',cor(predict(fit),Yhat[,1])^2))


K = Y - Z1 %*% as.matrix(A1) - b
write_delim(as.data.frame(K),'../results/K.tsv')
write_delim(as.data.frame(A2),'../results/A2.tsv')
write_delim(as.data.frame(Z1),'../results/Z1.tsv')
write_delim(as.data.frame(L1),'../results/L1.tsv')
write_delim(as.data.frame(L2),'../results/L2.tsv')

data_star <- cbind(L2,L1)
colnames(data_star) <- c(paste0('Wang_PC',seq(1,ncol(L2))), paste0('MPS_PC',seq(1,ncol(L1))))
cols <- colnames(data_star)[grepl('MPS',colnames(data_star))]
P <- matrix(nrow = ncol(Z1),ncol = ncol(Z2),0)
b <- NULL
i <- 1
for (var in cols){
  data_selected <- data_star[,c(colnames(Z2),var)]
  colnames(data_selected)[ncol(data_selected)] <- 'Y'
  data_selected <- as.data.frame(data_selected)
  fit_star <- lm(Y ~.,data_selected)
  P[i,] <- fit_star$coefficients[!grepl("(Intercept)",names(fit_star$coefficients))]
  b[i] <- fit_star$coefficients["(Intercept)"]
  i <- i+1
}

### evaluate------------------------------------- 
L1 = loadings_joint[,grepl('MPS',colnames(features_joint))]
L2 = loadings_joint[,!grepl('MPS',colnames(features_joint))]
L2_hat = L1 %*% P
# colnames(L2_hat) <- colnames(L2)
pcs_loads_corr <- cor(L2,L2_hat)
print(diag(pcs_loads_corr))
write_delim(as.data.frame(L2),'../results/L2.tsv')
write_delim(as.data.frame(L1),'../results/L1.tsv')
