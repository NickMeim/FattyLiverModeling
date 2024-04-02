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
source("distance_scores.R")
library(doFuture)
# parallel: set number of workers
cores <- 16
registerDoFuture()
plan(multisession,workers = cores)

### Load all raw data and re-do PCA---------------------------
# Load all data
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_mps <- data$data_B
data_wang <- read.delim2('../data/GSE166256_kallisto_counts.tsv') %>% column_to_rownames('Gene')
data_wang <- as.matrix(data_wang)
data_wang <- apply(data_wang,MARGIN = c(1,2),as.numeric)
genes_common <- reduce(list(rownames(data_mps),rownames(data_wang)),intersect)

# Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
data_mps <- apply(data_mps[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_wang <- apply(data_wang[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_mps <- log10(1 + data_mps)
data_mps <- t(data_mps - rowMeans(data_mps))
data_wang <- log10(1 + data_wang)
data_wang <- t(data_wang - rowMeans(data_wang))

# perform PCA in each dataset
# first for MPS
pca_mps <- prcomp(data_mps,center = TRUE,scale. = FALSE)

# secondly for Wang
pca_wang <- prcomp(data_wang,center = TRUE,scale. = FALSE)

### I want to project from relevent Wang to non-relevent MPS
load('../results/two_MPS_model.RData')
non_interesting_mps <- colnames(features_joint)[grep('MPS',colnames(features_joint))]
selected_wang <- colnames(features_joint)[grep('Wang',colnames(features_joint))]

L1 <- pca_mps$rotation
colnames(L1) <- paste0(colnames(L1),'_MPS')
# L1 <- L1[,which(!(colnames(L1) %in% non_interesting_mps))]
L2 <- pca_wang$rotation
colnames(L2) <- paste0(colnames(L2),'_Wang')
# L2 <- L2[,selected_wang]

### Use PRECISE approach: Optimal P = L2 * L1^T (if dim(L) = factors x genes-----------------------------
print(dim(L1))
print(dim(L2))
P  <-  t(L1)  %*% L2
max(P)
colnames(P) <- colnames(L2)
rownames(P) <- colnames(L1)
dim(P)

max_sim <- apply(P[,selected_wang],1,max)
max_sim <- max_sim[order(-max_sim)]
print(max_sim)
M <- t(P[,selected_wang])

### Find most similar loadings between two systems
df_similar <- as.data.frame(P) %>% rownames_to_column('MPS') %>%
  gather('Wang','similarity',-'MPS') %>% mutate(MPS=strsplit(MPS,'_')) %>% unnest(MPS) %>% filter(MPS!='MPS') %>%
  mutate(Wang=strsplit(Wang,'_')) %>% unnest(Wang) %>% filter(Wang!='Wang') #%>%
  # filter(!(MPS %in% str_split_fixed(non_interesting_mps,'_',2)[,1])) %>%
  # group_by(Wang) %>% mutate(max_sim = max(similarity)) %>% ungroup() %>%
  # filter(similarity==max_sim) %>% select(-max_sim) %>% unique() %>%
  # filter(Wang %in% str_split_fixed(selected_wang,'_',2)[,1])

df_corr <- as.data.frame(cor(L1,L2,method = 'pearson')) %>% rownames_to_column('MPS') %>% 
  gather('Wang','correlation',-'MPS') %>% mutate(MPS=strsplit(MPS,'_')) %>% unnest(MPS) %>% 
  filter(MPS!='MPS') %>%
  mutate(Wang=strsplit(Wang,'_')) %>% unnest(Wang) %>% filter(Wang!='Wang') #%>%
  # filter(!(MPS %in% str_split_fixed(non_interesting_mps,'_',2)[,1])) %>%
  # group_by(Wang) %>% mutate(max_corr = max(correlation)) %>% ungroup() %>%
  # filter(correlation==max_corr) %>% select(-max_corr) %>% unique() %>%
  # filter(Wang %in% str_split_fixed(selected_wang,'_',2)[,1])

# calculate GSEA distance
thresholds <- c(5,10,20,30,50,75,100,125,150,175,200)
# Initialize empty list for the results:
# Each element of the list (for each threshold)
# contains an NxN matrix with comparing all these
# samples. Each element of the matrix is the
# GSEA distance.
dist_gsea <- NULL
### SOS:
### RUN FIRST THE distance_scores.R
### SCRIPT TO LOAD THE FUNCTION!!!
### calculate distances: SEE distance_scores.R
# for information about the function inputs
dist_gsea <- foreach(thres = thresholds) %dopar% {
  distance_scores(num_table = cbind(L1,L2) ,
                  threshold_count = thres,names = colnames(cbind(L1,L2)))
}
# Transform list to array
distance <- do.call(cbind,dist_gsea)
distance <- array(distance,
                  c(dim=dim(dist_gsea[[1]]),length(dist_gsea)))
# Get the average distance across thresholds
mean_dist <- apply(distance, c(1,2), mean, na.rm = TRUE)
colnames(mean_dist) <- colnames(cbind(L1,L2))
rownames(mean_dist) <- colnames(cbind(L1,L2))
### Convert matrix into data frame
# Keep only unique (non-self) pairs
mean_dist[lower.tri(mean_dist,diag = T)] <- -100
df_GSEA <- reshape2::melt(mean_dist)
df_GSEA <- df_GSEA %>% filter(value != -100)
colnames(df_GSEA)[3] <- 'gsea_dist'
df_GSEA <- df_GSEA %>% mutate(keep = ifelse(grepl('MPS',Var1) & grepl('MPS',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
df_GSEA <- df_GSEA %>% mutate(keep = ifelse(grepl('Wang',Var1) & grepl('Wang',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
df_GSEA <- df_GSEA %>% mutate(gsea_dist=gsea_dist/2)
data_distribution <- df_GSEA$gsea_dist
hist(data_distribution,60)

### check std of calculating the distance
sd_dist <- apply(distance, c(1,2), sd, na.rm = TRUE)
colnames(sd_dist) <- colnames(cbind(L1,L2))
rownames(sd_dist) <- colnames(cbind(L1,L2))
sd_dist[lower.tri(sd_dist,diag = T)] <- -100
df_GSEA_sd <- reshape2::melt(sd_dist)
df_GSEA_sd <- df_GSEA_sd %>% filter(value != -100)
colnames(df_GSEA_sd)[3] <- 'gsea_dist'
df_GSEA_sd <- df_GSEA_sd %>% mutate(keep = ifelse(grepl('MPS',Var1) & grepl('MPS',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
df_GSEA_sd <- df_GSEA_sd %>% mutate(keep = ifelse(grepl('Wang',Var1) & grepl('Wang',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
df_GSEA_sd <- df_GSEA_sd %>% mutate(gsea_dist=gsea_dist/2)
colnames(df_GSEA_sd) <- c('MPS','Wang','std')
df_GSEA_sd$MPS <- as.character(df_GSEA_sd$MPS)
df_GSEA_sd$Wang <- as.character(df_GSEA_sd$Wang)

colnames(df_GSEA)[1:2] <- c('MPS','Wang')
df_GSEA$MPS <- as.character(df_GSEA$MPS)
df_GSEA$Wang <- as.character(df_GSEA$Wang)
df_GSEA <- left_join(df_GSEA,df_GSEA_sd)
df_GSEA <- df_GSEA %>% mutate(CV=std/gsea_dist)
hist(100*df_GSEA$CV,60,main='Coefficient of Variation',xlab = 'CV (%)')
df_GSEA <- df_GSEA %>% mutate(MPS=strsplit(MPS,'_')) %>% unnest(MPS) %>% filter(MPS!='MPS') %>%
  mutate(Wang=strsplit(Wang,'_')) %>% unnest(Wang) %>% filter(Wang!='Wang') #%>%
  # filter(!(MPS %in% str_split_fixed(non_interesting_mps,'_',2)[,1])) %>%
  # group_by(Wang) %>% mutate(min_dist = min(gsea_dist)) %>% ungroup() %>%
  # filter(gsea_dist==min_dist) %>% select(-min_dist) %>% unique() %>%
  # filter(Wang %in% str_split_fixed(selected_wang,'_',2)[,1])

## Combine all results
# df_combined <- left_join(df_GSEA,df_corr)
# df_combined <- left_join(df_combined,df_similar)
df_combined <- left_join(df_GSEA,df_similar)
# df_combined <- df_combined %>% mutate(score = (2*gsea_dist + (1 - correlation) +(1 - similarity))/3)
df_combined <- df_combined %>% mutate(score = (2*gsea_dist + (1 - similarity))*0.5)
data_score_distribution <- df_combined$score
hist(data_score_distribution,60)
df_combined <- df_combined %>% filter(!(MPS %in% str_split_fixed(non_interesting_mps,'_',2)[,1])) %>%
  group_by(Wang) %>% mutate(min_score = min(score)) %>% ungroup() %>%
  filter(score==min_score) %>% select(-min_score) %>% unique() %>%
  filter(Wang %in% str_split_fixed(selected_wang,'_',2)[,1])
df_combined <- df_combined %>% group_by(score) %>% 
  mutate(p.value = sum(score >= data_score_distribution)/length(data_score_distribution)) %>%
  ungroup()
df_combined <- df_combined %>% mutate(p.adj = p.value * n_distinct(Wang))

### Check how many genes on average are common for different distance thresholds-----------------------
common_genes <- NULL
NumberOfCommonElements <- function(x,th,names,name1,name2,type='up'){
  if (type=='bottom'){
    x_ranked <- apply(x,2,rank) 
  }else{
    x_ranked <- apply(-x,2,rank) 
  }
  # mask <- (x_ranked <= th) | (x_ranked >= (nrow(x_ranked) - th + 1))
  mask <- (x_ranked <= th)
  mask <- 1 * mask
  common <- t(mask) %*% mask 
  # common <- common/(2*th)
  common <- common/th
  colnames(common) <- names
  rownames(common) <- names
  common <- common[name1,name2]
  return(common)
}

common_top_genes <- foreach(thres = thresholds) %dopar% {
  NumberOfCommonElements(x = cbind(L1,L2) ,
                         th = thres,
                         names = colnames(cbind(L1,L2)),
                         name1 = colnames(L1),
                         name2= colnames(L2),
                         type='up')
}
common_bottom_genes <- foreach(thres = thresholds) %dopar% {
  NumberOfCommonElements(x = cbind(L1,L2) ,
                         th = thres,
                         names = colnames(cbind(L1,L2)),
                         name1 = colnames(L1),
                         name2= colnames(L2),
                         type='bottom')
}

# Transform list to array
common_top_matrix <- do.call(cbind,common_top_genes)
common_top_matrix <- t(common_top_matrix)
common_top_matrix <- as.data.frame(common_top_matrix) %>% rownames_to_column('Var2')
common_top_matrix$Var2 <- rep(colnames(L2),length(thresholds))
numbers <- NULL
for (i in 1:length(thresholds)){
  numbers <- c(numbers,rep(thresholds[i],length(colnames(L2))))
}
common_top_matrix$number <- numbers
common_top_matrix <- common_top_matrix %>% gather('Var1','num_genes',-Var2,-number)
common_top_matrix <- common_top_matrix %>% mutate(position='top')
### Repeat for bottom genes
common_bot_matrix <- do.call(cbind,common_bottom_genes)
common_bot_matrix <- t(common_bot_matrix)
common_bot_matrix <- as.data.frame(common_bot_matrix) %>% rownames_to_column('Var2')
common_bot_matrix$Var2 <- rep(colnames(L2),length(thresholds))
numbers <- NULL
for (i in 1:length(thresholds)){
  numbers <- c(numbers,rep(thresholds[i],length(colnames(L2))))
}
common_bot_matrix$number <- numbers
common_bot_matrix <- common_bot_matrix %>% gather('Var1','num_genes',-Var2,-number)
common_bot_matrix <- common_bot_matrix %>% mutate(position='bottom')

common_matrix <- rbind(common_top_matrix,common_bot_matrix)
df_common <- common_matrix %>% mutate(keep = ifelse(grepl('MPS',Var1) & grepl('MPS',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
df_common <- df_common %>% mutate(keep = ifelse(grepl('Wang',Var1) & grepl('Wang',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
common_distribution <- df_common$num_genes
hist(100*common_distribution,60)
colnames(df_common)[c(1,3)] <- c('Wang','MPS')
df_common$MPS <- as.character(df_common$MPS)
df_common$Wang <- as.character(df_common$Wang)
df_common <- df_common %>% mutate(MPS=strsplit(MPS,'_')) %>% unnest(MPS) %>% filter(MPS!='MPS') %>% filter(MPS!='Wang') %>%
  mutate(Wang=strsplit(Wang,'_')) %>% unnest(Wang) %>% filter(Wang!='Wang')%>% filter(Wang!='MPS')

### Retrieve all GSEA distances
new_dist_list <- NULL
for (i in 1:length(dist_gsea)){
  tmp <- dist_gsea[[i]]
  tmp <- tmp[colnames(L1),colnames(L2)]
  new_dist_list[[i]] <- tmp
}
dist_matrix <- do.call(cbind,new_dist_list)
dist_matrix <- t(dist_matrix)
dist_matrix <- as.data.frame(dist_matrix) %>% rownames_to_column('Var2')
dist_matrix$Var2 <- rep(colnames(L2),length(thresholds))
numbers <- NULL
for (i in 1:length(thresholds)){
  numbers <- c(numbers,rep(thresholds[i],length(colnames(L2))))
}
dist_matrix$number <- numbers
dist_matrix <- dist_matrix %>% gather('Var1','gsea_dist',-Var2,-number)
df_GSEA_all <- dist_matrix %>% mutate(keep = ifelse(grepl('MPS',Var1) & grepl('MPS',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
df_GSEA_all <- df_GSEA_all %>% mutate(keep = ifelse(grepl('Wang',Var1) & grepl('Wang',Var2),FALSE,TRUE)) %>% 
  filter(keep==TRUE) %>% select(-keep) %>% unique()
hist(df_GSEA_all$gsea_dist,60)
colnames(df_GSEA_all)[c(1,3)] <- c('Wang','MPS')
df_GSEA_all$MPS <- as.character(df_GSEA_all$MPS)
df_GSEA_all$Wang <- as.character(df_GSEA_all$Wang)
df_GSEA_all <- df_GSEA_all %>% mutate(MPS=strsplit(MPS,'_')) %>% unnest(MPS) %>% filter(MPS!='MPS') %>% filter(MPS!='Wang') %>%
  mutate(Wang=strsplit(Wang,'_')) %>% unnest(Wang) %>% filter(Wang!='Wang') %>% filter(Wang!='MPS')
df_GSEA_all <- left_join(df_GSEA_all,df_common)
# df_GSEA_all$number <- as.numeric(df_GSEA_all$number)

ggplot(df_GSEA_all,aes(x=gsea_dist/2,y=100 * num_genes,color=as.factor(number))) + geom_smooth()+
  geom_vline(xintercept = 0.33,linetype='dashed',color='black',linewidth=2)+
  xlab('GSEA-distance') + ylab('common (%)')+ labs(color='Number of genes')+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text= element_text(family = 'Arial',size = 24))+
  facet_wrap(~position)


ggplot(df_GSEA_all %>% filter(number>50),aes(x=gsea_dist/2,y=num_genes*number,color=as.factor(number))) + geom_smooth()+
  scale_y_continuous(breaks = c(0,5,10,15,20,25))+
  geom_vline(xintercept = 0.33,linetype='dashed',color='black',linewidth=2)+
  xlab('GSEA-distance') + ylab('number of common')+ labs(color='Number of genes')+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text= element_text(family = 'Arial',size = 24),
        panel.grid.major.y = element_line())+
  facet_wrap(~position)
