library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(pls)
library(pheatmap)
library(corrplot)
library(caret)
library(infotheo)
library(glmnet)
library(pls)
source("../modeling/aux_functions.R")
source("../modeling/functions_translation_nikos.R")

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

# combine with new data perturbed
data_lean <- data.table::fread('../results/optimized_mps/pertubed_control_lean.csv')
data_lean <- data_lean %>% column_to_rownames('V1')
data_fatty <- data.table::fread('../results/optimized_mps/pertubed_control_fatty.csv')
data_fatty <- data_fatty %>% column_to_rownames('V1')
data_fatty <- data_fatty[,paste0('perturbed',seq(0,11))]
data_lean <- data_lean[,paste0('perturbed',seq(0,11))]
data_lean <- t(data_lean)
data_fatty <- t(data_fatty)
print(all(colnames(data_fatty)==colnames(X_B)))
print(all(colnames(data_lean)==colnames(X_B)))
X_B <- rbind(X_B,data_fatty)
X_B <- rbind(X_B,data_lean)

# train translation model
df_all <- data.frame()
results_folds <- list()
valR2 <- NULL
valTau <- NULL
valMAE <- NULL
trainR2 <- NULL
trainTau <- NULL
trainMAE <- NULL
num_folds <- 5
for (i in 1:num_folds){
  train_inds <- readRDS(paste0('../preprocessing/TrainingValidationData/10fold_PCA_lasso/train_indices_',i,'.rds'))
  train_XA <- X_A[train_inds,]
  train_YA <- Y_A[train_inds]
  val_XA <- X_A[-train_inds,]
  val_YA <- Y_A[-train_inds]
  res <- translation_model(train_XA, X_B, train_YA,method = 'lmPCA',lasso=TRUE)
  trainR2[i] <- res$train_r2[1]
  trainTau[i] <- res$train_tau[1]
  trainMAE[i] <- res$tain_mae
  results_folds[[i]] <- res
  
  # validate
  model <- res$model_A_B
  y <- data.frame(out=val_YA)
  observed <- y$out
  X_A_B <- val_XA %*% res$PC_B$rotation
  predicted <- predict(model,newx = X_A_B)
  valMAE[i] <- mean(abs(observed - predicted))
  valTau[i] <- cor(predicted,observed,method = "kendall")
  valR2[i] <- cor(predicted,observed,method = 'pearson')^2
}
df <- data.frame(trainMAE,trainTau,trainR2,valMAE,valTau,valR2,fold = seq(1,num_folds))
df <- df %>% mutate(model = 'new model')
saveRDS(df,'../results/optimized_mps/df_all_with_perturbed.rds')

### Visualize performance comparison
df_all <- readRDS('../results/TransCompR_PCA/df_all_long.rds')
df <- df %>% gather('metric','value',-fold,-model) 
df <- df %>% mutate(set=ifelse(grepl('train',metric),'train','validation'))
df <- df %>% mutate(metric=ifelse(grepl('train',metric),substr(metric,6,nchar(metric)),substr(metric,4,nchar(metric))))

df_all <- rbind(df_all,df)

df_all$model <- factor(df_all$model,levels = c('lasso on genes','lasso on human PCs','model','new model','shuffled human genes'))


p1 <- ggboxplot(df_all %>% filter(model %in% c('model','new model')) %>% filter(set=='validation') %>% filter(metric=='MAE'),x='model',y='value',color = 'model',add='jitter') +
  ylab('MAE')+ xlab('')+
  stat_compare_means(comparisons = list(c('model','new model')),
                     method = 'wilcox.test',
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'))+
  coord_flip()
print(p1)

p2 <- ggboxplot(df_all %>% filter(model %in% c('model','new model')) %>% filter(set=='validation') %>% filter(metric=='Tau'),x='model',y='value',color = 'model',add='jitter') +
  ylab('kendall`s correlation')+xlab('')+
  stat_compare_means(comparisons = list(c('model','new model')),
                     method = 'wilcox.test',
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'))+
  coord_flip()
print(p2)

p3 <- ggboxplot(df_all %>% filter(model %in% c('model','new model')) %>% filter(set=='validation') %>% filter(metric=='R2'),x='model',y='value',color = 'model',add='jitter') +
  ylab('R2')+xlab('')+
  stat_compare_means(comparisons = list(c('model','new model')),
                     method = 'wilcox.test',
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'))+
  coord_flip()
print(p3)

p <- p1 +p3
print(p)

ggsave('../results/optimized_mps/allgenes_vs_translated_performance_perturbed_model.png',
       plot=p,
       width = 16,
       height = 8,
       units = 'in',
       dpi = 600)
