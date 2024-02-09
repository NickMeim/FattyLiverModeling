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
source("aux_functions.R")
source("functions_translation_nikos.R")
library(plotly)
library(sparsepca)

## Palatalization packages
library(doFuture)
# parallel set number of workers
# cores <- 8
# registerDoFuture()
# plan(multisession,workers = cores)

#### Load pre-processed data----------------------
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_A <- data$data_A
data_B <- data$data_B
Y_A <- data$Y_A

# Get mutliple phenotypes
pheno <- data.frame(NAS=Y_A)
rownames(pheno) <- colnames(data_A)
pheno <- pheno %>% rownames_to_column('sample')
patients_meta_data <- data.table::fread('../data/human_metadata_hoang.txt')
patients_meta_data <- patients_meta_data %>% 
  dplyr::select(c('sample'='GEO_Accession (exp)'),
                c('fibrosis'='Fibrosis_stage'),
                c('cytological_ballooning'='cytological_ballooning_grade'),
                #c('lobular_inflammation'='lobular_inflammation_grade'),
                c('steatosis'='steatosis_grade')) %>% unique()
pheno <- left_join(pheno,patients_meta_data)
pheno  <- pheno  %>% column_to_rownames('sample') #%>% select(-NAS)
Y_A <- as.matrix(pheno)

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

# # combine with new data perturbed
# data_lean <- data.table::fread('pertubed_control_lean.csv')
# data_lean <- data_lean %>% column_to_rownames('V1')
# data_fatty <- data.table::fread('pertubed_control_fatty.csv')
# data_fatty <- data_fatty %>% column_to_rownames('V1')
# data_fatty <- data_fatty[,paste0('perturbed',seq(0,11))]
# data_lean <- data_lean[,paste0('perturbed',seq(0,11))]
# data_lean <- t(data_lean)
# data_fatty <- t(data_fatty)
# print(all(colnames(data_fatty)==colnames(X_B)))
# print(all(colnames(data_lean)==colnames(X_B)))
# X_B <- rbind(X_B,data_fatty)
# X_B <- rbind(X_B,data_lean)

# save results to use in python
chip_meta_data <- read.table('../data/GSE168285/GSE168285_sample_meta_data.txt',header = TRUE, sep = "\t")
chip_meta_data <- chip_meta_data %>% filter(Number_of_cues==0) %>% filter(Background=='lean') %>% filter(NPC=='low') %>% 
  dplyr::select(sampleName,condition,treatment) %>% unique()
# data.table::fwrite(chip_meta_data,'../preprocessing/chip_lean_controls_indices.csv')
chip_meta_data <- read.table('../data/GSE168285/GSE168285_sample_meta_data.txt',header = TRUE, sep = "\t")
chip_meta_data <- chip_meta_data %>% filter(Number_of_cues==0) %>% filter(Background!='lean') %>% filter(NPC=='low') %>% 
  dplyr::select(sampleName,condition,treatment) %>% unique()
# data.table::fwrite(chip_meta_data,'chip_fatty_controls_indices.csv')
df_X_A = as.data.frame(X_A)
df_X_B = as.data.frame(X_B)
df_Y_A = as.data.frame(Y_A)
rownames(df_Y_A) <- rownames(X_A)
# data.table::fwrite(df_X_A,'../preprocessing/X_A.csv',row.names = T)
# data.table::fwrite(df_X_B,'../preprocessing/X_B.csv',row.names = T)
# data.table::fwrite(df_Y_A,'../preprocessing/Y_A.csv',row.names = T)

# ### Split in k fold cross-validation and save the data in folder----------------
num_folds <- 10
# min_outcomes <- 1
# min_steat <- 1
# thresh <-  3
# max_iter <- 10000
# iter <- 0
# while (min_outcomes<thresh | min_steat<thresh+1){
#   data_splits <- createMultiFolds(y = Y_A[,1], k = num_folds, times = 1)
#   i <- 1
#   for (ind in data_splits){
#     val_Y <- Y_A[-ind,]
#     if ((is.matrix(val_Y) | is.data.frame(val_Y))){
#       val_Y <- apply(val_Y,2,function(x) {return(length(unique(c(x))))})
#       min_outcomes <- min(val_Y)
#       min_steat <- val_Y[3]
#       saveRDS(ind,paste0('../preprocessing/TrainingValidationData/10fold_lasso/train_indices_',i,'.rds'))
#     }
#     # print(paste0('Finished fold ',i))
#     i <- i+1
#   }
#   iter <- iter +1
#   if (iter==max_iter){
#     iter <- 0
#     thresh <- thresh -1
#     print('Maximum iterations reached, decreased the threshold')
#   }
#   if (iter %% 1000 == 0 ){
#     message(paste0('Trials = ',iter))
#   }
# }


#### First perform tuning of SPARSE PCA-----------------------------------------------------------------------
alphas <- c(10^seq(-6,0),5*10^seq(-6,-1))
# Sparse PCA for multiple alphas in parallel
spca_results <- foreach(a = alphas) %dopar% {
  spca(X_B,alpha = a,verbose = FALSE)
}
# saveRDS(spca_results,'../results/TransCompR_sparsePCA/spca_results.rds')
captured_variance <- NULL
for (i in 1:length(alphas)){
  spca_B <- spca_results[[as.character(alphas[i])]]
  captured_variance[i] <- cumsum(spca_B$sdev^2/spca_B$var)[length(spca_B$sdev)]
}
df_spca_res <- data.frame(alpha = alphas,variance = captured_variance)
p <- ggplot(df_spca_res,aes(x=alphas,y=100*variance)) + 
  geom_point(size=3)+
  geom_line()+
  geom_hline(yintercept = 90,linetype='dashed',color='black',lwd = 2)+
  scale_x_continuous(breaks =10^seq(-6,0) ,
                     minor_breaks = do.call(c,lapply(10^seq(-6,-1),function(x) x*seq(2,9))),
                     trans = 'log10')+
  scale_y_continuous(n.breaks = 6)+
  xlab(expression(alpha))+
  ylab('% of captured variance')+
  ggtitle('Total captured variance by the sparse Principal Components')+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size = 24),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_line(linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.5))
print(p)
ggsave('../results/TransCompR_sparsePCA/alpha_optimize.png',
       plot=p,
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)
ggsave('../results/TransCompR_sparsePCA/alpha_optimize.eps',
       device = cairo_ps,
       plot=p,
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)
optimal_alpha <- max(df_spca_res$alpha[which(df_spca_res$variance>=0.9)])
print(optimal_alpha) # it is 1e-05
# saveRDS(spca_results[['1e-05']],'../results/TransCompR_sparsePCA/optimal_sPC_B.rds')

### Check if this optimal can also work for human data
sPC_A <- spca(X_A,alpha = 1e-05,verbose = TRUE,max_iter = 2000)
plot(log(sPC_A$objective), col='red', xlab='Number of iterations', ylab='Objective value')
print(paste0('Total variance captured : ',cumsum(sPC_A$sdev^2/sPC_A$var)[length(sPC_A$sdev)]))
print(paste0('Variance captured by sPCs 1-3: ',cumsum(sPC_A$sdev^2/sPC_A$var)[3]))
# saveRDS(sPC_A,'../results/TransCompR_sparsePCA/optimal_sPC_A.rds')

### First load saved pre-calculated sparse PCAs--------------------------------------------
spca_results <- readRDS('../results/TransCompR_sparsePCA/spca_results.rds')
# sPC_B <- spca_results[['1e-04']]
sPC_B <- readRDS('../results/TransCompR_sparsePCA/optimal_sPC_B.rds')
sPC_A <- readRDS('../results/TransCompR_sparsePCA/optimal_sPC_A.rds')

### Re-do comparing linear model with lasso and shuffled features---------------------------
df_all <- data.frame()
results_folds <- list()
valR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
all_laso_models <- list()
for (i in 1:num_folds){
  train_inds <- readRDS(paste0('train_indices_',i,'.rds'))
  train_XA <- X_A[train_inds,]
  train_YA <- Y_A[train_inds,]
  val_XA <- X_A[-train_inds,]
  val_YA <- Y_A[-train_inds,]
  res <- sparse_translation_model_multi(sPC_A, sPC_B, train_XA, train_YA,val_XA,val_YA)
  trainR2[i,] <- cbind(res$train_r^2)
  trainTau[i,] <- cbind(res$train_tau)
  trainMAE[i,] <- cbind(res$tain_mae)
  valMAE[i,] <- cbind(res$val_mae)
  valTau[i,] <- cbind(res$val_tau)
  valR2[i,] <- cbind(res$val_r^2)
  results_folds[[i]] <- res
  
  saveRDS(results_folds,'../results/TransCompR_sparsePCA/all_models.rds')
  
  print(paste0('Finished fold ',i))
}
colnames(trainMAE) <- colnames(train_YA)
colnames(trainTau) <- colnames(train_YA)
colnames(trainR2) <- colnames(train_YA)
colnames(valMAE) <- colnames(val_YA)
colnames(valTau) <- colnames(val_YA)
colnames(valR2) <- colnames(val_YA)
df <- rbind(as.data.frame(trainMAE) %>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valMAE)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold))
df <- df %>% mutate(model = 'model')
df_all <- rbind(df_all,df)
### shuffled genes-----------------------------------------
valR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
for (i in 1:num_folds){
  train_inds <- readRDS(paste0('train_indices_',i,'.rds'))
  train_XA <- X_A[train_inds,]
  train_XA <- train_XA[,sample(ncol(train_XA))]
  train_YA <- Y_A[train_inds,]
  val_XA <- X_A[-train_inds,]
  val_YA <- Y_A[-train_inds,]
  res <- sparse_translation_model_multi(sPC_A, sPC_B, train_XA, train_YA,val_XA,val_YA)
  trainR2[i,] <- cbind(res$train_r^2)
  trainTau[i,] <- cbind(res$train_tau)
  trainMAE[i,] <- cbind(res$tain_mae)
  valMAE[i,] <- cbind(res$val_mae)
  valTau[i,] <- cbind(res$val_tau)
  valR2[i,] <- cbind(res$val_r^2)
  print(paste0('Finished fold ',i))
}
colnames(trainMAE) <- colnames(train_YA)
colnames(trainTau) <- colnames(train_YA)
colnames(trainR2) <- colnames(train_YA)
colnames(valMAE) <- colnames(val_YA)
colnames(valTau) <- colnames(val_YA)
colnames(valR2) <- colnames(val_YA)
df <- rbind(as.data.frame(trainMAE) %>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valMAE)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold))
df <- df %>% mutate(model = 'shuffled human genes')
df_all <- rbind(df_all,df)
saveRDS(df_all,'../results/TransCompR_sparsePCA/df_all.rds')

### Train lasso model using X_A-------------------------------------------
valR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
lambda_grid <- 10^seq(-2, 1, length = 100)
# for (i in 1:num_folds){
#   train_inds <- readRDS(paste0('train_indices_',i,'.rds'))
#   train_XA <- X_A[train_inds,]
#   train_YA <- Y_A[train_inds,]
#   val_XA <- X_A[-train_inds,]
#   val_YA <- Y_A[-train_inds,]
#   gene_mapping <- data.frame(gene = colnames(train_XA))
#   gene_mapping$id <- paste0('gene_',seq(1:ncol(train_XA))) 
#   colnames(train_XA) <- gene_mapping$id
#   colnames(val_XA) <- gene_mapping$id
#   ctrl <- trainControl(method = "cv", number = 10)
#   # data_modeling <- cbind(as.matrix(train_XA),apply(as.matrix(train_YA[,which(colnames(train_YA)!='NAS')]),c(1,2),as.numeric))
#   # data_modeling_val <- cbind(as.matrix(val_XA),apply(as.matrix(val_YA[,which(colnames(val_YA)!='NAS')]),c(1,2),as.numeric))
#   data_modeling <- cbind(as.matrix(train_XA),apply(as.matrix(train_YA),c(1,2),as.numeric))
#   data_modeling_val <- cbind(as.matrix(val_XA),apply(as.matrix(val_YA),c(1,2),as.numeric))
#   j <- 1
#   for (outVar in colnames(train_YA)){
#     form <- as.formula(paste0(outVar,'~.'))
#     lasso_model <- train(form,
#                          data = data_modeling[,c(colnames(train_XA),outVar)],
#                          method = 'glmnet',
#                          trControl = ctrl,
#                          metric = "MAE",
#                          tuneGrid = expand.grid(alpha = 1, lambda = lambda_grid))
#     if (j == 1){
#       predicted <- predict(lasso_model,data_modeling[,c(colnames(train_XA),outVar)])
#       if (sd(predicted)==0){
#         predicted <- predicted + 1e-8*rnorm(length(predicted))
#       }
#       # validate
#       predicted_val <- predict(lasso_model,data_modeling_val[,c(colnames(val_XA),outVar)])
#       if (sd(predicted_val)==0){
#         predicted_val <- predicted_val + 1e-8*rnorm(length(predicted_val))
#       }
#     }else{
#       if (sd(predict(lasso_model,data_modeling[,c(colnames(train_XA),outVar)]))==0){
#         predicted <- cbind(predicted, predict(lasso_model,data_modeling[,c(colnames(train_XA),outVar)])+ 1e-8*rnorm(nrow(data_modeling)))
#       }else{
#         predicted <- cbind(predicted, predict(lasso_model,data_modeling[,c(colnames(train_XA),outVar)]))
#       }
#       if (sd(predict(lasso_model,data_modeling_val[,c(colnames(val_XA),outVar)]))==0){
#         predicted_val <- cbind(predicted_val, predict(lasso_model,data_modeling_val[,c(colnames(val_XA),outVar)])+ 1e-8*rnorm(nrow(data_modeling_val)))
#       }else{
#         predicted_val <- cbind(predicted_val, predict(lasso_model,data_modeling_val[,c(colnames(val_XA),outVar)]))
#       }
#     }
#     message(paste0('Finished output ',outVar))
#     j <- j+1
#   }
#   colnames(predicted) <- colnames(train_YA)
#   colnames(predicted_val) <- colnames(val_YA)
#   # nas <- predicted[,c('cytological_ballooning')] + predicted[,c('steatosis')] + predicted[,c('lobular_inflammation')]
#   # nas_val <- predicted_val[,c('cytological_ballooning')] + predicted_val[,c('steatosis')] + predicted_val[,c('lobular_inflammation')]
#   # predicted <- cbind(predicted,nas)
#   # colnames(predicted)[ncol(predicted)] <- 'NAS'
#   # predicted_val <- cbind(predicted_val,nas_val)
#   # colnames(predicted_val)[ncol(predicted_val)] <- 'NAS'
#   # evaluate 
#   mae <- apply(abs(train_YA - predicted),2,mean)
#   tau <- cor(predicted,train_YA,method = 'kendall')
#   tau <- diag(tau)
#   r <- cor(predicted,train_YA,method = 'pearson')^2
#   r <- diag(r)
#   trainR2[i,] <- r
#   trainTau[i,] <- tau
#   trainMAE[i,] <- mae
#   # validate
#   valMAE[i,] <- apply(abs(val_YA - predicted_val),2,mean)
#   valTau[i,] <- diag(cor(predicted_val,val_YA,method = 'kendall'))
#   valR2[i,] <- diag(cor(predicted_val,val_YA,method = 'pearson')^2)
#   print(paste0('Finished fold ',i))
# }
# colnames(trainMAE) <- colnames(predicted)
# colnames(trainR2) <- colnames(predicted)
# colnames(trainTau) <- colnames(predicted)
# colnames(valMAE) <- colnames(predicted_val)
# colnames(valTau) <- colnames(predicted_val)
# colnames(valR2) <- colnames(predicted_val)
# df <- rbind(as.data.frame(trainMAE) %>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='MAE') %>%
#               gather('output','value',-metric,-set,-fold),
#             as.data.frame(trainTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='Tau') %>%
#               gather('output','value',-metric,-set,-fold),
#             as.data.frame(trainR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='R2') %>%
#               gather('output','value',-metric,-set,-fold),
#             as.data.frame(valMAE)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='MAE') %>%
#               gather('output','value',-metric,-set,-fold),
#             as.data.frame(valTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='Tau') %>%
#               gather('output','value',-metric,-set,-fold),
#             as.data.frame(valR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='R2') %>%
#               gather('output','value',-metric,-set,-fold))
# df <- df %>% mutate(model = 'lasso on genes')
# # saveRDS(df,'../results/TransCompR_sparsePCA/df_genes_human_lasso.rds')
df <- readRDS('../results/TransCompR_sparsePCA/df_genes_human_lasso.rds')
df_all <- rbind(df_all,df)
saveRDS(df_all,'../results/TransCompR_sparsePCA/df_all_long.rds')

### Train in human sPCs-------------------------------------------
valR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
valMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainR2 <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainTau <- matrix(nrow= num_folds,ncol=ncol(Y_A))
trainMAE <- matrix(nrow= num_folds,ncol=ncol(Y_A))
for (i in 1:num_folds){
  train_inds <- readRDS(paste0('train_indices_',i,'.rds'))
  train_XA <- X_A[train_inds,]
  train_YA <- Y_A[train_inds,]
  val_XA <- X_A[-train_inds,]
  val_YA <- Y_A[-train_inds,]
  res <- sparse_translation_model_multi(sPC_A, sPC_A, train_XA, train_YA,val_XA,val_YA)
  trainR2[i,] <- cbind(res$train_r^2)
  trainTau[i,] <- cbind(res$train_tau)
  trainMAE[i,] <- cbind(res$tain_mae)
  valMAE[i,] <- cbind(res$val_mae)
  valTau[i,] <- cbind(res$val_tau)
  valR2[i,] <- cbind(res$val_r^2)
  print(paste0('Finished fold ',i))
}
colnames(trainMAE) <- colnames(train_YA)
colnames(trainTau) <- colnames(train_YA)
colnames(trainR2) <- colnames(train_YA)
colnames(valMAE) <- colnames(val_YA)
colnames(valTau) <- colnames(val_YA)
colnames(valR2) <- colnames(val_YA)
df <- rbind(as.data.frame(trainMAE) %>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='train') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valMAE)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valTau)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(valR2)%>% mutate(fold = seq(1,num_folds)) %>% mutate(set='validation') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold))
df <- df %>% mutate(model = 'model on human sPCs')
saveRDS(df,'../results/TransCompR_sparsePCA/df_performance_humanSPC.rds')

### Final visualize-----------------------------------------------------------
df_all <- readRDS('../results/TransCompR_sparsePCA/df_all_long.rds')
df_performance_human <- readRDS('../results/TransCompR_sparsePCA/df_performance_humanSPC.rds')
df_all <- rbind(df_all,df_performance_human)

df_all$model <- factor(df_all$model,levels = c('lasso on genes','model on human sPCs','model','shuffled human genes'))

p1 <- ggboxplot(df_all %>% filter(set=='validation') %>% filter(metric=='MAE'),x='model',y='value',color = 'model',add='jitter') +
  ylab('MAE')+ xlab('')+
  stat_compare_means(comparisons = list(c('model','model on human sPCs'),c('model','lasso on genes'),c('model','shuffled human genes')),
                     method = 'wilcox.test',
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'))+
  facet_wrap(~output)+
  coord_flip()
print(p1)

p2 <- ggboxplot(df_all %>% filter(set=='validation') %>% filter(metric=='Tau'),x='model',y='value',color = 'model',add='jitter') +
  ylab('kendall`s correlation')+xlab('')+
  stat_compare_means(comparisons = list(c('model','model on human sPCs'),c('model','lasso on genes'),c('model','shuffled human genes')),
                     method = 'wilcox.test',
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'))+
  facet_wrap(~output)+
  coord_flip()
print(p2)

p3 <- ggboxplot(df_all %>% filter(set=='validation') %>% filter(metric=='R2'),x='model',y='value',color = 'model',add='jitter') +
  ylab('R2')+xlab('')+
  stat_compare_means(comparisons = list(c('model','model on human sPCs'),c('model','lasso on genes'),c('model','shuffled human genes')),
                     method = 'wilcox.test',
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'))+
  facet_wrap(~output)+
  coord_flip()
print(p3)

p <- p1 / p3
p <- p +  plot_layout(guides = "collect") & theme(legend.position = 'top')
print(p)
ggsave('../results/TransCompR_sparsePCA/allgenes_vs_translated_performance.png',
       plot=p,
       width = 16,
       height = 8,
       units = 'in',
       dpi = 600)
ggsave('../results/TransCompR_sparsePCA/allgenes_vs_translated_mae.png',
       plot=p1,
       width = 16,
       height = 8,
       units = 'in',
       dpi = 600)

### PCs to predict NAS score--------------------
num_folds <- 10
results <- readRDS('../results/TransCompR_sparsePCA/all_models.rds')
df_res <- data.frame()
for (i in 1:length(results)){
  for (out in c("NAS","fibrosis","steatosis","cytological_ballooning")){
    pcs <- results[[i]][["lm_pcs"]][[out]]
    df_res <- rbind(df_res,data.frame(PC=pcs) %>% mutate(output=out))
  }
}
df_res <- df_res %>% group_by(output,PC) %>% mutate(counts = n())

## Check how many times can a PC appear by chance
# parallel set number of workers
cores <- 10
registerDoFuture()
plan(multisession,workers = cores)

df_res_random <- data.frame()
res_random <- foreach(i = seq(1,num_folds)) %dopar% {
  train_inds <- readRDS(paste0('train_indices_',i,'.rds'))
  train_XA <- X_A[train_inds,]
  train_YA <- Y_A[train_inds,]
  val_XA <- X_A[-train_inds,]
  val_YA <- Y_A[-train_inds,]
  X_B_star <- X_B[,sample(ncol(X_B))]
  sPC_B_star <- spca(X_B_star,alpha =  1e-05,verbose = F)
  sparse_translation_model_multi(sPC_A, sPC_B_star, train_XA,train_YA,val_XA,val_YA)
}
for (i in 1:num_folds){
  res <- res_random[[i]]
  for (out in c("NAS","fibrosis","steatosis","cytological_ballooning")){
    pcs <- res[["lm_pcs"]][[out]]
    df_res_random <- rbind(df_res_random,data.frame(PC=pcs) %>% mutate(output=out))
  }
  print(paste0('Finished fold ',i))
}
df_res_random <- df_res_random %>% group_by(output,PC) %>% mutate(counts = n())
saveRDS(df_res_random,'../results/TransCompR_sparsePCA/counts_random_pcs.rds')
df_res_random <- readRDS('../results/TransCompR_sparsePCA/counts_random_pcs.rds')

## Filter PCs
df_thresh <- data.frame()
for (out in c("NAS","fibrosis","steatosis","cytological_ballooning")){
  threshold <- apply(t(c(1,2,3,4,5,6,7,8,9,10)), 2, FUN = function(x) {
    return(sum(df_res_random$counts[which(df_res_random$output==out)]>=x)/length(which(df_res_random$output==out)))})
  df_thresh <- rbind(df_thresh,data.frame(threshold) %>% mutate(output=out))
}
# print(threshold)
df_thresh <- df_thresh %>% group_by(output) %>% 
  mutate(p.adj = p.adjust(threshold)) %>% ungroup()
# threshold <- which(threshold<0.05)[1]
df_thresh <- df_thresh %>% mutate(counts = rep(seq(1,10),4)) %>% filter(p.adj<=0.05)
df_thresh <- df_thresh %>% group_by(output) %>% mutate(min_count = min(counts)) %>% ungroup() %>% 
  filter(counts==min_count) %>% select(output,counts)%>% unique()
colnames(df_thresh)[2] <- 'threshold'
df_res <- left_join(df_res,df_thresh)
df_res <- df_res %>% filter(counts>=threshold)
colnames(df_res)[1] <- 'sPC'
print(unique(df_res$sPC))
Zh <- X_A %*% sPC_B$loadings
# Zfilt <- Zh[,pcs_to_keep]
varAB <- apply(Zh,2,sd)^2
varA <- sPC_A$var
perc <- 100*varAB/varA
colnames(Zh) <- paste0('sPC',seq(1:ncol(Zh)))
names(perc) <- colnames(Zh)
perc <- perc[order(-perc)][1:15]
perc <- as.data.frame(perc)
perc <- perc %>% rownames_to_column('sPC')
p1 <- ggbarplot(perc,x='sPC',y='perc',fill='purple') + ylab('% variance explained') + ggtitle('Percentage of variance explained of the human data')+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))

Zb <- sPC_B$scores
colnames(Zb) <- paste0('sPC',seq(1:ncol(Zb)))
varB <- sPC_B$var
# Zbfilt <- Zb[,pcs_to_keep]
varPCB <- apply(Zb,2,sd)^2
perc2 <- 100*varPCB/varB
perc2 <- perc2[perc$sPC]
perc2 <- as.data.frame(perc2)
perc2 <- perc2 %>% rownames_to_column('sPC')
p2 <- ggbarplot(perc2,x='sPC',y='perc2',fill='blue') + ylab('% variance explained') + ggtitle('Percentage of variance explained of the MPS data')+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
p <- p1/p2
print(p)

ggsave('../results/TransCompR_sparsePCA/allgenes_vs_translated_perc_variance.png',
       plot=p,
       width = 16,
       height = 8,
       units = 'in',
       dpi = 600)

### Correlation between PCs based on projected human samples
cor_mat <- cor(Zh)
cor_test_results <- apply(combn(ncol(Zh), 2), 2, function(pair) {
  cor_test_result <- cor.test(Zh[, pair[1]], Zh[, pair[2]])
  data.frame(
    variable1 = colnames(Zh)[pair[1]],
    variable2 = colnames(Zh)[pair[2]],
    correlation = cor_mat[pair[1], pair[2]],
    p_value = cor_test_result$p.value
  )
})
# Combine the results into a data frame
cor_test_results_df <- do.call(rbind, cor_test_results)
# cor_test_results_df$p.adj  <-  p.adjust(cor_test_results_df$p_value,method = 'fdr')
cor_test_results_df <- cor_test_results_df %>% 
  mutate(significance = ifelse(p_value<0.001,'***',
                               ifelse(p_value<0.01,'**',
                                      ifelse(p_value<0.05,'*',
                                             ifelse(p_value<0.1,'.','')))))
perc <- 100*varAB/varA
names(perc) <- colnames(Zh)
perc <- perc[order(-perc)][1:5]
spcs_to_see <- names(perc)
spcs_to_see <- union(unique(df_res$sPC),spcs_to_see)
inds <- sapply(spcs_to_see, substr,start = 4,stop = 5)
names(inds) <- NULL
inds <- as.numeric(inds)
spcs_to_see <- spcs_to_see[order(inds)]
cor_test_results_df <- cor_test_results_df %>% filter(variable1 %in% spcs_to_see) %>% filter(variable2 %in% spcs_to_see)
cor_test_results_df$variable1 <- factor(cor_test_results_df$variable1,levels = spcs_to_see)
cor_test_results_df$variable2 <- factor(cor_test_results_df$variable2,levels = spcs_to_see)
ggplot(cor_test_results_df, aes(variable1, variable2, fill = correlation))+
  geom_tile(color='black')+
  geom_text(aes(label=significance),size=3)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) + 
  ylab('') + xlab('')
ggsave('../results/TransCompR_sparsePCA/selected_pcs_correlation.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)

### visualize only PC12 vs NAS
projected <- as.data.frame(Zh) %>% rownames_to_column('sample') %>% gather('sPC','value',-sample)
Y_data <- as.data.frame(Y_A) %>% mutate(sample=rownames(Zh))
colnames(Y_data)[1] <- 'NAS'
projected <- left_join(projected,Y_data)
projected <- left_join(df_res %>% select(sPC,output) %>% unique(),projected)
keep_spc <- 'sPC12'
p1 <- ggscatter(projected %>% filter(sPC==keep_spc) %>% filter(output=='NAS') %>% unique(),x='NAS',y='value',
                cor.coef = T,cor.coef.size = 6) + 
  geom_smooth(se=TRUE,color='blue',method = 'lm') + 
  ylab(keep_spc)

p2 <- ggboxplot(projected %>% filter(sPC==keep_spc)%>% filter(output=='NAS') %>% unique(),x='NAS',y='value',color='NAS',add='jitter')  + 
  ylab(keep_spc) + theme(legend.position = 'none') +
  stat_compare_means()

p <- p1+p2
print(p)


keep_spcs <- unique(df_res$sPC)
p1_nas <- ggscatter(projected %>% filter(output=='NAS'),x='NAS',y='value',
                    cor.coef = T,cor.coef.size = 4) + 
  facet_wrap(~sPC)+
  geom_smooth(se=TRUE,color='blue',method = 'lm') 
p2_nas <- ggboxplot(projected%>% filter(output=='NAS'),x='NAS',y='value',color='NAS',add='jitter')  + 
  theme(legend.position = 'none') +
  facet_wrap(~sPC)+
  stat_compare_means(label.x = 2)
p_nas <- p1_nas+p2_nas
print(p_nas)
ggsave('../results/TransCompR_sparsePCA/important_pcs_vs_nas.png',
       plot = p_nas,
       width = 14,
       height = 9,
       units = 'in',
       dpi = 600)
### repeat fibrosis
p1_fibrosis <- ggscatter(projected %>% filter(output=='fibrosis'),x='fibrosis',y='value',
                         cor.coef = T,cor.coef.size = 4) + 
  facet_wrap(~sPC)+
  geom_smooth(se=TRUE,color='blue',method = 'lm') 
p2_fibrosis <- ggboxplot(projected%>% filter(output=='fibrosis'),x='fibrosis',y='value',color='fibrosis',add='jitter')  + 
  theme(legend.position = 'none') +
  facet_wrap(~sPC)+
  stat_compare_means(label.x = 2)
p_fibrosis <- p1_fibrosis+p2_fibrosis
print(p_fibrosis)
ggsave('../results/TransCompR_sparsePCA/important_pcs_vs_fibrosis.png',
       plot = p_fibrosis,
       width = 14,
       height = 9,
       units = 'in',
       dpi = 600)
### repeat steatosis
p1_steatosis <- ggscatter(projected %>% filter(output=='steatosis'),x='steatosis',y='value',
                          cor.coef = T,cor.coef.size = 4) + 
  facet_wrap(~sPC)+
  geom_smooth(se=TRUE,color='blue',method = 'lm') 
p2_steatosis <- ggboxplot(projected%>% filter(output=='steatosis'),x='steatosis',y='value',color='steatosis',add='jitter')  + 
  theme(legend.position = 'none') +
  facet_wrap(~sPC)+
  stat_compare_means(label.x = 2)
p_steatosis <- p1_steatosis+p2_steatosis
print(p_steatosis)
ggsave('../results/TransCompR_sparsePCA/important_pcs_vs_steatosis.png',
       plot = p_steatosis,
       width = 14,
       height = 9,
       units = 'in',
       dpi = 600)
### repeat cytological_ballooning
p1_balooning <- ggscatter(projected %>% filter(output=='cytological_ballooning'),x='cytological_ballooning',y='value',
                          cor.coef = T,cor.coef.size = 4) + 
  facet_wrap(~sPC)+
  geom_smooth(se=TRUE,color='blue',method = 'lm') 
p2_balooning <- ggboxplot(projected%>% filter(output=='cytological_ballooning'),x='cytological_ballooning',y='value',color='cytological_ballooning',add='jitter')  + 
  theme(legend.position = 'none') +
  facet_wrap(~sPC)+
  stat_compare_means(label.x = 1)
p_balooning <- p1_balooning+p2_balooning
print(p_balooning)
ggsave('../results/TransCompR_sparsePCA/important_pcs_vs_balooning.png',
       plot = p_balooning,
       width = 14,
       height = 9,
       units = 'in',
       dpi = 600)

## Visualize perturbations on PC12-----------------------------------------
# keep_spcs <- unique(df_res$sPC)
keep_spcs <- c('sPC1','sPC3','sPC5','sPC8','sPC11','sPC12','sPC14')
Zb <- sPC_B$scores
colnames(Zb) <- paste0('sPC',seq(1,ncol(Zb)))
# scores <- as.data.frame(Zb) %>% select(PC12) %>% rownames_to_column('sampleName')
scores <- as.data.frame(Zb) %>% select(all_of(keep_spcs)) %>% rownames_to_column('sampleName')
chip_meta_data <- read.table('../data/GSE168285/GSE168285_sample_meta_data.txt',header = TRUE, sep = "\t") %>% select(-X)
conditions <- data.table::fread('../data/GSE168285/sample_info_all.csv')
colnames(conditions)[1] <- 'sampleName'
colnames(conditions)[2] <- 'Background'
# scores <- left_join(scores,chip_meta_data)
scores <- left_join(scores,conditions)
scores <- left_join(scores,chip_meta_data %>% select(-Background,-NPC))
anov <- aov(cbind(sPC1,sPC3,sPC5,sPC8,sPC11,sPC12,sPC14)~Background+LPS+TGFb+Cholesterol+Fructose+NPC,data=scores)
summary(anov)
stats_results <- scores %>% select(all_of(keep_spcs),'Background','LPS','TGFb','Cholesterol','Fructose','NPC') %>%
  gather('PC','value',-Background,-LPS,-TGFb,-Cholesterol,-Fructose,-NPC) %>% group_by(PC) %>%
  rstatix::anova_test(value~Background+LPS+TGFb+Cholesterol+Fructose+NPC) %>% 
  rstatix::adjust_pvalue(method = 'BH') %>%
  ungroup()
coeffs <- coef(anov)
coeffs <- as.data.frame(coeffs)[2:nrow(coeffs),] %>% rownames_to_column('Effect') %>%
  gather('PC','coefficient',-Effect)
stats_results <- as.data.frame(stats_results)
stats_results <- stats_results %>% 
  mutate(significance = ifelse(p.adj<0.001,'***',
                               ifelse(p.adj<0.01,'**',
                                      ifelse(p.adj<0.05,'*',
                                             ifelse(p.adj<0.1,'.','')))))
stats_results <- left_join(stats_results,coeffs)

p_F <- ggplot(stats_results,aes(x=PC,y=Effect,fill=`F`)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient(low = 'white',high='red')+
  labs(fill = "F-statistic") + 
  ggtitle('ANOVA results')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
p_coeff <- ggplot(stats_results,aes(x=PC,y=Effect,fill=coefficient)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient2(low = 'blue',high='red',mid = 'white',midpoint = 0)+
  ggtitle('Linear model coefficients')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
p <- p_F + p_coeff
print(p)
ggsave('../results/TransCompR_sparsePCA/anova_stimuli.png',
       plot=p,
       height = 6,
       width = 12,
       units = 'in',
       dpi = 600)

### Print PC11-PC12-PC14--------------------------------------------------
projected <- as.data.frame(Zh) %>% rownames_to_column('sample') %>% gather('sPC','value',-sample)
projection_11_12_14 <- projected %>% filter(sPC %in% c('sPC11','sPC12','sPC14'))
patients_meta_data <- data.table::fread('../data/human_metadata_hoang.txt')
patients_meta_data <- patients_meta_data %>% 
  dplyr::select(c('sample'='GEO_Accession (exp)'),c('fibrosis'='Fibrosis_stage'),c('age'='age_at_biopsy'),sex) %>% unique()
projection_11_12_14 <- left_join(projection_11_12_14,patients_meta_data)
projection_11_12_14 <- projection_11_12_14 %>% spread('sPC','value')

p_11_12 <- ggplot(projection_11_12_14,aes(x=sPC11,y=sPC12,fill=age,shape=sex))+
  geom_point(colour="black", size=5)+
  scale_shape_manual(values = c(21,22))+
  scale_fill_gradient2(low = 'blue',mid = 'white',high='red',midpoint = 52)+
  theme_pubr(base_size = 24,base_family = 'Arial')+
  theme(text=element_text(size = 24,family = 'Arial'),
        legend.position = 'right')
p_11_13 <- ggplot(projection_11_12_14,aes(x=sPC11,y=sPC14,fill=age,shape=sex))+
  geom_point(colour="black", size=5)+
  scale_shape_manual(values = c(21,22))+
  scale_fill_gradient2(low = 'blue',mid = 'white',high='red',midpoint = 52)+
  theme_pubr(base_size = 24,base_family = 'Arial')+
  theme(text=element_text(size = 24,family = 'Arial'),
        legend.position = 'right')
p_12_13 <- ggplot(projection_11_12_14,aes(x=sPC12,y=sPC14,fill=age,shape=sex))+
  geom_point(colour="black", size=5)+
  scale_shape_manual(values = c(21,22))+
  scale_fill_gradient2(low = 'blue',mid = 'white',high='red',midpoint = 52)+
  theme_pubr(base_size = 24,base_family = 'Arial')+
  theme(text=element_text(size = 24,family = 'Arial'),
        legend.position = 'right')
p_all <- p_11_12 + p_12_13 + p_11_13 + plot_layout(guides = 'collect')
print(p_all)
ggsave('../results/TransCompR_sparsePCA/spc_11_12_14_vs_covariates.png',
       plot=p_all,
       height = 8,
       width = 16,
       units = 'in',
       dpi = 600)

#### Anova for hidden covariates
keep_spcs <- c('sPC1','sPC3','sPC5','sPC8','sPC11','sPC12','sPC14')
patients_meta_data <- data.table::fread('../data/human_metadata_hoang.txt')
patients_meta_data <- patients_meta_data %>% 
  dplyr::select(c('sample'='GEO_Accession (exp)'),
                c('fibrosis'='Fibrosis_stage'),
                c('age'='age_at_biopsy'),
                c('cytological_ballooning'='cytological_ballooning_grade'),
                c('lobular_inflammation'='lobular_inflammation_grade'),
                c('steatosis'='steatosis_grade'),
                c('NAS'='nafld_activity_score'),
                sex) %>% unique()
projected_for_anova <- left_join(projected %>% filter(sPC %in% keep_spcs),patients_meta_data)
projected_for_anova <- projected_for_anova %>% spread('sPC','value')
anov_projected <- aov(cbind(sPC1,sPC3,sPC5,sPC8,sPC11,sPC12,sPC14)~age+sex+cytological_ballooning+lobular_inflammation+steatosis+fibrosis,
                      data=projected_for_anova %>% mutate(sex = ifelse(sex=='male',1,0)))
summary(anov_projected)
stats_results <- projected_for_anova  %>% mutate(sex = ifelse(sex=='male',1,0)) %>% select(-sample) %>%
  gather('sPC','value',-age,-fibrosis,-sex,-NAS,-cytological_ballooning,-lobular_inflammation,-steatosis) %>% group_by(sPC) %>%
  rstatix::anova_test(value~age+sex+cytological_ballooning+lobular_inflammation+steatosis+fibrosis) %>% 
  rstatix::adjust_pvalue(method = 'BH') %>%
  ungroup()
coeffs <- coef(anov_projected)
coeffs <- as.data.frame(coeffs)[2:nrow(coeffs),] %>% rownames_to_column('Effect') %>%
  gather('sPC','coefficient',-Effect)
stats_results <- as.data.frame(stats_results)
# stats_results <- stats_results %>%
#   mutate(significance = ifelse(p.adj<0.001,'***',
#                                ifelse(p.adj<0.01,'**',
#                                       ifelse(p.adj<0.05,'*',
#                                              ifelse(p.adj<0.1,'.','')))))
stats_results <- stats_results %>%
  mutate(significance = ifelse(p<0.001,'***',
                               ifelse(p<0.01,'**',
                                      ifelse(p<0.05,'*',
                                             ifelse(p<0.1,'.','')))))
stats_results <- left_join(stats_results,coeffs)

p_F <- ggplot(stats_results,aes(x=sPC,y=Effect,fill=`F`)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient(low = 'white',high='red')+
  labs(fill = "F-statistic") + 
  ggtitle('ANOVA results')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
p_coeff <- ggplot(stats_results,aes(x=sPC,y=Effect,fill=coefficient)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient2(low = 'blue',high='red',mid = 'white',midpoint = 0)+
  ggtitle('Linear model coefficients')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
p <- p_F + p_coeff
print(p)
ggsave('../results/TransCompR_sparsePCA/uncorrected_anova_humas.png',
       plot=p,
       height = 6,
       width = 12,
       units = 'in',
       dpi = 600)

ggboxplot(projected_for_anova,x='sex',y='sPC12',add='jitter',color='sex') +
  theme(text=element_text(family = 'Arial',size=24),
        legend.position = 'none')+
  stat_compare_means(label.y = 7,size=8)
ggsave('../results/TransCompR_sparsePCA/sPC12_male_female.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)
