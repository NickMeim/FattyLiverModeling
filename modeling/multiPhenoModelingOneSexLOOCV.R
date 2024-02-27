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


#### Load pre-processed data----------------------
selected_sex <- 'male'
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
                # c('lobular_inflammation'='lobular_inflammation_grade'),
                sex,
                c('steatosis'='steatosis_grade')) %>% unique()
patients_meta_data <- patients_meta_data %>% filter(sex==selected_sex) %>% select(-sex)
pheno <- left_join(pheno,patients_meta_data)
pheno  <- pheno  %>% column_to_rownames('sample') #%>% select(-NAS)
pheno <- pheno %>% filter(!is.na(fibrosis))
Y_A <- as.matrix(pheno)
data_A <- data_A[,rownames(Y_A)]

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
# data.table::fwrite(chip_meta_data,'../preprocessing/chip_fatty_controls_indices.csv')
df_X_A = as.data.frame(X_A)
df_X_B = as.data.frame(X_B)
df_Y_A = as.data.frame(Y_A)
rownames(df_Y_A) <- rownames(X_A)

### Create list of indices for LOOCV
loocv_list <- NULL
all_inds <- rownames(Y_A)
for (i in 1:nrow(Y_A)){
  loocv_list[[i]] <- all_inds[which(all_inds!=rownames(Y_A)[i])]
}
### Re-do comparing linear model with lasso and shuffled features---------------------------
df_all <- data.frame()
results_folds <- list()
trainR2 <- matrix(nrow= nrow(Y_A),ncol=ncol(Y_A))
trainTau <- matrix(nrow= nrow(Y_A),ncol=ncol(Y_A))
trainMAE <- matrix(nrow= nrow(Y_A),ncol=ncol(Y_A))
all_laso_models <- list()
for (i in 1:nrow(Y_A)){
  train_inds <- loocv_list[[i]]
  train_XA <- X_A[which(rownames(X_A) %in% train_inds),]
  train_YA <- Y_A[which(rownames(Y_A) %in% train_inds),]
  val_XA <- X_A[which(!(rownames(X_A) %in% train_inds)),]
  val_YA <- Y_A[which(!(rownames(Y_A) %in% train_inds)),]
  res <- translation_model_multi_loocv(train_XA, X_B, train_YA)
  val_X_A_B <- val_XA %*% res$PC_B$rotation
  val_X_A_B <- val_X_A_B[,res$pcs_kept]
  predicted_val <- as.data.frame(predict(res$model_A_B,
                                         t(as.matrix(t(rbind(as.matrix(val_X_A_B),as.matrix(val_YA)))[,c(colnames(res$X_A_B),colnames(Y_A))]))))
  if (i==1){
    all_predicted_val <- predicted_val
  }else{
    all_predicted_val <- rbind(all_predicted_val,predicted_val)
  }
  trainR2[i,] <- cbind(res$train_r^2)
  trainTau[i,] <- cbind(res$train_tau)
  trainMAE[i,] <- cbind(res$tain_mae)
  results_folds[[i]] <- res
  
  saveRDS(results_folds,'../results/TransCompR_PCA_Males/all_laso_models_loocv.rds')
  
  print(paste0('Finished fold ',i))
}
colnames(trainMAE) <- colnames(train_YA)
colnames(trainTau) <- colnames(train_YA)
colnames(trainR2) <- colnames(train_YA)
df <- rbind(as.data.frame(trainMAE) %>% mutate(fold = seq(1,nrow(Y_A))) %>% mutate(set='train') %>% mutate(metric='MAE') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainTau)%>% mutate(fold = seq(1,nrow(Y_A))) %>% mutate(set='train') %>% mutate(metric='Tau') %>%
              gather('output','value',-metric,-set,-fold),
            as.data.frame(trainR2)%>% mutate(fold = seq(1,nrow(Y_A))) %>% mutate(set='train') %>% mutate(metric='R2') %>%
              gather('output','value',-metric,-set,-fold))
test_mae <- apply(abs(Y_A - all_predicted_val),2,mean)
test_tau <- cor(all_predicted_val,Y_A,method = 'kendall')
test_tau <- diag(test_tau)
test_r <- cor(all_predicted_val,Y_A,method = 'pearson')
test_r <- diag(test_r)
test_df <- cbind(test_mae,test_tau,test_r)
colnames(test_df) <- c('MAE','Tau','R2')
test_df <- as.data.frame(test_df) %>% rownames_to_column('output') %>% gather('metric','value',-output)
test_df <- test_df %>% mutate(set='test')
test_df <- test_df %>% mutate(fold = -1)
test_df <- test_df %>% select(all_of(colnames(df)))
df <- rbind(df,test_df)
df <- df %>% mutate(model = 'model')
rownames(all_predicted_val) <- rownames(Y_A)
predictions_df <- left_join(as.data.frame(Y_A) %>% rownames_to_column('sample') %>% gather('out','true',-sample),
                        as.data.frame(all_predicted_val) %>% rownames_to_column('sample') %>% gather('out','predicted',-sample))
df_all <- rbind(df_all,df)
ggscatter(predictions_df,x='true',y='predicted',cor.coef = T) + 
  facet_wrap(~out)
df <- df %>% group_by(set,metric,output) %>% mutate(mu=mean(value)) %>% mutate(se=sd(value)/sqrt(nrow(Y_A))) %>% ungroup()

(ggplot(df %>% filter(set=='train') %>% select(output,metric,mu,se) %>% unique(),
       aes(x=metric,y=mu,fill=metric)) +
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymax = mu + se,ymin=mu-se),size=1)+
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')+
  scale_y_continuous(breaks = c(-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25),limits = c(-0.9,1))+
  ylab('correlation')+
  theme_pubr(base_size = 24,base_family = 'Arial') +
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 24,family = 'Arial'),
        legend.title = element_blank()) +
  facet_wrap(~output)) +
  (ggplot(df %>% filter(set!='train') %>% select(output,metric,mu,se) %>% unique(),
          aes(x=metric,y=mu,fill=metric)) +
     geom_bar(stat = 'identity')+
     geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')+
     scale_y_continuous(breaks = c(-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25),limits = c(-0.9,1))+
     ylab('correlation')+
     theme_pubr(base_size = 24,base_family = 'Arial') +
     theme(text = element_text(size = 24,family = 'Arial'),
           axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           legend.text = element_text(size = 24,family = 'Arial'),
           legend.title = element_blank()) +
     facet_wrap(~output))
