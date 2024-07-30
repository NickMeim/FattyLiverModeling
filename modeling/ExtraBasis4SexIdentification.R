library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(caret)
library(ropls)
source("../utils/plotting_functions.R")
source("functions_translation.R")
source("CrossValidationUtilFunctions.R")
source('vector_space_interpretation.R')

### Load data--------------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Manually load also the other clinical datasets I have from ARCHS4
geo <- 'GSE162694' # only this one has multiple NAS and fibrosis scores
meta_data <- read.delim('../data/ARCHS4/FattyLiver_meta_data.tsv',row.names = 1)
meta_data <- meta_data %>% filter(series_id==geo)
old_cols <- colnames(meta_data)
meta_data <- meta_data %>% separate_rows(characteristics_ch1, sep = ",") %>%
  separate(characteristics_ch1, into = c("key", "value"), sep = ":") %>%
  mutate(key = str_trim(key), value = str_trim(value)) %>%
  spread(key, value)
new_cols <- colnames(meta_data)
new_cols <- new_cols[which(!(new_cols %in% old_cols))]
meta_data <- meta_data %>%
  mutate_at(vars(new_cols), ~ ifelse(grepl("\\d", .), as.numeric(.), .))
meta_data <- meta_data %>% filter(!(is.na(`nas score`) & is.na(`fibrosis stage`))) %>%
  mutate(`nas score`=ifelse(`nas score`=='NA',NA,`nas score`)) %>% 
  mutate(`fibrosis stage`=ifelse(`fibrosis stage`=='normal liver histology',0,`fibrosis stage`))
meta_data <- meta_data %>% filter(!is.na(`nas score`)) %>% filter(!is.na(`fibrosis stage`)) %>%
  filter(`nas score`!='NA') %>% filter(`fibrosis stage`!='NA')
expression_matrix <- readRDS('../data/ARCHS4/FattyLiver_expression_matrix.rds')
expression_matrix <- expression_matrix[,meta_data$sample]
data_list[['Pantano']] <- list(counts = expression_matrix,
                               metadata = meta_data,
                               genes = rownames(expression_matrix))
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Xh <- data_list[[ref_dataset]]$data_center %>% t()
sex_inferred <- apply(as.matrix(Xh[,c('RPS4Y1')]),2,sign) 
sex_inferred <- 1*(sex_inferred>0)
Yh <- as.factor(sex_inferred)
Yh <- as.matrix(Yh)
colnames(Yh) <- c('sex')
rownames(Yh) <- rownames(Xh)
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

## Compare MPS sex related genes expression with the clinical data------------------
cpm_fun <- function(x){return(log2(1 + x/sum(x)*1e6))}
sex_genes_expr <- rbind(data.frame(RPS4Y1=Xm[,c('RPS4Y1')]) %>% mutate(dataset = 'MPS'),
                        data.frame(RPS4Y1=Xh[,c('RPS4Y1')]) %>% mutate(dataset = ref_dataset))
ggviolin(sex_genes_expr,x='dataset',y='RPS4Y1',add = 'jitter')+
  theme(text = element_text(size=20,family = 'Arial'))
ggsave('../figures/sex_zscored_gene_comparisson_govaere_mps.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)
sex_genes_expr <- rbind(as.data.frame(t(data_list[[target_dataset]]$counts)[,c('RPS4Y1','XIST')]) %>% mutate(dataset = 'MPS'),
                        as.data.frame(t(data_list[[ref_dataset]]$counts)[,c('RPS4Y1','XIST')]) %>% mutate(dataset = ref_dataset)) %>%
  gather('gene','expression',-dataset)
ggboxplot(sex_genes_expr %>% filter(dataset=='MPS'),x='gene',y='expression',add = 'jitter')+
  scale_y_log10()+
  ylab('raw counts')+
  theme(text = element_text(size=20,family = 'Arial'))
ggsave('../figures/sex_gene_comparisson_mps_alone.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)

sex_genes_expr <- rbind(data.frame(RPS4Y1=Xm[,c('RPS4Y1')]) %>% mutate(dataset = 'MPS') %>% mutate(sex='unknown'),
                        as.data.frame(data_list[['Hoang']]$data_center %>% t()) %>% select(RPS4Y1) %>% 
                          mutate(dataset = 'Hoang') %>% mutate(sex = data_list[['Hoang']]$metadata$sex),
                        as.data.frame(data_list[['Pantano']]$data_center %>% t()) %>% select(RPS4Y1) %>% 
                          mutate(dataset = 'Pantano') %>% mutate(sex = data_list[['Pantano']]$metadata$Sex))
ggviolin(sex_genes_expr,x='dataset',y='RPS4Y1',color='sex',add = 'jitter')+
  scale_y_continuous(n.breaks = 10,limits = c(-4,NA))+
  theme(text = element_text(size=20,family = 'Arial'),
        panel.grid.major.y = element_line())
ggsave('../figures/sex_zscores_gene_comparisson_all.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)
# sex_genes_expr <- rbind(as.data.frame(t(data_list[[target_dataset]]$counts)[,c('RPS4Y1','XIST')]) %>% mutate(dataset = 'MPS') %>% mutate(sex='unknown'),
#                         as.data.frame(t(data_list[['Hoang']]$counts)) %>% select(RPS4Y1,XIST) %>% 
#                           mutate(dataset = 'Hoang') %>% mutate(sex = data_list[['Hoang']]$metadata$sex),
#                         as.data.frame(t(data_list[['Pantano']]$counts)) %>% select(RPS4Y1,XIST) %>% 
#                           mutate(dataset = 'Pantano') %>% mutate(sex = data_list[['Pantano']]$metadata$Sex)) %>%
#   gather('gene','expression',-dataset,-sex) %>% mutate(sex = tolower(sex))
# # sex_genes_expr <- sex_genes_expr %>% group_by(dataset,gene) %>% mutate(CPM=cpm_fun(expression))
# ggboxplot(sex_genes_expr,x='dataset',y='expression',color='sex',add = 'jitter')+
#   scale_y_continuous(n.breaks = 8,limits = c(0,NA))+
#   facet_wrap(~gene,scales = 'free_y')+
#   theme(text = element_text(size=20,family = 'Arial'),
#         panel.grid.major.y = element_line())
# ggsave('../figures/sex_gene_comparisson_all.png',
#        width = 9,
#        height = 9,
#        units = 'in',
#        dpi=600)

#### Train and evaluate PLS-DA model----------------
Wm_tot <- readRDS(paste0('../results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('../results/Wm_',tolower(target_dataset),'_extra.rds'))
Wm_combo <- readRDS(paste0('../results/Wm_',tolower(target_dataset),'_combo.rds'))
num_folds <- 10
loc <- '../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/'
num_LVS <- 2
train_acc<- NULL
train_f1<- NULL
val_f1<- NULL
val_acc<- NULL
val_f1_shuffle_x <- NULL
val_acc_shuffle_x <- NULL
val_f1_shuffle_y <- NULL
val_acc_shuffle_y <- NULL
# val_f1_random_x <- NULL
# val_acc_random_x <- NULL
train_acc_backproj<- NULL
val_f1_backproj<- NULL
val_f1_backproj<- NULL
val_acc_backproj<- NULL
train_acc_backproj_opt<- NULL
train_f1_backproj_opt<- NULL
val_f1_backproj_opt<- NULL
val_acc_backproj_opt<- NULL
val_acc_external <- matrix(0,nrow = 10,ncol = 2)
val_f1_external <- matrix(0,nrow = 10,ncol = 2)
val_acc_external_shuffle_y <- matrix(0,nrow = 10,ncol = 2)
val_f1_external_shuffle_y <- matrix(0,nrow = 10,ncol = 2)
val_acc_external_shuffle_x <- matrix(0,nrow = 10,ncol = 2)
val_f1_external_shuffle_x <- matrix(0,nrow = 10,ncol = 2)
external_clinical_datasets <- c("Hoang","Pantano")
for (j in 1:num_folds){
  message(paste0('Begun fold ',j))
  x_train <- readRDS(paste0(loc,'Xh_train',j,'.rds'))
  y_train <- Yh[rownames(x_train),]
  y_train <- factor(y_train,levels = c('0','1'))
  # y_train <- as.matrix(y_train)
  x_val <- readRDS(paste0(loc,'Xh_val',j,'.rds'))
  y_val<- Yh[rownames(x_val),]
  y_val <- factor(y_val,levels = c('0','1'))
  # y_val <- as.matrix(y_val)
  
  plsr_model <- opls(x = x_train, 
                     y = y_train,
                     predI = num_LVS,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  y_train_hat <- predict(plsr_model,x_train)
  cofusion_train <- caret::confusionMatrix(data = y_train_hat,reference = y_train,positive = '1')
  y_val_hat <- predict(plsr_model,x_val)
  cofusion_val <- caret::confusionMatrix(data = y_val_hat,reference = y_val,positive = '1')
  train_acc[j] <- as.numeric(cofusion_train$overall['Accuracy'])
  train_f1[j] <- as.numeric(cofusion_train$byClass['F1'])
  val_f1[j] <- as.numeric(cofusion_val$byClass['F1'])
  val_acc[j] <- as.numeric(cofusion_val$overall['Accuracy'])
  
  
  ### shuffled labels model
  y_train_shuffled <- y_train[sample.int(length(y_train))]
  names(y_train_shuffled) <- names(y_train)
  plsr_model_shuffle_y <- opls(x = x_train, 
                               y = y_train_shuffled,
                               predI = num_LVS,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  y_hat_val <- predict(plsr_model_shuffle_y,x_val)
  cofusion_val <- caret::confusionMatrix(data = y_hat_val,reference = y_val,positive = '1')
  val_f1_shuffle_y[j] <- as.numeric(cofusion_val$byClass['F1'])
  val_acc_shuffle_y[j] <- as.numeric(cofusion_val$overall['Accuracy'])
  
  ### shuffled features model
  x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
  colnames(x_train_shuffled) <- colnames(x_train)
  plsr_model_shuffle_x <- opls(x = x_train_shuffled, 
                               y = y_train,
                               predI = num_LVS,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  y_hat_val <- predict(plsr_model_shuffle_x,x_val)
  cofusion_val <- caret::confusionMatrix(data = y_hat_val,reference = y_val,positive = '1')
  val_f1_shuffle_x[j] <- as.numeric(cofusion_val$byClass['F1'])
  val_acc_shuffle_x[j] <- as.numeric(cofusion_val$overall['Accuracy'])
  
  # ### random features model
  # x_train_random <- matrix(rnorm(n = nrow(x_train)*ncol(x_train)),nrow = nrow(x_train))
  # colnames(x_train_random) <- colnames(x_train)
  # rownames(x_train_random) <- rownames(x_train)
  # plsr_model_random_x <- opls(x = x_train_random, 
  #                             y = y_train,
  #                             predI = num_LVS,
  #                             crossvalI = 1,
  #                             scaleC = "center",
  #                             fig.pdfC = "none",
  #                             info.txtC = "none")
  # y_hat_test <- predict(plsr_model_random_x,x_val)
  # cofusion_val <- caret::confusionMatrix(data = y_hat_test,reference = y_val,positive = '1')
  # val_f1_random_x[j] <- as.numeric(cofusion_val$byClass['F1'])
  # val_acc_random_x[j] <- as.numeric(cofusion_val$overall['Accuracy'])
  
  # Get Wh of PLSR
  Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
  rownames(Wh) <- colnames(x_train)
  colnames(Wh) <- colnames(plsr_model@weightMN)
  for (ii in 1:nrow(plsr_model@weightMN)){
    Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
  }
  # Get regression coefficients
  # Bh <- matrix(lm(y_train ~ x_train %*% Wh) %>% coef(), ncol = 2)
  Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
  
  print('Finished running initial PLSR for humans')
  
  y_hat_test <- predict(plsr_model,x_val %*% Wm_tot %*% t(Wm_tot))
  y_hat_train <- predict(plsr_model,x_train %*% Wm_tot %*% t(Wm_tot))
  cofusion_train <- caret::confusionMatrix(data = y_hat_train,reference = y_train,positive = '1')
  cofusion_val <- caret::confusionMatrix(data = y_hat_test,reference = y_val,positive = '1')
  train_acc_backproj[j] <- as.numeric(cofusion_train$overall['Accuracy'])
  train_f1_backproj[j] <- as.numeric(cofusion_train$byClass['F1'])
  val_f1_backproj[j] <- as.numeric(cofusion_val$byClass['F1'])
  val_acc_backproj[j] <- as.numeric(cofusion_val$overall['Accuracy'])
  print('Finished back-projection with NASH augmented in vitro space')
  
  ### find 3rd extra space and back-project
  phi <- Wh %*% Bh
  Wm_opt_sex <- analytical_solution_opt(y=as.matrix(y_train),
                                    W_invitro = Wm_tot,
                                    phi = phi)
  
  W_tot_2 <- cbind(Wm_tot,Wm_opt_sex)
  y_hat_test <- predict(plsr_model,x_val %*% W_tot_2 %*% t(W_tot_2))
  y_hat_train <- predict(plsr_model,x_train %*% W_tot_2 %*% t(W_tot_2))
  cofusion_train <- caret::confusionMatrix(data = y_hat_train,reference = y_train,positive = '1')
  cofusion_val <- caret::confusionMatrix(data = y_hat_test,reference = y_val,positive = '1')
  train_acc_backproj_opt[j] <- as.numeric(cofusion_train$overall['Accuracy'])
  train_f1_backproj_opt[j] <- as.numeric(cofusion_train$byClass['F1'])
  val_f1_backproj_opt[j] <- as.numeric(cofusion_val$byClass['F1'])
  val_acc_backproj_opt[j] <- as.numeric(cofusion_val$overall['Accuracy'])
  print('Finished back-projection with sex-extra in vitro space')
  
  k <- 1
  for (dataset in external_clinical_datasets){
    print(paste0('Begun ',dataset ,' external dataset'))
    X <- data_list[[dataset]]$data_center %>% t()
    if (dataset == 'Hoang'){
      Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(sex) %>% mutate(sex=ifelse(sex=='male',1,0)))
      
    }else if (dataset == 'Pantano'){
        Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(c('sex'='Sex')) %>% mutate(sex=ifelse(sex=='Male',1,0)))
    }
    Y <- factor(Y,levels = c(0,1))
    
    ## PLSR original
    y_val_hat <- predict(plsr_model,X)
    cofusion_val <- caret::confusionMatrix(data = y_val_hat,reference = Y,positive = '1')
    val_acc_external[j,k] <- as.numeric(cofusion_val$overall['Accuracy'])
    val_f1_external[j,k] <- as.numeric(cofusion_val$byClass['F1'])
    
    ## Shuffled labels
    y_hat_val <- predict(plsr_model_shuffle_y,X)
    cofusion_val <- caret::confusionMatrix(data = y_hat_val,reference = Y,positive = '1')
    val_acc_external_shuffle_y[j,k] <- as.numeric(cofusion_val$overall['Accuracy'])
    val_f1_external_shuffle_y[j,k] <- as.numeric(cofusion_val$byClass['F1'])
    
    ## Shuffled features
    y_hat_val <- predict(plsr_model_shuffle_x,X)
    cofusion_val <- caret::confusionMatrix(data = y_hat_val,reference = Y,positive = '1')
    val_acc_external_shuffle_x[j,k] <- as.numeric(cofusion_val$overall['Accuracy'])
    val_f1_external_shuffle_x[j,k] <- as.numeric(cofusion_val$byClass['F1'])
    k <- k +1
  }
}
res_human_genes <- rbind(data.frame(Accuracy = train_acc,F1=train_f1) %>% mutate(set = 'train'),
                         data.frame(Accuracy = val_acc,F1=val_f1) %>% mutate(set = 'test'),
                         data.frame(Accuracy = val_acc_shuffle_y,F1=val_f1_shuffle_y) %>% mutate(set = 'shuffle Y'),
                         data.frame(Accuracy = val_acc_shuffle_x,F1=val_f1_shuffle_x) %>% mutate(set = 'shuffle X')) %>%
  mutate(model = 'human genes')
# data.frame(Accuracy = val_acc_random_x,F1=val_f1_random_x) %>% mutate(set = 'random X')
saveRDS(res_human_genes,'../results/sex_analysis/plsda_sex_only_performance.rds')
res_train_all <- rbind(res_human_genes %>% filter(set =='train'),
                       data.frame(Accuracy = train_acc_backproj,F1=train_f1_backproj)%>% mutate(set = 'train') %>% 
                         mutate(model = 'back-projected'),
                       data.frame(Accuracy = train_acc_backproj_opt,F1=train_f1_backproj_opt) %>% mutate(set = 'train') %>% 
                         mutate(model = 'sex-optimized MPS'))
saveRDS(res_train_all,'../results/sex_analysis/plsda_sex_only_train_all_performance.rds')
res_test_all <- rbind(res_human_genes %>% filter(set =='test'),
                       data.frame(Accuracy = val_acc_backproj,F1=val_f1_backproj)%>% mutate(set = 'test') %>% 
                         mutate(model = 'back-projected'),
                       data.frame(Accuracy = val_acc_backproj_opt,F1=val_f1_backproj_opt) %>% mutate(set = 'test') %>% 
                         mutate(model = 'sex-optimized MPS'),
                      res_human_genes %>% filter(set =='shuffle X'))
saveRDS(res_test_all,'../results/sex_analysis/plsda_sex_only_test_all_performance.rds')

colnames(val_acc_external) <- external_clinical_datasets
colnames(val_acc_external_shuffle_x) <- external_clinical_datasets
colnames(val_acc_external_shuffle_y) <- external_clinical_datasets
colnames(val_f1_external) <- external_clinical_datasets
colnames(val_f1_external_shuffle_x) <- external_clinical_datasets
colnames(val_f1_external_shuffle_y) <- external_clinical_datasets
res_external <-  rbind(left_join(as.data.frame(val_acc_external) %>% mutate(fold=seq(1,num_folds)) %>%
                                   gather('dataset','Accuracy',-fold),
                             as.data.frame(val_f1_external)%>% mutate(fold=seq(1,num_folds)) %>%
                               gather('dataset','F1',-fold)) %>% select(-fold) %>% mutate(model= 'PLS-DA'),
                       left_join(as.data.frame(val_acc_external_shuffle_y) %>% mutate(fold=seq(1,num_folds)) %>%
                                   gather('dataset','Accuracy',-fold),
                                 as.data.frame(val_f1_external_shuffle_y)%>% mutate(fold=seq(1,num_folds)) %>%
                                   gather('dataset','F1',-fold)) %>% select(-fold)%>% mutate(model= 'shuffle Y'),
                       left_join(as.data.frame(val_acc_external_shuffle_x) %>% mutate(fold=seq(1,num_folds)) %>%
                                   gather('dataset','Accuracy',-fold),
                                 as.data.frame(val_f1_external_shuffle_x)%>% mutate(fold=seq(1,num_folds)) %>%
                                   gather('dataset','F1',-fold)) %>% select(-fold)%>% mutate(model= 'shuffle X'))
saveRDS(res_external,'../results/sex_analysis/plsda_sex_only_external_clinical_performance.rds')

min_val <- min(c(0.5,res_human_genes$F1,res_human_genes$Accuracy))
ggboxplot(res_human_genes %>% gather('metric','value',-set,-model),x='set',y='value',color='set',add='jitter')+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(min_val,1.1))+
  stat_compare_means(comparisons = list(c('test','shuffle Y'),
                                        c('test','shuffle X')),
                     step.increase = 0.05,
                     size=6)+
  facet_wrap(~metric)+
  theme(text = element_text(size=20,family = 'Arial'),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(),
        legend.position = 'none')
ggsave('../results/sex_analysis/plsda_sex_human_genes.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../results/sex_analysis/plsda_sex_human_genes.eps',
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
min_val <- min(c(0.5,res_external$F1,res_external$Accuracy))
ggboxplot(res_external %>% gather('metric','value',-model,-dataset),x='model',y='value',color='model',add='jitter')+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(min_val,1.1))+
  stat_compare_means(comparisons = list(c('PLS-DA','shuffle Y'),
                                        c('PLS-DA','shuffle X')),
                     step.increase = 0.05,
                     size=6)+
  facet_wrap(vars(dataset,metric))+
  theme(text = element_text(size=20,family = 'Arial'),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(),
        legend.position = 'none')
ggsave('../results/sex_analysis/plsda_sex_human_genes_external_clinical.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../results/sex_analysis/plsda_sex_human_genes_clinical.eps',
       device = cairo_ps,
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)

min_val <- min(c(0.5,res_train_all$F1,res_train_all$Accuracy))
ggboxplot(res_train_all %>% gather('metric','value',-set,-model),x='model',y='value',color='model',add='jitter')+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(min_val,1.1))+
  stat_compare_means(comparisons = list(c('human genes','back-projected'),
                                        c('sex-optimized MPS','back-projected'),
                                        c('human genes','sex-optimized MPS')),
                     label.y = c(1,1,1.05),
                     size=6)+
  facet_wrap(~metric,nrow = 2)+
  theme(text = element_text(size=20,family = 'Arial'),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(),
        legend.position = 'none')
ggsave('../results/sex_analysis/plsda_sex_comparison_10foldtrain.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../results/sex_analysis/plsda_sex_comparison_10foldtrain.eps',
       device = cairo_ps,
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)

min_val <- min(c(0.5,res_test_all$F1,res_test_all$Accuracy))
ggboxplot(res_test_all %>% mutate(model = ifelse(set %in% c('shuffle X','shuffle Y','random X'),set,model)) %>% 
            gather('metric','value',-set,-model),x='model',y='value',color='model',add='jitter')+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(min_val,1.1))+
  stat_compare_means(comparisons = list(c('human genes','back-projected'),
                                        c('sex-optimized MPS','back-projected'),
                                        c('human genes','sex-optimized MPS'),
                                        c('sex-optimized MPS','shuffle X')),
                     label.y = c(1,1,1.05,1),
                     size=6)+
  facet_wrap(~metric,nrow = 2)+
  theme(text = element_text(size=20,family = 'Arial'),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(),
        legend.position = 'none')
ggsave('../results/sex_analysis/plsda_sex_comparison_10foldtest.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../results/sex_analysis/plsda_sex_comparison_10foldtest.eps',
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
### Since we perform well, run PLS-DA with all the data and explore the extra LV-------------------
plsr_model <- opls(x = Xh, 
                   y = Yh,
                   predI = num_LVS,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
Zh_plsr <- as.data.frame(plsr_model@scoreMN)
Zh_plsr$sex <- Yh[,1]
Zh_plsr <- Zh_plsr %>% mutate(sex=ifelse(sex=='1','male','female'))
plt_scores <- ggplot(Zh_plsr, aes(p1, p2, fill = sex)) +
  geom_vline(xintercept = 0, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(color = "black",
                      size = 2.8,
                      stroke = 1.2,
                      shape = 21,
                      show.legend = TRUE) +
  labs(x = paste("LV1 (",
                          toString(round(plsr_model@modelDF$R2X[1] * 100,2)), "%)", sep = ""),
                y = paste("LV2 (",
                          toString(round(plsr_model@modelDF$R2X[2] * 100,2)), "%)", sep = "")) +
  theme_classic(base_size = 20,base_family = 'Arial') +
  theme(text = element_text(size = 20,family = 'Arial'),
        legend.position = "right",
                 aspect.ratio = 1,
                 axis.text = element_text(color = "black"))
plt_scores <- plt_scores +
  stat_ellipse(aes(color = sex),size=1.2) +
  scale_color_viridis_d()+
  scale_fill_viridis_d()
print(plt_scores)
ggsave('../results/sex_analysis/plsda_sex_visualize_scores.png',
       width = 6,
       height = 6,
       units = 'in',
       dpi = 600)
ggsave('../results/sex_analysis/plsda_sex_visualize_scores.eps',
       device = cairo_ps,
       width = 6,
       height = 6,
       units = 'in',
       dpi = 600)
# Get Wh of PLSR
Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
rownames(Wh) <- colnames(x_train)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)){
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
}
# Get regression coefficients
# Bh <- matrix(lm(y_train ~ x_train %*% Wh) %>% coef(), ncol = 2)
Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
### find 3rd extra space and back-project
phi <- Wh %*% Bh
Wm_opt_sex <- analytical_solution_opt(y=as.matrix(Yh),
                                      W_invitro = Wm_tot,
                                      phi = phi)
rownames(Wm_opt_sex) <- colnames(Xh)
Zh <- as.data.frame(Xh %*% Wm_opt_sex)
Zh$sex <- Yh[,1]
Zh <- Zh %>% mutate(sex=ifelse(sex=='1','male','female'))
Zh$sex <- factor(Zh$sex,levels = c('female','male'))
ggboxplot(Zh,x='sex',y='V1',color='sex',add='jitter')+
  scale_color_viridis_d()+
  ylab('sex-specific extra basis score')+
  stat_compare_means(comparisons = list(c('male','female')),
                     size=6)+
  theme(text = element_text(size=20,family = 'Arial'),
        # axis.title.x = element_blank(),
        panel.grid.major.y = element_line(),
        legend.position = 'none')
ggsave('../results/sex_analysis/plsda_sex_extra_basis_seperation.png',
       width = 6,
       height = 6,
       units = 'in',
       dpi = 600)
ggsave('../results/sex_analysis/plsda_sex_extra_basis_seperation.eps',
       device = cairo_ps,
       width = 6,
       height = 6,
       units = 'in',
       dpi = 600)

plot_extra_gene_loadings_lv1 <- plot_gene_loadings(loadings = Wm_opt_sex,
                                                   selection='V1',
                                                   y_lab = 'weight in sex-specific extra LV1',
                                                   top=20)
ggsave(paste0('../results/sex_analysis/gene_sec_woptimal_LV1_',
              tolower(target_dataset),
              '_loadings.png'),
       plot = plot_extra_gene_loadings_lv1,
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

### Analyze TFs------------------
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
extra_basis_TF_activity <- TF_activity_interpretation(Wm_opt_sex,
                                                      Wm_tot,
                                                      dorotheaData)
p <- extra_basis_TF_activity$figure[[1]] +
  geom_point(size=2,color = '#CD5C5C')  +
  xlab('Rank')+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        panel.grid.major =  element_line(),
        plot.title = element_blank(),
        legend.position = 'none')
ggsave(paste0('../results/sex_analysis/tfs_sex_plsda_',
              tolower(ref_dataset),
              '_barplot.png'),
       plot = p,
       width = 9,
       height = 9,
       dpi = 600)

### Analyze pathways------------------
extra_basis_pathway_activity <- pathway_activity_interpretation(Wm_opt_sex,
                                                                Wm_tot)
ggsave(paste0('../results/sex_analysis/progenies_sex_loadings_',
              tolower(ref_dataset),
              '_barplot.png'),
       plot = extra_basis_pathway_activity$figure[[1]] + theme(plot.title = element_blank()),
       width = 9,
       height = 9,
       dpi = 600)
setEPS()
postscript(paste0('../results/sex_analysis/progenies_sex_loadings_',
                  tolower(ref_dataset),
                  '_barplot.eps'),
           height = 9, width = 9)
print(extra_basis_pathway_activity$figure[[1]]+ theme(plot.title = element_blank()))
dev.off()
### Analyze hallmarks------------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Wm_opt_sex), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Wm_opt_sex[inds,])
rownames(meas) <- entrez_ids
msig <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes <- as.data.frame(msig$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark') 
colnames(msig_nes) <- c('Hallmark','NES')
msig_pval <- as.data.frame(msig$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')
colnames(msig_pval) <- c('Hallmark','padj')
df_msig <- left_join(msig_nes,msig_pval)
df_msig <- df_msig %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))
df_msig <- df_msig %>% mutate(Hallmark=str_replace_all(Hallmark,'_',' '))
df_msig <- df_msig %>% mutate(Hallmark = tolower(Hallmark)) %>% 
  mutate(Hallmark = paste0(toupper(substr(Hallmark, 1, 1)), tolower(substr(Hallmark, 2, nchar(Hallmark)))))
p1 <- (ggplot(df_msig %>% arrange(NES) %>%
                filter(padj<=0.1),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV extra 1')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_blank(),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p1)
ggsave(paste0('../results/sex_analysis/hallmark_sex_extra_lv_',
              tolower(target_dataset),
              '.png'),
       plot=p1,
       width=9,
       height=9,
       units = 'in',
       dpi = 600)
ggsave(paste0('../results/sex_analysis/hallmark_sex_extra_lv_',
              tolower(ref_dataset),
              '.eps'),
       device = cairo_ps,
       plot=p1,
       width=9,
       height=9,
       units = 'in',
       dpi = 600)
