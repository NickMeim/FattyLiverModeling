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

### Load in-vitro and in-vivo datasets-----------------------------------------
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
Yh <- as.matrix(data_list$Govaere$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# Load previously found extra basis
Wm_opt <- readRDS(paste0('../results/Wm_',target_dataset,'_extra.rds'))
Wm_tot <- readRDS(paste0('../results/Wm_',target_dataset,'_total.rds'))
  
### Load found extra basis from the Govaere and Kostrzewski datasets and check------------------------------
### if performance can be retrieved when backprojecting for other human datasets
external_clinical_datasets <- c("Hoang","Pantano")
df_results <- data.frame()
predictions_results <- data.frame()
for (dataset in external_clinical_datasets){
  message(paste0('Begun ',dataset ,' dataset'))
  X <- data_list[[dataset]]$data_center %>% t()
  if (dataset == 'Hoang'){
    Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(nafld_activity_score,Fibrosis_stage))

  }else if (dataset == 'Pantano'){
    Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(`nas score`,`fibrosis stage`))
    Y <- apply(Y,c(1,2),as.numeric)
  }
  colnames(Y) <- c('NAS','fibrosis')
  ### Run a small 5-fold cross validation procedure for 
  data_splits <- createMultiFolds(y = Y[,2], k = 10, times = 10)
  train_r <- NULL
  val_r <- NULL
  train_r_backproj <- NULL
  val_r_backproj <- NULL
  train_r_extra <- NULL
  val_r_extra <- NULL
  j <- 1
  for (ind in data_splits){
    name_fold_id <- names(data_splits)[j]
    y_val <- Y[-ind,]
    y_train <- Y[ind,]
    x_val <- X[-ind,]
    x_train <- X[ind,]
    plsr_model <- opls(x = x_train,
                       y = y_train,
                       predI = 8,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    y_train_hat <- predict(plsr_model,x_train)
    y_val_hat <- predict(plsr_model,x_val)
    train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
    
    # Get Wh of PLSR
    Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
    rownames(Wh) <- colnames(x_train)
    colnames(Wh) <- colnames(plsr_model@weightMN)
    for (ii in 1:nrow(plsr_model@weightMN)){
      Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    }
    # Get regression coefficients
    Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
    # Define projection matrices to make more readable
    Th_train <- x_train %*% Wh
    Thm_train <- x_train %*% Wm %*% t(Wm) %*% Wh
    Th_val<- x_val %*% Wh
    Thm_val <- x_val %*% Wm %*% t(Wm) %*% Wh
    
    y_train_hat <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    y_val_hat <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    
    train_r_extra[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r_extra[j] <- mean(diag(cor(y_val_hat,y_val)))
    
    ### simple backprojection 
    y_val_hat <- cbind(1, Thm_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    y_train_hat <- cbind(1, Thm_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    train_r_backproj[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r_backproj[j] <- mean(diag(cor(y_val_hat,y_val)))
    
    print(paste0('Finished ',j,' fold out of ',length(data_splits)))
    j <- j+1
  }
  tmp_train <- data.frame(PLSR = train_r,PC = train_r_backproj,'extra basis'=train_r_extra,fold_ids = names(data_splits)) %>% 
    mutate(dataset = dataset) %>% mutate(set = 'train')
  tmp_val <- data.frame(PLSR = val_r,PC =val_r_backproj,'extra basis'=val_r_extra,fold_ids = names(data_splits)) %>% 
    mutate(dataset = dataset) %>% mutate(set = 'test')
  df_results <- rbind(df_results,tmp_train,tmp_val)
}
# saveRDS(predictions_results,'../results/external_clinical_differences_results_of_extra_vector.rds')
# #df_results <- readRDS('../results/external_clinical_performance_of_extra_vector.rds')
df_results <- df_results %>% gather('model','r',-dataset,-set,-fold_ids)
# df_results <- df_results %>% mutate(model = ifelse(model=='extra.basis','extra basis',
#                                                    ifelse(model!='PLSR','in-vitro PCs','original PLSR')))
df_results <- df_results %>% mutate(model = ifelse(model=='extra.basis','optimized MPS',
                                                   ifelse(model!='PLSR','back-projected','human genes')))
# df_results$model <- factor(df_results$model,levels = c('original PLSR','extra basis','in-vitro PCs'))
df_results$model <- factor(df_results$model,levels = c('human genes','optimized MPS','back-projected'))
df_results <- df_results %>% separate('fold_ids',into = c('fold','repetition')) %>%
  group_by(dataset,model,set,repetition) %>% mutate(r_mu = mean(r))  %>% ungroup()%>%
  select(-repetition,-fold) %>% unique()

### Visualize 
# model_comparisons <- list(
#   c('original PLSR', 'extra basis'),
#   c('extra basis', 'in-vitro PCs')
# )
model_comparisons <- list(
  c('human genes', 'optimized MPS'),
  c('optimized MPS', 'back-projected')
)
df_results$set <- factor(df_results$set,levels=c('train','test')) 
df_results  <- left_join(df_results,
                 df_results %>% select(-r) %>% unique() %>% 
                   group_by(dataset,model,set) %>% mutate(total_se = sd(r_mu)) %>%
                   ungroup())
df_results <- df_results %>% group_by(dataset,model,set) %>%  mutate(all_mu = mean(r_mu)) %>% ungroup()
# p1 <- ggviolin(df_results %>% filter(dataset=='Hoang'),x='model',y='r',fill = 'model') +
#   geom_boxplot(data=df_results %>% select(dataset,model,set,r_mu) %>% unique() %>% filter(dataset=='Hoang'),
#                aes(x=model,y=r_mu,fill=model),
#                width = 0.3,
#                size=0.5,outlier.shape = NA)+
#   scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,NA))+
#   xlab('')+ ylab('pearson`s correlation')+
#   ggtitle('Hoang')+
#   theme(text = element_text(size = 20,family = 'Arial'),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_blank(),
#         panel.grid.major.y = element_line(linewidth = 1),
#         panel.grid.minor.y = element_line(linewidth = 0.5),
#         legend.position = 'bottom')+
#   facet_wrap(~set)+
#   stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
#                                         c('optimized MPS','back-projected'),
#                                         c('human genes','back-projected')),
#                      label.y = c(0.98,0.98,1.03),
#                      method = 'wilcox')
p1 <- ggboxplot(df_results %>% filter(dataset=='Hoang') %>% select(dataset,set,model,r_mu) %>% unique(),
                x='model',y='r_mu',color = 'model',add='jitter') +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,NA))+
  xlab('')+ ylab('pearson`s correlation')+
  ggtitle('Hoang')+
  theme(text = element_text(size = 20,family = 'Arial'),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(), 
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'bottom')+
  facet_wrap(~set)+
  stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
                                        c('optimized MPS','back-projected'),
                                        c('human genes','back-projected')),
                     label.y = c(0.98,0.98,1.03),
                     method = 'wilcox')
# print(p1)
p2 <- ggboxplot(df_results %>% filter(dataset=='Hoang') %>% select(dataset,set,model,r_mu) %>% unique(),
                x='model',y='r_mu',color = 'model',add='jitter') +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,NA))+
  xlab('')+ ylab('pearson`s correlation')+
  ggtitle('Pantano')+
  theme(text = element_text(size = 20,family = 'Arial'),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(), 
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'bottom')+
  facet_wrap(~set)+
  stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
                                        c('optimized MPS','back-projected'),
                                        c('human genes','back-projected')),
                     label.y = c(0.98,0.98,1.03),
                     method = 'wilcox')
# print(p2)
p <- (p1 + p2) + plot_layout(guides = "collect",axes = "collect") & theme(legend.position = 'bottom')
print(p)

ggsave('../figures/inter_clinical_dataset_performance.png',
       plot=p,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)
ggsave('../figures/inter_clinical_dataset_performance.eps',
       plot=p,
       device = cairo_ps,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)

### Check generalization of the PLSR model itself-----------------------------
external_clinical_datasets <- c("Hoang","Pantano")
df_results <- data.frame()
num_folds <- 10
num_LVs <- 8
val_r_shuffle_y <- matrix(0,nrow = 10,ncol = 2)
val_r_shuffle_x <- matrix(0,nrow = 10,ncol = 2)
val_r <- matrix(0,nrow = 10,ncol = 2)
predictions_results <- data.frame()
for (i in 1:num_folds){
  x_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/Xh_train',i,'.rds'))
  y_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/Yh_train',i,'.rds'))
  
  plsr_model <- opls(x = x_train,
                     y = y_train,
                     predI = num_LVs,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  
  ### shuffled labels model
  y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
  rownames(y_train_shuffled) <- rownames(y_train)
  plsr_model_shuffle_y <- opls(x = x_train, 
                               y = y_train_shuffled,
                               predI = num_LVs,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  
  ### shuffled features model
  x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
  colnames(x_train_shuffled) <- colnames(x_train)
  plsr_model_shuffle_x <- opls(x = x_train_shuffled, 
                               y = y_train,
                               predI = num_LVs,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  
  j <- 1
  for (dataset in external_clinical_datasets){
    message(paste0('Begun ',dataset ,' dataset'))
    X <- data_list[[dataset]]$data_center %>% t()
    if (dataset == 'Hoang'){
      Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(nafld_activity_score,Fibrosis_stage))
      
    }else if (dataset == 'Pantano'){
      Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(`nas score`,`fibrosis stage`))
      Y <- apply(Y,c(1,2),as.numeric)
    }
    colnames(Y) <- c('NAS','fibrosis')
    
    ## PLSR original
    y_val_hat <- predict(plsr_model,X)
    val_r[i,j] <- mean(diag(cor(y_val_hat,Y)))
    tmp_plsr <- rbind(data.frame(score= Y[,'NAS'],
                                 prediction = y_val_hat[,'NAS']) %>% mutate(phenotype='NAS'),
                      data.frame(score= Y[,'fibrosis'],
                                 prediction = y_val_hat[,'fibrosis']) %>% mutate(phenotype='fibrosis')) %>%
      mutate(fold = i) %>% mutate(dataset = dataset)
    predictions_results <- rbind(predictions_results,tmp_plsr)
    
    ## Shuffled labels
    y_hat_val <- predict(plsr_model_shuffle_y,X)
    val_r_shuffle_y[i,j]<- mean(diag(cor(y_hat_val,Y)))
    
    ## Shuffled features
    y_hat_val <- predict(plsr_model_shuffle_x,X)
    val_r_shuffle_x[i,j]<- mean(diag(cor(y_hat_val,Y)))
    
    j <- j+1
  }
  
  print(paste0('Finished ',i,' fold out of ',num_folds))
}
colnames(val_r) <- external_clinical_datasets
colnames(val_r_shuffle_x) <- external_clinical_datasets
colnames(val_r_shuffle_y) <- external_clinical_datasets
df_results <- rbind(as.data.frame(val_r) %>% mutate(model = 'PLSR'),
                    as.data.frame(val_r_shuffle_x) %>% mutate(model = 'shuffled features'),
                    as.data.frame(val_r_shuffle_y) %>% mutate(model = 'shuffled labels'))
df_results <- df_results %>% gather('dataset','r',-model)
df_results$model <- factor(df_results$model,levels = c('PLSR','shuffled labels','shuffled features'))

model_comparisons <- list(
  c('PLSR', 'shuffled labels'),
  c('PLSR', 'shuffled features')
)
p <- ggboxplot(df_results,x='model',y='r',color = 'model',add='jitter') +
  scale_y_continuous(breaks = seq(-0.2,1,0.1),limits = c(-0.25,1))+
  xlab('')+ ylab('pearson`s correlation')+
  # ggtitle('PLSR generalization in other in-vivo human datasets')+
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'none')+
  facet_wrap(~dataset)+
  stat_compare_means(comparisons = model_comparisons,
                     method = 'wilcox',
                     size=9)
print(p)
ggsave('../figures/inter_clinical_PLSR_performance.png',
       plot=p,
       width = 16,
       height = 10,
       units = 'in',
       dpi = 600)
ggsave('../figures/inter_clinical_PLSR_performance.eps',
       plot=p,
       device = cairo_ps,
       width = 16,
       height = 10,
       units = 'in',
       dpi = 600)

### Visualize predictions
predictions_results <- predictions_results %>% mutate(prediction = ifelse(prediction>8,8,
                                                                          ifelse(prediction<0,0,prediction)))
predictions_results <- predictions_results %>% mutate(prediction = ifelse(phenotype=='fibrosis',
                                                                          ifelse(prediction>4,4,prediction),
                                                                          prediction))
ggscatter(predictions_results %>% filter(phenotype=='NAS') %>% filter(dataset=='Hoang') %>%
                         mutate(fold = paste0('fold no:',fold)),
                       x='score',y='prediction',cor.coef = TRUE) +
  scale_y_continuous(limits = c(0,8))+
  xlab('Measured')+ ylab('Predicted')+
  geom_abline(slope=1)+
  ggtitle('NAS prediction in Hoang dataset')+
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'none')+
  facet_wrap(~fold)
ggsave('../figures/PLSR_predictions_nas_hoang.eps',
       device = cairo_ps,
       width = 16,
       height = 10,
       units = 'in',
       dpi = 600)
ggscatter(predictions_results %>% filter(phenotype=='fibrosis') %>% filter(dataset=='Hoang') %>%
            mutate(fold = paste0('fold no:',fold)),
          x='score',y='prediction',cor.coef = TRUE) +
  scale_y_continuous(limits = c(0,4))+
  xlab('Measured')+ ylab('Predicted')+
  geom_abline(slope=1)+
  ggtitle('Fibrosis stage prediction in Hoang dataset')+
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'none')+
  facet_wrap(~fold)
ggsave('../figures/PLSR_predictions_fibrosis_hoang.eps',
       device = cairo_ps,
       width = 16,
       height = 10,
       units = 'in',
       dpi = 600)
ggscatter(predictions_results %>% filter(phenotype=='NAS') %>% filter(dataset=='Pantano') %>%
            mutate(fold = paste0('fold no:',fold)),
          x='score',y='prediction',cor.coef = TRUE) +
  scale_y_continuous(limits = c(0,8))+
  xlab('Measured')+ ylab('Predicted')+
  geom_abline(slope=1)+
  ggtitle('NAS prediction in Pantano dataset')+
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'none')+
  facet_wrap(~fold)
ggsave('../figures/PLSR_predictions_nas_pantano.eps',
       device = cairo_ps,
       width = 16,
       height = 10,
       units = 'in',
       dpi = 600)
ggscatter(predictions_results %>% filter(phenotype=='fibrosis') %>% filter(dataset=='Pantano') %>%
            mutate(fold = paste0('fold no:',fold)),
          x='score',y='prediction',cor.coef = TRUE) +
  scale_y_continuous(limits = c(0,4))+
  xlab('Measured')+ ylab('Predicted')+
  geom_abline(slope=1)+
  ggtitle('Fibrosis stage prediction in Pantano dataset')+
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5),
        legend.position = 'none')+
  facet_wrap(~fold)
ggsave('../figures/PLSR_predictions_fibrosis_pantano.eps',
       device = cairo_ps,
       width = 16,
       height = 10,
       units = 'in',
       dpi = 600)
