library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(caret)
library(ropls)
source("../utils/plotting_functions.R")
source("functions_translation_jose.R")
source("CrossValidationUtilFunctions.R")

### Load in-vitro and in-vivo datasets and plit for Cross-Validation-----------------------------------------
dataset_names <- c("Hoang", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Hoang"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Hoang$metadata  %>% select(nafld_activity_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()


### Step 1: Determine number of LVs to keep for PLSR modeling--------------------------------------------
### Split into Train - Validation -Test
trials <- 10000
flag_nas <- FALSE
flag_fibr <- FALSE
i <- 1
best_mse <- 1e6
while((i < trials)){
  TrainValIndex <- createDataPartition(Yh[,1], p = .8, 
                                       list = FALSE, 
                                       times = 1)
  if(all(unique(Yh[TrainValIndex,1]) %in% unique(Yh[-TrainValIndex,1]))){
    flag_nas <- TRUE
    if(all(unique(Yh[TrainValIndex,2]) %in% unique(Yh[-TrainValIndex,2]))){
      flag_fibr <- TRUE
      hist_test_nas <- hist(Yh[-TrainValIndex,1],freq = T,breaks = seq(0,6,0.25),plot = F)
      hist_train_nas <- hist(Yh[TrainValIndex,1],freq = T,breaks = seq(0,6,0.25),plot = F)
      hist_test_fibr <- hist(Yh[-TrainValIndex,2],freq = T,breaks = seq(0,6,0.25),plot = F)
      hist_train_fibr <- hist(Yh[TrainValIndex,2],freq = T,breaks = seq(0,6,0.25),plot = F)
      mse <- mean(c(MLmetrics::MSE(hist_test_nas$density,hist_train_nas$density),
                    MLmetrics::MSE(hist_test_fibr$density,hist_train_fibr$density)))
      if (mse < best_mse) {
        best_mse <- mse
        best_sample <- TrainValIndex
        cat(paste0('Found better split in iteration ',i,' with mse=',best_mse,
                   '\nAll validation NAS are included in train NAS:',flag_nas,
                   '\nAll validation Fibrosis are included in train Fibrosis:',flag_fibr,"\n","\n"))
      }
    }
  }
  i <- i+1
}
hist(Yh[-best_sample,1],main='Test NAS',breaks = seq(0,6,0.25))
hist(Yh[best_sample,1],main='CrossVal NAS',breaks = seq(0,6,0.25))
hist(Yh[-best_sample,2],main='Test Fibrosis',breaks = seq(0,4,0.25))
hist(Yh[best_sample,2],main='CrossVal Fibrosis',breaks = seq(0,4,0.25))
Yh_test <- Yh[-best_sample,]
Xh_test <- Xh[-best_sample,]
Yh <- Yh[best_sample,]
Xh <- Xh[best_sample,]
# saveRDS(Xh,'../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh.rds')
# saveRDS(Yh,'../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh.rds')
# saveRDS(Xh_test,'../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_test.rds')
# saveRDS(Yh_test,'../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_test.rds')

### Perform 10-fold-cross-validation to tune the number of LVs of PLSR
### Take the best # of LVs on average across all validation sets
### Calculate performance of each of the 10 models for the selected number in validation and test
### Compare with randomized and null models
Xh <- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh.rds')
Yh <- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh.rds')
Xh_test<- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_test.rds')
Yh_test<- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_test.rds')

# First make the 10-fold-split
num_folds <- 10
trials <- 10000
flag_nas <- FALSE
flag_fibr <- FALSE
i <- 1
best_mse <- 1e6
### ATTENTION : I HAVE COMMENTED OUT SAVERDS SO THAT I DO NOT ACCIDENTALLY RUN THE LOOP AND OVERWRITE MY SAVED DATA
while((i < trials)){
  data_splits <- createMultiFolds(y = Yh[,2], k = num_folds, times = 1)
  j <- 1
  mse_all <- NULL
  flag_nas_all <- NULL
  flag_fibr_all <- rep(FALSE,num_folds)
  for (ind in data_splits){
    val_Y <- Yh[-ind,]
    if(all(unique(val_Y[,1]) %in% unique(Yh[ind,1]))){
      flag_nas_all[j] <- TRUE
    }else{
      flag_nas_all[j] <- FALSE
    }
    if(all(unique(val_Y[,2]) %in% unique(Yh[ind,2]))){
      flag_fibr_all[j] <- TRUE
    }
    if (sum(flag_nas_all==FALSE)>2){
      break
    }
    hist_val_nas <- hist(val_Y[,1],freq = T,breaks = seq(0,6,0.25),plot = F)
    hist_train_nas <- hist(Yh[ind,1],freq = T,breaks = seq(0,6,0.25),plot = F)
    hist_val_fibr <- hist(val_Y[,2],freq = T,breaks = seq(0,6,0.25),plot = F)
    hist_train_fibr <- hist(Yh[ind,2],freq = T,breaks = seq(0,6,0.25),plot = F)
    mse <- mean(c(MLmetrics::MSE(hist_val_nas$density,hist_train_nas$density),
                  MLmetrics::MSE(hist_val_fibr$density,hist_train_fibr$density)))
    mse_all[j] <- mse
    j <- j+1
  }
  if ((mean(mse_all) < best_mse) & (sum(flag_nas_all==FALSE)<=2)) {
    best_mse <- mean(mse_all)
    j <- 1
    for (ind in data_splits){
      # saveRDS(Xh[ind,],paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_train',j,'.rds'))
      # saveRDS(Yh[ind,],paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_train',j,'.rds'))
      # saveRDS(Xh[-ind,],paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_val',j,'.rds'))
      # saveRDS(Yh[-ind,],paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_val',j,'.rds'))
      j <- j+1
    }
    cat(paste0('Found better split in iteration ',i,' with mse=',best_mse,
               '\nAll validation NAS are included in train NAS:',sum(flag_nas_all),'/',num_folds,
               '\nAll validation Fibrosis are included in train Fibrosis:',sum(flag_fibr_all),'/',num_folds,
               "\n","\n"))
  }
  i <- i+1
}

# for (j in 1:10){
#   # x_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/Xh_train',j,'.rds'))
#   y_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/Yh_train',j,'.rds'))
#   # x_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/Xh_val',j,'.rds'))
#   y_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/Yh_val',j,'.rds'))
#   hist(y_train[,1],main=paste0('Train ',j,' NAS'),breaks = seq(0,6,0.25))
#   hist(y_val[,1],main=paste0('Validation ',j,' NAS'),breaks = seq(0,6,0.25))
#   # hist(y_train[,2],main=paste0('Train ',j,' Fibrosis'),breaks = seq(0,4,0.25))
#   # hist(y_val[,2],main=paste0('Validation ',j,' Fibrosis'),breaks = seq(0,4,0.25))
# }

### Begin PLSR modeling tuning
val_mae <- NULL
train_mae <- NULL
test_mae <- NULL
val_R2 <- NULL
train_R2 <- NULL
test_R2 <- NULL
val_r <- NULL
train_r <- NULL
test_r <- NULL
num_LVs <- seq(1,20)
num_folds <- 10
tuning_df <- data.frame()
all_models <- NULL
for (i in 1:length(num_LVs)){
  for (j in 1:num_folds){
    x_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_train',j,'.rds'))
    y_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_train',j,'.rds'))
    x_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_val',j,'.rds'))
    y_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_val',j,'.rds'))
    
    plsr_model <- opls(x = x_train, 
                       y = y_train,
                       predI = num_LVs[i],
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    y_train_hat <- predict(plsr_model,x_train)
    y_val_hat <- predict(plsr_model,x_val)
    # y_hat_test <- predict(plsr_model,Xh_test)
    
    train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
    # test_r[j]<- mean(diag(cor(y_hat_test,Yh_test)))
    train_R2[j] <- mean(diag(cor(y_train_hat,y_train)^2))
    val_R2[j] <- mean(diag(cor(y_val_hat,y_val)^2))
    # test_R2[j] <- mean(diag(cor(y_hat_test,Yh_test)^2))
    val_mae[j] <- mean(abs(y_val_hat-y_val))
    train_mae[j] <- mean(abs(y_train_hat-y_train))
    # test_mae[j] <- mean(abs(y_hat_test-Yh_test))
  }
  tuning_df <- rbind(tuning_df,
                     data.frame(set='train',r = train_r,MAE = train_mae,R2=train_R2,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)),
                     data.frame(set='validation',r = val_r,MAE = val_mae,R2=val_R2,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)))
  print(paste0('Finished fitting PLSR with ',num_LVs[i],' latent variables'))
}
# saveRDS(tuning_df,'../preprocessing/TrainingValidationData/WholePipeline/tuning_df.rds')

plotting_df <- tuning_df  %>% gather('metric','value',-set,-fold,-LVs) %>% 
  group_by(set,metric,LVs) %>% mutate(mu = mean(value)) %>% mutate(std = sd(value)) %>% ungroup() %>%
  select(set,metric,LVs,mu,std) %>% unique() %>%
  mutate(metric=ifelse(metric=='R2','R\u00B2',metric))
ggplot(plotting_df,aes(x=LVs,y=mu,color=set)) +
  geom_point() +
  geom_line(lwd=1)+
  geom_errorbar(aes(ymax = mu + std/num_folds, ymin = mu - std/num_folds))+
  scale_x_continuous(breaks = seq(1,20,2))+
  xlab('number of latent variables') + ylab('value') +
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family = 'Arial'),
        panel.grid.major = element_line(),
        strip.text = element_text(face = 'bold'))+
  facet_wrap(~metric,scales = 'free_y')

ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr.png',
       width = 16,
       height = 9,
       units = 'in',
       dpi=600)

### After selecting 9 LVs as the tuned parameter re-fit with only those
### But also calculate shuffled and null models performance
val_mae <- NULL
train_mae <- NULL
test_mae <- NULL
val_r <- NULL
train_r <- NULL
test_r <- NULL
num_folds <- 10
tuning_df <- data.frame()
all_models <- NULL
test_r_shuffle_y <- NULL
test_mae_shuffle_y <- NULL
test_r_shuffle_x <- NULL
test_mae_shuffle_x <- NULL
test_r_random_x <- NULL
test_mae_random_x <- NULL
for (j in 1:num_folds){
  message(paste0('Begun fold ',j))
  x_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_train',j,'.rds'))
  y_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_train',j,'.rds'))
  x_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_val',j,'.rds'))
  y_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_val',j,'.rds'))
  
  plsr_model <- opls(x = x_train, 
                     y = y_train,
                     predI = 9,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  y_train_hat <- predict(plsr_model,x_train)
  y_val_hat <- predict(plsr_model,x_val)
  y_hat_test <- predict(plsr_model,Xh_test)
  
  train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
  val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
  test_r[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  print(paste0('NAS = ',diag(cor(y_hat_test,Yh_test))['NAS'],' , Fibrosis = ',diag(cor(y_hat_test,Yh_test))['fibrosis']))
  # train_R2[j] <- mean(diag(cor(y_train_hat,y_train)^2))
  # val_R2[j] <- mean(diag(cor(y_val_hat,y_val)^2))
  # test_R2[j] <- mean(diag(cor(y_hat_test,Yh_test)^2))
  val_mae[j] <- mean(abs(y_val_hat-y_val))
  train_mae[j] <- mean(abs(y_train_hat-y_train))
  test_mae[j] <- mean(abs(y_hat_test-Yh_test))
  
  ### shuffled labels model
  y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
  rownames(y_train_shuffled) <- rownames(y_train)
  plsr_model_shuffle_y <- opls(x = x_train, 
                     y = y_train_shuffled,
                     predI = 9,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  y_hat_test <- predict(plsr_model_shuffle_y,Xh_test)
  test_r_shuffle_y[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  # test_R2_shuffle_y[j] <- mean(diag(cor(y_hat_test,Yh_test)^2))
  test_mae_shuffle_y[j] <- mean(abs(y_hat_test-Yh_test))
  
  print('Finished shuffled labels model')
  
  ### shuffled features model
  x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
  colnames(x_train_shuffled) <- colnames(x_train)
  plsr_model_shuffle_x <- opls(x = x_train_shuffled, 
                               y = y_train,
                               predI = 9,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  y_hat_test <- predict(plsr_model_shuffle_x,Xh_test)
  test_r_shuffle_x[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  # test_R2_shuffle_x[j] <- mean(diag(cor(y_hat_test,Yh_test)^2))
  test_mae_shuffle_x[j] <- mean(abs(y_hat_test-Yh_test))
  
  print('Finished shuffled features model')
  
  ### random features model
  x_train_random <- matrix(rnorm(n = nrow(x_train)*ncol(x_train)),nrow = nrow(x_train))
  colnames(x_train_random) <- colnames(x_train)
  rownames(x_train_random) <- rownames(x_train)
  plsr_model_random_x <- opls(x = x_train_random, 
                               y = y_train,
                               predI = 9,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  y_hat_test <- predict(plsr_model_random_x,Xh_test)
  test_r_random_x[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  # test_R2_random_x[j] <- mean(diag(cor(y_hat_test,Yh_test)^2))
  test_mae_random_x[j] <- mean(abs(y_hat_test-Yh_test))
  
  print('Finished random features model')
}
performance_df <- rbind(data.frame(set='train',r = train_r,MAE = train_mae,fold = seq(1,num_folds)),
                        data.frame(set='validation',r = val_r,MAE = val_mae,fold = seq(1,num_folds)),
                        data.frame(set='test',r = test_r,MAE = test_mae,fold = seq(1,num_folds)),
                        data.frame(set='shuffle Y',r = test_r_shuffle_y,MAE = test_mae_shuffle_y,fold = seq(1,num_folds)),
                        data.frame(set='shuffle X',r = test_r_shuffle_x,MAE = test_mae_shuffle_x,fold = seq(1,num_folds)),
                        data.frame(set='random X',r = test_r_random_x,MAE = test_mae_random_x,fold = seq(1,num_folds)))
saveRDS(performance_df,'../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned.rds')

avg_mae <- mean(abs(Yh_test-apply(Yh,2,mean)))
plotting_performance_df <- performance_df %>% 
  gather('metric','value',-set,-fold) %>%
  mutate(metric=ifelse(metric=='r','pearson`s r',metric))
plotting_performance_df$set <- factor(plotting_performance_df$set,
                                      levels = c('train','validation','test','shuffle Y','shuffle X','random X'))
(ggboxplot(plotting_performance_df %>% filter(metric=='MAE'),x='set',y='value',color = 'set',add='jitter')+
    scale_y_continuous(n.breaks = 15)+
    geom_hline(yintercept = avg_mae,linetype='dashed',color='black',lwd=1)+
    annotate('text',x=3,y=1.75,size=6,label = 'error from the mean of the data')+
    xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1))+
  stat_compare_means(comparisons = list(c('test','shuffle Y'),
                                        c('test','shuffle X'),
                                          c('test','random X')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01,
                     label.y = c(1.3,1.4,1.5)) +
  facet_wrap(~metric)) +
  (ggboxplot(plotting_performance_df %>% filter(metric!='MAE'),x='set',y='value',color = 'set',add='jitter')+
     scale_y_continuous(breaks = seq(-0.5,1,0.1))+
     # geom_hline(yintercept = 0,linetype='dashed',color='black',lwd=1)+
     xlab('')+
     theme(text = element_text(size=20,family = 'Arial'),
           legend.position = 'none',
           axis.text.x = element_text(size=16),
           strip.text = element_text(face = 'bold'),
           panel.grid.major.y = element_line(linewidth = 1))+
     stat_compare_means(comparisons = list(c('test','shuffle Y'),
                                           c('test','shuffle X'),
                                           c('test','random X')),
                        method = 'wilcox',label = 'p.signif',
                        tip.length = 0.01,
                        label.y = c(0.75,0.85,0.95)) +
     facet_wrap(~metric))

ggsave('../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned.png',
       width = 14,
       height = 9,
       units = 'in',
       dpi = 600)


### Since we do not overfit to validation now perform cross-validation with this dataset-------------
### Where you calculate validation and train performance of:
### a) human data PLSR to predict NAS,Fibrosis
### b) back-projected human data PLSR to predict NAS, Fibrosis
### c) back-projected human data PLSR to predict NAS, Fibrosis with translatable combination of MPS PCs
### d) c) back-projected human data PLSR to predict NAS, Fibrosis with optimal direction
num_folds <- 10
loc <- '../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/'
num_LVS <- 9
dataset_names <- c("Hoang", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Hoang"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Hoang$metadata  %>% select(nafld_activity_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()
# pca_humans <- prcomp(Xh,scale. = F)
# Wpch <- pca_humans$rotation
### Re-split without any test set since we so no drop in performance
trials <- 10000
flag_nas <- FALSE
flag_fibr <- FALSE
i <- 1
best_mse <- 1e6
### ATTENTION : I HAVE COMMENTED OUT SAVERDS SO THAT I DO NOT ACCIDENTALLY RUN THE LOOP AND OVERWRITE MY SAVED DATA
while((i < trials)){
  data_splits <- createMultiFolds(y = Yh[,2], k = num_folds, times = 1)
  j <- 1
  mse_all <- NULL
  flag_nas_all <- NULL
  flag_fibr_all <- NULL
  for (ind in data_splits){
    val_Y <- Yh[-ind,]
    if(all(unique(val_Y[,1]) %in% unique(Yh[ind,1]))){
      flag_nas_all[j] <- TRUE
    }else{
      flag_nas_all[j] <- FALSE
    }
    if(all(unique(val_Y[,2]) %in% unique(Yh[ind,2]))){
      flag_fibr_all[j] <- TRUE
    }else{
      flag_fibr_all[j] <- FALSE
    }
    if ((sum(flag_nas_all==FALSE)>2) | (sum(flag_fibr_all==FALSE)>2 & i>10)){
      break
    }
    hist_val_nas <- hist(val_Y[,1],freq = T,breaks = seq(0,6,0.25),plot = F)
    hist_train_nas <- hist(Yh[ind,1],freq = T,breaks = seq(0,6,0.25),plot = F)
    hist_val_fibr <- hist(val_Y[,2],freq = T,breaks = seq(0,6,0.25),plot = F)
    hist_train_fibr <- hist(Yh[ind,2],freq = T,breaks = seq(0,6,0.25),plot = F)
    mse <- mean(c(MLmetrics::MSE(hist_val_fibr$density,hist_train_nas$density),
                  MLmetrics::MSE(hist_val_fibr$density,hist_train_fibr$density)))
    mse_all[j] <- mse
    j <- j+1
  }
  if ((mean(mse_all) < best_mse) & (sum(flag_nas_all==FALSE)<=2) &(sum(flag_fibr_all==FALSE)<=2)) {
    best_mse <- mean(mse_all)
    j <- 1
    for (ind in data_splits){
      # saveRDS(Xh[ind,],paste0(loc,'Xh_train',j,'.rds'))
      # saveRDS(Yh[ind,],paste0(loc,'Yh_train',j,'.rds'))
      # saveRDS(Xh[-ind,],paste0(loc,'Xh_val',j,'.rds'))
      # saveRDS(Yh[-ind,],paste0(loc,'Yh_val',j,'.rds'))
      j <- j+1
    }
    cat(paste0('Found better split in iteration ',i,' with mse=',best_mse,
               '\nAll validation NAS are included in train NAS:',sum(flag_nas_all),'/',num_folds,
               '\nAll validation Fibrosis are included in train Fibrosis:',sum(flag_fibr_all),'/',num_folds,
               "\n","\n"))
  }
  i <- i+1
}
# Task a)
performance_1 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_plsr',
                                                    LVs = num_LVS)
# saveRDS(performance_1,'../results/performance_df_human_plsr.rds')
#Plot
performance_1 <- performance_1 %>% mutate(type = ifelse(type=='model',set,type))
performance_1$type <- factor(performance_1$type ,levels=c('train','test','shuffle Y','shuffle X','random X'))
ggboxplot(performance_1,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(-0.75,1,0.25))+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('test','shuffle Y'),
                                        c('test','shuffle X'),
                                        c('test','random X')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     label.y = c(0.82, 0.87, 0.92))
ggsave('../results/performance_df_just_human_plsr.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)

# Task b)
performance_2 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected',
                                                    LVs = num_LVS)
# saveRDS(performance_2,'../results/performance_df_human_backprojected.rds')
performance_2$type <- factor(performance_2$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_2$set <- factor(performance_2$set ,levels=c('train','test'))
ggboxplot(performance_2,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     step.increase = 0.04)
ggsave('../results/performance_df_just_human_backprojected.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)

# Task c)
performance_3 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_retrained',
                                                    LVs = num_LVS)
# saveRDS(performance_3,'../results/performance_df_human_backprojected_retrained.rds')
#Plot
performance_3$type <- factor(performance_3$type ,levels=c('model','shuffle W'))
performance_3$set <- factor(performance_3$set ,levels=c('train','test'))
ggboxplot(performance_3,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01)
ggsave('../results/performance_df_just_human_backprojected_retrained.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)

# Task d1)
performance_4 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_into_translatable_lvs',
                                                    LVs = num_LVS)
# saveRDS(performance_4,'../results/performance_df_translatable_lvs.rds')
#Plot
performance_4$type <- factor(performance_4$type ,levels=c('model','shuffle W','shuffle Bh'))
performance_4$set <- factor(performance_4$set ,levels=c('train','test'))
ggboxplot(performance_4,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     step.increase = 0.04)
ggsave('../results/performance_df_translatable_lvs.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)
# Task d2)
performance_5 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_into_optimal_lvs',
                                                    LVs = num_LVS)
saveRDS(performance_5,'../results/performance_df_optimal_direction.rds')
#plot
performance_5$type <- factor(performance_5$type ,levels=c('model','shuffle Wopt','shuffle Bh'))
performance_5$set <- factor(performance_5$set ,levels=c('train','test'))
ggboxplot(performance_5,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(round(min(performance_5$r),1),1,0.1),
                     limits = c(NA,1.15))+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle Wopt'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01,
                     step.increase = 0.04)

ggsave('../results/performance_df_optimal_direction.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)

performance_6 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'analytical_optimal',
                                                    LVs = 9)
saveRDS(performance_6,'../results/performance_df_analytical.rds')

performance_6$type <- factor(performance_6$type ,levels=c('model','shuffle Wopt','shuffle Bh'))
performance_6$set <- factor(performance_6$set ,levels=c('train','test'))
ggboxplot(performance_6,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(round(min(performance_6$r),1),1,0.1),
                     limits = c(NA,1.15))+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle Wopt'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01)# ,label.y = c(0.82, 0.87, 0.92)

ggsave('../results/performance_df_analytical.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)


### Combine all results and make a nice figure--------------------------------------------------------------------
performance_1 <- readRDS('../results/performance_df_human_plsr.rds')
performance_2 <- readRDS('../results/performance_df_human_backprojected.rds')
performance_3 <- readRDS('../results/performance_df_human_backprojected_retrained.rds')
performance_4 <- readRDS('../results/performance_df_translatable_lvs.rds')
performance_5 <- readRDS('../results/performance_df_optimal_direction.rds')
performance_6 <- readRDS('../results/performance_df_analytical.rds')

performance_all <- rbind(performance_1,
                         performance_2,
                         performance_3,
                         performance_4,
                         performance_5,
                         performance_6)
### Select what to show in the figure
performance_all_plot <- performance_all %>% 
  mutate(keep=ifelse(task=='human_plsr',ifelse(type %in% c('model','shuffle X'),TRUE,FALSE),
                     ifelse(type=='model',TRUE,FALSE))) %>% 
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(approach = ifelse(task=='human_plsr','PLSR',
                        ifelse(task=='human_backprojected','backprojected',
                               ifelse(task=='human_backprojected_retrained','backprojected retrained',
                                      ifelse(task=='human_backprojected_into_translatable_lvs','translatable LVs',
                                             ifelse(task=='analytical_optimal','analytical Wopt',
                                             'Wopt')))))) %>%
  mutate(approach = ifelse(approach=='PLSR',
                           ifelse(type=='model',approach,type),
                           approach)) %>%
  select(-type,-task)

performance_all_plot$approach <- factor(performance_all_plot$approach,
                                        levels = c('PLSR',
                                                   'backprojected',
                                                   'backprojected retrained',
                                                   'translatable LVs',
                                                   'Wopt',
                                                   'analytical Wopt',
                                                   'shuffle X'))
p_train <- ggboxplot(performance_all_plot %>% filter(set=='train'),x='approach',y='r',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(0.5,1,0.05),limits = c(0.5,NA))+
  xlab('')+
  ggtitle('10-fold Train')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=16,angle = 25),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('backprojected','translatable LVs'),
                                        c('backprojected','backprojected retrained'),
                                        c('PLSR','backprojected'),
                                        c('analytical Wopt','backprojected'),
                                        c('PLSR','Wopt'),
                                        c('PLSR','analytical Wopt'),
                                        c('Wopt','analytical Wopt')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     label.y = c(0.67,0.78,0.98,0.98,1,1.02,1.04))
print(p_train)  
ggsave('../results/approaches_comparison_training.png',
       plot = p_train,
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)

p_test <- ggboxplot(performance_all_plot %>% filter(set=='test'),x='approach',y='r',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(-0.5,1,0.1))+
  xlab('')+
  ggtitle('10-fold Test')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=16,angle = 25),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('backprojected','translatable LVs'),
                                        c('backprojected','backprojected retrained'),
                                        c('PLSR','backprojected'),
                                        c('analytical Wopt','backprojected'),
                                        c('PLSR','Wopt'),
                                        c('PLSR','analytical Wopt'),
                                        c('Wopt','analytical Wopt')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     label.y = c(0.75,0.8,0.85,0.85,0.9,0.95,0.95))
print(p_test)  
ggsave('../results/approaches_comparison_10foldtest.png',
       plot = p_test,
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)


# Explore how the evolutionary algorithm VS analytical approach works-------------------------------------

#Load data
dataset_names <- c("Hoang", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Hoang"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Hoang$metadata  %>% select(nafld_activity_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
rownames(Yh) <- rownames(Xh)
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# Parameters of PLSR
num_LVS <- 9

# Begin iterations
partitions <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2)
iterations <- 10
analytical_performance <- data.frame()
evol_performance <- data.frame()
plsr_performance <- data.frame()
number_of_optimal_LVs <- data.frame()
analytical_performance_val <- data.frame()
evol_performance_val <- data.frame()
plsr_performance_val <- data.frame()
for (p in partitions){
  message(paste0('Begun partition:',100*p,'%'))
  for (i in 1:iterations){
    Xh_parted <- as.matrix(as.data.frame(Xh) %>% sample_frac(size=p))
    Yh_parted <- Yh[rownames(Xh_parted),]
    
    ### Keep the rest randomly for validation
    if (p <1){
      Xh_val <- Xh[which(!(rownames(Xh) %in% rownames(Xh_parted))),]
      Yh_val <- Yh[rownames(Xh_val),]
    }
    
    #### First run human PLSR
    plsr_model <- suppressMessages(opls(x = Xh_parted, 
                       y = Yh_parted,
                       predI = num_LVS,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none"))
    y_hat_plsr <- predict(plsr_model,Xh_parted)
    if (p <1){
      y_hat_val_plsr <- predict(plsr_model,Xh_val)
    }
    # Get Wh of PLSR
    Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh_parted))
    rownames(Wh) <- colnames(Xh_parted)
    colnames(Wh) <- colnames(plsr_model@weightMN)
    for (ii in 1:nrow(plsr_model@weightMN)){
      Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    }
    # Get regression coefficients
    Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
    # Define projection matrices to make more readable
    Th <- Xh_parted %*% Wh
    Thm <- Xh_parted %*% Wm %*% t(Wm) %*% Wh
    
    ### Run evolutionary algorithm
    Wm_opt <- get_translatable_LV(Xh_parted, Yh_parted, Wh, Wm,
                                  rbind(apply(Yh_parted,2,mean),Bh),
                                  find_extra = TRUE,
                                  verbose = FALSE)

    Wm_opt <- Wm_opt$Wm_new
    number_of_optimal_LVs <- rbind(number_of_optimal_LVs,
                                   data.frame(model = 'evolutionary',
                                              partition = p,
                                              sample=i,
                                              num_lvs=ncol(Wm_opt)))
    colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
    # Extend latent variables
    Wm_tot <- cbind(Wm, Wm_opt)
    # predict
    y_hat_evol <- cbind(1, Xh_parted %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(Yh_parted,2,mean),Bh)
    if (p <1){
      y_hat_val_evol <- cbind(1, Xh_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(Yh_parted,2,mean),Bh)
    }
    
    ### Run analytical approach
    phi <- Wh %*% Bh
    Wm_opt <- analytical_solution_opt(y=Yh_parted,
                                      W_invitro = Wm,
                                      phi = phi)
    # Extend latent variables
    Wm_tot <- cbind(Wm, Wm_opt)
    # predict
    y_hat_analytical <- cbind(1, Xh_parted %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(Yh_parted,2,mean),Bh)
    if (p <1){
      y_hat_val_analytical <- cbind(1, Xh_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(Yh_parted,2,mean),Bh)
    }
    
    ### Evaluate pearson correlation
    r_analytical <- diag(cor(y_hat_analytical,Yh_parted))
    analytical_performance <- rbind(analytical_performance,
                                    data.frame(model = 'analytical',
                                               partition = p,
                                               sample=i,
                                               NAS = r_analytical['NAS'],
                                               fibrosis = r_analytical['fibrosis']))
    r_evol <- diag(cor(y_hat_evol,Yh_parted))
    evol_performance <- rbind(evol_performance,
                                    data.frame(model = 'evolutionary',
                                               partition = p,
                                               sample=i,
                                               NAS = r_evol['NAS'],
                                               fibrosis = r_evol['fibrosis']))
    r_plsr <- diag(cor(y_hat_plsr,Yh_parted))
    plsr_performance <- rbind(plsr_performance,
                                    data.frame(model = 'PLSR',
                                               partition = p,
                                               sample=i,
                                               NAS = r_plsr['NAS'],
                                               fibrosis = r_plsr['fibrosis']))
    
    if (p <1){
      r_analytical <- diag(cor(y_hat_val_analytical,Yh_val))
      analytical_performance_val <- rbind(analytical_performance_val,
                                      data.frame(model = 'analytical',
                                                 partition = p,
                                                 sample=i,
                                                 NAS = r_analytical['NAS'],
                                                 fibrosis = r_analytical['fibrosis']))
      r_evol <- diag(cor(y_hat_val_evol,Yh_val))
      evol_performance_val <- rbind(evol_performance_val,
                                data.frame(model = 'evolutionary',
                                           partition = p,
                                           sample=i,
                                           NAS = r_evol['NAS'],
                                           fibrosis = r_evol['fibrosis']))
      r_plsr <- diag(cor(y_hat_val_plsr,Yh_val))
      plsr_performance_val <- rbind(plsr_performance_val,
                                data.frame(model = 'PLSR',
                                           partition = p,
                                           sample=i,
                                           NAS = r_plsr['NAS'],
                                           fibrosis = r_plsr['fibrosis']))
    }
    
    message(paste0('Finished iteration ',i))

  }
  cat(paste0('Performance for partition ',100*p,'%',
             '\nPLSR: NAS:',mean(plsr_performance_val$NAS),' , Fibrosis:',mean(plsr_performance_val$fibrosis),
             '\nEvolutionary algorithm: NAS:',mean(evol_performance_val$NAS),' , Fibrosis:',mean(evol_performance_val$fibrosis),
             '\nAnalytical approach: NAS:',mean(analytical_performance_val$NAS),' , Fibrosis:',mean(analytical_performance_val$fibrosis),
             '\n'))
}
all_performance_train <- rbind(plsr_performance,analytical_performance,evol_performance)
all_performance_val <- rbind(plsr_performance_val,analytical_performance_val,evol_performance_val)
# saveRDS(all_performance_train,'../results/all_performance_train_analytical_vs_evol_parted.rds')
# saveRDS(all_performance_val,'../results/all_performance_val_analytical_vs_evol_parted.rds')
# saveRDS(number_of_optimal_LVs,'../results/number_of_optimal_LVs_parted.rds')

all_performance_train_4_plotting <- all_performance_train %>% gather('phenotype','r',-model,-partition,-sample) %>% 
  group_by(model,partition,phenotype) %>% mutate(mu = mean(r)) %>% mutate(std=sd(r)) %>% ungroup()
all_performance_val_4_plotting <- all_performance_val %>% gather('phenotype','r',-model,-partition,-sample) %>% 
  group_by(model,partition,phenotype) %>% mutate(mu = mean(r)) %>% mutate(std=sd(r)) %>% ungroup()
p <- (ggplot(all_performance_train_4_plotting %>% select(model,partition,phenotype,mu,std) %>% unique(),
       aes(x=100*partition,y=mu,color=model,fill=model)) +
  geom_errorbar(aes(ymax = mu + std/sqrt(iterations), ymin = mu - std/sqrt(iterations)),width = 2)+
  geom_point(aes(y=mu))+
  geom_line() +
  xlab('') + ylab('pearson`s r') +
  ggtitle('Training performance')+
  ylim(c(0.5,1))+
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~phenotype)) / 
  (ggplot(all_performance_val_4_plotting %>% select(model,partition,phenotype,mu,std) %>% unique(),
          aes(x=100*partition,y=mu,color=model,fill=model)) +
     geom_errorbar(aes(ymax = mu + std/sqrt(iterations), ymin = mu - std/sqrt(iterations)),width = 2)+
     geom_point(aes(y=mu))+
     geom_line() +
     xlab('partition (%)') + ylab('pearson`s r') +
     ggtitle('Test performance')+
     ylim(c(0.5,1))+
     theme_pubr(base_size = 20,base_family = 'Arial')+
     theme(text = element_text(size=20,family = 'Arial'),
           legend.position = 'right',
           plot.title = element_text(hjust = 0.5))+
     facet_wrap(~phenotype))
p <- p + plot_layout(guides = "collect")
print(p)
ggsave(filename = '../results/compare_analytical_vs_evol_across_partitions.png',
       plot = p,
       height = 9,
       width=12,
       unit = 'in',
       dpi=600)

ggplot(number_of_optimal_LVs %>% mutate(partition=100*partition) %>% mutate(partition=as.factor(partition)),
       aes(x=partition,y=num_lvs)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0.1) +
  xlab('partition (%)') + ylab('number of extra latent vectors used') +
  scale_y_continuous(breaks = c(1,2,3,4,5),limits = (c(1,5.5)))+
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))+ 
  stat_compare_means(method = 'kruskal',size=5)
ggsave(filename = '../results/evol_number_lvs_inferred.png',
       height = 9,
       width=12,
       unit = 'in',
       dpi=600)

### Check how similar are the the LVs found by the evolutionary algorithm
### with the analytical solution
### with the averaged LV across multiple runs of the evolutionary algorithm
### with the result of the evolutionary algorithm from shuffled data
### with the result of the analytical approach from shuffled data
iterations <- 100
similarity_df <- data.frame()
for (i in 1:iterations){
  Xh_parted <- as.matrix(as.data.frame(Xh) %>% sample_frac(size=1))
  Yh_parted <- Yh[rownames(Xh_parted),]
  #### First run human PLSR
  plsr_model <- suppressMessages(opls(x = Xh_parted, 
                                      y = Yh_parted,
                                      predI = num_LVS,
                                      crossvalI = 1,
                                      scaleC = "center",
                                      fig.pdfC = "none",
                                      info.txtC = "none"))
  # Get Wh of PLSR
  Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh_parted))
  rownames(Wh) <- colnames(Xh_parted)
  colnames(Wh) <- colnames(plsr_model@weightMN)
  for (ii in 1:nrow(plsr_model@weightMN)){
    Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
  }
  # Get regression coefficients
  Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
  # Define projection matrices to make more readable
  Th <- Xh_parted %*% Wh
  Thm <- Xh_parted %*% Wm %*% t(Wm) %*% Wh
  
  ### Run evolutionary algorithm
  Wm_opt_evol <- get_translatable_LV(Xh_parted, Yh_parted, Wh, Wm,
                                     rbind(apply(Yh_parted,2,mean),Bh),
                                     find_extra = TRUE,
                                     verbose = FALSE)
  Wm_opt_evol <- Wm_opt_evol$Wm_new
  colnames(Wm_opt_evol) <- paste0('evol_LV_',seq(1,ncol(Wm_opt_evol)))
  # Extend latent variables
  Wm_tot_evol <- cbind(Wm, Wm_opt_evol)
  
  ### Run analytical approach
  phi <- Wh %*% Bh
  # div <- sqrt(apply(phi^2,2,sum))
  # for (jj in 1:length(div)){
  #   phi[,jj] <- phi[,jj]/div[jj]
  # }
  Wm_opt_analyt <- analytical_solution_opt(y=Yh_parted,
                                    W_invitro = Wm,
                                    phi = phi)
  colnames(Wm_opt_analyt) <- paste0('analyt_LV_',seq(1,ncol(Wm_opt_analyt)))
  # Extend latent variables
  Wm_tot_analyt <- cbind(Wm, Wm_opt_analyt)
  
  ### Evaluate similarity between evol and analytical vectors
  cosine_sim <- lsa::cosine(cbind(Wm_opt_analyt,Wm_opt_evol))
  # cosine_sim[lower.tri(cosine_sim,diag = T)] <- -100
  cosine_sim <- as.data.frame(cosine_sim) %>% rownames_to_column('vector_1') %>% gather('vector_2','sim',-vector_1)
  # cosine_sim <- cosine_sim %>% filter(sim!=-100)
  # ## Find the two maximum similar things
  # cosine_sim_max1 <- cosine_sim %>% group_by(vector_1) %>% mutate(max_sim = max(abs_sim)) %>% ungroup() %>%
  #   filter(abs_sim==max_sim) %>% select(-max_sim) %>% unique()
  # cosine_sim_max2 <- anti_join(cosine_sim,cosine_sim_max1)
  # cosine_sim_max2 <- cosine_sim_max2 %>% group_by(vector_1) %>% mutate(max_sim = max(abs_sim)) %>% ungroup() %>%
  #   filter(abs_sim==max_sim) %>% select(-max_sim) %>% unique()
  # cosine_sim <- rbind(cosine_sim_max1,cosine_sim_max2)
  similarity_df <- rbind(similarity_df,
                          cosine_sim %>% mutate(iteration=i))
  if ((i %% 10 ==0) |(i==1)){
    print(paste0('Finished iteration ',i))
  }
}
# saveRDS(similarity_df,'../result/cosine_similarity_df_analyt_evol.rds')

corrr::network_plot(tt[1:4,1:4],colors = c("red",'white', "blue"),min_cor = 0.01)
ggsave('../results/analyt_evol_cosine_sim.png',
       height = 9,
       width = 9,
       units = 'in',
       dpi=600)
