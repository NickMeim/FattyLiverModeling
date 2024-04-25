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
saveRDS(performance_1,'../results/performance_df_human_plsr.rds')
#Plot
performance_1$type <- factor(performance_1$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_1$set <- factor(performance_1$set ,levels=c('train','test'))
ggboxplot(performance_1,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 6)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01)
# Task b)
performance_2 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected',
                                                    LVs = num_LVS)
saveRDS(performance_2,'../results/performance_df_human_backprojected.rds')
performance_2$type <- factor(performance_2$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_2$set <- factor(performance_2$set ,levels=c('train','test'))
ggboxplot(performance_2,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 6)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01)

# Task c)
performance_3 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_retrained',
                                                    LVs = num_LVS)
saveRDS(performance_3,'../results/performance_df_human_backprojected_retrained.rds')
#Plot
performance_3$type <- factor(performance_3$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_3$set <- factor(performance_3$set ,levels=c('train','test'))
ggboxplot(performance_3,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 6)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01)

# Task d1)
performance_4 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_into_translatable_lvs',
                                                    LVs = num_LVS)
saveRDS(performance_4,'../results/performance_df_translatable_lvs.rds')
#Plot
performance_4$type <- factor(performance_4$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_4$set <- factor(performance_4$set ,levels=c('train','test'))
ggboxplot(performance_4,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 6)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01)
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
performance_5$type <- factor(performance_5$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_5$set <- factor(performance_5$set ,levels=c('train','test'))
ggboxplot(performance_5,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 6)+
  xlab('')+
  theme(text = element_text(size=20,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=16),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',label = 'p.signif',
                     tip.length = 0.01)

ggsave('../results/performance_df_optimal_direction.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)








