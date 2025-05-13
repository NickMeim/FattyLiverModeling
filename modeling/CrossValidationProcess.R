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

### Load in-vitro and in-vivo datasets and split for Cross-Validation-----------------------------------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver","Govaere")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
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
      hist_test_nas <- hist(Yh[-TrainValIndex,1],freq = T,breaks = seq(0,8,0.25),plot = F)
      hist_train_nas <- hist(Yh[TrainValIndex,1],freq = T,breaks = seq(0,8,0.25),plot = F)
      hist_test_fibr <- hist(Yh[-TrainValIndex,2],freq = T,breaks = seq(0,8,0.25),plot = F)
      hist_train_fibr <- hist(Yh[TrainValIndex,2],freq = T,breaks = seq(0,8,0.25),plot = F)
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
hist(Yh[-best_sample,1],main='Test NAS',breaks = seq(0,8,0.25))
hist(Yh[best_sample,1],main='CrossVal NAS',breaks = seq(0,8,0.25))
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
rownames(Yh) <- rownames(Xh)
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
    hist_val_nas <- hist(val_Y[,1],freq = T,breaks = seq(0,8,0.25),plot = F)
    hist_train_nas <- hist(Yh[ind,1],freq = T,breaks = seq(0,8,0.25),plot = F)
    hist_val_fibr <- hist(val_Y[,2],freq = T,breaks = seq(0,8,0.25),plot = F)
    hist_train_fibr <- hist(Yh[ind,2],freq = T,breaks = seq(0,8,0.25),plot = F)
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
#   hist(y_train[,1],main=paste0('Train ',j,' NAS'),breaks = seq(0,8,0.25))
#   hist(y_val[,1],main=paste0('Validation ',j,' NAS'),breaks = seq(0,8,0.25))
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
train_rho <- NULL
val_rho <- NULL
test_rho <- NULL
val_R2_shuffled <- NULL
val_r_shuffled <- NULL
val_mae_shuffled <- NULL
val_rho_shuffled <- NULL
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
    
    ### Shuffled X model
    x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
    colnames(x_train_shuffled) <- colnames(x_train)
    plsr_model_shuffled <- opls(x = x_train_shuffled, 
                       y = y_train,
                       predI = num_LVs[i],
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    y_val_hat_shuffled <- predict(plsr_model_shuffled,x_val)

    train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
    val_r_shuffled[j] <- mean(diag(cor(y_val_hat_shuffled,y_val)))
    train_R2[j] <- mean(diag(cor(y_train_hat,y_train)^2))
    val_R2[j] <- mean(diag(cor(y_val_hat,y_val)^2))
    val_R2_shuffled[j] <- mean(diag(cor(y_val_hat_shuffled,y_val)^2))
    val_mae[j] <- mean(abs(y_val_hat-y_val))
    val_mae_shuffled[j] <- mean(abs(y_val_hat_shuffled-y_val))
    train_mae[j] <- mean(abs(y_train_hat-y_train))
    train_rho[j] <- mean(diag(cor(y_train_hat,y_train, method = "spearman")))
    val_rho[j] <- mean(diag(cor(y_val_hat,y_val, method = "spearman")))
    val_rho_shuffled[j] <- mean(diag(cor(y_val_hat_shuffled,y_val, method = "spearman")))
  }
  tuning_df <- rbind(tuning_df,
                     data.frame(set='train',rho=train_rho,r = train_r,MAE = train_mae,R2=train_R2,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)),
                     data.frame(set='validation',rho=val_rho,r = val_r,MAE = val_mae,R2=val_R2,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)),
                     data.frame(set='shuffled',rho=val_rho_shuffled,r = val_r_shuffled,MAE = val_mae_shuffled,R2=val_R2_shuffled,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)))
  print(paste0('Finished fitting PLSR with ',num_LVs[i],' latent variables'))
}

saveRDS(tuning_df,'../preprocessing/TrainingValidationData/WholePipeline/tuning_df_with_spearman.rds')

plotting_df <- tuning_df  %>% gather('metric','value',-set,-fold,-LVs) %>% 
  group_by(set,metric,LVs) %>% mutate(mu = mean(value)) %>% mutate(std = sd(value)) %>% ungroup() %>%
  select(set,metric,LVs,mu,std) %>% unique() %>%
  mutate(metric=ifelse(metric=='R2','R\u00B2',metric))
ggplot(plotting_df,aes(x=LVs,y=mu,color=set)) +
  geom_point() +
  geom_line(lwd=1)+
  geom_errorbar(aes(ymax = mu + std/sqrt(num_folds), ymin = mu - std/sqrt(num_folds)))+
  scale_x_continuous(breaks = seq(1,20,2))+
  scale_y_continuous(n.breaks = 12)+
  xlab('number of latent variables') + ylab('value') +
  theme_pubr(base_size = 22,base_family = 'Arial')+
  theme(text = element_text(size=22,family = 'Arial'),
        panel.grid.major = element_line(),
        strip.text = element_text(face = 'bold'))+
  facet_wrap(~metric,scales = 'free_y')

ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr.png',
       width = 10,
       height = 10,
       units = 'in',
       dpi=600)

ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr.eps',
       device = cairo_ps,
       width = 10,
       height = 10,
       units = 'in',
       dpi=600)

### Check different random partitions
num_lvs <- 20
iterations <- 10
partitions <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2)
performance_train <- data.frame() 
performance_val <- data.frame()
for (p in partitions){
  message(paste0('Begun partition ',100*p,' %'))
  if (round(p*nrow(Xh))<num_lvs){
    num_lvs <- round(p*nrow(Xh))-1
  }
  for (LV in 1:num_lvs){
    for (i in 1:iterations){
      Xh_parted <- as.matrix(as.data.frame(Xh) %>% sample_frac(size=p))
      Yh_parted <- Yh[rownames(Xh_parted),]
      Xh_val <- Xh[which(!(rownames(Xh) %in% rownames(Xh_parted))),]
      Yh_val <- Yh[rownames(Xh_val),]
      #### First run human PLSR
      plsr_model <- suppressMessages(opls(x = Xh_parted, 
                                          y = Yh_parted,
                                          predI = LV,
                                          crossvalI = 1,
                                          scaleC = "center",
                                          fig.pdfC = "none",
                                          info.txtC = "none"))
      y_hat_plsr <- predict(plsr_model,Xh_parted)
      y_hat_val_plsr <- predict(plsr_model,Xh_val)
      
      rho_plsr <- diag(cor(y_hat_plsr,Yh_parted,method = "spearman"))
      performance_train <- rbind(performance_train,
                                 data.frame(model = 'PLSR',
                                            num_LVs = LV,
                                            partition = p,
                                            sample=i,
                                            NAS = rho_plsr['NAS'],
                                            fibrosis = rho_plsr['fibrosis']))
      
      rho_plsr <- diag(cor(y_hat_val_plsr,Yh_val,method = "spearman"))
      performance_val <- rbind(performance_val,
                               data.frame(model = 'PLSR',
                                          num_LVs = LV,
                                          partition = p,
                                          sample=i,
                                          NAS = rho_plsr['NAS'],
                                          fibrosis = rho_plsr['fibrosis']))
      
    }
    print(paste0('Finished PLSR with ',LV,' latent variables'))
  }
}
# saveRDS(performance_val,'../results/performance_val_random_partitions_many_lvs_spearman.rds')
# saveRDS(performance_train,'../results/performance_train_random_partitions_many_lvs_spearman.rds')

val_plot <- performance_val %>% gather('phenotype','rho',-model,-num_LVs,-partition,-sample)  %>% 
  group_by(model,num_LVs,partition,phenotype) %>% mutate(mu = mean(rho)) %>% 
  mutate(std = sd(rho)) %>% ungroup() %>%
  unique() 
train_plot <- performance_train %>% gather('phenotype','rho',-model,-num_LVs,-partition,-sample)  %>% 
  group_by(model,num_LVs,partition,phenotype) %>% mutate(mu = mean(rho)) %>% 
  mutate(std = sd(rho)) %>% ungroup() %>%
  unique()
performance_plot <- rbind(train_plot %>% mutate(set='train'),
                          val_plot %>% mutate(set='test'))
performance_plot$partition <- factor(performance_plot$partition,levels = partitions)
performance_plot <- performance_plot %>% mutate(partition = paste0('train size = ',partition,'%'))

ggplot(performance_plot,
       aes(x=num_LVs,y=rho,colour=set,shape=phenotype))+
  geom_smooth(aes(linetype=phenotype),lwd=1)+ 
  # geom_point()+ 
  # geom_errorbar(aes(ymax = mu + std/sqrt(iterations),ymin=mu - std/sqrt(iterations)))+
  scale_x_continuous(breaks = seq(1,20,2))+
  scale_y_continuous(breaks = seq(0.5,1,0.1))+
  xlab('number of latent variables') +
  ylab('spearman`s rank correlation coefficient')+
  theme_pubr(base_family = 'Arial',base_size = 18)+
  theme(text = element_text(size=18,family = 'Arial'),
        legend.position = 'right',
        panel.grid.major = element_line(),
        strip.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~partition,scale='free_x')

ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr_multiple_partitions.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)

ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr_multiple_partitions.eps',
       device = cairo_ps,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)

### After selecting 8 LVs as the tuned parameter re-fit with only those
### But also calculate shuffled and null models performance
val_mae <- NULL
train_mae <- NULL
test_mae <- NULL
val_r <- NULL
train_r <- NULL
test_r <- NULL
train_rho <- NULL
val_rho <- NULL
test_rho <- NULL
num_folds <- 10
tuning_df <- data.frame()
all_models <- NULL
test_r_shuffle_y <- NULL
test_rho_shuffle_y <- NULL
test_mae_shuffle_y <- NULL
test_r_shuffle_x <- NULL
test_mae_shuffle_x <- NULL
test_rho_shuffle_x <- NULL
test_r_random_x <- NULL
test_rho_random_x <- NULL
test_mae_random_x <- NULL
for (j in 1:num_folds){
  message(paste0('Begun fold ',j))
  x_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_train',j,'.rds'))
  y_train <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_train',j,'.rds'))
  x_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_val',j,'.rds'))
  y_val <- readRDS(paste0('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_val',j,'.rds'))
  
  plsr_model <- opls(x = x_train, 
                     y = y_train,
                     predI = 8,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  y_train_hat <- predict(plsr_model,x_train)
  y_val_hat <- predict(plsr_model,x_val)
  y_hat_test <- predict(plsr_model,Xh_test)
  
  train_rho[j] <- mean(diag(cor(y_train_hat,y_train,method = 'spearman')))
  val_rho[j] <- mean(diag(cor(y_val_hat,y_val,method = 'spearman')))
  test_rho[j]<- mean(diag(cor(y_hat_test,Yh_test,method = 'spearman')))
  print(paste0('NAS = ',diag(cor(y_hat_test,Yh_test,method = 'spearman'))['NAS'],' , Fibrosis = ',diag(cor(y_hat_test,Yh_test,method = 'spearman'))['fibrosis']))
  train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
  val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
  test_r[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  val_mae[j] <- mean(abs(y_val_hat-y_val))
  train_mae[j] <- mean(abs(y_train_hat-y_train))
  test_mae[j] <- mean(abs(y_hat_test-Yh_test))
  
  ### shuffled labels model
  y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
  rownames(y_train_shuffled) <- rownames(y_train)
  plsr_model_shuffle_y <- opls(x = x_train, 
                               y = y_train_shuffled,
                               predI = 8,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  y_hat_test <- predict(plsr_model_shuffle_y,Xh_test)
  test_r_shuffle_y[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  test_rho_shuffle_y[j] <- mean(diag(cor(y_hat_test,Yh_test,method = 'spearman')))
  test_mae_shuffle_y[j] <- mean(abs(y_hat_test-Yh_test))
  
  print('Finished shuffled labels model')
  
  ### shuffled features model
  x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
  colnames(x_train_shuffled) <- colnames(x_train)
  plsr_model_shuffle_x <- opls(x = x_train_shuffled, 
                               y = y_train,
                               predI = 8,
                               crossvalI = 1,
                               scaleC = "center",
                               fig.pdfC = "none",
                               info.txtC = "none")
  y_hat_test <- predict(plsr_model_shuffle_x,Xh_test)
  test_r_shuffle_x[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  test_rho_shuffle_x[j] <- mean(diag(cor(y_hat_test,Yh_test,method = 'spearman')))
  test_mae_shuffle_x[j] <- mean(abs(y_hat_test-Yh_test))
  
  print('Finished shuffled features model')
  
  ### random features model
  x_train_random <- matrix(rnorm(n = nrow(x_train)*ncol(x_train)),nrow = nrow(x_train))
  colnames(x_train_random) <- colnames(x_train)
  rownames(x_train_random) <- rownames(x_train)
  plsr_model_random_x <- opls(x = x_train_random, 
                              y = y_train,
                              predI = 8,
                              crossvalI = 1,
                              scaleC = "center",
                              fig.pdfC = "none",
                              info.txtC = "none")
  y_hat_test <- predict(plsr_model_random_x,Xh_test)
  test_r_random_x[j]<- mean(diag(cor(y_hat_test,Yh_test)))
  test_rho_random_x[j] <- mean(diag(cor(y_hat_test,Yh_test,method = 'spearman')))
  test_mae_random_x[j] <- mean(abs(y_hat_test-Yh_test))
  
  print('Finished random features model')
}
performance_df <- rbind(data.frame(set='train',rho=train_rho,r = train_r,MAE = train_mae,fold = seq(1,num_folds)),
                        data.frame(set='validation',rho=val_rho,r = val_r,MAE = val_mae,fold = seq(1,num_folds)),
                        data.frame(set='test',rho=test_rho,r = test_r,MAE = test_mae,fold = seq(1,num_folds)),
                        data.frame(set='shuffle Y',rho = test_rho_shuffle_y,r = test_r_shuffle_y,MAE = test_mae_shuffle_y,fold = seq(1,num_folds)),
                        data.frame(set='shuffle X',rho = test_rho_shuffle_x,r = test_r_shuffle_x,MAE = test_mae_shuffle_x,fold = seq(1,num_folds)),
                        data.frame(set='random X',rho = test_rho_random_x,r = test_r_random_x,MAE = test_mae_random_x,fold = seq(1,num_folds)))
# saveRDS(performance_df,'../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned_spearman.rds')

avg_mae <- mean(abs(Yh_test-apply(Yh,2,mean)))
plotting_performance_df <- performance_df %>% 
  gather('metric','value',-set,-fold) %>%
  mutate(metric=ifelse(metric=='r','pearson`s r',ifelse(metric=='rho','spearman`s rank correlation',metric)))
plotting_performance_df$set <- factor(plotting_performance_df$set,
                                      levels = c('train','validation','test','shuffle Y','shuffle X','random X'))
p <- (ggboxplot(plotting_performance_df %>% filter(metric=='MAE'),x='set',y='value',color = 'set',add='jitter')+
    scale_y_continuous(n.breaks = 15)+
    geom_hline(yintercept = avg_mae,linetype='dashed',color='black',lwd=1)+
    annotate('text',x=2,y=1.95,size=6,label = 'error from the mean of the data')+
    xlab('')+
    theme(text = element_text(size=20,family = 'Arial'),
          legend.position = 'none',
          axis.text.x = element_text(size=20),
          strip.text = element_text(face = 'bold'),
          panel.grid.major.y = element_line(linewidth = 1))+
    stat_compare_means(comparisons = list(c('validation','test'),
                                          c('test','shuffle Y'),
                                          c('test','shuffle X'),
                                          c('test','random X')),
                       method = 'wilcox',label = 'p.signif',
                       tip.length = 0.01,
                       label.y = c(1,1.7,1.8,1.9)) +
    facet_wrap(~metric)) /
  (ggboxplot(plotting_performance_df %>% filter(metric!='MAE'),x='set',y='value',color = 'set',add='jitter')+
     scale_y_continuous(breaks = seq(-0.5,1,0.1))+
     # geom_hline(yintercept = 0,linetype='dashed',color='black',lwd=1)+
     xlab('')+
     theme(text = element_text(size=20,family = 'Arial'),
           legend.position = 'none',
           axis.text.x = element_text(size=20),
           strip.text = element_text(face = 'bold'),
           panel.grid.major.y = element_line(linewidth = 1))+
     stat_compare_means(comparisons = list(c('validation','test'),
                                           c('test','shuffle Y'),
                                           c('test','shuffle X'),
                                           c('test','random X')),
                        method = 'wilcox',label = 'p.signif',
                        tip.length = 0.01,
                        label.y = c(0.75,0.75,0.85,0.95)) +
     facet_wrap(~metric))
print(p)
ggsave('../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned.png',
       plot = p,
       width = 15,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned.eps',
       plot = p,
       device = cairo_ps,
       width = 15,
       height = 9,
       units = 'in',
       dpi = 600)
### Since we do not overfit to validation now perform cross-validation with this dataset-------------
### Where you calculate validation and train performance of:
### a) human data PLSR to predict NAS,Fibrosis
### b) truncated human data PLSR to predict NAS, Fibrosis
### c) truncated human data PLSR to predict NAS, Fibrosis with translatable combination of MPS PCs
### d) c) truncated human data PLSR to predict NAS, Fibrosis with optimal direction
num_folds <- 10
loc <- '../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/'
num_LVS <- 8
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
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
    hist_val_nas <- hist(val_Y[,1],freq = T,breaks = seq(0,8,0.25),plot = F)
    hist_train_nas <- hist(Yh[ind,1],freq = T,breaks = seq(0,8,0.25),plot = F)
    hist_val_fibr <- hist(val_Y[,2],freq = T,breaks = seq(0,8,0.25),plot = F)
    hist_train_fibr <- hist(Yh[ind,2],freq = T,breaks = seq(0,8,0.25),plot = F)
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
# saveRDS(performance_1,'../results/performance_df_human_plsr_spearman.rds')


# Task b)
performance_2 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected',
                                                    LVs = num_LVS)
# saveRDS(performance_2,'../results/performance_df_human_backprojected_spearman.rds')
performance_2$type <- factor(performance_2$type ,levels=c('model','shuffle W','shuffle X','shuffle Y','shuffle Bh'))
performance_2$set <- factor(performance_2$set ,levels=c('train','test'))
performance_2 <- performance_2 %>% filter(metric=='rho') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(rho = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_2,x='type',y='rho',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+ylab('spearman`s rank correlation')+
  theme(text = element_text(size=22,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=18),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     step.increase = 0.04,
                     size=6)

# Task c)
performance_3 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_retrained',
                                                    LVs = num_LVS)
saveRDS(performance_3,'../results/performance_df_human_backprojected_retrained_spearman.rds')
#Plot
performance_3$type <- factor(performance_3$type ,levels=c('model','shuffle W'))
performance_3$set <- factor(performance_3$set ,levels=c('train','test'))
performance_3 <- performance_3 %>% filter(metric=='rho') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(rho = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_3,x='type',y='rho',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+ylab('spearman`s rank correlation')+
  theme(text = element_text(size=22,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=18),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     size=6)

# Task d1)
performance_4 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'human_backprojected_into_translatable_lvs',
                                                    LVs = num_LVS)
# saveRDS(performance_4,'../results/performance_df_translatable_lvs_spearman.rds')
#Plot
performance_4$type <- factor(performance_4$type ,levels=c('model','shuffle W','shuffle Bh'))
performance_4$set <- factor(performance_4$set ,levels=c('train','test'))
performance_4 <- performance_4 %>% filter(metric=='rho') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(rho = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_4,x='type',y='rho',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+ylab('spearman`s rank correlation')+
  theme(text = element_text(size=22,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=18),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle W'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     step.increase = 0.04,
                     size=6)

performance_6 <- cross_validation_complete_pipeline(Wm,
                                                    folds = num_folds,
                                                    file_loc = loc ,
                                                    target_dataset = target_dataset,
                                                    task = 'analytical_optimal',
                                                    LVs = num_LVS)
saveRDS(performance_6,'../results/performance_df_analytical_spearman.rds')

performance_6$type <- factor(performance_6$type ,levels=c('model','shuffle Wopt','shuffle Bh'))
performance_6$set <- factor(performance_6$set ,levels=c('train','test'))
performance_6 <- performance_6 %>% filter(metric=='rho') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(rho = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_6,x='type',y='rho',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(round(min(performance_6$rho),1),1,0.1),
                     limits = c(NA,1.15))+
  xlab('')+ylab('pspearman`s rank correlation')+
  theme(text = element_text(size=22,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=18),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle Wopt'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',
                     tip.length = 0.01,size=6)# ,label.y = c(0.82, 0.87, 0.92)

ggsave('../results/performance_df_analytical_spearman.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)



### Combine all results and make a nice figure--------------------------------------------------------------------
performance_1 <- readRDS('../results/performance_df_human_plsr_spearman.rds')
performance_2 <- readRDS('../results/performance_df_human_backprojected_spearman.rds')
performance_3 <- readRDS('../results/performance_df_human_backprojected_retrained_spearman.rds')
performance_4 <- readRDS('../results/performance_df_translatable_lvs_spearman.rds')
performance_6 <- readRDS('../results/performance_df_analytical_spearman.rds')

performance_all <- rbind(performance_1,
                         performance_2,
                         performance_3,
                         performance_4,
                         performance_6)

### Select what to show in the figure
performance_all_plot <- performance_all %>%  filter(metric=='rho') %>% select(-metric) %>% mutate(rho=value) %>% select(-value)%>%
  mutate(keep=ifelse(task=='human_plsr',ifelse(type %in% c('model','shuffle X'),TRUE,FALSE),
                     ifelse(type=='model',TRUE,FALSE))) %>% 
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(approach = ifelse(task=='human_plsr','human genes',
                           ifelse(task=='human_backprojected','truncated',
                                  ifelse(task=='human_backprojected_retrained','truncated re-trained',
                                         ifelse(task=='human_backprojected_into_translatable_lvs','translatable LVs',
                                                ifelse(task=='analytical_optimal','optimized MPS',
                                                       'Wopt')))))) %>%
  mutate(approach = ifelse(approach=='human genes',
                           ifelse(type=='model',approach,type),
                           approach)) %>%
  select(-type,-task)

performance_all_truncated <- performance_all_plot %>%
  filter(approach %in% c('truncated','truncated re-trained','human genes','shuffle X'))
performance_all_truncated$approach <- factor(performance_all_truncated$approach,
                                        levels = c('human genes',
                                                   'truncated',
                                                   'truncated re-trained',
                                                   'shuffle X'))
performance_all_truncated <- performance_all_truncated %>% mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage','MAS'))
performance_all_truncated$phenotype <- factor(performance_all_truncated$phenotype)
performance_all_truncated$phenotype <- factor(performance_all_truncated$phenotype,
                                         levels = rev(levels(performance_all_truncated$phenotype)))


performance_all_plot <- performance_all_plot %>%
  filter(approach %in% c('human genes','truncated','optimized MPS','shuffle X'))
performance_all_plot$approach <- factor(performance_all_plot$approach,
                                        levels = c('human genes',
                                                   'truncated',
                                                   'optimized MPS',
                                                   'shuffle X'))
performance_all_plot <- performance_all_plot %>% mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
# performance_all_plot$phenotype <- factor(performance_all_plot$phenotype,levels = c('NAS','Fibrosis stage'))
performance_all_plot$phenotype <- factor(performance_all_plot$phenotype)
performance_all_plot$phenotype <- factor(performance_all_plot$phenotype,
                                         levels = rev(levels(performance_all_plot$phenotype)))
# Save data frame
saveRDS(performance_all_plot,'../results/performanceall_plot_spearman.rds')

p_train <- ggboxplot(performance_all_plot %>% filter(set=='train') %>% filter(approach!='Wopt'),
                     x='approach',y='rho',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(0.4,1,0.05),limits = c(NA,1.02))+
  xlab('')+ ylab('speaman`s rank correlation') +
  # ggtitle('10-fold train performance in predicting phenotype')+
  theme(text = element_text(size=26,family = 'Arial'),
        legend.position = 'none',
        # plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=26), #angle = 25
        # strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
                                        c('human genes','truncated'),
                                        c('optimized MPS','truncated')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     label.y = c(0.98,0.94,0.94),
                     size = 6)+
  facet_wrap(~phenotype,nrow = 2)
print(p_train)  
ggsave('../figures/approaches_comparison_training.png',
       plot = p_train,
       height = 8,
       width = 10,
       units = 'in',
       dpi=600)

p_test <- ggboxplot(performance_all_plot %>% filter(set=='test')%>% filter(approach!='Wopt'),
                    x='approach',y='rho',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(-0.5,1,0.15),limits = c(NA,1.05))+
  xlab('')+ylab('spearman`s rank correlation') +
  # ggtitle('10-fold test performance in predicting phenotype')+
  theme(text = element_text(size=26,family = 'Arial'),
        legend.position = 'none',
        axis.text.y = element_text(size=22),
        # plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=24), #angle = 25
        # strip.text = element_text(face = 'bold'),
        panel.spacing = unit(0, "lines"),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
                                        c('human genes','truncated'),
                                        c('optimized MPS','truncated'),
                                        c('optimized MPS','shuffle X')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     label.y = c(0.9,0.8,0.8,0.8),
                     size = 8)+
  facet_wrap(~phenotype,nrow = 2)
print(p_test)  
ggsave('../figures/approaches_comparison_10foldtest.png',
       plot = p_test,
       height = 9,
       width = 11,
       units = 'in',
       dpi=600)
ggsave('../figures/approaches_comparison_10foldtest.eps',
       plot = p_test,
       device = cairo_ps,
       height = 9,
       width = 11,
       units = 'in',
       dpi=600)

### Compare the truncated versions only
p_train_truncated <- ggboxplot(performance_all_truncated %>% filter(set=='train') %>%
                                 filter(approach!='Wopt'),
                     x='approach',y='rho',color='approach',add='jitter') +
  scale_y_continuous(n.breaks = 10,limits = c(NA,NA))+
  xlab('')+ ylab('spearman`s rank correlation') +
  theme(text = element_text(size=26,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=26), #angle = 25
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('truncated','truncated re-trained'),
                                        c('truncated','human genes'),
                                        c('human genes','truncated re-trained')),
                     method = 'wilcox',
                     paired = TRUE,
                     tip.length = 0.01,
                     label.y = c(0.7,0.95,1),
                     size = 6)+
  facet_wrap(~phenotype,nrow = 2)
print(p_train_truncated)  
ggsave('../figures/truncated_all_genes_comparison_training.png',
       plot = p_train_truncated,
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)
ggsave('../figures/truncated_all_genes_comparison_training.eps',
       device = cairo_ps,
       plot = p_train_truncated,
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)

p_test_truncated <- ggboxplot(performance_all_truncated %>% filter(set=='test')%>% filter(approach!='Wopt'),
                    x='approach',y='rho',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(-0.4,1,0.2),limits = c(NA,1.05))+
  xlab('')+ylab('spearman`s rank correlation') +
  # ggtitle('10-fold test performance in predicting phenotype')+
  theme(text = element_text(size=26,family = 'Arial'),
        legend.position = 'none',
        axis.text.y = element_text(size=22),
        # plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=24), #angle = 25
        # strip.text = element_text(face = 'bold'),
        panel.spacing = unit(0, "lines"),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('truncated','truncated re-trained'),
                                        c('truncated','human genes'),
                                        c('truncated','shuffle X'),
                                        c('human genes','truncated re-trained')),
                     method = 'wilcox',
                     paired=TRUE,
                     tip.length = 0.01,
                     label.y = c(0.7,0.8,0.8,0.9),
                     size = 6)+
  facet_wrap(~phenotype,nrow = 2)
print(p_test_truncated)  
ggsave('../figures/truncated_all_genes_comparison_10foldtest.png',
       plot = p_test_truncated,
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)
ggsave('../figures/truncated_all_genes_comparison_10foldtest.eps',
       device = cairo_ps,
       plot = p_test_truncated,
       height = 12,
       width = 12,
       units = 'in',
       dpi=600)

### See how training performance converges
### when including incrementally the extra basis--------------------------------------------------------------------
train_rho <- NULL
train_rho_translatables <- NULL
train_rho_basis_1 <- NULL
train_rho_all_extra_bases <- NULL

test_rho <- NULL
test_rho_shuffled <- NULL
test_rho_translatables <- NULL
test_rho_basis_1 <- NULL
test_rho_all_extra_bases <- NULL

df_scatterPlot <- data.frame()
df_scatterPlot_backproj <- data.frame()

for (j in 1:num_folds){
  print(paste0('Begun fold ',j))
  x_train <- readRDS(paste0(loc,'Xh_train',j,'.rds'))
  y_train <- readRDS(paste0(loc,'Yh_train',j,'.rds'))
  x_val <- readRDS(paste0(loc,'Xh_val',j,'.rds'))
  y_val <- readRDS(paste0(loc,'Yh_val',j,'.rds'))
  
  plsr_model <- opls(x = x_train, 
                     y = y_train,
                     predI = num_LVS,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  y_train_hat <- predict(plsr_model,x_train)
  y_val_hat <- predict(plsr_model,x_val)
  train_rho[j] <- mean(diag(cor(y_train_hat,y_train,method='spearman')))
  test_rho[j] <- mean(diag(cor(y_val_hat,y_val,method='spearman')))
  
  df_scatterPlot <- rbind(df_scatterPlot,
                          left_join(data.frame(y_val) %>% 
                                      mutate(id = seq(1,nrow(y_val))) %>% 
                                      gather('phenotype','true',-id),
                                    data.frame(y_val_hat) %>% 
                                      mutate(id = seq(1,nrow(y_val_hat))) %>% 
                                      gather('phenotype','prediction',-id)) %>%
                            select(-id))

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
  
  # Define projection matrices to make more readable
  Th_train <- x_train %*% Wh
  Thm_train <- x_train %*% Wm %*% t(Wm) %*% Wh
  Th_val<- x_val %*% Wh
  Thm_val <- x_val %*% Wm %*% t(Wm) %*% Wh
  
  ## backprojecting
  y_hat_test_backproj <- cbind(1, Thm_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
  df_scatterPlot_backproj <- rbind(df_scatterPlot_backproj,
                                   left_join(data.frame(y_val) %>% 
                                               mutate(id = seq(1,nrow(y_val))) %>% 
                                               gather('phenotype','true',-id),
                                             data.frame(y_hat_test_backproj) %>% 
                                               mutate(id = seq(1,nrow(y_hat_test_backproj))) %>% 
                                               gather('phenotype','prediction',-id)) %>%
                                     select(-id))
  
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
  test_rho_shuffled[j]<- mean(diag(cor(y_hat_val,y_val,method='spearman')))
  
  message('Finished running initial PLSR for humans')
  
  ### Get linear combinaiton of translatable LVs
  Wm_new <- get_translatable_LV_2phenotype(x_train, y_train, Wh, Wm,Bh)
  Wm_new <- Wm_new$Wm_TC
  colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
  y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  train_rho_translatables[j] <- mean(diag(cor(y_hat_train,y_train,method='spearman')))
  y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  test_rho_translatables[j] <- mean(diag(cor(y_hat_test,y_val,method='spearman')))
  
  message('Finished linear combination of translatable components')
  
  ## Start analytical solution for optimal extra basis
  phi <- Wh %*% Bh
  # phi <- phi/sqrt(apply(phi^2,2,sum))
  Wm_opt <- analytical_solution_opt(y=y_train,
                                    W_invitro = Wm,
                                    phi = phi)
  # Extend latent variables
  Wm_tot_1 <- cbind(Wm, Wm_opt[,1])
  Wm_tot <- cbind(Wm, Wm_opt)
  # predict with one extra basis vectors
  y_hat_train <- cbind(1, x_train %*% Wm_tot_1 %*% t(Wm_tot_1) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  train_rho_basis_1[j] <- mean(diag(cor(y_hat_train,y_train,method='spearman')))
  y_hat_test <- cbind(1, x_val %*% Wm_tot_1 %*% t(Wm_tot_1) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  test_rho_basis_1[j] <- mean(diag(cor(y_hat_test,y_val,method='spearman')))
  message('Finished analytical solution performance with 1 extra vector')
  # predict with one extra basis vectors
  y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  train_rho_all_extra_bases[j] <- mean(diag(cor(y_hat_train,y_train,method='spearman')))
  y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  test_rho_all_extra_bases[j] <- mean(diag(cor(y_hat_test,y_val,method='spearman')))
  message('Finished analytical solution performance with 2 extra vectors')
}

## First plot scatter plots of prediction VS ground truth
# Calculate correlation coefficients
cor_results <- df_scatterPlot %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction,method='spearman'))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))
ggplot(df_scatterPlot, aes(x = true, y = prediction)) +
  geom_jitter(width = 0.05, color = '#4682B4') + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black', linewidth = 1.5) +
  geom_text(data = cor_results, 
            aes(x = -Inf, y = Inf, 
                label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
            hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
  facet_wrap(~phenotype, scales = 'free') +
  theme_pubr(base_family = 'Arial', base_size = 25) +
  theme(text = element_text(family = 'Arial', size = 25),
        panel.grid.major = element_line())
ggsave('../figures/10foldTest_Scatterplot_human_plsr.png',
          height = 6,
          width=9,
          units = 'in',
          dpi=600)

# repreat for backprojection
cor_resultst_backproj <- df_scatterPlot_backproj %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction,method = 'spearman'))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))
ggplot(df_scatterPlot_backproj,aes(x = true,y=prediction)) +
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
 geom_text(data = cor_resultst_backproj, 
            aes(x = -Inf, y = Inf, 
                label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
            hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
  facet_wrap(~phenotype,scales = 'free')+
  theme_pubr(base_family = 'Arial',base_size=25)+
  theme(text = element_text(family = 'Arial',size=25),
        panel.grid.major = element_line())
ggsave('../figures/10foldTest_Scatterplot_backprojecting.png',
       height = 6,
       width=9,
       units = 'in',
       dpi=600)

## scatter plot for using all the data
plsr_model <- opls(x = Xh, 
                   y = Yh,
                   predI = num_LVS,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
Y_pred <- predict(plsr_model,Xh)
all_scatter_plot <- left_join(data.frame(Yh) %>% 
                                mutate(id = seq(1,nrow(Yh))) %>% 
                                gather('phenotype','true',-id),
                              data.frame(Y_pred) %>% 
                                mutate(id = seq(1,nrow(Y_pred))) %>% 
                                gather('phenotype','prediction',-id)) %>%
  select(-id)%>% 
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
all_cor_results <- all_scatter_plot %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction,method = 'spearman'))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
ggplot(all_scatter_plot,aes(x = true,y=prediction)) +
  ylab('Predicted') + xlab('Measured')+
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = all_cor_results, 
            aes(x = -Inf, y = Inf, 
                label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
            hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
  facet_wrap(~phenotype,scales = 'free')+
  theme_pubr(base_family = 'Arial',base_size=28)+
  theme(text = element_text(family = 'Arial',size=28),
        panel.grid.major = element_line())
ggsave('../figures/AllData_Scatterplot_human_plsr.png',
       height = 6,
       width=9,
       units = 'in',
       dpi=600)
ggplot(all_scatter_plot,aes(x = true,y=prediction)) +
  ylab('Predicted') + xlab('Measured')+
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = all_cor_results, 
            aes(x = -Inf, y = Inf, 
                label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
            hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
  facet_wrap(~phenotype,scales = 'free')+
  theme_pubr(base_family = 'Arial',base_size=28)+
  theme(text = element_text(family = 'Arial',size=28),
        plot.title = element_text(hjust = 0.5,face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face='bold'),
        panel.grid.major = element_line())
ggsave('../figures/AllData_Scatterplot_human_plsr.eps',
       device = cairo_ps,
       height = 6,
       width=12,
       units = 'in',
       dpi=600)

# repeat for backprojection
# Get Wh of PLSR
Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
rownames(Wh) <- colnames(Xh)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)){
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
}
Th <- Xh %*% Wm %*% t(Wm) %*% Wh
Y_pred_backproj <- cbind(1, Th)  %*% rbind(apply(Yh,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
Y_pred_backproj_retrained <- matrix(0,nrow = nrow(Yh),ncol = ncol(Yh))
for (i in 1:ncol(Yh)){
  train_data <- cbind(Yh[,i], Th)
  colnames(train_data)[1] <- 'V1'
  mod_retrain <- lm(V1 ~., data = as.data.frame(train_data))
  Y_pred_backproj_retrained[,i] <- predict(mod_retrain)
}
colnames(Y_pred_backproj_retrained) <- colnames(Y_pred_backproj)
all_scatter_plot_backproj <- left_join(data.frame(Yh) %>% 
                                mutate(id = seq(1,nrow(Yh))) %>% 
                                gather('phenotype','true',-id),
                              data.frame(Y_pred_backproj) %>% 
                                mutate(id = seq(1,nrow(Y_pred_backproj))) %>% 
                                gather('phenotype','prediction',-id)) %>%
  select(-id) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
all_cor_results_backproj <- all_scatter_plot_backproj %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction,method = 'spearman'))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
(ggplot(all_scatter_plot_backproj %>% filter(phenotype=='Fibrosis stage'),aes(x = true,y=prediction)) +
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = all_cor_results_backproj%>% filter(phenotype=='Fibrosis stage'), 
              aes(x = -Inf, y = Inf, 
                  label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
              hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
  ylab('Predicted') + xlab('Measured')+
  ylim(c(0,4))+
  facet_wrap(~phenotype,scales = 'free')+
  theme_pubr(base_family = 'Arial',base_size=25)+
  theme(text = element_text(family = 'Arial',size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face='bold'),
        panel.grid.major = element_line()))+
  (ggplot(all_scatter_plot_backproj %>% filter(phenotype!='Fibrosis stage'),aes(x = true,y=prediction)) +
     geom_jitter(width = 0.05,color='#4682B4') + 
     geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
     geom_text(data = all_cor_results_backproj%>% filter(phenotype!='Fibrosis stage'), 
               aes(x = -Inf, y = Inf, 
                   label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
               hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
     ylab('Predicted') + xlab('Measured')+
     ylim(c(0,8))+
     facet_wrap(~phenotype,scales = 'free')+
     theme_pubr(base_family = 'Arial',base_size=25)+
     theme(text = element_text(family = 'Arial',size=25),
           axis.title.x = element_text(hjust = -0.6,face='bold'),
           axis.title.y = element_blank(),
           panel.grid.major = element_line())) + 
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(size = 25, family = "Arial", hjust = 0.5,face='bold'))
  )
ggsave('../figures/AllData_backproj_Scatterplot_human_plsr.png',
       height = 6,
       width=9,
       units = 'in',
       dpi=600)
(ggplot(all_scatter_plot_backproj %>% filter(phenotype=='Fibrosis stage'),aes(x = true,y=prediction)) +
    geom_jitter(width = 0.05,color='#4682B4') + 
    geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
    geom_text(data = all_cor_results_backproj%>% filter(phenotype=='Fibrosis stage'), 
              aes(x = -Inf, y = Inf, 
                  label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
              hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
    ylab('Predicted') + xlab('Measured')+
    ylim(c(0,4))+
    facet_wrap(~phenotype,scales = 'free')+
    theme_pubr(base_family = 'Arial',base_size=28)+
    theme(text = element_text(family = 'Arial',size=28),
          plot.title = element_text(hjust = 0.5,face = 'bold'),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face='bold'),
          # axis.title.y = element_blank(),
          panel.grid.major = element_line()))+
  (ggplot(all_scatter_plot_backproj %>% filter(phenotype!='Fibrosis stage'),aes(x = true,y=prediction)) +
     geom_jitter(width = 0.05,color='#4682B4') + 
     geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
     geom_text(data = all_cor_results_backproj%>% filter(phenotype!='Fibrosis stage'), 
               aes(x = -Inf, y = Inf, 
                   label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
               hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
     ylab('Predicted') + xlab('Measured')+
     ylim(c(0,8))+
     facet_wrap(~phenotype,scales = 'free')+
     theme_pubr(base_family = 'Arial',base_size=28)+
     theme(text = element_text(family = 'Arial',size=28),
           # axis.title.x = element_text(hjust = -0.4,face='bold'),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.grid.major = element_line())) + 
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(size = 28, family = "Arial", hjust = 0.5,face='bold'))
  )
ggsave('../figures/AllData_backproj_Scatterplot_human_plsr.eps',
       device = cairo_ps,
       height = 6,
       width=12,
       units = 'in',
       dpi=600)

### plot scatterplot when re-training Bh
all_scatter_plot_backproj_retrained <- left_join(data.frame(Yh) %>% 
                                         mutate(id = seq(1,nrow(Yh))) %>% 
                                         gather('phenotype','true',-id),
                                       data.frame(Y_pred_backproj_retrained) %>% 
                                         mutate(id = seq(1,nrow(Y_pred_backproj_retrained))) %>% 
                                         gather('phenotype','prediction',-id)) %>%
  select(-id) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
all_cor_results_backproj_retrained <- all_scatter_plot_backproj_retrained %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction,method='spearman'))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
(ggplot(all_scatter_plot_backproj_retrained %>% filter(phenotype=='Fibrosis stage'),aes(x = true,y=prediction)) +
    geom_jitter(width = 0.05,color='#4682B4') + 
    geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
    geom_text(data = all_cor_results_backproj_retrained%>% filter(phenotype=='Fibrosis stage'), 
              aes(x = -Inf, y = Inf, 
                  label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
              hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
    ylab('Predicted') + xlab('Measured')+
    ylim(c(0,4))+
    facet_wrap(~phenotype,scales = 'free')+
    theme_pubr(base_family = 'Arial',base_size=25)+
    theme(text = element_text(family = 'Arial',size=25),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face='bold'),
          panel.grid.major = element_line()))+
  (ggplot(all_scatter_plot_backproj_retrained %>% filter(phenotype!='Fibrosis stage'),aes(x = true,y=prediction)) +
     geom_jitter(width = 0.05,color='#4682B4') + 
     geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
     geom_text(data = all_cor_results_backproj_retrained%>% filter(phenotype!='Fibrosis stage'), 
               aes(x = -Inf, y = Inf, 
                   label = sprintf("rho~'='~%.2f*','~p~'='~%.2g", cor_coef, p_value)),
               hjust = -0.1, vjust = 1.5, size = 8, family = 'Arial', parse = TRUE) +
     ylab('Predicted') + xlab('Measured')+
     ylim(c(0,8))+
     facet_wrap(~phenotype,scales = 'free')+
     theme_pubr(base_family = 'Arial',base_size=25)+
     theme(text = element_text(family = 'Arial',size=25),
           axis.title.x = element_text(hjust = -0.6,face='bold'),
           axis.title.y = element_blank(),
           panel.grid.major = element_line())) + 
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(size = 25, family = "Arial", hjust = 0.5,face='bold'))
  )
ggsave('../figures/AllData_backproj_retrained_Scatterplot_human_plsr.png',
       height = 6,
       width=9,
       units = 'in',
       dpi=600)
### Now see convergence of performance
train_performance_res <- rbind(data.frame(rho = train_rho,
                                          fold = seq(1,length(train_rho)),
                                          model = rep('PLSR',length(train_rho))),
                               data.frame(rho = train_rho_translatables,
                                          fold = seq(1,length(train_rho_translatables)),
                                          model = rep('translatable PCs',length(train_rho_translatables))),
                               data.frame(rho = train_rho_basis_1,
                                          fold = seq(1,length(train_rho_basis_1)),
                                          model = rep('PCs + extra LV1',length(train_rho_basis_1))),
                               data.frame(rho = train_rho_all_extra_bases,
                                          fold = seq(1,length(train_rho_all_extra_bases)),
                                          model = rep('PCs + extra LV1 + extra LV2',length(train_rho_all_extra_bases))))
train_performance_res <- train_performance_res %>% group_by(model) %>% 
  mutate(mu = mean(rho)) %>% mutate(std = sd(rho)) %>%
  mutate(min_rho = min(rho)) %>% mutate(max_rho = max(rho))%>%
  ungroup() %>% mutate(set='train')
mu_plsr_train <- mean(train_rho)
sd_plsr_train <- sd(train_rho)
test_performance_res <- rbind(data.frame(rho = test_rho,
                                          fold = seq(1,length(test_rho)),
                                          model = rep('PLSR',length(test_rho))),
                               data.frame(rho = test_rho_translatables,
                                          fold = seq(1,length(test_rho_translatables)),
                                          model = rep('translatable PCs',length(test_rho_translatables))),
                               data.frame(rho = test_rho_basis_1,
                                          fold = seq(1,length(test_rho_basis_1)),
                                          model = rep('PCs + extra LV1',length(test_rho_basis_1))),
                               data.frame(rho = test_rho_all_extra_bases,
                                          fold = seq(1,length(test_rho_all_extra_bases)),
                                          model = rep('PCs + extra LV1 + extra LV2',length(test_rho_all_extra_bases))))
test_performance_res <- test_performance_res %>% group_by(model) %>% 
  mutate(mu = mean(rho)) %>% mutate(std = sd(rho)) %>%
  mutate(min_rho = min(rho)) %>% mutate(max_rho = max(rho))%>%
  ungroup() %>% mutate(set='test')
mu_plsr_test <- mean(test_rho)
sd_plsr_test <- sd(test_rho)
all_performance_res <- rbind(train_performance_res,
                             test_performance_res)
all_performance_res$model <- factor(all_performance_res$model,
                                    levels = c('PLSR',
                                               'translatable PCs',
                                               'PCs + extra LV1',
                                               'PCs + extra LV1 + extra LV2'))

# Save
saveRDS(all_performance_res,'../results/all_performance_res_spearman.rds')
saveRDS(test_rho_shuffled,'../results/test_r_shuffled_cumulative_lvs_spearman.rds')

plt_required_extra_basis <- 
  ggplot(all_performance_res %>% filter(model!='PLSR') %>% select(model,set,mu,std),
       aes(x=model,y=mu,color = set , group =set))+
  geom_point(size=1.5)+
  geom_line(lwd = 0.75)+
  geom_errorbar(aes(ymax = mu + std/sqrt(num_folds),ymin = mu - std/sqrt(num_folds)),
                width = 0.05,size=0.75)+
  ylim(NA,1) +
  ylab('Spearman`s rank correlation')+
  ### train performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mu_plsr_train - sd_plsr_train/sqrt(num_folds), 
                  ymax = mu_plsr_train + sd_plsr_train/sqrt(num_folds)), 
              fill = "#01B8BB", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mu_plsr_train, 
           yend = mu_plsr_train, 
           color = "#01B8BB", linewidth = 1,linetype='dashed') +
  annotate('text',x=1,
           y=mu_plsr_train + 0.03,
           label="human PLSR train performance", 
           hjust = 0 , 
           size=size_annotation)+
  ### test performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mu_plsr_test - sd_plsr_test/sqrt(num_folds), 
                  ymax = mu_plsr_test + sd_plsr_test/sqrt(num_folds)), 
              fill = "#E0766D", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mu_plsr_test, 
           yend = mu_plsr_test, 
           color = "#E0766D", linewidth = 1,linetype='dashed') +
  annotate('text',x=1,
           y=mu_plsr_test + 0.03,
           label="human PLSR test performance", 
           hjust = 0 , 
           size=size_annotation)+
  ### random performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mean(test_rho_shuffled) - sd(test_rho_shuffled)/sqrt(num_folds), 
                  ymax = mean(test_rho_shuffled) + sd(test_rho_shuffled)/sqrt(num_folds)), 
              fill = "#F564E3", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mean(test_rho_shuffled), 
           yend = mean(test_rho_shuffled), 
           color = "#F564E3", linewidth = 1,linetype='dashed') +
  annotate('text',x=1,
           y=mean(test_rho_shuffled) + 0.03,
           label="shuffled model performance", 
           hjust = 0 , 
           size=size_annotation)+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line())

plt_required_extra_basis <- add_theme(plt_required_extra_basis)

ggsave('../figures/required_extra_basis_to_retrieve_performance.png',
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)
ggsave('../figures/required_extra_basis_to_retrieve_performance.eps',
       device = cairo_ps,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)

### Explore how analytical LVs change by: ----------------------
### a) Human data partition
num_folds <- 10
loc <- '../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/'
num_LVS <- 8
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Govaere$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
rownames(Yh) <- rownames(Xh)
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()
# a) Human data partition
# Begin iterations
partitions <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2)
iterations <- 10
cos_sim_all <- data.frame()
Wm_opt_all <- NULL
for (p in partitions){
  message(paste0('Begun partition:',100*p,'%'))
  for (i in 1:iterations){
    Xh_parted <- as.matrix(as.data.frame(Xh) %>% sample_frac(size=p))
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
    
    ### Run analytical approach
    phi <- Wh %*% Bh
    Wm_opt <- analytical_solution_opt(y=Yh_parted,
                                      W_invitro = Wm,
                                      phi = phi)
    colnames(Wm_opt) <- paste0('LVstar',seq(1,ncol(Wm_opt)))
    # Wm_opt_all <- rbind(Wm_opt_all,
                        # as.data.frame(Wm_opt) %>% mutate(partition=p) %>% mutate(iteration=i))
    # Extend latent variables
    Wm_tot <- cbind(Wm, Wm_opt)
    
    # Calculate cosine similarity
    cos_sim <- lsa::cosine(Wm_tot)
    cos_sim <- cos_sim[colnames(Wm_opt),]
    cos_sim <- as.data.frame(cos_sim) %>% rownames_to_column('extra_basis') %>% gather('basis','sim',-extra_basis)
    cos_sim_all <- rbind(cos_sim_all,
                         cos_sim %>% mutate(partition = p) %>% mutate(iteration=i))
    colnames(Wm_opt) <- paste0(colnames(Wm_opt),'_part',p,'_iter',i)
    Wm_opt_all <- cbind(Wm_opt_all,Wm_opt)
    print(paste0('Finished iteration ',i,' / ',iterations))
  }
  
  cat(paste0('Average absolute cosine of found LVs for partition ',100*p,'%',
             '\nLV1-LV2 : ',mean(as.matrix(cos_sim_all %>% filter(partition==p) %>% 
                                           filter(extra_basis=='LVstar1') %>% filter(basis=='LVstar2') %>% 
                                           unique() %>% select(sim) %>% mutate(sim = abs(sim)))),
             '\nLV1 to all PCs : ',mean(as.matrix(cos_sim_all %>% filter(partition==p) %>% 
                                                  filter(extra_basis=='LVstar1') %>% 
                                                  filter(!(basis %in% c('LVstar1','LVstar2'))) %>% 
                                                  unique() %>% select(sim) %>% mutate(sim = abs(sim)))),
             '\nLV2 to all PCs : ',mean(as.matrix(cos_sim_all %>% filter(partition==p) %>% 
                                                  filter(extra_basis=='LVstar2') %>% 
                                                  filter(!(basis %in% c('LVstar1','LVstar2'))) %>% 
                                                  unique() %>% select(sim) %>% mutate(sim = abs(sim)))),
             '\n'))
}
### Infer pathway activity with progeny
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
rownames(Wm_opt_all) <- rownames(Wm)
WPaths_opt_all <- decoupleR::run_viper(Wm_opt_all, net_prog,minsize = 1,verbose = TRUE)
WPaths_opt_all <- WPaths_opt_all %>% dplyr::select(source,condition,score) %>%
  spread(condition,score)
WPaths_opt_all <- as.matrix(WPaths_opt_all %>% column_to_rownames('source'))
# saveRDS(Wm_opt_all,'../results/Wm_opt_all_different_partitions.rds')
# saveRDS(WPaths_opt_all,'../results/WPaths_opt_all_different_partitions.rds')
# saveRDS(cos_sim_all,'../results/cosine_sim_all_different_partitions.rds')

# Wm_opt_all <- readRDS('../results/Wm_opt_all_different_partitions.rds')
# cos_sim_all <- readRDS('../results/cosine_sim_all_different_partitions.rds')

cos_optimal <- lsa::cosine(Wm_opt_all)
dend <- hclust(dist(cos_optimal))
ordered_matrix <- cos_optimal[dend$order, dend$order]

lv_groups <- ifelse(grepl('LVstar1',colnames(ordered_matrix)),'LV1','LV2')
partition_groups <- substr(colnames(ordered_matrix),13,15)
partition_groups[grep('_i',partition_groups)] <- '1'
partition_groups <- as.numeric(partition_groups)
groups_col <- data.frame(LV = lv_groups,partition=partition_groups)
iteration_groups_cols <- str_split_fixed(colnames(ordered_matrix),'iter',2)
iteration_groups_cols <- as.numeric(iteration_groups_cols[,2])

partition_groups <- substr(rownames(ordered_matrix),13,15)
partition_groups[grep('_i',partition_groups)] <- '1'
partition_groups <- as.numeric(partition_groups)
lv_groups <- ifelse(grepl('LVstar1',rownames(ordered_matrix)),'LV1','LV2')
groups_row <- data.frame(LV = lv_groups,partition=partition_groups)

rownames(groups_col) <- colnames(ordered_matrix)
rownames(groups_row) <- rownames(ordered_matrix)
# png('../results/similarity_of_optimal_analytical_lvs_for_many_partitions.png',
#     height = 14,width = 16,units = 'in',res=600)
pheatmap::pheatmap(ordered_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         breaks = seq(-1, 1, length.out = 100),
         main = "Clustered Cosine Similarity Heatmap",
         annotation_col = groups_col,
         annotation_row = groups_row,
         fontsize = 18,
         show_colnames = FALSE,  
         show_rownames = FALSE,
         filename = '../results/similarity_of_optimal_analytical_lvs_for_many_partitions.png',
         height = 12,
         width = 14,
         units = 'in',
         res=600
         )
# dev.off()

cos_sim_all <- cos_sim_all %>% filter(extra_basis!=basis)
cos_sim_all <- cos_sim_all %>% filter(grepl('PC',basis))
ggplot(cos_sim_all,aes(x=sim,fill=extra_basis)) + geom_histogram() +xlab('cosine similarity with in vitro PCs')
ggsave('../results/cosine_similarity_PC_extra_basis.png',
       width = 9,
       height = 9,
       units = 'in',
       dpi=600)

groups_col <- groups_col %>% mutate(iteration_groups = iteration_groups_cols)
groups_col <- groups_col %>% rownames_to_column('var') %>% select(var,partition,iteration_groups,LV)
iteration_groups <- str_split_fixed(rownames(ordered_matrix),'iter',2)
iteration_groups <- as.numeric(iteration_groups[,2])
plot_df <- as.data.frame(ordered_matrix) %>%
  mutate(partition = partition_groups) %>% mutate(LV = lv_groups) %>% mutate(iteration_groups=iteration_groups) %>%
  gather('var','sim',-partition,-LV,-iteration_groups)
plot_df <- left_join(plot_df,groups_col,by='var')
plot_df <- plot_df %>% 
  mutate(keep = ifelse(partition.x==partition.y & iteration_groups.x==iteration_groups.y & LV.x==LV.y,
                       FALSE,TRUE)) %>% filter(keep == TRUE) %>%
  select(-keep) %>% mutate(comparison = paste0(LV.x,'/',LV.y)) %>%
  mutate(keep = ifelse(LV.x!=LV.y & partition.x==partition.y & iteration_groups.x==iteration_groups.y,
                       ifelse(LV.x=='LV1',TRUE,FALSE),
                       TRUE)) %>% filter(keep == TRUE) %>%
  select(-keep) %>% filter(partition.x==partition.y) %>%
  select(comparison,partition.x,iteration_groups.x,sim) %>%
  mutate(comparison = ifelse(comparison=="LV2/LV1","LV1/LV2",comparison)) %>%
  unique()
colnames(plot_df) <- c('comparison','partition','iteration','sim')
# cos_sim_all <- cos_sim_all %>% 
#   mutate(comparison = ifelse(grepl('star1',extra_basis),'LV1/PCs','LV2/PCs')) %>%
#   select(comparison,partition,iteration,sim) 
# plot_df <- rbind(plot_df,cos_sim_all)
plot_df <- plot_df %>% group_by(comparison,partition) %>% 
  mutate(abs_mu = mean(abs(sim))) %>% mutate(abs_se = sd(abs(sim))/sqrt(max(iteration_groups))) %>%
  mutate(mu = mean(sim)) %>% mutate(se = sd(sim)/sqrt(max(iteration_groups))) %>%
  ungroup()
ggplot(plot_df %>% select(comparison,partition,mu,se) %>% unique(),
       aes(x = partition,y = mu,color = comparison)) +
  scale_x_continuous(breaks = seq(0,1,0.2))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  geom_point(size=1.5)+
  geom_line(lwd = 0.75)+
  geom_errorbar(aes(ymax = mu + se,ymin = mu - se),
                width = 0.05,size=0.75)+
  xlab('data partition (%)')+
  ylab('cosine similarity')+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line())
ggsave('../figures/cosine_similarity_PC_extra_basis_convergence_plot.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../figures/cosine_similarity_PC_extra_basis_convergence_plot.eps',
       device = cairo_ps,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)

### Cumulative comparison
plot_df <- as.data.frame(ordered_matrix) %>%
  mutate(partition = partition_groups) %>% mutate(LV = lv_groups) %>% mutate(iteration_groups=iteration_groups)
row_inds <- which(grepl('LVstar1',rownames(plot_df)) & grepl('_part1_',rownames(plot_df)))
col_inds <- c(which(grepl('LVstar1',colnames(plot_df))))
              # which(colnames(plot_df)=='partition'),
              # which(colnames(plot_df)=='LV'),
              # which(colnames(plot_df)=='iteration_groups'))
LV1s_similarities <- plot_df[row_inds,col_inds]
# LV1s_similarities[upper.tri(LV1s_similarities,diag = TRUE)] <- -100
LV1s_similarities <- LV1s_similarities  %>% rownames_to_column('var.x') %>%
  gather('var.y','sim',-var.x)
# LV1s_similarities <- LV1s_similarities %>% filter(sim!=-100)
LV1s_similarities <- LV1s_similarities %>% filter(var.x!=var.y)
LV1s_similarities <- left_join(LV1s_similarities,
                               groups_col,
                               by=c('var.y'='var'))
## Repeat for LV2
row_inds <- which(grepl('LVstar2',rownames(plot_df)) & grepl('_part1_',rownames(plot_df)))
col_inds <- c(which(grepl('LVstar2',colnames(plot_df))))
LV2s_similarities <- plot_df[row_inds,col_inds]
LV2s_similarities <- LV2s_similarities  %>% rownames_to_column('var.x') %>%
  gather('var.y','sim',-var.x)
LV2s_similarities <- LV2s_similarities %>% filter(var.x!=var.y)
LV2s_similarities <- left_join(LV2s_similarities,
                               groups_col,
                               by=c('var.y'='var'))

all_cumulative_similarity <- rbind(LV1s_similarities,LV2s_similarities)
all_cumulative_similarity <- all_cumulative_similarity %>%
  group_by(LV,partition) %>% mutate(mu = mean(sim))  %>%
  mutate(std = sd(sim)) %>% ungroup()
ggplot(all_cumulative_similarity %>% select(var.x,var.y,LV,partition,mu,std) %>% unique(),
       aes(x=partition,y=mu,color = LV))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  geom_point(size=1.5)+
  geom_line(lwd = 0.75)+
  geom_errorbar(aes(ymax = mu + 3*std,ymin = mu - 3*std),
                width = 0.05,size=0.75)+
  xlab('data partition (%)')+
  ylab('cosine similarity with 100% data partition')+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line())
ggsave('../figures/all_cumulative_similarity.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../figures/all_cumulative_similarity.eps',
       device = cairo_ps,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)



#### Repeat for pathway activity similarity-----------------------------
cos_optimal <- lsa::cosine(WPaths_opt_all)
dend <- hclust(dist(cos_optimal))
ordered_matrix <- cos_optimal[dend$order, dend$order]

lv_groups <- ifelse(grepl('LVstar1',colnames(ordered_matrix)),'LV1','LV2')
partition_groups <- substr(colnames(ordered_matrix),13,15)
partition_groups[grep('_i',partition_groups)] <- '1'
partition_groups <- as.numeric(partition_groups)
groups_col <- data.frame(LV = lv_groups,partition=partition_groups)
iteration_groups_cols <- str_split_fixed(colnames(ordered_matrix),'iter',2)
iteration_groups_cols <- as.numeric(iteration_groups_cols[,2])

partition_groups <- substr(rownames(ordered_matrix),13,15)
partition_groups[grep('_i',partition_groups)] <- '1'
partition_groups <- as.numeric(partition_groups)
lv_groups <- ifelse(grepl('LVstar1',rownames(ordered_matrix)),'LV1','LV2')
groups_row <- data.frame(LV = lv_groups,partition=partition_groups)

rownames(groups_col) <- colnames(ordered_matrix)
rownames(groups_row) <- rownames(ordered_matrix)

pheatmap::pheatmap(ordered_matrix, 
                   color = colorRampPalette(c("blue", "white", "red"))(100), 
                   breaks = seq(-1, 1, length.out = 100),
                   main = "Clustered Cosine Similarity Heatmap",
                   annotation_col = groups_col,
                   annotation_row = groups_row,
                   fontsize = 18,
                   show_colnames = FALSE,  
                   show_rownames = FALSE,
                   filename = '../results/pathway_activity_similarity.png',
                   height = 12,
                   width = 14,
                   units = 'in',
                   res=600
)

groups_col <- groups_col %>% mutate(iteration_groups = iteration_groups_cols)
groups_col <- groups_col %>% rownames_to_column('var') %>% select(var,partition,iteration_groups,LV)
iteration_groups <- str_split_fixed(rownames(ordered_matrix),'iter',2)
iteration_groups <- as.numeric(iteration_groups[,2])
plot_df <- as.data.frame(ordered_matrix) %>%
  mutate(partition = partition_groups) %>% mutate(LV = lv_groups) %>% mutate(iteration_groups=iteration_groups) %>%
  gather('var','sim',-partition,-LV,-iteration_groups)
plot_df <- left_join(plot_df,groups_col,by='var')
plot_df <- plot_df %>% 
  mutate(keep = ifelse(partition.x==partition.y & iteration_groups.x==iteration_groups.y & LV.x==LV.y,
                       FALSE,TRUE)) %>% filter(keep == TRUE) %>%
  select(-keep) %>% mutate(comparison = paste0(LV.x,'/',LV.y)) %>%
  mutate(keep = ifelse(LV.x!=LV.y & partition.x==partition.y & iteration_groups.x==iteration_groups.y,
                       ifelse(LV.x=='LV1',TRUE,FALSE),
                       TRUE)) %>% filter(keep == TRUE) %>%
  select(-keep) %>% filter(partition.x==partition.y) %>%
  select(comparison,partition.x,iteration_groups.x,sim) %>%
  mutate(comparison = ifelse(comparison=="LV2/LV1","LV1/LV2",comparison)) %>%
  unique()
colnames(plot_df) <- c('comparison','partition','iteration','sim')
plot_df <- plot_df %>% group_by(comparison,partition) %>% 
  mutate(abs_mu = mean(abs(sim))) %>% mutate(abs_se = sd(abs(sim))/sqrt(max(iteration_groups))) %>%
  mutate(mu = mean(sim)) %>% mutate(se = sd(sim)/sqrt(max(iteration_groups))) %>%
  mutate(std = sd(sim)) %>%
  ungroup()
ggplot(plot_df %>% select(comparison,partition,mu,se) %>% unique(),
       aes(x = partition,y = mu,color = comparison)) +
  scale_x_continuous(breaks = seq(0,1,0.2))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  geom_point(size=1.5)+
  geom_line(lwd = 0.75)+
  geom_errorbar(aes(ymax = mu + se,ymin = mu - se),
                width = 0.05,size=0.75)+
  xlab('data partition (%)')+
  ylab('cosine similarity')+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line())
ggsave('../figures/cosine_similarity_pathways_convergence_plot.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../figures/cosine_similarity_pathways_convergence_plot.eps',
       device = cairo_ps,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)

### Cumulative comparison
plot_df <- as.data.frame(ordered_matrix) %>%
  mutate(partition = partition_groups) %>% mutate(LV = lv_groups) %>% mutate(iteration_groups=iteration_groups)
row_inds <- which(grepl('LVstar1',rownames(plot_df)) & grepl('_part1_',rownames(plot_df)))
col_inds <- c(which(grepl('LVstar1',colnames(plot_df))))
LV1s_similarities <- plot_df[row_inds,col_inds]
LV1s_similarities <- LV1s_similarities  %>% rownames_to_column('var.x') %>%
  gather('var.y','sim',-var.x)
LV1s_similarities <- LV1s_similarities %>% filter(var.x!=var.y)
LV1s_similarities <- left_join(LV1s_similarities,
                               groups_col,
                               by=c('var.y'='var'))
## Repeat for LV2
row_inds <- which(grepl('LVstar2',rownames(plot_df)) & grepl('_part1_',rownames(plot_df)))
col_inds <- c(which(grepl('LVstar2',colnames(plot_df))))
LV2s_similarities <- plot_df[row_inds,col_inds]
LV2s_similarities <- LV2s_similarities  %>% rownames_to_column('var.x') %>%
  gather('var.y','sim',-var.x)
LV2s_similarities <- LV2s_similarities %>% filter(var.x!=var.y)
LV2s_similarities <- left_join(LV2s_similarities,
                               groups_col,
                               by=c('var.y'='var'))

all_cumulative_similarity <- rbind(LV1s_similarities,LV2s_similarities)
all_cumulative_similarity <- all_cumulative_similarity %>%
  group_by(LV,partition) %>% mutate(mu = mean(sim))  %>%
  mutate(std = sd(sim)) %>% ungroup()
ggplot(all_cumulative_similarity,
       aes(x=partition,y=mu,color = LV))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  geom_point(size=1.5)+
  geom_line(lwd = 0.75)+
  geom_errorbar(aes(ymax = mu + 3*std,ymin = mu - 3*std),
                width = 0.05,size=0.75)+
  xlab('data partition (%)')+
  ylab('cosine similarity with 100% data partition')+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line())
ggsave('../figures/all_cumulative_similarity_pathways.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../figures/all_cumulative_similarity_pathways.eps',
       device = cairo_ps,
       width = 14,
       height = 10.5,
       units = 'in',
       dpi = 600)

### Cumulative pathway activity
activity_df <- as.data.frame(t(WPaths_opt_all))
activity_df <- activity_df %>% rownames_to_column('var') %>%
  gather('pathway','activity',-var)
activity_df <- left_join(activity_df,groups_col)
ggplot(activity_df %>% filter(LV=='LV1'),
       aes(x=as.factor(partition),y=activity))+
  ggtitle('extra LV 1')+
  scale_y_continuous(breaks = seq(-10,10,2))+
  geom_violin(fill = '#F8766D')+
  xlab('data partition (%)')+
  ylab('activity')+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line()) +
  facet_wrap(~pathway,ncol = 3,scales = 'free_y')
ggsave('../figures/LV1_cumulative_singed_activity_pathways.png',
       width = 18,
       height = 18,
       units = 'in',
       dpi = 600)
ggsave('../figures/LV1_cumulative_singed_activity_pathways.eps',
       device = cairo_ps,
       width = 16,
       height = 16,
       units = 'in',
       dpi = 600)

ggplot(activity_df %>% filter(LV=='LV2'),
       aes(x=as.factor(partition),y=activity))+
  ggtitle('extra LV 2')+
  scale_y_continuous(breaks = seq(-10,10,2))+
  geom_violin(fill = '#00BFC4')+
  xlab('data partition (%)')+
  ylab('activity')+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line()) +
  facet_wrap(~pathway,ncol = 3,scales = 'free_y')
ggsave('../figures/LV2_cumulative_singed_activity_pathways.png',
       width = 18,
       height = 18,
       units = 'in',
       dpi = 600)
ggsave('../figures/LV2_cumulative_singed_activity_pathways.eps',
       device = cairo_ps,
       width = 16,
       height = 16,
       units = 'in',
       dpi = 600)
