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
val_R2_shuffled <- NULL
val_r_shuffled <- NULL
val_mae_shuffled <- NULL
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
  }
  tuning_df <- rbind(tuning_df,
                     data.frame(set='train',r = train_r,MAE = train_mae,R2=train_R2,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)),
                     data.frame(set='validation',r = val_r,MAE = val_mae,R2=val_R2,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)),
                     data.frame(set='shuffled',r = val_r_shuffled,MAE = val_mae_shuffled,R2=val_R2_shuffled,fold = seq(1,num_folds),LVs=rep(num_LVs[i],num_folds)))
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
  scale_y_continuous(n.breaks = 12)+
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

ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr.eps',
       device = cairo_ps,
       width = 16,
       height = 9,
       units = 'in',
       dpi=600)

### Check different partitions of data based on the average similarity of training and test set
### First find similarity distribution of all data
Xh <- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh.rds')
Yh <- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh.rds')
rownames(Yh) <- rownames(Xh)
Xh_test<- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Xh_test.rds')
Yh_test<- readRDS('../preprocessing/TrainingValidationData/WholePipeline/TrainTestValPLSR/Yh_test.rds')
all_correlations <- cor(t(Xh))
# all_correlations <- all_correlations[upper.tri(all_correlations,diag = FALSE)]
diag(all_correlations) <- NA
all_correlations <- apply(all_correlations,1,max,na.rm=TRUE)
# min_lvl <- sign(min(all_correlations)) * ceiling(10*abs(min(all_correlations)))/10
# max_lvl <- ceiling(10*max(all_correlations))/10
# # binned_cor <- cut(all_correlations,breaks = seq(min_lvl,max_lvl,0.2))
# bins_thresh_corr <- seq(min_lvl,max_lvl,0.2)[2:length(seq(min_lvl,max_lvl,0.2))]
# hist(all_correlations,40,main='Pearson`s r of human samples')
bins_thresh_corr <- c(0.2,0.3,0.4,0.5,0.7)

partition_indices <- list()
for (i in 1:length(bins_thresh_corr)){
  # if (i==1){
  #   [[i]] <- names(which(all_correlations<=bins_thresh_corr[i]))
  # }else{
  #   partition_indices[[i]] <- names(which(all_correlations<=bins_thresh_corr[i] & all_correlations>bins_thresh_corr[i-1]))
  # }
  partition_indices[[i]] <- names(which(all_correlations<=bins_thresh_corr[i]))
}

# Begin iterations
performance_train <- data.frame()
performance_val <- data.frame()
iterations <- 10
number_samp <- ceiling(0.2 * nrow(Xh))
num_lvs <- 20
for (j in 1:length(partition_indices)){
  message(paste0('Begun withhelding samples with maximum correlation of ',bins_thresh_corr[j],' with all the data'))
  for (LV in 1:num_lvs){
    if (length(partition_indices[[j]])>number_samp){
      kk <- iterations
    }else{
      kk <- 1
    }
    for (i in 1:kk){
      Xh_val <- Xh[partition_indices[[j]],]
      if (kk>1){
        Xh_val <- as.matrix(as.data.frame(Xh_val) %>% sample_n(number_samp))
      }
      Xh_parted <- Xh[which(!(rownames(Xh) %in% rownames(Xh_val))),]
      Yh_parted <- Yh[rownames(Xh_parted),]
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
      
      r_plsr <- diag(cor(y_hat_plsr,Yh_parted))
      performance_train <- rbind(performance_train,
                                 data.frame(model = 'PLSR',
                                            num_LVs = LV,
                                            max_corr = bins_thresh_corr[j],
                                            sample=i,
                                            NAS = r_plsr['NAS'],
                                            fibrosis = r_plsr['fibrosis']))
      
      r_plsr <- diag(cor(y_hat_val_plsr,Yh_val))
      performance_val <- rbind(performance_val,
                               data.frame(model = 'PLSR',
                                          num_LVs = LV,
                                          max_corr = bins_thresh_corr[j],
                                          sample=i,
                                          NAS = r_plsr['NAS'],
                                          fibrosis = r_plsr['fibrosis']))
      
      print(paste0('Finished iteration ',i))
    }
    print(paste0('Finished PLSR with ',LV,' latent variables'))
  }
}
# saveRDS(performance_val,'../results/performance_val_difficult_partitions.rds')
# saveRDS(performance_train,'../results/performance_train_difficult_partitions.rds')

val_plot <- performance_val %>% gather('phenotype','r',-model,-num_LVs,-max_corr,-sample)  %>% 
  # group_by(model,num_LVs,max_corr,sample) %>% mutate(r = mean(r)) %>% select(-phenotype) %>% ungroup() %>%
  group_by(model,num_LVs,max_corr,phenotype) %>% mutate(mu = mean(r)) %>% 
  mutate(std = sd(r)) %>% ungroup() %>%
  select(-sample,-r,) %>% unique() %>%
  mutate(std = ifelse(is.na(std),0,std)) 
val_plot$max_corr <- factor(val_plot$max_corr,levels = bins_thresh_corr)
val_plot <- val_plot %>% mutate(max_corr = paste0('max r = ',max_corr))
val_plot <- val_plot %>% filter(!is.na(mu))

train_plot <- performance_train %>% gather('phenotype','r',-model,-num_LVs,-max_corr,-sample)  %>% 
  # group_by(model,num_LVs,max_corr,sample) %>% mutate(r = mean(r)) %>% select(-phenotype) %>% ungroup() %>%
  group_by(model,num_LVs,max_corr,phenotype) %>% mutate(mu = mean(r)) %>% 
  mutate(std = sd(r)) %>% ungroup() %>%
  select(-sample,-r) %>% unique() %>%
  mutate(std = ifelse(is.na(std),0,std)) 
train_plot$max_corr <- factor(train_plot$max_corr,levels = bins_thresh_corr)
train_plot <- train_plot %>% mutate(max_corr = paste0('max r = ',max_corr))
train_plot <- train_plot %>% filter(max_corr %in% val_plot$max_corr)
performance_plot <- rbind(train_plot %>% mutate(set='train'),
                          val_plot %>% mutate(set='test'))

ggplot(performance_plot,
       aes(x=num_LVs,y=mu,colour=phenotype,shape=set))+
  geom_line(aes(linetype = set),lwd=1)+ 
  geom_point()+ 
  geom_errorbar(aes(ymax = mu + std/sqrt(iterations),ymin=mu - std/sqrt(iterations)))+
  scale_x_continuous(breaks = seq(1,20,2))+
  xlab('number of latent variables') +
  ylab('pearson`s r')+
  ggtitle('Performance for test sets of increasing difficulty')+
  theme_pubr(base_family = 'Arial',base_size = 18)+
  theme(text = element_text(size=18,family = 'Arial'),
        legend.position = 'right',
        panel.grid.major = element_line(),
        strip.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~max_corr)
ggsave('../preprocessing/TrainingValidationData/WholePipeline/tunning_plsr_difficult.png',
       width = 12,
       height = 9,
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
      
      r_plsr <- diag(cor(y_hat_plsr,Yh_parted))
      performance_train <- rbind(performance_train,
                                 data.frame(model = 'PLSR',
                                            num_LVs = LV,
                                            partition = p,
                                            sample=i,
                                            NAS = r_plsr['NAS'],
                                            fibrosis = r_plsr['fibrosis']))
      
      r_plsr <- diag(cor(y_hat_val_plsr,Yh_val))
      performance_val <- rbind(performance_val,
                               data.frame(model = 'PLSR',
                                          num_LVs = LV,
                                          partition = p,
                                          sample=i,
                                          NAS = r_plsr['NAS'],
                                          fibrosis = r_plsr['fibrosis']))
      
    }
    print(paste0('Finished PLSR with ',LV,' latent variables'))
  }
}
# saveRDS(performance_val,'../results/performance_val_random_partitions_many_lvs.rds')
# saveRDS(performance_train,'../results/performance_train_random_partitions_many_lvs.rds')

val_plot <- performance_val %>% gather('phenotype','r',-model,-num_LVs,-partition,-sample)  %>% 
  group_by(model,num_LVs,partition,phenotype) %>% mutate(mu = mean(r)) %>% 
  mutate(std = sd(r)) %>% ungroup() %>%
  unique() 
train_plot <- performance_train %>% gather('phenotype','r',-model,-num_LVs,-partition,-sample)  %>% 
  group_by(model,num_LVs,partition,phenotype) %>% mutate(mu = mean(r)) %>% 
  mutate(std = sd(r)) %>% ungroup() %>%
  unique()
performance_plot <- rbind(train_plot %>% mutate(set='train'),
                          val_plot %>% mutate(set='test'))
performance_plot$partition <- factor(performance_plot$partition,levels = partitions)
performance_plot <- performance_plot %>% mutate(partition = paste0('train size = ',partition,'%'))

ggplot(performance_plot,
       aes(x=num_LVs,y=r,colour=set,shape=phenotype))+
  geom_smooth(aes(linetype=phenotype),lwd=1)+ 
  # geom_point()+ 
  # geom_errorbar(aes(ymax = mu + std/sqrt(iterations),ymin=mu - std/sqrt(iterations)))+
  scale_x_continuous(breaks = seq(1,20,2))+
  scale_y_continuous(breaks = seq(0.5,1,0.1))+
  xlab('number of latent variables') +
  ylab('pearson`s r')+
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
                     predI = 8,
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
                               predI = 8,
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
                               predI = 8,
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
                              predI = 8,
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
# saveRDS(performance_df,'../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned.rds')

avg_mae <- mean(abs(Yh_test-apply(Yh,2,mean)))
plotting_performance_df <- performance_df %>% 
  gather('metric','value',-set,-fold) %>%
  mutate(metric=ifelse(metric=='r','pearson`s r',metric))
plotting_performance_df$set <- factor(plotting_performance_df$set,
                                      levels = c('train','validation','test','shuffle Y','shuffle X','random X'))
p <- (ggboxplot(plotting_performance_df %>% filter(metric=='MAE'),x='set',y='value',color = 'set',add='jitter')+
    scale_y_continuous(n.breaks = 15)+
    geom_hline(yintercept = avg_mae,linetype='dashed',color='black',lwd=1)+
    annotate('text',x=3,y=1.95,size=6,label = 'error from the mean of the data')+
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
                       label.y = c(1,1.5,1.6,1.7)) +
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
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
ggsave('../preprocessing/TrainingValidationData/WholePipeline/performance_df_tuned.eps',
       plot = p,
       device = cairo_ps,
       width = 12,
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
# saveRDS(performance_1,'../results/performance_df_human_plsr.rds')
#Plot
performance_1 <- performance_1 %>% mutate(type = ifelse(type=='model',set,type))
performance_1 <- performance_1 %>% filter(metric=='r') %>% select(-metric)
performance_1$type <- factor(performance_1$type ,levels=c('train','test','shuffle Y','shuffle X','random X'))
performance_1 <- performance_1 %>% mutate(phenotype = ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
ggboxplot(performance_1,x='type',y='value',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(-0.75,1,0.25),limits = c(NA,1))+
  xlab('')+ylab('pearson`s correlation')+
  theme(text = element_text(size=28,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=22),
        # strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('train','test'),
                                        c('test','shuffle Y'),
                                        c('test','shuffle X'),
                                        c('test','random X')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.01,
                     label.y = c(0.9,0.75, 0.8, 0.85),
                     size=7)+
  facet_wrap(~phenotype,nrow = 2)
ggsave('../figures/performance_df_just_human_plsr.png',
       height = 9,
       width = 12,
       units = 'in',
       dpi=600)
ggsave('../figures/performance_df_just_human_plsr.eps',
       device = cairo_ps,
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
performance_2 <- performance_2 %>% filter(metric=='r') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(r = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_2,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+ylab('pearson`s correlation')+
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
performance_3 <- performance_3 %>% filter(metric=='r') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(r = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_3,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+ylab('pearson`s correlation')+
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
performance_4 <- performance_4 %>% filter(metric=='r') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(r = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_4,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(n.breaks = 10)+
  xlab('')+ylab('pearson`s correlation')+
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
# saveRDS(performance_5,'../results/performance_df_optimal_direction.rds')
# performance_5_tosave <- performance_5
#plot
performance_5$type <- factor(performance_5$type ,levels=c('model','shuffle Wopt','shuffle Bh'))
performance_5$set <- factor(performance_5$set ,levels=c('train','test'))
performance_5 <- performance_5 %>% filter(metric=='r') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(r = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_5,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(round(min(performance_5$r),1),1,0.1),
                     limits = c(NA,1.15))+
  xlab('')+ylab('pearson`s correlation')+
  theme(text = element_text(size=22,family = 'Arial'),
        legend.position = 'none',
        axis.text.x = element_text(size=18),
        strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  facet_wrap(~set,scales = 'free')+
  stat_compare_means(comparisons = list(c('model','shuffle Wopt'),
                                        c('model','shuffle Bh')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     step.increase = 0.04,
                     size=6)

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
                                                    LVs = num_LVS)
# saveRDS(performance_6,'../results/performance_df_analytical.rds')

performance_6$type <- factor(performance_6$type ,levels=c('model','shuffle Wopt','shuffle Bh'))
performance_6$set <- factor(performance_6$set ,levels=c('train','test'))
performance_6 <- performance_6 %>% filter(metric=='r') %>% select(-metric) %>% group_by(set,type,task,fold) %>% mutate(r = mean(value)) %>%
  ungroup() %>% select(-value,-phenotype) %>% unique()
ggboxplot(performance_6,x='type',y='r',color='type',add='jitter') +
  scale_y_continuous(breaks = seq(round(min(performance_6$r),1),1,0.1),
                     limits = c(NA,1.15))+
  xlab('')+ylab('pearson`s correlation')+
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
# performance_all_plot <- performance_all %>% 
#   mutate(keep=ifelse(task=='human_plsr',ifelse(type %in% c('model','shuffle X'),TRUE,FALSE),
#                      ifelse(type=='model',TRUE,FALSE))) %>% 
#   filter(keep==TRUE) %>% select(-keep) %>%
#   mutate(approach = ifelse(task=='human_plsr','PLSR',
#                         ifelse(task=='human_backprojected','backprojected',
#                                ifelse(task=='human_backprojected_retrained','backprojected retrained',
#                                       ifelse(task=='human_backprojected_into_translatable_lvs','translatable LVs',
#                                              ifelse(task=='analytical_optimal','analytical Wopt',
#                                              'Wopt')))))) %>%
#   mutate(approach = ifelse(approach=='PLSR',
#                            ifelse(type=='model',approach,type),
#                            approach)) %>%
#   select(-type,-task)
performance_all_plot <- performance_all %>%  filter(metric=='r') %>% select(-metric) %>% mutate(r=value) %>% select(-value)%>%
  mutate(keep=ifelse(task=='human_plsr',ifelse(type %in% c('model','shuffle X'),TRUE,FALSE),
                     ifelse(type=='model',TRUE,FALSE))) %>% 
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(approach = ifelse(task=='human_plsr','human genes',
                           ifelse(task=='human_backprojected','back-projected',
                                  ifelse(task=='human_backprojected_retrained','back-projected re-trained',
                                         ifelse(task=='human_backprojected_into_translatable_lvs','translatable LVs',
                                                ifelse(task=='analytical_optimal','optimized MPS',
                                                       'Wopt')))))) %>%
  mutate(approach = ifelse(approach=='human genes',
                           ifelse(type=='model',approach,type),
                           approach)) %>%
  select(-type,-task)
performance_all_plot <- performance_all_plot %>%
  filter(approach %in% c('human genes','back-projected','optimized MPS','shuffle X'))

# performance_all_plot$approach <- factor(performance_all_plot$approach,
#                                         levels = c('PLSR',
#                                                    'backprojected',
#                                                    'backprojected retrained',
#                                                    'translatable LVs',
#                                                    'Wopt',
#                                                    'analytical Wopt',
#                                                    'shuffle X'))
performance_all_plot$approach <- factor(performance_all_plot$approach,
                                        levels = c('human genes',
                                                   'back-projected',
                                                   'optimized MPS',
                                                   'shuffle X'))
performance_all_plot <- performance_all_plot %>% mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
# performance_all_plot$phenotype <- factor(performance_all_plot$phenotype,levels = c('NAS','Fibrosis stage'))
performance_all_plot$phenotype <- factor(performance_all_plot$phenotype)
performance_all_plot$phenotype <- factor(performance_all_plot$phenotype,
                                         levels = rev(levels(performance_all_plot$phenotype)))
p_train <- ggboxplot(performance_all_plot %>% filter(set=='train') %>% filter(approach!='Wopt'),
                     x='approach',y='r',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(0.4,1,0.05),limits = c(0.4,NA))+
  xlab('')+ ylab('pearson`s correlation') +
  ggtitle('10-fold train performance in predicting phenotype')+
  theme(text = element_text(size=26,family = 'Arial'),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=26), #angle = 25
        # strip.text = element_text(face = 'bold'),
        panel.grid.major.y = element_line(linewidth = 1)) +
  stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
                                        c('human genes','back-projected'),
                                        c('optimized MPS','back-projected')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     label.y = c(1,0.95,0.95),
                     size = 8)+
  facet_wrap(~phenotype,nrow = 2)
print(p_train)  
ggsave('../figures/approaches_comparison_training.png',
       plot = p_train,
       height = 9,
       width = 11,
       units = 'in',
       dpi=600)

p_test <- ggboxplot(performance_all_plot %>% filter(set=='test')%>% filter(approach!='Wopt'),
                    x='approach',y='r',color='approach',add='jitter') +
  scale_y_continuous(breaks = seq(-0.5,1,0.15),limits = c(NA,1.05))+
  xlab('')+ylab('pearson`s correlation') +
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
                                        c('human genes','back-projected'),
                                        c('optimized MPS','back-projected'),
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



### See how training performance converges
### when including incrementally the extra basis--------------------------------------------------------------------
train_r <- NULL
train_r_translatables <- NULL
train_r_basis_1 <- NULL
train_r_all_extra_bases <- NULL

test_r <- NULL
test_r_shuffled <- NULL
test_r_translatables <- NULL
test_r_basis_1 <- NULL
test_r_all_extra_bases <- NULL

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
  train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
  test_r[j] <- mean(diag(cor(y_val_hat,y_val)))
  
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
  test_r_shuffled[j]<- mean(diag(cor(y_hat_val,y_val)))
  
  message('Finished running initial PLSR for humans')
  
  ### Get linear combinaiton of translatable LVs
  Wm_new <- get_translatable_LV(x_train, y_train, Wh, Wm,
                                rbind(apply(y_train,2,mean),Bh),
                                find_extra = FALSE)
  Wm_new <- Wm_new$Wm_new
  colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
  y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  train_r_translatables[j] <- mean(diag(cor(y_hat_train,y_train)))
  y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  test_r_translatables[j] <- mean(diag(cor(y_hat_test,y_val)))
  
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
  train_r_basis_1[j] <- mean(diag(cor(y_hat_train,y_train)))
  y_hat_test <- cbind(1, x_val %*% Wm_tot_1 %*% t(Wm_tot_1) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  test_r_basis_1[j] <- mean(diag(cor(y_hat_test,y_val)))
  message('Finished analytical solution performance with 1 extra vector')
  # predict with one extra basis vectors
  y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  train_r_all_extra_bases[j] <- mean(diag(cor(y_hat_train,y_train)))
  y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
  test_r_all_extra_bases[j] <- mean(diag(cor(y_hat_test,y_val)))
  message('Finished analytical solution performance with 2 extra vectors')
}

## First plot scatter plots of prediction VS ground truth
# Calculate correlation coefficients
cor_results <- df_scatterPlot %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))
ggplot(df_scatterPlot,aes(x = true,y=prediction)) +
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = cor_results, 
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
              hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
  facet_wrap(~phenotype,scales = 'free')+
  theme_pubr(base_family = 'Arial',base_size=25)+
  theme(text = element_text(family = 'Arial',size=25),
        panel.grid.major = element_line())
ggsave('../figures/10foldTest_Scatterplot_human_plsr.png',
          height = 6,
          width=9,
          units = 'in',
          dpi=600)

# repreat for backprojection
cor_resultst_backproj <- df_scatterPlot_backproj %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))
ggplot(df_scatterPlot_backproj,aes(x = true,y=prediction)) +
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = cor_resultst_backproj, 
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
            hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
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
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
ggplot(all_scatter_plot,aes(x = true,y=prediction)) +
  ylab('Predicted') + xlab('Measured')+
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = all_cor_results, 
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
            hjust = 0, vjust =  1.5, size = 9, family = 'Arial') +
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
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
            hjust = 0, vjust =  1.5, size = 9, family = 'Arial') +
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
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))
(ggplot(all_scatter_plot_backproj %>% filter(phenotype=='Fibrosis stage'),aes(x = true,y=prediction)) +
  geom_jitter(width = 0.05,color='#4682B4') + 
  geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
  geom_text(data = all_cor_results_backproj%>% filter(phenotype=='Fibrosis stage'), 
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
            hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
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
               aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
               hjust = 0, vjust =  1.5, size = 8, family = 'Arial') +
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
              aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
              hjust = 0, vjust =  1.5, size = 9, family = 'Arial') +
    ylab('Predicted') + xlab('Measured')+
    ylim(c(0,4))+
    facet_wrap(~phenotype,scales = 'free')+
    theme_pubr(base_family = 'Arial',base_size=28)+
    theme(text = element_text(family = 'Arial',size=28),
          plot.title = element_text(hjust = 0.5,face = 'bold'),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face='bold'),
          panel.grid.major = element_line()))+
  (ggplot(all_scatter_plot_backproj %>% filter(phenotype!='Fibrosis stage'),aes(x = true,y=prediction)) +
     geom_jitter(width = 0.05,color='#4682B4') + 
     geom_abline(slope=1,intercept = 0,linetype = 'dashed',color='black',linewidth = 1.5)+
     geom_text(data = all_cor_results_backproj%>% filter(phenotype!='Fibrosis stage'), 
               aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
               hjust = 0, vjust =  1.5, size = 9, family = 'Arial') +
     ylab('Predicted') + xlab('Measured')+
     ylim(c(0,8))+
     facet_wrap(~phenotype,scales = 'free')+
     theme_pubr(base_family = 'Arial',base_size=28)+
     theme(text = element_text(family = 'Arial',size=28),
           axis.title.x = element_text(hjust = -0.4,face='bold'),
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

### Now see convergence of performance
train_performance_res <- rbind(data.frame(r = train_r,
                                          fold = seq(1,length(train_r)),
                                          model = rep('PLSR',length(train_r))),
                               data.frame(r = train_r_translatables,
                                          fold = seq(1,length(train_r_translatables)),
                                          model = rep('translatable PCs',length(train_r_translatables))),
                               data.frame(r = train_r_basis_1,
                                          fold = seq(1,length(train_r_basis_1)),
                                          model = rep('PCs + extra LV1',length(train_r_basis_1))),
                               data.frame(r = train_r_all_extra_bases,
                                          fold = seq(1,length(train_r_all_extra_bases)),
                                          model = rep('PCs + extra LV1 + extra LV2',length(train_r_all_extra_bases))))
train_performance_res <- train_performance_res %>% group_by(model) %>% 
  mutate(mu = mean(r)) %>% mutate(std = sd(r)) %>%
  mutate(min_r = min(r)) %>% mutate(max_r = max(r))%>%
  ungroup() %>% mutate(set='train')
mu_plsr_train <- mean(train_r)
sd_plsr_train <- sd(train_r)
test_performance_res <- rbind(data.frame(r = test_r,
                                          fold = seq(1,length(test_r)),
                                          model = rep('PLSR',length(test_r))),
                               data.frame(r = test_r_translatables,
                                          fold = seq(1,length(test_r_translatables)),
                                          model = rep('translatable PCs',length(test_r_translatables))),
                               data.frame(r = test_r_basis_1,
                                          fold = seq(1,length(test_r_basis_1)),
                                          model = rep('PCs + extra LV1',length(test_r_basis_1))),
                               data.frame(r = test_r_all_extra_bases,
                                          fold = seq(1,length(test_r_all_extra_bases)),
                                          model = rep('PCs + extra LV1 + extra LV2',length(test_r_all_extra_bases))))
test_performance_res <- test_performance_res %>% group_by(model) %>% 
  mutate(mu = mean(r)) %>% mutate(std = sd(r)) %>%
  mutate(min_r = min(r)) %>% mutate(max_r = max(r))%>%
  ungroup() %>% mutate(set='test')
mu_plsr_test <- mean(test_r)
sd_plsr_test <- sd(test_r)
all_performance_res <- rbind(train_performance_res,
                             test_performance_res)
all_performance_res$model <- factor(all_performance_res$model,
                                    levels = c('PLSR',
                                               'translatable PCs',
                                               'PCs + extra LV1',
                                               'PCs + extra LV1 + extra LV2'))

ggplot(all_performance_res %>% filter(model!='PLSR') %>% select(model,set,mu,std),
       aes(x=model,y=mu,color = set , group =set))+
  geom_point(size=1.5)+
  geom_line(lwd = 0.75)+
  geom_errorbar(aes(ymax = mu + std/sqrt(num_folds),ymin = mu - std/sqrt(num_folds)),
                width = 0.05,size=0.75)+
  ylim(NA,1) +
  ylab('pearson`s correlation')+
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
           size=9)+
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
           size=9)+
  ### random performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mean(test_r_shuffled) - sd(test_r_shuffled)/sqrt(num_folds), 
                  ymax = mean(test_r_shuffled) + sd(test_r_shuffled)/sqrt(num_folds)), 
              fill = "#F564E3", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mean(test_r_shuffled), 
           yend = mean(test_r_shuffled), 
           color = "#F564E3", linewidth = 1,linetype='dashed') +
  annotate('text',x=1,
           y=mean(test_r_shuffled) + 0.03,
           label="shuffled model performance", 
           hjust = 0 , 
           size=9)+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75)+
  theme_pubr(base_family = 'Arial',base_size = 27)+
  theme(text = element_text(family = 'Arial',size = 27),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.75),
        panel.grid.major = element_line())
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

### Explore how the evolutionary algorithm VS analytical approach works-------------------------------------

#Load data
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

# Parameters of PLSR
num_LVS <- 8

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
# saveRDS(similarity_df,'../results/cosine_similarity_df_analyt_evol.rds')
similarity_2_plot <- similarity_df %>% group_by(vector_1,vector_2) %>%
  mutate(mean_abs = mean(abs(sim))) %>% ungroup() %>%
  select(vector_1,vector_2,mean_abs) %>% unique() %>%
  spread('vector_2','mean_abs') %>% column_to_rownames('vector_1')
png('../results/analyt_evol_cosine_sim.png',width = 9,height = 9,units = 'in',res = 600)
corrplot::corrplot(as.matrix(similarity_2_plot), is.corr = F,addCoef.col = T,
                   diag = F,method = 'shade',tl.col = 'black',col = brewer.pal(n = 9, name = "Reds"),
                   type="upper")
dev.off()
# corrr::network_plot(similarity_2_plot,colors = c('white',"red"),min_cor = 0.01)
# ggsave('../results/analyt_evol_cosine_sim_net.png',
#        height = 9,
#        width = 9,
#        units = 'in',
#        dpi=600)


### Explore how analytical LVs change by: ----------------------
### a) Human data partition
### b) number of LVs used in PLSR of human genes
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
# saveRDS(Wm_opt_all,'../results/Wm_opt_all_different_partitions.rds')
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
### b) number of LVs used in PLSR of human genes
