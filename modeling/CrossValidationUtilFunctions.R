library(tidyverse)
library(caret)
library(ropls)
cross_validation_complete_pipeline <- function(W_invitro,
                                               folds,
                                               file_loc,
                                               target_dataset,
                                               task = c('human_plsr',
                                                        'human_backprojected',
                                                        'human_backprojected_retrained',
                                                        'human_backprojected_into_translatable_lvs',
                                                        'human_backprojected_into_optimal_lvs'),
                                               LVs=9){
  val_mae <- NULL
  train_mae <- NULL
  val_r <- NULL
  train_r <- NULL
  val_r <- NULL
  all_models <- NULL
  val_r_shuffle_y <- NULL
  val_mae_shuffle_y <- NULL
  val_r_shuffle_x <- NULL
  val_mae_shuffle_x <- NULL
  val_r_random_x <- NULL
  val_mae_random_x <- NULL
  performance_df <- NULL
  train_r_shuffle_w <- NULL
  val_r_shuffle_w <- NULL
  val_mae_shuffle_w <- NULL
  train_mae_shuffle_w <- NULL
  train_r_shuffle_bh <- NULL
  val_r_shuffle_bh <- NULL
  val_mae_shuffle_bh <- NULL
  train_mae_shuffle_bh <- NULL
  for (j in 1:folds){
    message(paste0('Begun fold ',j))
    x_train <- readRDS(paste0(file_loc,'Xh_train',j,'.rds'))
    y_train <- readRDS(paste0(file_loc,'Yh_train',j,'.rds'))
    x_val <- readRDS(paste0(file_loc,'Xh_val',j,'.rds'))
    y_val <- readRDS(paste0(file_loc,'Yh_val',j,'.rds'))
    
    plsr_model <- opls(x = x_train, 
                       y = y_train,
                       predI = LVs,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    y_train_hat <- predict(plsr_model,x_train)
    y_val_hat <- predict(plsr_model,x_val)
    train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
    val_mae[j] <- mean(abs(y_val_hat-y_val))
    train_mae[j] <- mean(abs(y_train_hat-y_train))

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
    y_hat_val <- predict(plsr_model_shuffle_y,x_val)
    val_r_shuffle_y[j]<- mean(diag(cor(y_hat_val,y_val)))
    val_mae_shuffle_y[j] <- mean(abs(y_hat_val-y_val))
    
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
    y_hat_val <- predict(plsr_model_shuffle_x,x_val)
    val_r_shuffle_x[j]<- mean(diag(cor(y_hat_val,y_val)))
    val_mae_shuffle_x[j] <- mean(abs(y_hat_val-y_val))
    
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
    y_hat_test <- predict(plsr_model_random_x,x_val)
    val_r_random_x[j]<- mean(diag(cor(y_hat_test,y_val)))
    val_mae_random_x[j] <- mean(abs(y_hat_test-y_val))
    
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
    Thm_train <- x_train %*% W_invitro %*% t(W_invitro) %*% Wh
    Th_val<- x_val %*% Wh
    Thm_val <- x_val %*% W_invitro %*% t(W_invitro) %*% Wh
    
    print('Finished running initial PLSR for humans')
    
    if (task=='human_backprojected'){
      y_hat_test <- cbind(1, Thm_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      y_hat_train <- cbind(1, Thm_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      train_r[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae[j] <- mean(abs(y_hat_test-y_val))
      train_mae[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle W_invitro
      W_invitro_suffled <- W_invitro[sample.int(nrow(W_invitro)),]
      rownames(W_invitro_suffled) <- rownames(W_invitro)
      colnames(W_invitro_suffled) <- colnames(W_invitro)
      Thm_train_shuffled <- x_train %*% W_invitro_suffled %*% t(W_invitro_suffled) %*% Wh
      Thm_val_shuffled <- x_val %*% W_invitro_suffled %*% t(W_invitro_suffled) %*% Wh
      y_hat_test <- cbind(1, Thm_val_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      y_hat_train <- cbind(1, Thm_train_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      train_r_shuffle_w[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_w[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_w[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_w[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      y_hat_test <- cbind(1, Thm_val)  %*% Bshuffled
      y_hat_train <- cbind(1, Thm_train)  %*% Bshuffled
      train_r_shuffle_bh[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_bh[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_bh[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_bh[j] <- mean(abs(y_hat_train-y_train))
      
    }else if (task=='human_backprojected_retrained'){
      y_hat_test <- cbind(1, Thm_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      y_hat_train <- cbind(1, Thm_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      train_r[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae[j] <- mean(abs(y_hat_test-y_val))
      train_mae[j] <- mean(abs(y_hat_train-y_train))
      
      ### Shuffle W_invitro
      W_invitro_suffled <- W_invitro[sample.int(nrow(W_invitro)),]
      rownames(W_invitro_suffled) <- rownames(W_invitro)
      colnames(W_invitro_suffled) <- colnames(W_invitro)
      Thm_train_shuffled <- x_train %*% W_invitro_suffled %*% t(W_invitro_suffled) %*% Wh
      Thm_val_shuffled <- x_val %*% W_invitro_suffled %*% t(W_invitro_suffled) %*% Wh
      y_hat_test <- cbind(1, Thm_val_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      y_hat_train <- cbind(1, Thm_train_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      train_r_shuffle_w[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_w[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_w[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_w[j] <- mean(abs(y_hat_train-y_train))
      
    }else if (task=='human_backprojected_into_translatable_lvs'){
      Wm_new <- get_translatable_LV(x_train, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = FALSE)
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae[j] <- mean(abs(y_hat_test-y_val))
      train_mae[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle W_invitro
      W_invitro_suffled <- W_invitro[sample.int(nrow(W_invitro)),]
      rownames(W_invitro_suffled) <- rownames(W_invitro)
      colnames(W_invitro_suffled) <- colnames(W_invitro)
      Wm_new <- get_translatable_LV(x_train, y_train, Wh, W_invitro_suffled,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = FALSE)
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_w[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_w[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_w[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_w[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle x_train
      x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
      colnames(x_train_shuffled) <- colnames(x_train)
      Wm_new <- get_translatable_LV(x_train_shuffled, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = FALSE)
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      val_r_shuffle_x[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_x[j] <- mean(abs(y_hat_test-y_val))
      
      ### shuffled labels model
      y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
      rownames(y_train_shuffled) <- rownames(y_train)
      Wm_new <- get_translatable_LV(x_train, y_train_shuffled, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = FALSE)
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      val_r_shuffle_y[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_y[j] <- mean(abs(y_hat_test-y_val))
      
      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      Wm_new <- get_translatable_LV(x_train, y_train, Wh, W_invitro,
                                    Bshuffled,
                                    find_extra = FALSE)
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% Bshuffled
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% Bshuffled
      train_r_shuffle_bh[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_bh[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_bh[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_bh[j] <- mean(abs(y_hat_train-y_train))
      
    }else if (task=='human_backprojected_into_optimal_lvs'){
      Wm_opt <- get_translatable_LV(x_train, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = TRUE)
      Wm_opt <- Wm_opt$Wm_new
      colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae[j] <- mean(abs(y_hat_test-y_val))
      train_mae[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle W_invitro
      W_invitro_suffled <- W_invitro[sample.int(nrow(W_invitro)),]
      rownames(W_invitro_suffled) <- rownames(W_invitro)
      colnames(W_invitro_suffled) <- colnames(W_invitro)
      Wm_opt <- get_translatable_LV(x_train, y_train, Wh, W_invitro_suffled,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = TRUE)
      Wm_opt <- Wm_opt$Wm_new
      colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
      # Extend latent variables
      Wm_tot <- cbind(W_invitro_suffled, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_w[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_w[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_w[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_w[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle x_train
      x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
      colnames(x_train_shuffled) <- colnames(x_train)
      Wm_opt <- get_translatable_LV(x_train_shuffled, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = TRUE)
      Wm_opt <- Wm_opt$Wm_new
      colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      val_r_shuffle_x[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_x[j] <- mean(abs(y_hat_test-y_val))
      
      ### shuffled labels model
      y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
      rownames(y_train_shuffled) <- rownames(y_train)
      Wm_opt <- get_translatable_LV(x_train, y_train_shuffled, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = TRUE)
      Wm_opt <- Wm_opt$Wm_new
      colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      val_r_shuffle_y[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_y[j] <- mean(abs(y_hat_test-y_val))

      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      Wm_opt <- get_translatable_LV(x_train, y_train, Wh, W_invitro,
                                    Bshuffled,
                                    find_extra = TRUE)
      Wm_opt <- Wm_opt$Wm_new
      colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% Bshuffled
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% Bshuffled
      train_r_shuffle_bh[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_bh[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_bh[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_bh[j] <- mean(abs(y_hat_train-y_train))
    }
    
  }
  if (task=='human_plsr'){
    performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Y',r = val_r_shuffle_y,MAE = val_mae_shuffle_y,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle X',r = val_r_shuffle_x,MAE = val_mae_shuffle_x,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='random X',r = val_r_random_x,MAE = val_mae_random_x,fold = seq(1,num_folds),task=task))
  }else if (task=='human_backprojected'){
    performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle W',r = train_r_shuffle_w,MAE = train_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle W',r = val_r_shuffle_w,MAE = val_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle Bh',r = train_r_shuffle_bh,MAE = train_mae_shuffle_bh,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Bh',r = val_r_shuffle_bh,MAE = val_mae_shuffle_bh,fold = seq(1,num_folds),task=task))
  }else if (task=='human_backprojected_retrained'){
    performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle W',r = train_r_shuffle_w,MAE = train_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle W',r = val_r_shuffle_w,MAE = val_mae_shuffle_w,fold = seq(1,num_folds),task=task))
  }else{
    #### also Add shuffle X in this approaches
    performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle W',r = train_r_shuffle_w,MAE = train_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle W',r = val_r_shuffle_w,MAE = val_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle Bh',r = train_r_shuffle_bh,MAE = train_mae_shuffle_bh,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Bh',r = val_r_shuffle_bh,MAE = val_mae_shuffle_bh,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Y',r = val_r_shuffle_y,MAE = val_mae_shuffle_y,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle X',r = val_r_shuffle_x,MAE = val_mae_shuffle_x,fold = seq(1,num_folds),task=task))
  }
  return(performance_df)
}