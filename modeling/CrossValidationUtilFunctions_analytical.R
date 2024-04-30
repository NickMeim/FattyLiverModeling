library(tidyverse)
library(caret)
library(ropls)
cross_validation_complete_pipeline_2 <- function(W_invitro,
                                               folds,
                                               file_loc,
                                               target_dataset,
                                               task = c('human_plsr',
                                                        'human_backprojected',
                                                        'human_backprojected_retrained',
                                                        'human_backprojected_into_translatable_lvs',
                                                        'human_backprojected_into_optimal_lvs',
                                                        'analytical'),
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
  # train_r_random_w <- NULL
  # val_r_random_w <- NULL
  # val_mae_random_w <- NULL
  # train_mae_random_w <- NULL
  train_r_shuffle_wopt <- NULL
  val_r_shuffle_wopt <- NULL
  val_mae_shuffle_wopt <- NULL
  train_mae_shuffle_wopt <- NULL
  val_r_scramble_x <- NULL
  train_r_scramble_x <- NULL
  val_mae_scramble_x <- NULL
  train_mae_scramble_x <- NULL
  for (j in 1:folds){
    message(paste0('Begun fold ',j))
    x_train <- readRDS(paste0(file_loc,'Xh_train',j,'.rds'))
    y_train <- readRDS(paste0(file_loc,'Yh_train',j,'.rds'))
    x_val <- readRDS(paste0(file_loc,'Xh_val',j,'.rds'))
    y_val <- readRDS(paste0(file_loc,'Yh_val',j,'.rds'))
    
    # plsr_model_og <- opls(x = x_train, 
    #                    y = y_train,
    #                    predI = LVs,
    #                    crossvalI = 1,
    #                    scaleC = "center",
    #                    fig.pdfC = "none",
    #                    info.txtC = "none")
    # # Get Wh of PLSR
    # Wh <- matrix(data = 0, ncol = ncol(plsr_model_og@weightMN), nrow = ncol(Xh))
    # rownames(Wh) <- colnames(x_train)
    # colnames(Wh) <- colnames(plsr_model_og@weightMN)
    # for (ii in 1:nrow(plsr_model_og@weightMN)){
    #   Wh[rownames(plsr_model_og@weightMN)[ii], ] <- plsr_model_og@weightMN[ii,]
    # }
    # # Get regression coefficients
    # # Bh <- matrix(lm(y_train ~ x_train %*% Wh) %>% coef(), ncol = 2)
    # Bh <- t(plsr_model_og@weightMN) %*% plsr_model_og@coefficientMN
    # 
    # # Define projection matrices to make more readable
    # Th_train <- x_train %*% Wh
    # Thm_train <- x_train %*% W_invitro %*% t(W_invitro) %*% Wh
    # Th_val<- x_val %*% Wh
    # Thm_val <- x_val %*% W_invitro %*% t(W_invitro) %*% Wh
    # 
    # 
    # ## Backprojected invitro
    # y_train_back <- cbind(1, Thm_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model_og@weightMN) %*% plsr_model_og@coefficientMN)
    
    pca_humans <- prcomp(x_train,scale. = FALSE,center=TRUE)
    Wpch <- pca_humans$rotation
    M <- t(W_invitro) %*% Wpch
    Anull <- MASS::Null(t(M))
    Nh <- Wpch %*% Anull
    P <- Nh %*%  matlib::inv(t(Nh) %*% Nh) %*% t(Nh)
    z_train <- x_train %*% P
    z_val <-  x_val %*% P
    
    plsr_model <- opls(x = z_train, 
                       y = y_train,
                       predI = LVs,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    # # Get Wh of PLSR
    # Wopt <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
    # rownames(Wopt) <- colnames(x_train)
    # colnames(Wopt) <- colnames(plsr_model@weightMN)
    # for (ii in 1:nrow(plsr_model@weightMN)){
    #   Wopt[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    # }
    # # Get regression coefficients
    # # Bh <- matrix(lm(y_train ~ x_train %*% Wh) %>% coef(), ncol = 2)
    # Bh_null <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
    # 
    # Wm_tot <- cbind(W_invitro, Wopt)
    # # predict
    # y_train_hat <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    # y_val_hat <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    y_train_hat <- predict(plsr_model,z_train)
    y_val_hat <- predict(plsr_model,z_val)
    train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
    val_mae[j] <- mean(abs(y_val_hat-y_val))
    train_mae[j] <- mean(abs(y_train_hat-y_train))
    
    ### shuffled labels model
    y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
    rownames(y_train_shuffled) <- rownames(y_train)
    plsr_model_shuffle_y <- opls(x = z_train, 
                                 y = y_train_shuffled,
                                 predI = LVs,
                                 crossvalI = 1,
                                 scaleC = "center",
                                 fig.pdfC = "none",
                                 info.txtC = "none")
    # # Get Wh of PLSR
    # Wopt <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
    # rownames(Wopt) <- colnames(x_train)
    # colnames(Wopt) <- colnames(plsr_model@weightMN)
    # for (ii in 1:nrow(plsr_model@weightMN)){
    #   Wopt[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    # }
    # # Get regression coefficients
    # # Bh <- matrix(lm(y_train ~ x_train %*% Wh) %>% coef(), ncol = 2)
    # Bh_null <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
    # Wm_tot <- cbind(W_invitro, Wopt)
    # # predict
    # y_hat_val <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    y_val_hat <- predict(plsr_model_shuffle_y,z_val)
    val_r_shuffle_y[j]<- mean(diag(cor(y_val_hat,y_val)))
    val_mae_shuffle_y[j] <- mean(abs(y_val_hat-y_val))
    
    ### shuffled features model
    x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
    colnames(x_train_shuffled) <- colnames(x_train)
    pca_humans_shuffled <- prcomp(x_train_shuffled,scale. = FALSE,center=TRUE)
    Wpch_shuffled <- pca_humans_shuffled$rotation
    M <- t(W_invitro) %*% Wpch_shuffled
    Anull <- MASS::Null(t(M))
    Nh <- Wpch_shuffled %*% Anull
    P <- Nh %*%  matlib::inv(t(Nh) %*% Nh) %*% t(Nh)
    z_train <- x_train_shuffled %*% P
    z_val <- x_val %*% P
    plsr_model_shuffle_x <- opls(x = z_train, 
                                 y = y_train,
                                 predI = 9,
                                 crossvalI = 1,
                                 scaleC = "center",
                                 fig.pdfC = "none",
                                 info.txtC = "none")
    # # Get Wh of PLSR
    # Wopt <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
    # rownames(Wopt) <- colnames(x_train)
    # colnames(Wopt) <- colnames(plsr_model@weightMN)
    # for (ii in 1:nrow(plsr_model@weightMN)){
    #   Wopt[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    # }
    # # Get regression coefficients
    # # Bh <- matrix(lm(y_train ~ x_train %*% Wh) %>% coef(), ncol = 2)
    # Bh_null <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
    # Wm_tot <- cbind(W_invitro, Wopt)
    # predict
    # y_hat_val <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    y_val_hat <- predict(plsr_model_shuffle_x,z_val)
    val_r_shuffle_x[j]<- mean(diag(cor(y_val_hat,y_val)))
    val_mae_shuffle_x[j] <- mean(abs(y_val_hat-y_val))
    
    
    print(paste0('Finished RREF approach ',j))
    
  }
  performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                          data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                          data.frame(set='test',type='shuffle Y',r = val_r_shuffle_y,MAE = val_mae_shuffle_y,fold = seq(1,num_folds),task=task),
                          data.frame(set='test',type='shuffle X',r = val_r_shuffle_x,MAE = val_mae_shuffle_x,fold = seq(1,num_folds),task=task))
  return(performance_df)
}