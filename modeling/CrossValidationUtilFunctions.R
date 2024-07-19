library(tidyverse)
library(caret)
library(ropls)

analytical_solution_opt <- function(y,W_invitro,phi){
  Wopt <- matrix(0,nrow = nrow(W_invitro),ncol=ncol(y))
  for (i in 1:ncol(y)){
    if (i == 1){
      alpha <- t(W_invitro) %*% phi[,i]
      Wopt[,i] <- (phi[,i] - W_invitro %*% alpha)/sqrt(sum(phi[,i]^2) -sum(alpha^2))
    }else{
      Wnew <- cbind(W_invitro,Wopt[,1:(i-1)])
      alpha <- t(Wnew) %*% phi[,i]
      Wopt[,i] <- (phi[,i] - Wnew %*% alpha)/sqrt(sum(phi[,i]^2) -sum(alpha^2))
    }
  }
  return(Wopt)
}

cross_validation_complete_pipeline <- function(W_invitro,
                                               folds,
                                               file_loc,
                                               target_dataset,
                                               task = c('human_plsr',
                                                        'human_backprojected',
                                                        'human_backprojected_retrained',
                                                        'human_backprojected_into_translatable_lvs',
                                                        'human_backprojected_into_optimal_lvs',
                                                        'analytical_optimal'),
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
                                 predI = LVs,
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
                                 predI = LVs,
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
                                predI = LVs,
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
      # Retrain model to overfit projected data
      y_hat_train <- matrix(0,nrow = nrow(y_train),ncol = ncol(y_train))
      y_hat_test <- matrix(0,nrow = nrow(y_val),ncol = ncol(y_val))
      for (i in 1:ncol(y_train)){
        train_data <- cbind(y_train[,i], Thm_train)
        colnames(train_data)[1] <- 'V1'
        mod_retrain <- lm(V1 ~., data = as.data.frame(train_data))
        y_hat_train[,i] <- predict(mod_retrain)
        y_hat_test[,i] <- predict(mod_retrain,newdata = as.data.frame(Thm_val))
      }
      # y_hat_test <- cbind(1, Thm_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      # y_hat_train <- cbind(1, Thm_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
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
      # Retrain model to overfit projected data
      y_hat_train <- matrix(0,nrow = nrow(y_train),ncol = ncol(y_train))
      y_hat_test <- matrix(0,nrow = nrow(y_val),ncol = ncol(y_val))
      for (i in 1:ncol(y_train)){
        train_data <- cbind(y_train[,i], Thm_train_shuffled)
        colnames(train_data)[1] <- 'V1'
        mod_retrain <- lm(V1 ~., data = as.data.frame(train_data))
        y_hat_train[,i] <- predict(mod_retrain)
        y_hat_test[,i] <- predict(mod_retrain,newdata = as.data.frame(Thm_val))
      }
      # y_hat_test <- cbind(1, Thm_val_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      # y_hat_train <- cbind(1, Thm_train_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
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
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_bh[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_bh[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_bh[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_bh[j] <- mean(abs(y_hat_train-y_train))
      
    }else if (task=='human_backprojected_into_optimal_lvs'){
      Wm_opt <- get_translatable_LV(x_train, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh),
                                    find_extra = TRUE,
                                    verbose = FALSE)
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
      
      ##shuffle W_opt
      Wm_opt_suffled <- Wm_opt[sample.int(nrow(Wm_opt)),]
      rownames(Wm_opt_suffled) <- rownames(Wm_opt_suffled)
      colnames(Wm_opt_suffled) <- colnames(Wm_opt_suffled)
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt_suffled) 
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_wopt[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_wopt[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_wopt[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_wopt[j] <- mean(abs(y_hat_train-y_train))

      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      Wm_opt <- get_translatable_LV(x_train, y_train, Wh, W_invitro,
                                    Bshuffled,
                                    find_extra = TRUE,
                                    verbose = FALSE)
      Wm_opt <- Wm_opt$Wm_new
      colnames(Wm_opt) <- paste0(target_dataset,"_LVopt",1:ncol(Wm_opt))
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_bh[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_bh[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_bh[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_bh[j] <- mean(abs(y_hat_train-y_train))
    } else if (task=='analytical_optimal'){
      # phi <- cbind(1, Wh)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      phi <- Wh %*% Bh
      # phi <- phi/sqrt(apply(phi^2,2,sum))
      Wm_opt <- analytical_solution_opt(y=y_train,
                                        W_invitro = W_invitro,
                                        phi = phi)
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae[j] <- mean(abs(y_hat_test-y_val))
      train_mae[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle W_opt
      Wm_opt_suffled <- Wm_opt[sample.int(nrow(Wm_opt)),]
      rownames(Wm_opt_suffled) <- rownames(Wm_opt_suffled)
      colnames(Wm_opt_suffled) <- colnames(Wm_opt_suffled)
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt_suffled) 
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_wopt[j] <- mean(diag(cor(y_hat_train,y_train)))
      val_r_shuffle_wopt[j] <- mean(diag(cor(y_hat_test,y_val)))
      val_mae_shuffle_wopt[j] <- mean(abs(y_hat_test-y_val))
      train_mae_shuffle_wopt[j] <- mean(abs(y_hat_train-y_train))
      
      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      # Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      phi <- Wh %*% Bshuffled
      # phi <- phi/sqrt(apply(phi^2,2,sum))
      Wm_opt <- analytical_solution_opt(y=y_train,
                                        W_invitro = W_invitro,
                                        phi = phi)
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
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
  }else if (task=='human_backprojected_into_translatable_lvs'){
    #### also Add shuffle X in this approaches
    performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle W',r = train_r_shuffle_w,MAE = train_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle W',r = val_r_shuffle_w,MAE = val_mae_shuffle_w,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle Bh',r = train_r_shuffle_bh,MAE = train_mae_shuffle_bh,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Bh',r = val_r_shuffle_bh,MAE = val_mae_shuffle_bh,fold = seq(1,num_folds),task=task))
                            # data.frame(set='test',type='shuffle Y',r = val_r_shuffle_y,MAE = val_mae_shuffle_y,fold = seq(1,num_folds),task=task),
                            # data.frame(set='test',type='shuffle X',r = val_r_shuffle_x,MAE = val_mae_shuffle_x,fold = seq(1,num_folds),task=task)
  }else{
    performance_df <- rbind(data.frame(set='train',type='model',r = train_r,MAE = train_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='model',r = val_r,MAE = val_mae,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle Wopt',r = train_r_shuffle_wopt,MAE = train_mae_shuffle_wopt,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Wopt',r = val_r_shuffle_wopt,MAE = val_mae_shuffle_wopt,fold = seq(1,num_folds),task=task),
                            data.frame(set='train',type='shuffle Bh',r = train_r_shuffle_bh,MAE = train_mae_shuffle_bh,fold = seq(1,num_folds),task=task),
                            data.frame(set='test',type='shuffle Bh',r = val_r_shuffle_bh,MAE = val_mae_shuffle_bh,fold = seq(1,num_folds),task=task))
                            # data.frame(set='train',type='random W',r = train_r_random_w,MAE = train_mae_random_w,fold = seq(1,num_folds),task=task),
                            # data.frame(set='test',type='random W',r = val_r_random_w,MAE = val_mae_random_w,fold = seq(1,num_folds),task=task),
  }
  return(performance_df)
}


GeneralMLCompletePipeLineCrossVal <- function(model,cv_location,external_datasets,data_list_all,
                                              Wm,Wtot,
                                              pheno,
                                              num_folds=10,
                                              dim_reduction=TRUE,
                                              shuffling=NULL){
  # Initialize some further parameters
  if (model=='lasso'){
    lambda_grid <- c(0,10^seq(-5, -1, length = 1000),1)
    alpha_grid <- 1
  }
  if (model=='ridge'){
    lambda_grid <- 10^seq(-3, 1, length = 1000)
    alpha_grid <- 0
  }
  if (model=='elasticnet'){
    lambda_grid <- 10^seq(-3, 1, length = 100)
    alpha_grid <- seq(0.01,0.99,0.05)
  }
  all_models <- NULL
  val_r_extrl <- matrix(0,nrow = num_folds,ncol = length(external_datasets))
  val_r_extrl_backproj <- matrix(0,nrow = num_folds,ncol = length(external_datasets))
  val_r_extrl_extra <- matrix(0,nrow = num_folds,ncol = length(external_datasets))
  val_r <- NULL
  train_r <- NULL
  val_r_backproj <- NULL
  train_r_backproj <- NULL
  val_r_extra <- NULL
  train_r_extra <- NULL
  for (i in 1:num_folds){
    x_train <- readRDS(paste0(cv_location,'Xh_train',i,'.rds'))
    y_train <- readRDS(paste0(cv_location,'Yh_train',i,'.rds'))
    x_val <- readRDS(paste0(cv_location,'Xh_val',i,'.rds'))
    y_val <- readRDS(paste0(cv_location,'Yh_val',i,'.rds'))
    
    if (dim_reduction==TRUE){
      x_train_og <- x_train
      x_val_og <- x_val
      pca_train <- prcomp(x_train,center = TRUE,scale. = FALSE)
      x_train <- pca_train$x
      x_val <- x_val %*% pca_train$rotation
      
    }
    
    if (!is.null(shuffling)){
      if (shuffling=='X'){
        message('A randomized X model will be trained')
        x_train <-x_train[,sample.int(ncol(x_train))]
        colnames(x_train) <- colnames(x_val)
      }else{
        message('A randomized Y model will be trained')
        y_train <-y_train[sample.int(nrow(y_train)),]
        rownames(y_train) <- rownames(x_train)
      }
    }
    train_data_new <- as.data.frame(cbind(x_train,y_train))
    train_data_new <- train_data_new[,c(colnames(x_train),pheno)]
    # colnames(train_data_new)[ncol(train_data_new)] <- 'out'
    # build new model
    ctrl <- trainControl(method = "cv", number = 10)
    if (model %in% c('lasso','ridge','elasticnet')){
      form <- as.formula(paste0(pheno,'~.'))
      mdl <- train(form,
                   data = train_data_new,
                   method = 'glmnet', 
                   metric = "MAE",
                   trControl = ctrl,
                   tuneGrid = expand.grid(alpha = alpha_grid, lambda = lambda_grid))
      
    }else{
      form <- as.formula(paste0(pheno,'~.'))
      invisible(capture.output(mdl <- train(form, 
                                            data = train_data_new, 
                                            method = model, 
                                            metric = "MAE",
                                            trControl = ctrl)))
    }
    all_models[[i]] <- mdl
    message('Finished training model')
    
    j <- 1
    for (dataset in external_datasets){
      # message(paste0('Begun ',dataset ,' dataset'))
      X <- data_list_all[[dataset]]$data_center %>% t()
      if (dataset == 'Hoang'){
        Y <- as.matrix(data_list_all[[dataset]]$metadata  %>% select(nafld_activity_score,Fibrosis_stage))
        
      }else if (dataset == 'Pantano'){
        Y <- as.matrix(data_list_all[[dataset]]$metadata  %>% select(`nas score`,`fibrosis stage`))
        Y <- apply(Y,c(1,2),as.numeric)
      }
      colnames(Y) <- c('NAS','fibrosis')
      Y <- Y[,pheno]
      
      ## PLSR original
      if (dim_reduction==TRUE){
        X2 <- X %*% pca_train$rotation
        y_val_hat_extrl <- predict(mdl,newdata = X2)
      }else{
        y_val_hat_extrl <- predict(mdl,newdata = X)
      }
      val_r_extrl[i,j] <- mean(diag(as.matrix(cor(y_val_hat_extrl,Y))))
      
      ## PLSR back-projected
      if (dim_reduction==TRUE){
        X2 <- X %*% Wm %*% t(Wm) %*% pca_train$rotation
      }else{
        X2 <- X %*% Wm %*% t(Wm)
      }
      y_val_hat_extrl <- predict(mdl,newdata = X2)
      val_r_extrl_backproj[i,j] <- mean(diag(as.matrix(cor(y_val_hat_extrl,Y))))
      
      ## PLSR back-projected with extra basis
      if (dim_reduction==TRUE){
        X2 <- X %*% Wtot %*% t(Wtot) %*% pca_train$rotation
      }else{
        X2 <- X %*% Wtot %*% t(Wtot)
      }
      y_val_hat_extrl <- predict(mdl,newdata = X2)
      val_r_extrl_extra[i,j] <- mean(diag(as.matrix(cor(y_val_hat_extrl,Y))))
      
      j <- j + 1
    }
    message('Finished evaluating model in external datasets')
    
    y_train <- y_train[,pheno]
    y_val <- y_val[,pheno]
    
    y_hat_val <- predict(mdl,newdata = x_val)
    y_hat_train <- predict(mdl,newdata = x_train)
    val_r[i] <- mean(diag(as.matrix(cor(y_hat_val,y_val))))
    train_r[i] <- mean(diag(as.matrix(cor(y_hat_train,y_train))))
    
    ## evaluate what the model does with back-projecting
    if (dim_reduction==TRUE){
      x_train_2 <- x_train_og %*% Wm %*% t(Wm) %*% pca_train$rotation
      x_val_2 <- x_val_og %*% Wm %*% t(Wm) %*% pca_train$rotation
    }else{
      x_train_2 <- x_train %*% Wm %*% t(Wm)
      x_val_2 <- x_val %*% Wm %*% t(Wm)
    }
    y_hat_val <- predict(mdl,newdata = x_val_2)
    y_hat_train <- predict(mdl,newdata = x_train_2)
    val_r_backproj[i] <- mean(diag(as.matrix(cor(y_hat_val,y_val))))
    train_r_backproj[i] <- mean(diag(as.matrix(cor(y_hat_train,y_train))))
    
    ## evaluate what the model does with back-projecting with extra basis
    if (dim_reduction==TRUE){
      x_train_2 <- x_train_og %*% Wtot %*% t(Wtot) %*% pca_train$rotation
      x_val_2 <- x_val_og %*% Wtot %*% t(Wtot) %*% pca_train$rotation
    }else{
      x_train_2 <- x_train %*% Wtot %*% t(Wtot)
      x_val_2 <- x_val %*% Wtot %*% t(Wtot)
    }
    y_hat_val <- predict(mdl,newdata = x_val_2)
    y_hat_train <- predict(mdl,newdata = x_train_2)
    val_r_extra[i] <- mean(diag(as.matrix(cor(y_hat_val,y_val))))
    train_r_extra[i] <- mean(diag(as.matrix(cor(y_hat_train,y_train))))
    
    message('Finished evaluating model in CV')
    print(paste0('Finished fold ',i,' out of ',num_folds))
  }
  res_val <- data.frame(model=val_r,'back-projected'=val_r_backproj,'optimized MPS'=val_r_extra) %>% mutate(set = 'test')
  colnames(res_val)[1] <- model
  res_train <- data.frame(model=train_r,'back-projected'=train_r_backproj,'optimized MPS'=train_r_extra) %>% mutate(set = 'train')
  colnames(res_train)[1] <- model
  res_cv <- rbind(res_train,res_val) %>% gather('model','r',-set)
  colnames(val_r_extrl) <- external_datasets
  colnames(val_r_extrl_backproj) <- external_datasets
  colnames(val_r_extrl_extra) <- external_datasets
  val_r_extrl <- as.data.frame(val_r_extrl) %>% mutate(model = model)
  val_r_extrl_backproj <- as.data.frame(val_r_extrl_backproj) %>% mutate(model='back-projected')
  val_r_extrl_extra <- as.data.frame(val_r_extrl_extra) %>% mutate(model = 'optimized MPS')
  res_external <- rbind(val_r_extrl,val_r_extrl_backproj,val_r_extrl_extra) %>% gather('dataset','r',-model)
  results <- list(all_models=all_models,cv_resuls=res_cv,external_results = res_external)
  return(results)
}
