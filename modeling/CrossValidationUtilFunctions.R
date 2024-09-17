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
                                                        'analytical_optimal'),
                                               LVs=9){
  val_mae <- matrix(0,nrow = folds,ncol = 2)
  train_mae <- matrix(0,nrow = folds,ncol = 2)
  val_r <- matrix(0,nrow = folds,ncol = 2)
  train_r <- matrix(0,nrow = folds,ncol = 2)
  val_r <- matrix(0,nrow = folds,ncol = 2)
  all_models <- NULL
  val_r_shuffle_y <- matrix(0,nrow = folds,ncol = 2)
  val_mae_shuffle_y <- matrix(0,nrow = folds,ncol = 2)
  val_r_shuffle_x <- matrix(0,nrow = folds,ncol = 2)
  val_mae_shuffle_x <- matrix(0,nrow = folds,ncol = 2)
  val_r_random_x <- matrix(0,nrow = folds,ncol = 2)
  val_mae_random_x <- matrix(0,nrow = folds,ncol = 2)
  performance_df <- NULL
  train_r_shuffle_w <- matrix(0,nrow = folds,ncol = 2)
  val_r_shuffle_w <- matrix(0,nrow = folds,ncol = 2)
  val_mae_shuffle_w <- matrix(0,nrow = folds,ncol = 2)
  train_mae_shuffle_w <- matrix(0,nrow = folds,ncol = 2)
  train_r_shuffle_bh <- matrix(0,nrow = folds,ncol = 2)
  val_r_shuffle_bh <- matrix(0,nrow = folds,ncol = 2)
  val_mae_shuffle_bh <- matrix(0,nrow = folds,ncol = 2)
  train_mae_shuffle_bh <- matrix(0,nrow = folds,ncol = 2)
  # train_r_random_w <- matrix(0,nrow = folds,ncol = 2)
  # val_r_random_w <- matrix(0,nrow = folds,ncol = 2)
  # val_mae_random_w <- matrix(0,nrow = folds,ncol = 2)
  # train_mae_random_w <- matrix(0,nrow = folds,ncol = 2)
  train_r_shuffle_wopt <- matrix(0,nrow = folds,ncol = 2)
  val_r_shuffle_wopt <- matrix(0,nrow = folds,ncol = 2)
  val_mae_shuffle_wopt <- matrix(0,nrow = folds,ncol = 2)
  train_mae_shuffle_wopt <- matrix(0,nrow = folds,ncol = 2)
  val_r_scramble_x <- matrix(0,nrow = folds,ncol = 2)
  train_r_scramble_x <- matrix(0,nrow = folds,ncol = 2)
  val_mae_scramble_x <- matrix(0,nrow = folds,ncol = 2)
  train_mae_scramble_x <- matrix(0,nrow = folds,ncol = 2)
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
    train_r[j,] <- diag(cor(y_train_hat,y_train))
    val_r[j,] <- diag(cor(y_val_hat,y_val))
    val_mae[j,] <- apply(abs(y_val_hat-y_val),2,mean)
    train_mae[j,] <- apply(abs(y_train_hat-y_train),2,mean)

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
    val_r_shuffle_y[j,]<- diag(cor(y_hat_val,y_val))
    val_mae_shuffle_y[j,] <- apply(abs(y_hat_val-y_val),2,mean)
    
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
    val_r_shuffle_x[j,]<- diag(cor(y_hat_val,y_val))
    val_mae_shuffle_x[j,] <- apply(abs(y_hat_val-y_val),2,mean)
    
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
    val_r_random_x[j,]<- diag(cor(y_hat_test,y_val))
    val_mae_random_x[j,] <- apply(abs(y_hat_test-y_val),2,mean)
    
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
      train_r[j,] <- diag(cor(y_hat_train,y_train))
      val_r[j,] <- diag(cor(y_hat_test,y_val))
      val_mae[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
      ##shuffle W_invitro
      W_invitro_suffled <- W_invitro[sample.int(nrow(W_invitro)),]
      rownames(W_invitro_suffled) <- rownames(W_invitro)
      colnames(W_invitro_suffled) <- colnames(W_invitro)
      Thm_train_shuffled <- x_train %*% W_invitro_suffled %*% t(W_invitro_suffled) %*% Wh
      Thm_val_shuffled <- x_val %*% W_invitro_suffled %*% t(W_invitro_suffled) %*% Wh
      y_hat_test <- cbind(1, Thm_val_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      y_hat_train <- cbind(1, Thm_train_shuffled)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
      train_r_shuffle_w[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_w[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_w[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_w[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      y_hat_test <- cbind(1, Thm_val)  %*% Bshuffled
      y_hat_train <- cbind(1, Thm_train)  %*% Bshuffled
      train_r_shuffle_bh[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_bh[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_bh[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_bh[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
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
      train_r[j,] <- diag(cor(y_hat_train,y_train))
      val_r[j,] <- diag(cor(y_hat_test,y_val))
      val_mae[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
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
      train_r_shuffle_w[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_w[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_w[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_w[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
    }else if (task=='human_backprojected_into_translatable_lvs'){
      Wm_new <- get_translatable_LV_2phenotype(x_train, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh))
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r[j,] <- diag(cor(y_hat_train,y_train))
      val_r[j,] <- diag(cor(y_hat_test,y_val))
      val_mae[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
      ##shuffle W_invitro
      W_invitro_suffled <- W_invitro[sample.int(nrow(W_invitro)),]
      rownames(W_invitro_suffled) <- rownames(W_invitro)
      colnames(W_invitro_suffled) <- colnames(W_invitro)
      Wm_new <- get_translatable_LV_2phenotype(x_train, y_train, Wh, W_invitro_suffled,
                                    rbind(apply(y_train,2,mean),Bh))
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_w[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_w[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_w[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_w[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
      ##shuffle x_train
      x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
      colnames(x_train_shuffled) <- colnames(x_train)
      Wm_new <- get_translatable_LV_2phenotype(x_train_shuffled, y_train, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh))
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      val_r_shuffle_x[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_x[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      
      ### shuffled labels model
      y_train_shuffled <- y_train[sample.int(nrow(y_train)),]
      rownames(y_train_shuffled) <- rownames(y_train)
      Wm_new <- get_translatable_LV_2phenotype(x_train, y_train_shuffled, Wh, W_invitro,
                                    rbind(apply(y_train,2,mean),Bh))
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      val_r_shuffle_y[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_y[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      
      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      Bshuffled <- rbind(apply(y_train,2,mean),Bshuffled)
      Wm_new <- get_translatable_LV_2phenotype(x_train, y_train, Wh, W_invitro,
                                    Bshuffled)
      Wm_new <- Wm_new$Wm_new
      colnames(Wm_new) <- paste0(target_dataset,"_LVdata",1:ncol(Wm_new))
      y_hat_train <- cbind(1, x_train %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_new %*% t(Wm_new) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_bh[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_bh[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_bh[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_bh[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
    }else if (task=='analytical_optimal'){
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
      train_r[j,] <- diag(cor(y_hat_train,y_train))
      val_r[j,] <- diag(cor(y_hat_test,y_val))
      val_mae[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
      ##shuffle W_opt
      Wm_opt_suffled <- Wm_opt[sample.int(nrow(Wm_opt)),]
      rownames(Wm_opt_suffled) <- rownames(Wm_opt_suffled)
      colnames(Wm_opt_suffled) <- colnames(Wm_opt_suffled)
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt_suffled) 
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_wopt[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_wopt[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_wopt[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_wopt[j,] <- apply(abs(y_hat_train-y_train),2,mean)
      
      ##shuffle Bh
      Bshuffled <- Bh[sample.int(nrow(Bh)),]
      rownames(Bshuffled) <- rownames(Bh)
      colnames(Bshuffled) <- colnames(Bh)
      phi <- Wh %*% Bshuffled
      Wm_opt <- analytical_solution_opt(y=y_train,
                                        W_invitro = W_invitro,
                                        phi = phi)
      # Extend latent variables
      Wm_tot <- cbind(W_invitro, Wm_opt)
      # predict
      y_hat_train <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      y_hat_test <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
      train_r_shuffle_bh[j,] <- diag(cor(y_hat_train,y_train))
      val_r_shuffle_bh[j,] <- diag(cor(y_hat_test,y_val))
      val_mae_shuffle_bh[j,] <- apply(abs(y_hat_test-y_val),2,mean)
      train_mae_shuffle_bh[j,] <- apply(abs(y_hat_train-y_train),2,mean)
    }
    
  }
  if (task=='human_plsr'){
    performance_df_r <- rbind(data.frame(set='train',type='model',NAS = train_r[,1],fibrosis = train_r[,2],fold = seq(1,num_folds),task=task) %>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                            data.frame(set='test',type='model',NAS = val_r[,1],fibrosis = val_r[,2],fold = seq(1,num_folds),task=task)%>%
                              gather('phenotype','value',-set,-type,-fold,-task),
                            data.frame(set='test',type='shuffle Y',NAS = val_r_shuffle_y[,1],fibrosis = val_r_shuffle_y[,2],fold = seq(1,num_folds),task=task)%>%
                              gather('phenotype','value',-set,-type,-fold,-task),
                            data.frame(set='test',type='shuffle X',NAS = val_r_shuffle_x[,1],fibrosis = val_r_shuffle_x[,2],fold = seq(1,num_folds),task=task)%>%
                              gather('phenotype','value',-set,-type,-fold,-task),
                            data.frame(set='test',type='random X',NAS = val_r_random_x[,1],fibrosis = val_r_random_x[,2],fold = seq(1,num_folds),task=task)%>%
                              gather('phenotype','value',-set,-type,-fold,-task))
    performance_df_mae <- rbind(data.frame(set='train',type='model',NAS = train_mae[,1],fibrosis = train_mae[,2],fold = seq(1,num_folds),task=task) %>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='model',NAS = val_mae[,1],fibrosis = val_mae[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle Y',NAS = val_mae_shuffle_y[,1],fibrosis = val_mae_shuffle_y[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle X',NAS = val_mae_shuffle_x[,1],fibrosis = val_mae_shuffle_x[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='random X',NAS = val_mae_random_x[,1],fibrosis = val_mae_random_x[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task))
    
    performance_df <- rbind(performance_df_r %>% mutate(metric = 'r'),performance_df_mae %>% mutate(metric = 'MAE'))
  }else if (task=='human_backprojected'){
    performance_df_r <- rbind(data.frame(set='train',type='model',NAS = train_r[,1],fibrosis = train_r[,2],fold = seq(1,num_folds),task=task) %>%
                               gather('phenotype','value',-set,-type,-fold,-task),
                             data.frame(set='test',type='model',NAS = val_r[,1],fibrosis = val_r[,2],fold = seq(1,num_folds),task=task)%>%
                               gather('phenotype','value',-set,-type,-fold,-task),
                             data.frame(set='train',type='shuffle W',NAS = train_r_shuffle_w[,1],fibrosis = train_r_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                               gather('phenotype','value',-set,-type,-fold,-task),
                             data.frame(set='test',type='shuffle W',NAS = val_r_shuffle_w[,1],fibrosis = val_r_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                               gather('phenotype','value',-set,-type,-fold,-task),
                             data.frame(set='train',type='shuffle Bh',NAS = train_r_shuffle_bh[,1],fibrosis = train_r_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                               gather('phenotype','value',-set,-type,-fold,-task),
                             data.frame(set='test',type='shuffle Bh',NAS = val_r_shuffle_bh[,1],fibrosis = val_r_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                               gather('phenotype','value',-set,-type,-fold,-task))
    performance_df_mae <- rbind(data.frame(set='train',type='model',NAS = train_mae[,1],fibrosis = train_mae[,2],fold = seq(1,num_folds),task=task) %>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='model',NAS = val_mae[,1],fibrosis = val_mae[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle W',NAS = train_mae_shuffle_w[,1],fibrosis = train_mae_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle W',NAS = val_mae_shuffle_w[,1],fibrosis = val_mae_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle Bh',NAS = train_mae_shuffle_bh[,1],fibrosis = train_mae_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle Bh',NAS = val_mae_shuffle_bh[,1],fibrosis = val_mae_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task))
    performance_df <- rbind(performance_df_r %>% mutate(metric = 'r'),performance_df_mae %>% mutate(metric = 'MAE'))
  }else if (task=='human_backprojected_retrained'){
    performance_df_r <- rbind(data.frame(set='train',type='model',NAS = train_r[,1],fibrosis = train_r[,2],fold = seq(1,num_folds),task=task) %>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='model',NAS = val_r[,1],fibrosis = val_r[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='train',type='shuffle W',NAS = train_r_shuffle_w[,1],fibrosis = train_r_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='shuffle W',NAS = val_r_shuffle_w[,1],fibrosis = val_r_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task))
    performance_df_mae <- rbind(data.frame(set='train',type='model',NAS = train_mae[,1],fibrosis = train_mae[,2],fold = seq(1,num_folds),task=task) %>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='model',NAS = val_mae[,1],fibrosis = val_mae[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle W',NAS = train_mae_shuffle_w[,1],fibrosis = train_mae_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle W',NAS = val_mae_shuffle_w[,1],fibrosis = val_mae_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task))
    performance_df <- rbind(performance_df_r %>% mutate(metric = 'r'),performance_df_mae %>% mutate(metric = 'MAE'))
  }else if (task=='human_backprojected_into_translatable_lvs'){
    performance_df_r <- rbind(data.frame(set='train',type='model',NAS = train_r[,1],fibrosis = train_r[,2],fold = seq(1,num_folds),task=task) %>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='model',NAS = val_r[,1],fibrosis = val_r[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='train',type='shuffle W',NAS = train_r_shuffle_w[,1],fibrosis = train_r_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='shuffle W',NAS = val_r_shuffle_w[,1],fibrosis = val_r_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='train',type='shuffle Bh',NAS = train_r_shuffle_bh[,1],fibrosis = train_r_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='shuffle Bh',NAS = val_r_shuffle_bh[,1],fibrosis = val_r_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task))
    performance_df_mae <- rbind(data.frame(set='train',type='model',NAS = train_mae[,1],fibrosis = train_mae[,2],fold = seq(1,num_folds),task=task) %>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='model',NAS = val_mae[,1],fibrosis = val_mae[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle W',NAS = train_mae_shuffle_w[,1],fibrosis = train_mae_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle W',NAS = val_mae_shuffle_w[,1],fibrosis = val_mae_shuffle_w[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle Bh',NAS = train_mae_shuffle_bh[,1],fibrosis = train_mae_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle Bh',NAS = val_mae_shuffle_bh[,1],fibrosis = val_mae_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task))
    performance_df <- rbind(performance_df_r %>% mutate(metric = 'r'),performance_df_mae %>% mutate(metric = 'MAE'))
  }else{
    performance_df_r <- rbind(data.frame(set='train',type='model',NAS = train_r[,1],fibrosis = train_r[,2],fold = seq(1,num_folds),task=task) %>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='model',NAS = val_r[,1],fibrosis = val_r[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='train',type='shuffle Wopt',NAS = train_r_shuffle_wopt[,1],fibrosis = train_r_shuffle_wopt[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='shuffle Wopt',NAS = val_r_shuffle_wopt[,1],fibrosis = val_r_shuffle_wopt[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='train',type='shuffle Bh',NAS = train_r_shuffle_bh[,1],fibrosis = train_r_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task),
                              data.frame(set='test',type='shuffle Bh',NAS = val_r_shuffle_bh[,1],fibrosis = val_r_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                gather('phenotype','value',-set,-type,-fold,-task))
    performance_df_mae <- rbind(data.frame(set='train',type='model',NAS = train_mae[,1],fibrosis = train_mae[,2],fold = seq(1,num_folds),task=task) %>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='model',NAS = val_mae[,1],fibrosis = val_mae[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle Wopt',NAS = train_mae_shuffle_wopt[,1],fibrosis = train_mae_shuffle_wopt[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle Wopt',NAS = val_mae_shuffle_wopt[,1],fibrosis = val_mae_shuffle_wopt[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='train',type='shuffle Bh',NAS = train_mae_shuffle_bh[,1],fibrosis = train_mae_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task),
                                data.frame(set='test',type='shuffle Bh',NAS = val_mae_shuffle_bh[,1],fibrosis = val_mae_shuffle_bh[,2],fold = seq(1,num_folds),task=task)%>%
                                  gather('phenotype','value',-set,-type,-fold,-task))
    performance_df <- rbind(performance_df_r %>% mutate(metric = 'r'),performance_df_mae %>% mutate(metric = 'MAE'))
  }
  return(performance_df)
}
