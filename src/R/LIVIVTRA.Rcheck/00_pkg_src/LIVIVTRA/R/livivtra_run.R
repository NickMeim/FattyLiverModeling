library(tidyverse)
library(ropls)
source('R/functions_for_translation.R')

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

livivtra_run <- function(X_invivo,Y_invivo,X_invitro,W_invitro,
                         LV_number=8){

  message(paste0('Fit PLSR model with ',LV_number,' latent variables...'))
  plsr_model <- opls(x = X_invivo,
                     y = Y_invivo,
                     predI = LV_number,
                     crossvalI = 1,
                     scaleC = "center",
                     fig.pdfC = "none",
                     info.txtC = "none")
  # Get Wh of PLSR
  Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(X_invivo))
  rownames(Wh) <- colnames(X_invivo)
  colnames(Wh) <- colnames(plsr_model@weightMN)
  for (ii in 1:nrow(plsr_model@weightMN)){
    Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
  }
  # Get regression coefficients
  Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN

  message('Calculating extra basis to incease the clinical predictive power...')
  # Find extra basis
  phi <- Wh %*% Bh
  Wm_opt <- analytical_solution_opt(y=Y_invivo,
                                    W_invitro = W_invitro,
                                    phi = phi)
  Wm_tot <- cbind(W_invitro, Wm_opt)

  message('Calculating translatable components that contain most of the current predictive power...')
  ### Find translatable LV of the in vitro system
  ### Run evolutionary algorithm
  Wm_combo <- get_translatable_LV(X_invivo, Y_invivo, Wh, W_invitro,
                                  rbind(apply(Y_invivo,2,mean),Bh),
                                  find_extra = FALSE,
                                  verbose = TRUE)
  Wm_combo <- Wm_combo$Wm_new

  rownames(Wm_tot) <- rownames(W_invitro)
  rownames(Wm_opt) <- rownames(W_invitro)
  rownames(Wm_combo) <- rownames(W_invitro)
  rownames(Wh) <- rownames(W_invitro)

  # evaluate performance of plsr model, and truncation with W_invitro , with
  # Wm_combo, and W_opt
  Th <- X_invivo %*% Wm_tot %*% t(Wm_tot) %*% Wh
  Th0 <- X_invivo %*% W_invitro %*% t(W_invitro) %*% Wh
  yhat <- predict(plsr_model,X_invivo)
  yhat_extra_basis <- cbind(1, Th)  %*% rbind(apply(Y_invivo,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
  yhat_truncated <-  cbind(1, Th0)  %*% rbind(apply(Y_invivo,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
  if (ncol(Y_invivo)>1){
    r_model <- mean(diag(cor(yhat,Y_invivo)))
    r_extra_basis <- mean(diag(cor(yhat_extra_basis,Y_invivo)))
    r_truncated <- mean(diag(cor(yhat_truncated,Y_invivo)))
    mae_model <- apply(abs(yhat-Y_invivo),2,mean)
    mae_extra_basis <- apply(abs(yhat_extra_basis-Y_invivo),2,mean)
    mae_truncated <- apply(abs(yhat_truncated-Y_invivo),2,mean)
  }else{
    r_model <- mean(cor(yhat,Y_invivo))
    r_extra_basis <- mean(cor(yhat_extra_basis,Y_invivo))
    r_truncated <- mean(cor(yhat_truncated,Y_invivo))
    mae_model <- mean(abs(yhat-Y_invivo))
    mae_extra_basis <- mean(abs(yhat_extra_basis-Y_invivo))
    mae_truncated <- mean(abs(yhat_truncated-Y_invivo))
  }

  res_df = data.frame(r = c(r_model,r_extra_basis,r_truncated),
                      MAE = c(mae_model,mae_extra_basis,mae_truncated),
                      input = c('human features',
                                'human features truncated via the original in vitro model',
                                'human features truncated via the optimized in vitro model'),
                      model = rep('PLSR',3))

  message('Done!')
  # put all results into a list object
  model_list = list(model = plsr_model,
                    W_invitro=W_invitro,Wh=Wh,
                    W_opt=Wm_opt,W_translatable = Wm_combo,W_tot = Wm_tot,
                    performance = res_df)
  return(model_list)
}

livivtra_run_CV <- function(X_invivo,Y_invivo,X_invitro,W_invitro,
                         num_folds,LV_number=8){

  data_splits <- createMultiFolds(y = Y_invivo[,1], k = num_folds, times = 1)
  train_r_model <- NULL
  train_r_extra_basis <- NULL
  train_r_truncated <- NULL
  train_mae_model <- NULL
  train_mae_extra_basis <- NULL
  train_mae_truncated <- NULL
  val_r_model <- NULL
  val_r_extra_basis <- NULL
  val_r_truncated <- NULL
  val_r_shuffled <- NULL
  val_mae_model <- NULL
  val_mae_extra_basis <- NULL
  val_mae_truncated <- NULL
  val_mae_shuffled <- NULL

  j <- 1
  for (ind in data_splits){
    y_val <- as.matrix(Y_invivo[-ind,])
    y_train <- as.matrix(Y_invivo[ind,])
    x_val <- X_invivo[-ind,]
    x_train <- X_invivo[ind,]

    print(paste0('Begun cross-validation for fold  ',j,' out of ',num_folds,
                 ' using ',LV_number,' latent variables'))
    plsr_model <- opls(x = x_train,
                       y = y_train,
                       predI = LV_number,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    ## train also a random model
    x_train_shuffled <- x_train[,sample.int(ncol(x_train))]
    colnames(x_train_shuffled) <- colnames(x_train)
    plsr_model_random <- opls(x = x_train_shuffled,
                       y = y_train,
                       predI = LV_number,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    # Get Wh of PLSR
    Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(X_invivo))
    rownames(Wh) <- colnames(X_invivo)
    colnames(Wh) <- colnames(plsr_model@weightMN)
    for (ii in 1:nrow(plsr_model@weightMN)){
      Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    }
    # Get regression coefficients
    Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN

    message('Calculating extra basis to incease the clinical predictive power...')
    # Find extra basis
    phi <- Wh %*% Bh
    Wm_opt <- analytical_solution_opt(y=y_train,
                                      W_invitro = W_invitro,
                                      phi = phi)
    Wm_tot <- cbind(W_invitro, Wm_opt)

    rownames(Wm_tot) <- rownames(W_invitro)
    rownames(Wm_opt) <- rownames(W_invitro)
    rownames(Wh) <- rownames(W_invitro)

    # evaluate performance of plsr model, and truncation with W_invitro , with
    # Wm_combo, and W_opt
    Th_train <- x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh
    Th0_train <- x_train %*% W_invitro %*% t(W_invitro) %*% Wh
    Th_val <- x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh
    Th0_val <- x_val %*% W_invitro %*% t(W_invitro) %*% Wh
    # predict train data
    yhat <- predict(plsr_model,x_train)
    yhat_extra_basis <- cbind(1, Th_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    yhat_truncated <-  cbind(1, Th0_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    # predict validation data
    yhat_val <- predict(plsr_model,x_val)
    yhat_val_shuffled <- predict(plsr_model_random,x_val)
    yhat_extra_basis_val <- cbind(1, Th_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    yhat_truncated_val <-  cbind(1, Th0_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    if (ncol(Y_invivo)>1){
      train_r_model[j] <- mean(diag(cor(yhat,y_train)))
      train_r_extra_basis[j] <- mean(diag(cor(yhat_extra_basis,y_train)))
      train_r_truncated[j] <- mean(diag(cor(yhat_truncated,y_train)))
      train_mae_model[j] <- apply(abs(yhat-y_train),2,mean)
      train_mae_extra_basis[j] <- apply(abs(yhat_extra_basis-y_train),2,mean)
      train_mae_truncated[j] <- apply(abs(yhat_truncated-y_train),2,mean)

      val_r_model[j] <- mean(diag(cor(yhat_val,y_val)))
      val_r_extra_basis[j] <- mean(diag(cor(yhat_extra_basis_val,y_val)))
      val_r_truncated[j] <- mean(diag(cor(yhat_truncated_val,y_val)))
      val_r_shuffled[j] <- mean(diag(cor(yhat_val_shuffled,y_val)))
      val_mae_model[j] <- apply(abs(yhat_val-y_val),2,mean)
      val_mae_extra_basis[j] <- apply(abs(yhat_extra_basis_val-y_val),2,mean)
      val_mae_truncated[j] <- apply(abs(yhat_truncated_val-y_val),2,mean)
      val_mae_shuffled[j] <- apply(abs(yhat_val_shuffled-y_val),2,mean)
    }else{
      train_r_model[j] <- mean(cor(yhat,y_train))
      train_r_extra_basis[j] <- mean(cor(yhat_extra_basis,y_train))
      train_r_truncated[j] <- mean(cor(yhat_truncated,y_train))
      train_mae_model[j] <- mean(abs(yhat-y_train))
      train_mae_extra_basis[j] <- mean(abs(yhat_extra_basis-y_train))
      train_mae_truncated[j] <- mean(abs(yhat_truncated-y_train))

      val_r_model[j] <- mean(cor(yhat_val,y_val))
      val_r_extra_basis[j] <- mean(cor(yhat_extra_basis_val,y_val))
      val_r_truncated[j] <- mean(cor(yhat_truncated_val,y_val))
      val_r_shuffled[j] <- mean(cor(yhat_val_shuffled,y_val))
      val_mae_model[j] <- mean(abs(yhat_val-y_val))
      val_mae_extra_basis[j] <- mean(abs(yhat_extra_basis_val-y_val))
      val_mae_truncated[j] <- mean(abs(yhat_truncated_val-y_val))
      val_mae_shuffled[j] <- mean(abs(yhat_val_shuffled-y_val))
    }
    j <- j+1
  }
  res_train = rbind(data.frame(type = rep('PLSR',num_folds),
                               r = train_r_model,
                               MAE = train_mae_model,
                               fold = seq(1,num_folds),
                               model = rep('human features',num_folds),
                               set= rep('train',num_folds)),
                    data.frame(type = rep('PLSR',num_folds),
                               r = train_r_truncated,
                               MAE = train_mae_truncated,
                               fold = seq(1,num_folds),
                               input = rep('human features truncated via the original in vitro model',num_folds),
                               set= rep('train',num_folds)),
                    data.frame(type = rep('PLSR',num_folds),
                               r = train_r_extra_basis,
                               MAE = train_mae_extra_basis,
                               fold = seq(1,num_folds),
                               input = rep('human features truncated via the optimized in vitro model',num_folds),
                               set= rep('train',num_folds)))

  res_val = rbind(data.frame(type = rep('PLSR',num_folds),
                               r = val_r_model,
                               MAE = val_mae_model,
                               fold = seq(1,num_folds),
                               model = rep('human features',num_folds),
                               set= rep('validation',num_folds)),
                    data.frame(type = rep('PLSR',num_folds),
                               r = val_r_truncated,
                               MAE = val_mae_truncated,
                               fold = seq(1,num_folds),
                               input = rep('human features truncated via the original in vitro model',num_folds),
                               set= rep('validation',num_folds)),
                    data.frame(type = rep('PLSR',num_folds),
                               r = val_r_extra_basis,
                               MAE = val_mae_extra_basis,
                               fold = seq(1,num_folds),
                               input = rep('human features truncated via the optimized in vitro model',num_folds),
                               set= rep('validation',num_folds)),
                  data.frame(type = rep('PLSR',num_folds),
                             r = val_r_shuffled,
                             MAE = val_mae_shuffled,
                             fold = seq(1,num_folds),
                             input = rep('human features truncated via the optimized in vitro model',num_folds),
                             set= rep('validation',num_folds)))
  res_df <- rbind(res_train,res_val)
  message('Done!')
  # put all results into a list object
  return(res_df)
}
