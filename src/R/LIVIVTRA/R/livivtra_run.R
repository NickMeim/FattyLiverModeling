library(tidyverse)
library(ropls)
source('functions_translation.R')

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
