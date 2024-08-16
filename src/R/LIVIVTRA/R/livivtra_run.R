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
  Wm_tot <- cbind(Wm, Wm_opt)

  message('Calculating translatable components that contain most of the current predictive power...')
  ### Find translatable LV of the in vitro system
  ### Run evolutionary algorithm
  Wm_combo <- get_translatable_LV(X_invivo, Y_invivo, Wh, W_invitro,
                                  rbind(apply(Y_invivo,2,mean),Bh),
                                  find_extra = FALSE,
                                  verbose = TRUE)
  Wm_combo <- Wm_combo$Wm_new

  rownames(Wm_tot) <- rownames(Wm)
  rownames(Wm_opt) <- rownames(Wm)
  rownames(Wm_combo) <- rownames(Wm)
  rownames(Wh) <- rownames(Wm)

  # evaluate performance of plsr model, and truncation with W_invitro , with
  # Wm_combo, and W_opt

  message('Done!')
  # put all results into a list object
  model_list = list(model = plsr_model,
                    W_invitro=W_invitro,Wh=Wh,
                    W_opt=Wm_opt,W_translatable = Wm_combo,W_tot = Wm_tot,
                    performance = res_df)
  return(model_list)
}
