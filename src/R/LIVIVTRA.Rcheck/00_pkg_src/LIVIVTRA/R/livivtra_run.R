library(tidyverse)
library(ropls)

### Utility functions-----------------------------------
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

# Fast version with explicit solution to least squares and simplification
# of matrix operations to reduce redundancies
find_extra_basis <- function(Wh, Wm, Xh, ncomp){

  # Define useful matrices
  WhWm <- Wh - Wm %*% (t(Wm) %*% Wh)
  E_tot_vec <- Xh %*% WhWm
  Wnew <- matrix(data = 0, nrow = nrow(Wm), ncol = ncomp)
  # Loop for columns of Wh
  E_tot <- rep(0, ncomp)
  for (ii in 1:ncomp){
    # Regress without intercept - direct solution
    res <- WhWm[,ii] - Wnew %*% (t(Wnew) %*% Wh[,ii])
    # Normalize to unit length
    res <- res/sqrt(sum(res^2))
    Wnew[,ii] <- res
    # Calculate projection error - reduce successively
    E_tot_vec <- E_tot_vec - Xh %*% Wnew[,ii] %*% (t(Wnew[,ii]) %*% Wh)
    E_tot[ii] <- E_tot_vec^2 %>% sum()
  }

  # Get new vectors
  colnames(Wnew) <- paste0("LV_extra", 1:ncomp)
  return(list(theta = Wnew,
              E_tot = E_tot))
}


# Target objective function to find missing components

# The function to optimize is  || T - X*theta*alpha*alpha'*theta'*Wh|| which will find a weight
# vector alpha (with norm 1) that minimizes the remaining error of projection. It could be further weighted
# to max the correlation between NAFLD predicted and measured using the coefficients of a linear model (exlucing intercept)

# There is no trivial analytical solution because alpha is quite nested in the equation and it becomes non linear
eval_f <- function(x, T, X_hat, W_hat, mod_coefs = NULL){
  if (is.null(mod_coefs)){
    E <- (T - (X_hat %*% x %*% t(x) %*% W_hat))^2 %>% sum()
  } else {
    E <- ((T - (X_hat %*% x %*% t(x) %*% W_hat)) %*% mod_coefs)^2 %>% sum()
  }
  return(E)
}

# Evolutionary-like algorithm to optimize the objective function - minimize error of correlation
optim_function_evolutionary <- function(Xh, Wh, Wm = NULL, theta, mod_coefs, num_it = 100, pop_size = 5){
  # Population size - a multiplier times the size of alpha (number of columns of theta)
  # Guarantee at least 40 points so that there can be four top points
  nP <- max(pop_size*ncol(theta), 40)
  # Number of objects that are top 10%
  nTop <- round(0.1*nP)
  # Initial population - random unit norm vectors
  alpha0 <- matrix(rnorm(ncol(theta)*nP), ncol = nP)
  if (nrow(alpha0)>1){
    alpha0 <- apply(alpha0, MARGIN = 2, FUN = function(x){x/sqrt(sum(x^2))})
  }else{
    alpha0 <- alpha0/sqrt(sum(alpha0^2))
  }
  if (is.null(nrow(alpha0))){
    alpha0 <- t(alpha0)
  }
  # Initialize best objective function
  best_obj <- rep(0,num_it)
  var_obj <- best_obj
  median_obj <- var_obj
  # Calculate matrices
  X_hat <- Xh %*% theta
  W_hat <- t(theta) %*% Wh

  # Calculate remainder
  if (!is.null(Wm)){
    T <- (Xh %*% Wh) - (Xh %*% Wm %*% t(Wm) %*% Wh)
  } else {
    T <- Xh %*% Wh
  }


  # Iterate - successive rounds of testing objective function across a population
  # of vectors, followed by selection of top 10% and combination/mutation.
  # Doing for 100 iterations but could stop prematurely based on changes in the best result
  for (kk in 1:num_it){
    # Evaluate function
    obj <- rep(0, nP)
    for (ii in 1:nP){
      obj[ii] <- eval_f(alpha0[,ii], T, X_hat, W_hat, mod_coefs = mod_coefs[2:nrow(mod_coefs),])
    }
    # Choose best performers and store
    best_obj[kk] <- min(obj)
    var_obj[kk] <- var(obj)
    median_obj[kk] <- median(obj)
    alpha_opt <- alpha0[, best_obj[kk] == obj]
    # Pick best 10%
    idx <- which(rank(obj) <= nTop)

    # Combine them with unit-norm weights
    w <- matrix(rnorm(nTop*nP), nrow = nTop)
    w <- apply(w, MARGIN = 2, FUN = function(x){x/sqrt(sum(x^2))})
    alpha0 <- alpha0[,idx] %*% w

    # Mutate slightly - reduce mutation range as we go
    alpha_mut <- matrix(rnorm(ncol(theta)*nP), ncol = nP)*exp(-kk/10)
    alpha0 <- alpha0 + alpha_mut
    if (nrow(alpha0)>1){
      alpha0 <- apply(alpha0, MARGIN = 2, FUN = function(x){x/sqrt(sum(x^2))})
    }else{
      alpha0 <- alpha0/sqrt(sum(alpha0^2))
    }

  }

  # Return new optimal PC and performance metrics
  return(list(best_obj = best_obj,
              var_obj = var_obj,
              median_obj = median_obj,
              Wopt = theta %*% alpha_opt))
}

### Rank LVs of MPS based on how their successive removal drops translation performance
rank_LV_MPS <- function(Xh, Yh, Wh, Bm, Wm, retrain = FALSE){
  # Get error without removing LVs as a reference

  # If retrain, we not only remove the component but also retrain the model
  if (retrain){
    Ypred <- predict(lm(Yh ~ Xh %*% Wm %*% t(Wm) %*% Wh))
  } else {
    Ypred <- Bm[1] + ((Xh %*% Wm %*% t(Wm) %*% Wh) %*% Bm[-1])
  }
  E_ref <- 1 - sum((Yh - Ypred)^2)/sum((Yh - mean(Yh))^2)

  # Loop until all LVs have been removed
  nLV <- ncol(Wm)
  vec_rem <- NULL
  E_top <- rep(0, nLV)
  E_rem_out <- E_top
  idx_check <- 1:nLV
  for (nn in 1:nLV){
    E_rem <- rep(0, nLV)
    E_rem[!(1:nLV %in% idx_check)] <- 1000 #Do this to make sure the same component is not chosen twice
    for (ii in idx_check){
      # Calculate when removing component ii
      if (retrain){
        Ypred <- predict(lm(Yh ~ Xh %*% Wm[,-c(vec_rem,ii)] %*% t(Wm[,-c(vec_rem,ii)]) %*% Wh))
      } else {
        Ypred <- Bm[1] + ((Xh %*% Wm[,-c(vec_rem,ii)] %*% t(Wm[,-c(vec_rem,ii)]) %*% Wh) %*% Bm[-1])
      }

      E_rem[ii] <- 1 - sum((Yh - Ypred)^2)/sum((Yh - mean(Yh))^2) #sum((Yh - Ypred)^2)
    }
    # Find E_rem that leads to higher error upon removal
    E_top[nn] <- which(E_rem == min(E_rem))[1]
    E_rem_out[nn] <- min(E_rem)
    # Update idx to remove
    vec_rem <- c(vec_rem, E_top[nn])
    idx_check <- idx_check[idx_check != E_top[nn]]
  }
  return(list(E_top = E_top,
              E_rem_out = E_rem_out,
              E_ref = E_ref))

}

# Function to get representative latent variables that are translatable linear combinations of
# the PCs of the data
get_translatable_LV <- function(Xh, Yh, Wh, Wm, Bh, find_extra = FALSE,verbose=TRUE){
  # If there is a single candidate vector, then we simply return it
  if (ncol(Wm) == 1){
    return(list(Wm_new = Wm))
  }
  # Define a vector basis depending on whether we want an extra LV or a spanned LV
  if (find_extra){
    if(verbose){
      print("Finding extra LVs outside of MPS data")
    }
    theta <- find_extra_basis(Wh, Wm, Xh, ncomp = ncol(Wh))
    theta <- theta$theta
    LV_lab <- "LV_extra"
  } else {
    if(verbose){
      print("Finding translatable LVs in data")
    }
    # Set the PCs of the MPS (Wm) as a set of basis vectors
    theta <- Wm
    LV_lab <- "LV_data"
  }
  # Initialize a NULL matrix to store LVs
  Wm_new <- NULL
  # Store sum of squared error between predicted phenotype and prediction when filtering
  ErrorY <- NULL
  dErrorY <- 100
  ii <- 0
  # flag <- TRUE
  # thresh <- 0.01
  # th <- 0.25
  # counter <- 1
  base_err <- 1
  # Iterate and add latent variables until the error metric does not improve by more than 1%
  while ((dErrorY > base_err) & (ii<ncol(Bh)+2)){
    ii <- ii +1
    # Find optimal weights to combine basis vectors - reduce population size multiplier
    if(verbose){
      print(paste0("Finding LV #",ii,"..."))
    }
    res_opt <- optim_function_evolutionary(Xh = Xh, Wh = Wh, Wm = Wm_new,
                                           theta = theta, mod_coefs = Bh, pop_size = 5)
    if(verbose){
      print("Done!")
    }

    # Extract new latent variable
    Wopt <- res_opt$Wopt
    colnames(Wopt) <- paste0("LV_opt",ii)

    Wm_new <- cbind(Wm_new, Wopt)
    # Calculate error
    Ypred <- cbind(1, Xh %*% Wm_new %*% t(Wm_new) %*% Wh) %*% Bh
    ErrorY[ii] <- sum((Ypred - Yh)^2)
    dErrorY <- ifelse(ii == 1, ErrorY[ii], 100*abs(ErrorY[ii] - ErrorY[ii-1])/ErrorY[ii])

    # Find new basis for next iteration - find vectors orthogonal to the current new basis
    # that will still span the MPS space
    if(verbose){
      print("Finding new vector basis...")
    }
    if (find_extra){
      theta <- find_extra_basis(Wh, Wm, Xh, ncomp = ncol(Wh) - ii)
    } else {
      theta <- find_extra_basis(Wm, Wm_new, Xh, ncomp = ncol(Wm) - ii)
    }
    theta <- theta$theta
    # Update index
    # ii <- ii + 1
    if(verbose){
      print("Done!")
    }
  }
  # After exiting it means the last component was unnecessary, so we exclude it and convert to matrix
  if((dErrorY <= 1)){
    Wm_new <- matrix(data = Wm_new[,-ii], ncol = ii - 1)
  }
  colnames(Wm_new) <- paste0(LV_lab,1:ncol(Wm_new))
  rownames(Wm_new) <- rownames(Wm)

  return(list(Wm_new = Wm_new,
              ErrorY = ErrorY))



}
#### Main run functions-----------------------------------------
livivtra_run <- function(X_invivo,Y_invivo,X_invitro,W_invitro,
                         LV_number=8){

  library(tidyverse)
  library(ropls)

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
    mae_model <- mean(apply(abs(yhat-Y_invivo),2,mean))
    mae_extra_basis <- mean(apply(abs(yhat_extra_basis-Y_invivo),2,mean))
    mae_truncated <- mean(apply(abs(yhat_truncated-Y_invivo),2,mean))
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
                                'human features truncated via the optimized in vitro model',
                                'human features truncated via the original in vitro model'),
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

  library(tidyverse)
  library(ropls)
  library(caret)

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
      train_mae_model[j] <- mean(apply(abs(yhat-y_train),2,mean))
      train_mae_extra_basis[j] <- mean(apply(abs(yhat_extra_basis-y_train),2,mean))
      train_mae_truncated[j] <- mean(apply(abs(yhat_truncated-y_train),2,mean))

      val_r_model[j] <- mean(diag(cor(yhat_val,y_val)))
      val_r_extra_basis[j] <- mean(diag(cor(yhat_extra_basis_val,y_val)))
      val_r_truncated[j] <- mean(diag(cor(yhat_truncated_val,y_val)))
      val_r_shuffled[j] <- mean(diag(cor(yhat_val_shuffled,y_val)))
      val_mae_model[j] <- mean(apply(abs(yhat_val-y_val),2,mean))
      val_mae_extra_basis[j] <- mean(apply(abs(yhat_extra_basis_val-y_val),2,mean))
      val_mae_truncated[j] <- mean(apply(abs(yhat_truncated_val-y_val),2,mean))
      val_mae_shuffled[j] <- mean(apply(abs(yhat_val_shuffled-y_val),2,mean))
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
                               input = rep('human features',num_folds),
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
                               input = rep('human features',num_folds),
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
                             input = rep('shuffled features',num_folds),
                             set= rep('validation',num_folds)))
  res_df <- rbind(res_train,res_val)
  message('Done!')
  # put all results into a list object
  return(res_df)
}
