# analytical solution to get extra basis
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
