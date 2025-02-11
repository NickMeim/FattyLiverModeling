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

# Scale fector to unit norm
norm_v <- function(V){
  return(V/sqrt(sum(V^2)))
}

# Function for loading datasets and appending them to a list
load_datasets <- function(dataset_names, dir_data){
  # Look at files in data directory
  files <- dir(dir_data)
  # Trim to datasets
  files <- files[grepl("dataset.RData", files)]
  # Loop for dataset names given and load one by one, appending to a list.
  # dataset names don't have to match perfectly
  data_list <- list(NULL)
  if (length(files) > 0){ # Loop if any match was found
    for (name in dataset_names){
      if (any(grepl(name, files, ignore.case = TRUE))){ # Avoid loading bad matches
        load(file = paste0(dir_data, files[grep(name, files, ignore.case = TRUE)]))
        print(paste0("Loading file: ", paste0(dir_data, files[grep(name, files, ignore.case = TRUE)])))
        data_list[[name]]$counts <- data
        data_list[[name]]$metadata <- metadata
        data_list[[name]]$genes <- rownames(data)
        if (exists("exp_factors")){
          data_list[[name]]$exp_factors <- exp_factors
        }
      }
    }
    # Remove redundant first element
    data_list[[1]] <- NULL
  }
  
  return(data_list)

} # End function

# Function for processing datasets - intercept gene lists, CPM normalize, and log2 + 1
# Save some plots along the way
process_datasets <- function(data_list, plt_list = NULL, filter_variance = FALSE){
  
  # Find common features
  genes_common <- NULL

  for (ff in names(data_list)){
    if (is.null(genes_common)){
      genes_common <- data_list[[ff]]$genes
    } else {
      genes_common <- intersect(genes_common, data_list[[ff]]$genes)
    }
  }
  
  # Filter based on variance (a -5 threshold seems good) cut more genes
  # using OR instead of AND since a gene might be variable in a dataset and be important in it
  if (filter_variance){
    keep_genes <- rep(F, length(genes_common))
    #keep_genes <- rep(T, length(genes_common))
    # Transform counts to log2 + 1 CPM and get variances
    for (ff in names(data_list)){
      cpm <- apply(data_list[[ff]]$counts[genes_common,], 
                   MARGIN = 2, FUN = function(x){log2(1 + x/sum(x)*1e6)})
      # Apply threshold
      keep_genes <- keep_genes | (log2(rowVars(cpm)) > -5)
      #keep_genes <- keep_genes & (log2(rowVars(cpm)) > -5)
    }
    genes_common <- genes_common[keep_genes]

  }

  # Calculate percentage of counts in each dataset that are retained in the common features as a way of QC
  # Also run PCA in each dataset
  for (ff in names(data_list)){
    data_list[[ff]]$per_counts_common <- 100*colSums(data_list[[ff]]$counts[genes_common,])/colSums(data_list[[ff]]$counts)
    # PCA using only counts of common features
    pca <- pca_center(data_list[[ff]]$counts[genes_common,])
    # Append PCA info - PCA and centered data
    data_list[[ff]]$pca <- pca$pca
    data_list[[ff]]$data_center <- pca$data
    # Get normalized variance
    data_list[[ff]]$pca$var_per <- 100*data_list[[ff]]$pca$sdev^2/sum(data_list[[ff]]$pca$sdev^2)
    
    # Plot explained variance
    plt_list[[ff]]$variance_explained <- plot_bars(data_bar = data.frame(PC = 1:20, 
                                                                         Var = data_list[[ff]]$pca$var_per[1:20]),
                                                    mapping_bar = aes(x = PC, y = Var),
                                                    plt_labs = list(xaxis = paste0("Principal component (", ff, ")"),
                                                                   yaxis = "Explained variance (%)"),
                                                    extra_args = list(fill = "steelblue")) +
                                          geom_vline(xintercept = 0.5 + data_list[[ff]]$pca$nPC, 
                                                     linetype="dotted", color = "black", size = 1)
    
    
    
    # Plot PC1 and PC2
    plt_list[[ff]]$pca <- plot_scatter(x = data_list[[ff]]$pca$x[,1], 
                                       y = data_list[[ff]]$pca$x[,2], 
                                       xlab = paste0("PC1 (", round(data_list[[ff]]$pca$var_per[1],2),"%)"), 
                                       ylab = paste0("PC2 (", round(data_list[[ff]]$pca$var_per[2],2),"%)"))
    
    # If the dataset contains an argument called exp_factor that is not "NA", it corresponds to an in vitro model
    # and we find a reduced PC space by averaging the tech reps
    if (!is.null(data_list[[ff]]$exp_factors)){
      exp_factors <- data_list[[ff]]$exp_factors
      groups <- apply(data_list[[ff]]$metadata, MARGIN = 1, FUN = function(x){paste(x[exp_factors], collapse = "_")})
      groups_unique <- unique(groups)
      data_grouped <- matrix(0, nrow = nrow(data_list[[ff]]$data_center), ncol = length(groups_unique))
      rownames(data_grouped) <- rownames(data_list[[ff]]$data_center)
      colnames(data_grouped) <- groups_unique
      for (ii in 1:length(groups_unique)){
        data_grouped[,ii] <- rowMeans(data_list[[ff]]$data_center[, groups == groups_unique[ii]])
      }
      # Re-center and PCA (no log transform)
      pca_group <- pca_center(data_grouped, log_cpm = FALSE, center = TRUE)
      # We only need to return the rotation matrix (minus the final PC that has variance zero)
      data_list[[ff]]$Wm_group <- pca_group$pca$rotation[,-length(groups_unique)]
      data_list[[ff]]$Xm_grouped <- data_grouped
    }
    
  }

  return(list(data_list = data_list,
         plt_list = plt_list))
}

# Aux function to get log2 + 1 CPM, center dataset and PCA
pca_center <- function(data, log_cpm = TRUE, center = TRUE){
  # log2 + 1 cpm
  if (log_cpm){
    data <- apply(data, MARGIN = 2, FUN = function(x){log2(1 + x/sum(x)*1e6)})
  }
  # Center
  if (center){
    data <- data - rowMeans(data)
  }
  #print(dim(data))
  # PCA - no features removed 
  pca <- prcomp(t(data), scale. = FALSE)
  # Calculate percentage variance and cap PCs to 95% total variance or 1% in PC (tunable)
  per_var <- 100*pca$sdev^2/sum(pca$sdev^2)
  # Get number of PCs as Meelim did: more strict of either number of PCs with at least 1% variance
  # or total number with cumulative variance 95%
  nPC <- min(max(which(per_var >= 1)), min(which(cumsum(per_var) > 95)))
  # Append
  pca$nPC <- nPC
  
  return(list(pca = pca, data = data))
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

# It seems that a single pass is enough since contribution is additive, but check carefully
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

# This function projects and backprojects a source dataset onto a target dataset
# to assess how much predictive information is lost based on a regression model
get_info_loss <- function(Xh, Yh, Wh, Wm, Bh = NULL){
  # Define projection matrices to make more readable
  Th <- Xh %*% Wh
  Thm <- Xh %*% Wm %*% t(Wm) %*% Wh
  
  # Estimate regression coefficients if not given and predict Y
  # variable for original space and filtered (backprojected) space
  if (is.null(Bh)){
    Ypred <- matrix(0,nrow = nrow(Yh),ncol = ncol(Yh))
    Ypred_loss <- matrix(0,nrow = nrow(Yh),ncol = ncol(Yh))
    for (i in 1:ncol(Yh)){
      train_data <- cbind(Yh[,i], Th)
      mod <- lm(V1 ~. ,data = as.data.frame(train_data))
      Ypred[,i] <- predict(mod)
      Ypred_loss[,i] <- predict(mod, newdata = data.frame(Thm))
    }
  } else {
    Ypred <- cbind(1, Th) %*% Bh
    Ypred_loss <- cbind(1, Thm) %*% Bh
  }

  # Retrain model to overfit projected data
  Ypred_refit <- matrix(0,nrow = nrow(Yh),ncol = ncol(Yh))
  for (i in 1:ncol(Yh)){
    train_data <- cbind(Yh[,i], Thm)
    mod_retrain <- lm(V1 ~., data = as.data.frame(train_data))
    Ypred_refit[,i] <- predict(mod_retrain)
  }
  
  # Generate plot
  plt_mod_project <- rbind(data.frame(x1 = Yh[,1], 
                                      y1 = Ypred[,1], 
                                      source = "Original"),
                           data.frame(x1 = Yh[,1], 
                                      y1 = Ypred_loss[,1], 
                                      source = "Filtered"),
                           data.frame(x1 = Yh[,1], 
                                      y1 = Ypred_refit[,1], 
                                      source = "Filtered retrain")) %>%
    mutate(source = factor(source, levels = c("Original", "Filtered", "Filtered retrain"))) %>%
    ggplot(aes(x = x1, y = y1, color = factor(source))) +
    geom_point(show.legend = F, size = size_dot) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "indianred", linewidth = size_line) +
    facet_wrap(facets = vars(source), nrow = 1) +
    labs(x = "Y measured", y = "Y predicted", color = "Reference space") +
    stat_cor(aes(label = after_stat(r.label)), show.legend = FALSE, color = "black", label.y = max(Yh[,1]), size = size_annotation)
  
  plt_mod_project <- add_theme(plt_mod_project) + scale_color_brewer(palette = color_palette)
  
  plt_mod_project_fibrosis <- rbind(data.frame(x2 = Yh[,2], 
                                      y2 = Ypred[,2], 
                                      source = "Original"),
                           data.frame(x2 = Yh[,2], 
                                      y2 = Ypred_loss[,2], 
                                      source = "Filtered"),
                           data.frame(x2 = Yh[,2], 
                                      y2 = Ypred_refit[,2], 
                                      source = "Filtered retrain")) %>%
    mutate(source = factor(source, levels = c("Original", "Filtered", "Filtered retrain"))) %>%
    ggplot(aes(x = x2, y = y2, color = factor(source))) +
    geom_point(show.legend = FALSE, size = size_dot) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "indianred", linewidth = size_line) +
    facet_wrap(facets = vars(source), nrow = 1) +
    labs(x = "Y measured", y = "Y predicted", color = "Reference space") +
    stat_cor(aes(label = after_stat(r.label)), show.legend = FALSE, color = "black", label.y = max(Yh[,2]), size = size_annotation)
  
  plt_mod_project_fibrosis <- add_theme(plt_mod_project_fibrosis) + scale_color_brewer(palette = color_palette)
  
  return(list(plt_mod_project,plt_mod_project_fibrosis))
}

# Function to get translatable components for two phenotypes - semi analytical solution
get_translatable_LV_2phenotype <- function(Xh, Yh, Wh, Wm, Bh){
  # Get mean of phenotypes
  Yh_mean <- colMeans(Yh)
  # Center phenotypes
  for (ii in 1:ncol(Yh)){
    Yh[,ii] <- Yh[,ii] - Yh_mean[ii]
  }
  
  # Get predicted phenotypes (centered) with all PCs
  # Find translatable components
  Yh_pred <- Xh %*% Wm %*% t(Wm) %*% Wh %*% Bh
  # Auxiliary matrices
  Xhat <- Xh %*% Wm
  phi <- t(Wm) %*% Wh %*% Bh
  
  # Initialize optimal weights alpha
  alpha_opt <- NULL
  
  # The TCs seemed to be spanned by the columns of phi, so we use those as a first estimate.
  # If phenotypes are fully independent this will be the final solution with no particular ordering
  alpha_0 <- apply(phi, MARGIN = 2, FUN = function(x){x/sqrt(sum(x^2))})
  
  # Find a single vector that provides a good tradeoff in predicting all phenotypes
  # NOTE: Line search for the case of two phenotypes to find an optimal tradeoff.
  w <- seq(0,1,0.01)
  error_TC <- matrix(0, nrow = length(w), ncol = 2)
  for (ii in 1:length(w)){
    v <- w[ii]*alpha_0[,1] + (1-w[ii])*alpha_0[,2]
    v <- norm_v(v) # Normalize combination
    error_TC[ii,1] <- sum( (Yh_pred[,1] - Xhat %*% v %*% (t(v) %*% phi[,1]))^2 )
    error_TC[ii,2] <- sum( (Yh_pred[,2] - Xhat %*% v %*% (t(v) %*% phi[,2]))^2 )
  }
  # Normalize each column by the variance of each phenotype
  error_TC[,1] <- error_TC[,1]/var(Yh_pred[,1])
  error_TC[,2] <- error_TC[,2]/var(Yh_pred[,2])
  # Find optimal combo (for more phenotypes a grid search might not be the best solution)
  w_min <- w[min(rowSums(error_TC)) == rowSums(error_TC)]
  # Append optimal solution as the first TC
  alpha_opt <- cbind(alpha_opt, norm_v(w_min*alpha_0[,1] + (1-w_min)*alpha_0[,2]))
  
  # For two phenotypes a second TC is enough and can be calculated directly to be orthogonal to
  # the first one but still be spanned by phi1 and phi2. For more phenotypes (i.e. components)
  # we can find a new set of basis by doing regression on each deflated phenotype
  # i.e. lm(Xhat %*% phi - Xhat %*% alpha_opt %*% t(alpha_opt) %*% phi ~ Xhat - 1)
  
  # For now, simply get the second TC as an orthogonal vector. Analytical solution
  c1 <- w_min + (1-w_min)*as.numeric(t(alpha_0[,1]) %*% alpha_0[,2])
  c2 <- w_min*as.numeric(t(alpha_0[,2]) %*% alpha_0[,1]) + (1-w_min)
  
  alpha_opt <- cbind(alpha_opt,
                     norm_v(c2*alpha_0[,1] - c1*alpha_0[,2]))
  
  # Convert to TCs and return
  Wm_TC <- Wm %*% alpha_opt
  colnames(Wm_TC) <- paste0("TC", 1:ncol(Wm_TC))
  return(list(Wm_TC = Wm_TC, alpha0 = alpha_0))
  
}
