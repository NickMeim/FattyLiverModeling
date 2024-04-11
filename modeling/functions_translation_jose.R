

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
      if (any(grepl(name, files, ignore.case = T))){ # Avoid loading bad matches
        load(file = paste0(dir_data, files[grep(name, files, ignore.case = T)]))
        print(paste0("Loading file: ", paste0(dir_data, files[grep(name, files, ignore.case = T)])))
        data_list[[name]]$counts <- data
        data_list[[name]]$metadata <- metadata
        data_list[[name]]$genes <- rownames(data)
      }
    }
    # Remove redundant first element
    data_list[[1]] <- NULL
  }
  
  return(data_list)

} # End function

# Function for processing datasets - intercept gene lists, CPM normalize, and log2 + 1
# Save some plots along the way
process_datasets <- function(data_list, plt_list = NULL, filter_variance = F){
  
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
      
  }
  
  return(list(data_list = data_list,
         plt_list = plt_list))
}

# Aux function for LASSO based model selection - implementing this to learn, but can copy from Nikos
# Assuming Y is a vector. If it were multiple variables, we could run this for each predicted variable
get_lasso_model <- function(X, Y, idx_train_list = NULL){
  
  
  # If no training indices are given (should be a list), create 10 folds
  if (is.null(idx_train_list)){
    idx_train_list <- createFolds(Y, k = 10, returnTrain = T)
  }
  
  # For each fold, run lasso
  ctrl <- trainControl(method = "cv", number = 10)
  lasso_features_list <- list()
  lasso_MAE <- rep(0, length(idx_train_list))
  lasso_pearson <- lasso_MAE
  # Initialize vector to tally how often a feature appears
  feature_tally <- rep(0, ncol(X))
  # Loop
  for (ff in 1:length(idx_train_list)){
    # Split data
    idx_train <- idx_train_list[[ff]]
    Xtrain <- X[idx_train, ]
    Ytrain <- Y[idx_train]
    Xtest <- X[-idx_train, ]
    Ytest <- Y[-idx_train]
    
    # Run LASSO with a grid
    lasso_model <- train(y = Ytrain,
                         x = Xtrain,
                         method = 'glmnet',
                         trControl = ctrl,
                         metric = "MAE",
                         tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-2, 1, length = 100)))
    
    # Get the coefficients of the final model with the optimal lambda. Exclude intercept
    lasso_features <- coef(lasso_model$finalModel,s=lasso_model$finalModel$lambdaOpt)
    lasso_features <- lasso_features@Dimnames[[1]][lasso_features@i]
    # Store selected features in a list
    lasso_features_list[[ff]] <- lasso_features[which(lasso_features!="(Intercept)")]
    # Calculate MAE in the test dataset
    lasso_MAE[ff] <- (Ytest - predict(lasso_model, Xtest)) %>% abs() %>% mean()
    lasso_pearson[ff] <- cor(Ytest, predict(lasso_model, Xtest))
    # Update feature tally
    feature_tally[colnames(X) %in% lasso_features] <- 1 + feature_tally[colnames(X) %in% lasso_features]
  }
  
  # Create a finalized model with features that appear more than 60% of the times (double-check threshold)
  mod_final <- lm(Y ~ X[, feature_tally > 8])
  
  # Return
  return(list(feature_tally = feature_tally, 
              lasso_MAE = lasso_MAE, 
              lasso_pearson = lasso_pearson,
              lasso_features_list = lasso_features_list,
              final_features = colnames(X)[feature_tally > 6],
              final_model = mod_final))
  
}


# Model cross validation - simple linear regression - could lump with previous function with an additional argument
get_model_CV <- function(X, Y, idx_train_list = NULL){
  
  
  # If no training indices are given (should be a list), create 10 folds
  if (is.null(idx_train_list)){
    idx_train_list <- createFolds(Y, k = 10, returnTrain = T)
  }
  
  
  mod_MAE <- rep(0, length(idx_train_list))
  mod_pearson <- mod_MAE
  mod_R2Y <- mod_MAE 
  mod_Q2Y <- mod_MAE
  
  # Loop
  for (ff in 1:length(idx_train_list)){
    # Split data
    idx_train <- idx_train_list[[ff]]
    Xtrain <- X[idx_train, ]
    Ytrain <- Y[idx_train]
    Xtest <- X[-idx_train, ]
    Ytest <- Y[-idx_train]
    
    df_train <- data.frame(y = Ytrain, X = Xtrain)
    df_test <- data.frame(y = Ytest, X = Xtest)
    # Train model
    mod <- lm(y ~., data = df_train)
    # Get R2Y in train
    Ypred <- predict(mod, df_train)
    mod_R2Y[ff] <- 1 - sum((Ytrain - Ypred)^2)/sum((Ytrain - mean(Ytrain))^2)
    # Get performance in test dataset
    Ypred <- predict(mod, df_test)
    mod_Q2Y[ff] <- 1 - sum((Ytest - Ypred)^2)/sum((Ytest - mean(Ytest))^2)
    mod_MAE[ff] <- (Ytest - Ypred) %>% abs() %>% mean()
    mod_pearson[ff] <- cor(Ytest, Ypred)
    
  }
  
  
  # Return
  return(list(mod_MAE = mod_MAE, 
              mod_pearson = mod_pearson,
              mod_R2Y = mod_R2Y,
              mod_Q2Y = mod_Q2Y))
  
}


# Aux function to get log2 + 1 CPM, center dataset and PCA
pca_center <- function(data){
  # log2 + 1 cpm
  data <- apply(data, MARGIN = 2, FUN = function(x){log2(1 + x/sum(x)*1e6)})
  # Center
  data <- data - rowMeans(data)
  #print(dim(data))
  # PCA - no features removed 
  pca <- prcomp(t(data), scale. = F)
  # Calculate percentage variance and cap PCs to 95% total variance or 1% in PC (tunable)
  per_var <- 100*pca$sdev^2/sum(pca$sdev^2)
  # Get number of PCs as Meelim did: more strict of either number of PCs with at least 1% variance
  # or total number with cumulative variance 95%
  nPC <- min(max(which(per_var >= 1)), min(which(cumsum(per_var) > 95)))
  # Append
  pca$nPC <- nPC
  
  return(list(pca = pca, data = data))
}


# Perform PCA and return number of PCs based on % variance. Could be defined by
# boostrapping to find significant PC. This is simply a wrapper function

# Data is a matrix with observations in rows and features in columns (transposed
# from the common format)

get_PC_dimensions <- function(data){
  # Do PCA - do nos scale, but center
  PCA <- prcomp(data, scale. = F, center = T)
  # Calculate % variance
  per_var <- 100*PCA$sdev^2/sum(PCA$sdev^2)
  # Get number of PCs as Meelim did: more strict of either number of PCs with at least 1% variance
  # or total number with cumulative variance 95%
  n_PC <- min(max(which(per_var >= 1)), min(which(cumsum(per_var) > 95)))
  # Return full scores, rotation matrix and number of PCs
  return(list(scores = PCA$x, rotation = PCA$rotation, n_PC = n_PC))
}

# Train a linear model to predict a phenotypic response based on PC scores.
# Since PCs are orthogonal, a sequential approach  where we add PCs one by one
# based on their individual performance works well. This needs to be refined

# X is the scores (columns are PCs, rows are observations)
# Y is the phenotype to model
# n_PC_max is the max number of PCs to consider
# Default 5-fold CV repeated 10 times

# Note: This does not work so well for the projected data since now the PC scores are
# not orthogonal

get_linear_model_PC <- function(X, Y, n_PC_max, n_CV = 5, repeats = 10){
  # Predict response from each PC
  R2_single <- rep(0, n_PC_max)
  # Get R2 for each PC (full dataset)
  for (ii in 1:n_PC_max){
    mod <- lm(Y ~ X[,ii])
    R2_single[ii] <- summary(mod)[["r.squared"]]
  }
  
  # Rank PCs in descending order of R2 (mult by -1 since rank function does ascending)
  PC_rank <- rank(-R2_single)
  
  # Run linear models in CV with different number of PCs
  set.seed(123)
  train.control <- trainControl(method = "cv", number = n_CV)
  # Repeat CV and get the number of PCs that 
  # minimize RMSE. We then take the smallest number of out these runs
  n_opt <- rep(0,repeats)
  for (rep in 1:repeats){
    RMSE_vec <- rep(0, n_PC_max)
    # Train with different numbers of features chosen sequentally
    for (n_PC in 1:n_PC_max){
      # Select PCs to regress
      X_sub <- X[, which(PC_rank <= n_PC)] %>% as.matrix()
      colnames(X_sub) <- paste0("X",1:n_PC)
      # Train in cross validation
      mod_CV <- train(x = X_sub, 
                      y = Y, 
                      method = "lm",
                      trControl = train.control)
      # Get mean RMSE and store
      RMSE_vec[n_PC] <- mod_CV[["results"]][["RMSE"]]
      
    }
    # Find number of PCs that minimize RMSE
    n_opt[rep] <- which(RMSE_vec == min(RMSE_vec))
  }
  # Use the least number of PCs out of all the repeats and choose the PCs
  PC_opts <- which(PC_rank <= min(n_opt))
  # Fit the full model and return
  mod_out <- lm(Y ~ X[,PC_opts])
  
  return(list(PC_opts = PC_opts, mod = mod_out))
  
  
  #mod_CV_df %>% ggplot(aes(x = n_PC, y = RMSE)) + geom_line()
  
 

}


# Alternative approach - fit the full model and only keep statistically significant features

get_linear_model_PC_V2 <- function(X, Y, n_PC_max, n_CV = 5){
  X <- X[,1:n_PC_max]
  mod <- lm(Y ~ X)
  coefs_mod <- summary(mod)[["coefficients"]]
  # Index of significant PCs (exlucing intercept and subtracting 1 because of intercept)
  PC_opts <- which(coefs_mod[,"Pr(>|t|)"] < 0.05)[-1] - 1
  print(PC_opts)
  mod_out <- lm(Y ~ X[, PC_opts])
  
  # Cross-validate reduced and full models
  set.seed(123)
  train.control <- trainControl(method = "cv", number = n_CV)
  mod_CV_full <- train(x = X, 
                       y = Y, 
                       method = "lm",
                       trControl = train.control)
  # Do this colnames because the package gets annoying if only a single vector is subset
  X_sub <- X[, PC_opts] %>% as.matrix()
  colnames(X_sub) <- paste0("X", PC_opts)
  
  mod_CV_trim <- train(x = X_sub, 
                       y = Y, 
                       method = "lm",
                       trControl = train.control)
  
  
  return(list(PC_opts = PC_opts, 
              mod = mod_out, 
              R2_CV_full = mod_CV_full[["resample"]][["Rsquared"]],
              R2_CV_trim = mod_CV_trim[["resample"]][["Rsquared"]]
              )
         )
}


# Alternative approach - recursive feature elimination with caret package

get_linear_model_PC_V3 <- function(X, Y, n_PC_max, n_CV = 5, repeats = 10, sizes_model = NULL){
  # Preset random seed
  set.seed(123)
  # If not model size ranges are given, we accept ALL sizes
  if (is.null(sizes_model)){
    sizes_model = 1:n_PC_max
  }
  # Control options
  control <- rfeControl(functions = lmFuncs, # random forest
                       method = "repeatedcv", # repeated cv
                       repeats = repeats, # number of repeats
                       number = n_CV)
  # Run RFE
  results_rfe <- rfe(x = X[,1:n_PC_max],
                     y = Y,
                     sizes = sizes_model,
                     rfeControl = control)
  
  return(results_rfe)
  
  
}


# Master function for translation - TO BE CLEANED
# This function takes two datasets (data_A, data_B) and a phenotype (Y_A) and translates from A to B to model Y_A
# It keeps common features, does log2 + 1 transform after getting CPM (arguments are count matrices but could use TPM) and uses linear modeling to predict Y_A
# (will be an argument depending on phenotype class).
# Rownames of datasets are gene names (at least the same naming systems to get matching features)
# This function is directional, so it only projects A to B

# TO DO: Add parameters for filtering samples that lost many common features??

get_model_translation <- function(data_A, data_B, Y_A){
  
  # Intersect datasets - keep common genes
  genes_common <- intersect(rownames(data_A), rownames(data_B))
  # Check whether most counts are represented in the common gene set
  #hist(100*colSums(data_B[genes_common,])/ colSums(data_B), main = "Data system B", xlab = "% counts in common genes", ylab = "Count")
  #hist(100*colSums(data_A[genes_common,])/ colSums(data_A), main = "Data system A", xlab = "% counts in common genes", ylab = "Count")
  # Hoang dataset has a couple of mismatched samples. Check them later
  
  # Get percentage of mito counts per sample
  #per_mito_Hoang <- 100*colSums(data_Hoang[grepl("MT-", rownames(data_Hoang)), ])/ colSums(data_Hoang)

  
  # Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
  data_A <- apply(data_A[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
  data_B <- apply(data_B[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
  # Remove features absent from at least 10% in each samples
  keep_gene <- (rowSums(data_A >= 1) >= 0.1*ncol(data_A)) &
               (rowSums(data_B >= 1) >= 0.1*ncol(data_B))

  # Run PCA on each dataset. Not scaling PCA
  # Log and center dataset A and run PCA
  X_A <- log2(1 + data_A[keep_gene,])
  X_A <- t(X_A - rowMeans(X_A))
  PC_A <- get_PC_dimensions(X_A)
  # Log and center dataset B and run PCA
  X_B <- log2(1 + data_B[keep_gene,])
  X_B <- t(X_B - rowMeans(X_B))
  PC_B <- get_PC_dimensions(X_B)
  # Project dataset A onto the PCs of dataset B
  X_A_B <- X_A %*% PC_B$rotation

  
  # Make plots of PCA
  plt_PCA_A <- data.frame(x = PC_A$scores[,1], y = PC_A$scores[,2]) %>%
                ggplot(aes(x = x, y = y)) +
                geom_point(color = "steelblue") +
                labs(x = "ref PC1", y = "ref PC2", title = "Reference data")
  
  plt_PCA_B <- data.frame(x = PC_B$scores[,1], y = PC_B$scores[,2]) %>%
                ggplot(aes(x = x, y = y)) +
                geom_point(color = "red") +
                labs(x = "mod PC1", y = "mod PC2", title = "Model data")
  
  plt_PCA_A_B <- data.frame(x = X_A_B[,1], y = X_A_B[,2]) %>%
                ggplot(aes(x = x, y = y)) +
                geom_point(color = "steelblue") +
                labs(x = "mod PC1", y = "mod PC2", title = "Projected reference")
  
  
  # Use linear models to predict Y_A from PCs - using recursive feature elimination and cross validation
  # to do feature selection
  # Predict Y_A from PC scores in space A
  #mod_A <- get_linear_model_PC_V3(PC_A$scores, Y_A, PC_A$n_PC, sizes_model = 1:min(10, PC_A$n_PC))
  # Using greedy approach
  mod_A <- get_linear_model_PC(PC_A$scores, Y_A, PC_A$n_PC, n_CV = 5, repeats = 10)
  # Predict Y_A from PC scores in space B (translated)
  mod_A_B <- get_linear_model_PC_V3(X_A_B, Y_A, PC_B$n_PC, sizes_model = 1:min(10, PC_B$n_PC))
  
  # Plot the results of the linear models
  dummy <- data.frame(Y_A = Y_A, Reference = predict(mod_A$mod), Translated = predict(mod_A_B$fit))
  
  plt_pred1 <- dummy %>% 
               ggplot(aes(x = Reference, y = Translated)) +
               geom_point(color = "steelblue") +
               geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
               labs(x = "Predicted phenotype (ref)", y = "Predicted phenotype (trans)")
  
  plt_pred2 <- dummy %>% 
              pivot_longer(cols = 2:3, names_to = "source", values_to = "Y_A_pred") %>%
              ggplot(aes(x = Y_A, y = Y_A_pred, color = source)) +
                geom_point() +
                geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
                facet_wrap(facets = vars(source), ncol = 2) +
                labs(x = "Measured phenotype", y = "Predicted phenotype", color = "Model")
  
  # Plot tvalues of the coefficients in the models - clean later and make prettier
  dummy <- summary(mod_A[["mod"]])[["coefficients"]] %>% as.data.frame()
  colnames(dummy) <- c("Estimate", "Std_error", "tval", "pval")
  dummy$PC <- gsub("X[, PC_opts]","",rownames(dummy), fixed = T)
  dummy <- dummy[-1,] %>% arrange(desc(abs(tval))) %>% mutate(PC = factor(PC, levels = PC))
  plt_tval_mod_A <- dummy %>% 
                    ggplot(aes(y = tval, x = PC, fill = pval < 0.05)) +
                    geom_col() +
                    labs(y = "t value", x = NULL)
  
  dummy <- summary(mod_A_B[["fit"]])[["coefficients"]] %>% as.data.frame()
  colnames(dummy) <- c("Estimate", "Std_error", "tval", "pval")
  dummy$PC <- rownames(dummy)
  dummy <- dummy[-1,] %>% arrange(desc(abs(tval))) %>% mutate(PC = factor(PC, levels = PC))
  plt_tval_mod_A_B <- dummy %>% 
                      ggplot(aes(y = tval, x = PC, fill = pval < 0.05)) +
                        geom_col() +
                        labs(y = "t value", x = NULL)
  
  
  # Return PC and translation models
  return(list(PC_A = PC_A, 
              PC_B = PC_B, 
              X_A_B = X_A_B,
              mod_A = mod_A, 
              mod_A_B = mod_A_B, 
              plots = list(plt_PCA_A = plt_PCA_A,
                           plt_PCA_B = plt_PCA_B,
                           plt_PCA_A_B = plt_PCA_A_B,
                           plt_pred1 = plt_pred1,
                           plt_pred2 = plt_pred2,
                           plt_tval_mod_A = plt_tval_mod_A,
                           plt_tval_mod_A_B = plt_tval_mod_A_B)))
}





# This could be evaluated in CV to determine how many more components to add

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
# TO DO - determine number of iterations with a stop criterion 
optim_function_evolutionary <- function(Xh, Wh, Wm = NULL, theta, mod_coefs, num_it = 100, pop_size = 5){
  # Population size - a multiplier times the size of alpha (number of columns of theta)
  # Guarantee at least 40 points so that there can be four top points
  nP <- max(pop_size*ncol(theta), 40)
  # Number of objects that are top 10%
  nTop <- round(0.1*nP)
  # Initial population - random unit norm vectors
  alpha0 <- matrix(rnorm(ncol(theta)*nP), ncol = nP)
  alpha0 <- apply(alpha0, MARGIN = 2, FUN = function(x){x/sqrt(sum(x^2))})
  
  # Initialize best objective function
  best_obj <- rep(0,100)
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
      obj[ii] <- eval_f(alpha0[,ii], T, X_hat, W_hat, mod_coefs = mod_coefs[-1])
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
    alpha0 <- apply(alpha0, MARGIN = 2, FUN = function(x){x/sqrt(sum(x^2))})
    
  }
  
  # Return new optimal PC and performance metrics
  return(list(best_obj = best_obj,
              var_obj = var_obj,
              median_obj = median_obj,
              Wopt = theta %*% alpha_opt))
}


# Function for building PLSR model and returning the rotation matrix, linear weights, cross validation
# and plots. Just a wrapper function
# X is gene exp matrix (samples x gene)
# Y is phenotype matrix (samples x readouts)

# TO DO: Extend to multicolumn Y

get_plsr_mod <- function(X, Y, nCV = 10, idx_train_CV = NULL, ncomp_PLSR = NA){
  # Run PLSR
  plsr_mod <- opls(x = X, 
                   y = Y,
                   scaleC = "center",
                   crossvalI = nCV,
                   predI = ncomp_PLSR)
  
  # Get model predictions
  Y_pred <- predict(plsr_mod)
  
  # The ROPLS package removes features with low variance so we re-include them
  # in the rotation matrix with zero weights. Make sure the features are in the correct order
  # there might be a faster way, but the loop works
  W_plsr <- plsr_mod@weightMN
  W_plsr_full <- matrix(data = 0, ncol = ncol(W_plsr), nrow = ncol(X))
  rownames(W_plsr_full) <- colnames(X)
  colnames(W_plsr_full) <- colnames(W_plsr)
  for (ii in 1:nrow(W_plsr)){
    W_plsr_full[rownames(W_plsr)[ii], ] <- W_plsr[ii,]
  }
  
  # Get regression model based on the rotation matrix (recast the PLSR matrices)
  S <- X %*% W_plsr_full
  colnames(S) <- colnames(W_plsr_full)
  mod_lm <- lm(Y ~ S)
  # Get cross validation in regression model
  mod_CV <- get_model_CV(X = S, Y = Y, idx_train_list = idx_train_CV)
  
  # Get data variance in PLSR components
  var_plsr <- colVars(S)
  
  # Plot predicted Y vs measured Y for full model
  plt <- plot_scatter(x = Y, y = Y_pred, xlab = "Measured", ylab = "Predicted", xjitter = T) +
          geom_abline(slope = 1, intercept = 0, color = "indianred",
                      linetype="dotted", size = 1)
  
  # Return stuff - we only need the regression matrices for now - but return full model for diagnostics
  return(list(W = W_plsr_full,
              mod_lm = mod_lm,
              mod_CV = mod_CV,
              plot = plt,
              var_plsr = var_plsr,
              Y_pred = Y_pred,
              plsr_full = plsr_mod))
          
  
}

# Fixed training of PLSR model with proper CV - double check with Nikos
train_plsr <- function(X, Y, idx_train_CV, nFolds_internal = 10, nLV_max = 10){
  
  ### First pass - create new training/test split to choose optimal number of LVs
  
  # 10-fold CV by default
  nFolds <- nFolds_internal
  idx_train_list <- createFolds(Y, k = nFolds, returnTrain = T)
  # Initialize metrics
  nLV <- nLV_max
  Q2Y <- matrix(0, ncol = nLV, nrow = nFolds)
  R2Y <- Q2Y
  
  # For each fold, train with maximum number of LV and get CV performance in sequential
  # number of LVs (because PLSR is built sequentially, this works faster than repeating PLSR
  # for multiple numbers of LVs and it's the same)
  for (ff in 1:nFolds){
    idx_train <- idx_train_list[[ff]]
    Xtrain <- X[idx_train, ]
    Ytrain <- Y[idx_train]
    Xtest <- X[-idx_train, ]
    Ytest <- Y[-idx_train]
    
    # PLSR components are added sequentially, so train once for a high number of LVs
    # and do prediction with increasing number of those LVs
    plsr_mod <- opls(x = Xtrain, 
                     y = Ytrain,
                     scaleC = "center",
                     crossvalI = 2, # keeping internal CVs to low since we are overriding this
                     predI = nLV,
                     fig.pdfC = "none",
                     info.txtC = "none")
    
    # Extract weights and convert them to a matrix that does not exclude features
    W <- plsr_mod@weightMN
    Wf <- matrix(data = 0, ncol = ncol(W), nrow = ncol(Xtrain))
    rownames(Wf) <- colnames(Xtrain)
    colnames(Wf) <- colnames(W)
    for (ii in 1:nrow(W)){
      Wf[rownames(W)[ii], ] <- W[ii,]
    }
    # Predict - fit linear model
    for (j in 1:nLV){
      B <- coef(lm(Ytrain ~ Xtrain %*% Wf[,1:j]))
      Ypred_train <- cbind(1, Xtrain %*% Wf[,1:j]) %*% B
      Ypred <- cbind(1, Xtest %*% Wf[,1:j]) %*% B
      Q2Y[ff, j] <- 1 - sum((Ytest - Ypred)^2)/sum((Ytest - mean(Ytest))^2)
      R2Y[ff, j] <- 1 - sum((Ytrain - Ypred_train)^2)/sum((Ytrain - mean(Ytrain))^2)
    }
    
  } # Next fold
  
  
  
  # Plot performance metrics
  plt_Q2Y_CV <- t(Q2Y) %>% as.data.frame() %>% mutate(nLV = 1:10) %>%
    pivot_longer(cols = 1:nFolds, names_to = "fold", values_to = "Q2Y") %>%
    ggplot(aes(x = factor(nLV), y = Q2Y)) + geom_boxplot()
  
  plt_R2Y_CV <- t(R2Y) %>% as.data.frame() %>% mutate(nLV = 1:10) %>%
    pivot_longer(cols = 1:nFolds, names_to = "fold", values_to = "R2Y") %>%
    ggplot(aes(x = factor(nLV), y = R2Y)) + geom_boxplot()
  
  # Choose nLVs such that Q2Y peaks
  nLV_opt <- which(colMeans(Q2Y) == max(colMeans(Q2Y)))
  
  ### Second pass: Now perform CV with external data splits and also compare with a shuffled model
  ### using the optimal number of LVs. The data splits should be different than those used to get
  ### the optimal number of LVs (double-checked with Nikos)
  
  # Make new fold split
  idx_train_list <- idx_train_CV
  Q2Y_CV <- rep(0,10)
  R2Y_CV <- Q2Y_CV
  MAE_CV <- Q2Y_CV
  Q2Y_CV_shuffle <- Q2Y_CV
  R2Y_CV_shuffle <- Q2Y_CV
  MAE_CV_shuffle <- Q2Y_CV
  
  for (ff in 1:length(idx_train_list)){
    idx_train <- idx_train_list[[ff]]
    Xtrain <- X[idx_train, ]
    Ytrain <- Y[idx_train]
    Xtest <- X[-idx_train, ]
    Ytest <- Y[-idx_train]
    
    # Train PLSR
    plsr_mod <- opls(x = Xtrain, 
                     y = Ytrain,
                     scaleC = "center",
                     crossvalI = 2,
                     predI = nLV_opt,
                     fig.pdfC = "none",
                     info.txtC = "none")
    
    # Predict performance
    Ypred <- predict(plsr_mod, Xtest)
    Q2Y_CV[ff] <- 1 - sum((Ytest - Ypred)^2)/sum((Ytest - mean(Ytest))^2)
    MAE_CV[ff] <- (Ytest - Ypred) %>% abs() %>% mean()
    Ypred <- predict(plsr_mod)
    R2Y_CV[ff] <- 1 - sum((Ytrain - Ypred)^2)/sum((Ytrain - mean(Ytrain))^2)
    
    # Now retrain with shuffled Y in the training data
    Ytrain_shuffle <- sample(Ytrain)
    plsr_mod <- opls(x = Xtrain, 
                     y = Ytrain_shuffle,
                     scaleC = "center",
                     crossvalI = 2,
                     predI = nLV_opt,
                     fig.pdfC = "none",
                     info.txtC = "none")
    
    # Predict performance
    Ypred <- predict(plsr_mod, Xtest)
    Q2Y_CV_shuffle[ff] <- 1 - sum((Ytest - Ypred)^2)/sum((Ytest - mean(Ytest))^2)
    MAE_CV_shuffle[ff] <- (Ytest - Ypred) %>% abs() %>% mean()
    Ypred <- predict(plsr_mod)
    R2Y_CV_shuffle[ff] <- 1 - sum((Ytrain_shuffle - Ypred)^2)/sum((Ytrain_shuffle - mean(Ytrain_shuffle))^2)
  }
  
  ### Final pass: Retrain the final model with all of the data
  plsr_mod_full <- opls(x = X, 
                        y = Y,
                        scaleC = "center",
                        crossvalI = 2,
                        predI = nLV_opt,
                        fig.pdfC = "none",
                        info.txtC = "none")
  # Extract weights and convert them to a matrix that does not exclude features
  W <- plsr_mod_full@weightMN
  Wf <- matrix(data = 0, ncol = ncol(W), nrow = ncol(X))
  rownames(Wf) <- colnames(X)
  colnames(Wf) <- colnames(W)
  for (ii in 1:nrow(W)){
    Wf[rownames(W)[ii], ] <- W[ii,]
  }
  
  ### Return
  return(list(plsr_mod_full = plsr_mod_full,
              Wf = Wf,
              Q2Y_CV = Q2Y_CV,
              R2Y_CV = R2Y_CV,
              MAE_CV = MAE_CV,
              MAE_CV_shuffle = MAE_CV_shuffle,
              Q2Y_CV_shuffle = Q2Y_CV_shuffle,
              R2Y_CV_shuffle = R2Y_CV_shuffle,
              LV_selection = list(plt_Q2Y_CV = plt_Q2Y_CV,
                                  plt_R2Y_CV = plt_R2Y_CV,
                                  Q2Y = Q2Y,
                                  R2Y = R2Y)
              )
         )
}

### Rank LVs of MPS based on how their successive removal drops translation performance

# It seems that a single pass is enough since contribution is additive, but check carefully
# Talk to Nikos in case we need to do this with CV
rank_LV_MPS <- function(Xh, Yh, Wh, Bm, Wm, retrain = F){
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
# To do: Extend with CV
get_info_loss <- function(Xh, Yh, Wh, Wm, Bh = NULL){
  # Define projection matrices to make more readable
  Th <- Xh %*% Wh
  Thm <- Xh %*% Wm %*% t(Wm) %*% Wh
  
  # Estimate regression coefficients if not given and predict Y
  # variable for original space and filtered (backprojected) space
  # TO DO: Extend to multiblock Y
  if (is.null(Bh)){
    mod <- lm(Yh ~. ,data = data.frame(Yh, Th))
    Ypred <- predict(mod)
    Ypred_loss <- predict(mod, newdata = data.frame(Thm))
  } else {
    Ypred <- cbind(1, Th) %*% Bh
    Ypred_loss <- cbind(1, Thm) %*% Bh
  }

  # Retrain model to overfit projected data
  mod_retrain <- lm(Yh ~., data = data.frame(Yh, Thm))
  Ypred_refit <- predict(mod_retrain)
  
  # Generate plot
  plt_mod_project <- rbind(data.frame(x = Yh, 
                                      y = Ypred, 
                                      source = "Original"),
                           data.frame(x = Yh, 
                                      y = Ypred_loss, 
                                      source = "Filtered"),
                           data.frame(x = Yh, 
                                      y = Ypred_refit, 
                                      source = "Filtered retrain")) %>%
    mutate(source = factor(source, levels = c("Original", "Filtered", "Filtered retrain"))) %>%
    ggplot(aes(x = x, y = y, color = factor(source))) +
    geom_point(show.legend = F, size = size_dot) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "indianred", linewidth = size_line) +
    facet_wrap(facets = vars(source), nrow = 1) +
    labs(x = "Y measured", y = "Y predicted", color = "Reference space") +
    stat_cor(aes(label = after_stat(r.label)), show.legend = F, color = "black", label.y = max(Yh), size = size_annotation)
  
  plt_mod_project <- add_theme(plt_mod_project) + scale_color_brewer(palette = color_palette)
  
  # TO DO: Add CV and some more statistics beyond the plot
  return(plt_mod_project)
}





# Function to get representative latent variables that are translatable linear combinations of 
# the PCs of the data
# TO DO: can be extended to use CV to determine optimal number of translatable LVs
get_translatable_LV <- function(Xh, Yh, Wh, Wm, Bh, find_extra = F){
  # Define a vector basis depending on whether we want an extra LV or a spanned LV
  if (find_extra){
    print("Finding extra LVs outside of MPS data")
    theta <- find_extra_basis(Wh, Wm, Xh, ncomp = ncol(Wh))
    theta <- theta$theta
    LV_lab <- "LV extra"
  } else {
    print("Finding translatable LVs in data")
    # Set the PCs of the MPS (Wm) as a set of basis vectors
    theta <- Wm
    LV_lab <- "LV data"
  }
  # Initialize a NULL matrix to store LVs
  Wm_new <- NULL
  # Store sum of squared error between predicted phenotype and prediction when filtering
  ErrorY <- NULL
  dErrorY <- 100
  ii <- 0
  # Iterate and add latent variables until the error metric does not improve by more than 1%
  while (dErrorY > 1){
    # Update index
    ii <- ii + 1
    # Find optimal weights to combine basis vectors - reduce population size multiplier
    print(paste0("Finding LV #",ii,"..."))
    res_opt <- optim_function_evolutionary(X = Xh, Wh = Wh, Wm = Wm_new, 
                                           theta = theta, mod_coefs = Bh, pop_size = 5)
    
    

    
    print("Done!")
    print(res_opt$best_obj[90:100])
    
    
    # Extract new latent variable
    Wopt <- res_opt$Wopt
    colnames(Wopt) <- paste0("LV_opt",ii)
    
    # Augment new basis 
    Wm_new <- cbind(Wm_new, Wopt)
    
    # Calculate error
    Ypred <- cbind(1, Xh %*% Wm_new %*% t(Wm_new) %*% Wh) %*% Bh
    ErrorY[ii] <- sum((Ypred - Yh)^2)
    dErrorY <- ifelse(ii == 1, ErrorY[ii], 100*abs(ErrorY[ii] - ErrorY[ii-1])/ErrorY[ii])
    
    # Find new basis for next iteration - find vectors orthogonal to the current new basis
    # that will still span the MPS space
    print("Finding new vector basis...")
    if (find_extra){
      theta <- find_extra_basis(Wh, Wm, Xh, ncomp = ncol(Wh) - ii)
    } else {
      theta <- find_extra_basis(Wm, Wm_new, Xh, ncomp = ncol(Wm) - ii)
    }
    theta <- theta$theta
    print("Done!")
  }
  # After exiting it means the last component was unnecessary, so we exclude it and convert to matrix
  Wm_new <- matrix(data = Wm_new[,-ii], ncol = ii - 1)
  colnames(Wm_new) <- paste0(LV_lab,1:ncol(Wm_new))
  rownames(Wm_new) <- rownames(Wm)
  
  return(list(Wm_new = Wm_new,
              ErrorY = ErrorY))
 

  
}
