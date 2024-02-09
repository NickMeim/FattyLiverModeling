# Auxiliary functions for linear model translation (similar to the ones used for TransCompR)

# Using sum for now, but could include an arg for a function. 
# row_ID_col is a column within data that includes the presumed row names
process_duplicated_rows <- function(data, row_ID_col){
  row_ID <- data[, row_ID_col]
  cols_data <- which(colnames(data) != row_ID_col)
  # Identify duplicated IDs
  duplicated_ID <- row_ID[which(duplicated(row_ID))] %>% unique()
  # Change column name of data
  colnames(data)[colnames(data) == row_ID_col] <- "X"
  # Split data into duplicated and non-duplicated
  keep_row <- row_ID %in% duplicated_ID
  data_duplicated <- data[keep_row,]
  # New dataframe to dump the collapsed rows
  data_duplicated_unique <- matrix(0, nrow = length(duplicated_ID), ncol = ncol(data) - 1)
  for (ii in 1:length(duplicated_ID)){
    data_duplicated_unique[ii,] <- colSums(data_duplicated[data_duplicated$X == duplicated_ID[ii], cols_data])
  }
  colnames(data_duplicated_unique) <- colnames(data)[cols_data]
  data_duplicated_unique <- cbind(data.frame(X = duplicated_ID), data_duplicated_unique)
  # Combine with unique data and move row ID column to row names
  data <- rbind(data[!keep_row,], data_duplicated_unique) 
  rownames(data) <- data$X
  data <- data[, cols_data]
  return(data)
}

# Wrapper to zscore (can use function scale with transposition also)
zscore <- function(X, by.rows = F){
  # Convert to matrix
  X <- as.matrix(X)
  # Transpose if necessary
  if (by.rows){
    X <- t(X)
  }
  
  X <- apply(X, 2, function(x) (x - mean(x))/sd(x))
  
  if (by.rows){
    X <- t(X)
  }
  
  return(X)
}

# 



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

#genes_NAFLD <- c('ARPC5', 'YWHAH', 'ARF4', 'TNFRSF12A', 'ADHFE1', 'USP33', 'CD52','ACVR2B',
                # 'ING5', 'ASB3', 'IFI30', 'ERVW-1', 'ERBB3', 'KPNA2', 'COQ10B', 'MAGI1', 'MAPRE1', 'ABCA6')
#dummy <- lm(Y_A ~ t(data_A[genes_NAFLD,]))
