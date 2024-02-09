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

get_PC_dimensions <- function(data,keep_variance = 100){
  # Do PCA - do nos scale, but center
  PCA <- prcomp(data, scale. = F, center = T)
  # Calculate % variance
  per_var <- 100*PCA$sdev^2/sum(PCA$sdev^2)
  if (keep_variance<100){
    # Get number of PCs as Meelim did: more strict of either number of PCs with at least 1% variance
    # or total number with cumulative variance 95%
    n_PC <- min(max(which(per_var >= 1)), min(which(cumsum(per_var) > keep_variance)))
  }else{
    n_PC <- ncol(PCA$x)
  }
  # Keep PCs that at least explain 1% variance
  pc_inds <- paste0('PC',seq(1:n_PC))
  inds <- which(per_var[1:n_PC]>=1)
  pc_inds <- pc_inds[inds]
  
  # Return full scores, rotation matrix and number of PCs
  return(list(scores = PCA$x, rotation = PCA$rotation, keep_PCs = pc_inds,n_PC=n_PC))
}

### Nikos function ###
### method: It can be 1) 'lmPCA' using linear regression in a TransCompR-based space, 
### 2) 'plsrPCA' using PLSR in a TransCompR-based space. The default is lmPCA.
### lasso: The default is FALSE It specifies whether to perform LASSO to drop features in the latent/PCA space.
custom_mae <- function(data, validation) {
  observed <- data$y
  predicted <- fitted(validation)
  mae <- mean(abs(observed - predicted))
  return(mae)
}

translation_model <- function(X_A, X_B, Y_A,method = 'lmPCA',lasso=FALSE,ncomp = 4){
  library(glmnet)
  library(pls)
  # varA <-  sum(apply(X_A,2,sd)^2)
  PC_A <- get_PC_dimensions(X_A,keep_variance=90)
  PC_B <- get_PC_dimensions(X_B,keep_variance=90)
  # Project dataset A onto the PCs of dataset B
  X_A_B <- X_A %*% PC_B$rotation
  # Keep only PCs that explain at least 1% of the MPS variance or 1% of the human variance
  varAB <- apply(X_A_B,2,sd)^2
  ## Select PCs that capture 90% of the variance of human that the MPS can capture
  percAB <- 100*varAB/sum(varAB)
  pc_inds_AB <- paste0('PC',seq(1:ncol(X_A_B)))
  inds_AB <- which(cumsum(percAB)>=95)[1]
  # inds_AB <- which(percAB>=1)
  # pc_inds_AB <- pc_inds_AB[inds_AB]
  pc_inds_AB <- pc_inds_AB[1:inds_AB]
  pcs_2_keep  <- unique(c(PC_B$keep_PCs,pc_inds_AB))
  # pcs_2_keep  <- Reduce(x = list(PC_B$keep_PCs,pc_inds_AB),intersect)
  X_A_B <- X_A_B[,pcs_2_keep]

  
  ### Begin modeling
  if (lasso==TRUE){
    ctrl <- trainControl(method = "cv", number = 10)
    data_modeling <- cbind(as.matrix(X_A_B),as.matrix(Y_A))
    colnames(data_modeling)[ncol(data_modeling)] <- 'NAS' 
    lasso_model <- train(NAS~.,data=data_modeling, method = 'glmnet', trControl = ctrl,metric = "MAE",
                 tuneGrid = expand.grid(alpha = 1, lambda = lambda_grid <- 10^seq(-5, 1, length = 1000)))
    predicted <- predict(lasso_model,data_modeling)
    #produce plot of test MSE by lambda value
    lasso_pcs <- coef(lasso_model$finalModel,s=lasso_model$finalModel$lambdaOpt)
    lasso_pcs <- lasso_pcs@Dimnames[[1]][lasso_pcs@i]
    lasso_pcs <- lasso_pcs[which(lasso_pcs!="(Intercept)")]
    X_A_B_filtered <- X_A_B[,lasso_pcs]
    if (method=='plsrPCA'){
      y <- data.frame(out=Y_A)
      #min(c(nrow(X_A_B_filtered),ncol(X_A_B_filtered)))-1
      model <- plsr(out~.,ncomp =ncomp, data=cbind(X_A_B_filtered,y), validation ="LOO")
      # ncomp.permut <- selectNcomp(model,method="randomization", plot=TRUE,ncomp =length(lasso_pcs),alpha = 0.05)
      # if (ncomp.permut<2){
      #   ncomp <- 2
      # }else{
      #   ncomp <- ncomp.permut
      # }
      observed <- y$out
      predicted <- predict(model)
      predicted <- predicted[,,paste0(ncomp,' comps')]
      mae <- mean(abs(observed - predicted))
      tau <- cor(predicted,observed,method = "kendall")
      r2 <- cor(predicted,observed,method = 'pearson')^2
    }else{
      y <- data.frame(out=Y_A)
      model <- lasso_model
      observed <- y$out
      # predicted <- predict(model,newx = X_A_B)
      mae <- mean(abs(observed - predicted))
      tau <- cor(predicted,observed,method = "kendall")
      r2 <- cor(predicted,observed,method = 'pearson')^2
    }
  }else{
    lasso_pcs <- NULL
    if (method=='lmPCA'){
      y <- data.frame(out=Y_A)
      ctrl <- trainControl(method = "cv", number = 10)
      model <- train(out ~ ., data = cbind(X_A_B,y), method = "glm", trControl = ctrl,trace=F)
      observed <- y$out
      predicted <-predict(model,newdata =X_A_B)
      mae <- mean(abs(observed - predicted))
      tau <- cor(predicted,observed,method = "kendall")
      r2 <- cor(predicted,observed,method = 'pearson')^2
    } else{
      y <- data.frame(out=Y_A)
      #min(c(nrow(X_A_B),ncol(X_A_B)))-1
      model <- plsr(out~.,ncomp =ncomp, data=cbind(X_A_B,y), validation ="LOO")
      # ncomp.permut <- selectNcomp(model,method="randomization", plot=TRUE,ncomp =ncol(X_A_B),alpha = 0.05)
      # if (ncomp.permut<2){
      #   ncomp <- 2
      # }else{
      #   ncomp <- ncomp.permut
      # }
      observed <- y$out
      predicted <- predict(model,ncomp =ncomp,newdata = X_A_B)
      predicted <- predicted[,,paste0(ncomp,' comps')]
      mae <- mean(abs(observed - predicted))
      tau <- cor(predicted,observed,method = "kendall")
      r2 <- cor(predicted,observed,method = 'pearson')^2
    }
  }
  
  # Return PC and translation models
  return(list(PC_A = PC_A, 
              PC_B = PC_B, 
              X_A_B = X_A_B,
              model_A_B = model,
              lasso_pcs = lasso_pcs,
              pcs_kept = pcs_2_keep,
              tain_mae = mae,
              train_r2= r2,
              train_tau = tau))
}


#### Modeling function for multi phenotype
translation_model_multi <- function(X_A, X_B, Y_A,valX,valY){
  library(glmnet)
  library(pls)
  # varA <-  sum(apply(X_A,2,sd)^2)
  PC_A <- get_PC_dimensions(X_A,keep_variance=90)
  PC_B <- get_PC_dimensions(X_B,keep_variance=90)
  # Project dataset A onto the PCs of dataset B
  X_A_B <- X_A %*% PC_B$rotation
  # Keep only PCs that explain at least 1% of the MPS variance or 1% of the human variance
  varAB <- apply(X_A_B,2,sd)^2
  ## Select PCs that capture 90% of the variance of human that the MPS can capture
  percAB <- 100*varAB/sum(varAB)
  pc_inds_AB <- paste0('PC',seq(1:ncol(X_A_B)))
  inds_AB <- which(cumsum(percAB)>=95)[1]
  # inds_AB <- which(percAB>=1)
  # pc_inds_AB <- pc_inds_AB[inds_AB]
  pc_inds_AB <- pc_inds_AB[1:inds_AB]
  pcs_2_keep  <- unique(c(PC_B$keep_PCs,pc_inds_AB))
  # pcs_2_keep  <- Reduce(x = list(PC_B$keep_PCs,pc_inds_AB),intersect)
  X_A_B <- X_A_B[,pcs_2_keep]
  
  ## Do not use lobular_inflammation
  # Y_A <- Y_A[,which(colnames(Y_A)!='NAS')]
  # valY <- valY[,which(colnames(valY)!='NAS')]
  
  ### Begin modeling
  ctrl <- trainControl(method = "cv", number = 10)
  data_modeling <- cbind(as.matrix(X_A_B),apply(as.matrix(Y_A),c(1,2),as.numeric))
  val_X_A_B <- valX %*% PC_B$rotation
  val_X_A_B <- val_X_A_B[,pcs_2_keep]
  data_modeling_validation <- cbind(as.matrix(val_X_A_B),apply(as.matrix(valY),c(1,2),as.numeric))
  lasso_pcs <- NULL
  model <- NULL
  k <- 1
  for (outVar in colnames(Y_A)){
    form <- as.formula(paste0(outVar,'~.'))
    #which(colnames(data_modeling) %in% c(outVar,colnames(X_A_B)))
    lasso_model <- train(form,
                         data = data_modeling[,c(colnames(X_A_B),outVar)],
                         method = 'glmnet',
                         trControl = ctrl,
                         metric = "MAE",
                         tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-2, 1, length = 100)))
    if (k == 1){
      predicted <- predict(lasso_model,data_modeling[,c(colnames(X_A_B),outVar)])
      if (sd(predicted)==0){
        predicted <- predicted + 1e-8*rnorm(length(predicted))
      }
      # validate
      predicted_val <- predict(lasso_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)])
      if (sd(predicted_val)==0){
        predicted_val <- predicted_val + 1e-8*rnorm(length(predicted_val))
      }
    }else{
      if (sd(predict(lasso_model,data_modeling[,c(colnames(X_A_B),outVar)]))==0){
        predicted <- cbind(predicted, predict(lasso_model,data_modeling[,c(colnames(X_A_B),outVar)])+ 1e-8*rnorm(nrow(data_modeling)))
      }else{
        predicted <- cbind(predicted, predict(lasso_model,data_modeling[,c(colnames(X_A_B),outVar)]))
      }
      if (sd(predict(lasso_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)]))==0){
        predicted_val <- cbind(predicted_val, predict(lasso_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)])+ 1e-8*rnorm(nrow(data_modeling_validation)))
      }else{
        predicted_val <- cbind(predicted_val, predict(lasso_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)]))
      }
    }
    pcs <- coef(lasso_model$finalModel,s=lasso_model$finalModel$lambdaOpt)
    pcs <- pcs@Dimnames[[1]][pcs@i]
    lasso_pcs[[outVar]] <- pcs[which(pcs!="(Intercept)")]
    model[[outVar]] <- lasso_model 
    message(paste0('Finished output ',outVar))
    k <- k+1
  }
  colnames(predicted) <- colnames(Y_A)
  colnames(predicted_val) <- colnames(valY)
  # evaluate 
  mae <- apply(abs(Y_A - predicted),2,mean)
  tau <- cor(predicted,Y_A,method = 'kendall')
  tau <- diag(tau)
  r <- cor(predicted,Y_A,method = 'pearson')
  r <- diag(r)
  # validation
  valmae <- apply(abs(valY - predicted_val),2,mean)
  valtau <- cor(predicted_val,valY,method = 'kendall')
  valtau <- diag(valtau)
  valr <- cor(predicted_val,valY,method = 'pearson')
  valr <- diag(valr)
  
  # Return PC and translation models
  return(list(PC_A = PC_A, 
              PC_B = PC_B, 
              X_A_B = X_A_B,
              model_A_B = model,
              lasso_pcs = lasso_pcs,
              pcs_kept = pcs_2_keep,
              tain_mae = mae,
              train_r= r,
              train_tau = tau,
              val_mae = valmae,
              val_r= valr,
              val_tau = valtau))
}


#### Modeling function for multi phenotype with Sparse PCA
get_sparePCs <- function(sPCA,min_variance=1){
  per_var <- 100*sPCA$sdev^2/sPCA$var
  z <- sPCA$scores
  pc_inds <- which(apply(z,2,sd)>0)
  pc_inds <- paste0('sPC',pc_inds)
  inds <- which(per_var>=min_variance)
  inds <- paste0('sPC',inds)
  pc_inds <- reduce(list(pc_inds,inds),intersect)
  
  # Return full scores, rotation matrix and number of PCs
  return(list(sparse_scores = sPCA$scores, rotation = sPCA$loadings, keep_PCs = pc_inds))
}
sparse_translation_model_multi <- function(sPC_A, sPC_B, X_A,Y_A,valX,valY){
  library(glmnet)
  library(pls)
  # varA <-  sum(apply(X_A,2,sd)^2)
  PC_A <- get_sparePCs(sPC_A,min_variance=0.5)
  PC_B <- get_sparePCs(sPC_B,min_variance=0.5)
  # Project dataset A onto the PCs of dataset B
  X_A_B <- X_A %*% PC_B$rotation
  colnames(X_A_B) <- paste0('sPC',seq(1:ncol(X_A_B)))
  varAB <- apply(X_A_B,2,sd)^2
  # Keep only sPCs that explain at least 1% of the MPS variance or 1% of the human variance
  percAB <- 100*varAB/sum(varAB)
  pc_inds_AB <- paste0('sPC',seq(1:ncol(X_A_B)))
  inds_AB <- which(percAB>=1)
  pc_inds_AB <- pc_inds_AB[inds_AB]
  # pcs_2_keep  <- Reduce(x = list(PC_B$keep_PCs,pc_inds_AB),intersect)
  pcs_2_keep  <- unique(c(PC_B$keep_PCs,pc_inds_AB))
  X_A_B <- X_A_B[,pcs_2_keep]
  
  ### Begin modeling
  # ctrl <- trainControl(method = "cv", number = 10)
  ctrl <- trainControl(method = "LOOCV")
  data_modeling <- cbind(as.matrix(X_A_B),apply(as.matrix(Y_A),c(1,2),as.numeric))
  val_X_A_B <- valX %*% PC_B$rotation
  colnames(val_X_A_B) <- paste0('sPC',seq(1:ncol(val_X_A_B)))
  val_X_A_B <- val_X_A_B[,pcs_2_keep]
  data_modeling_validation <- cbind(as.matrix(val_X_A_B),apply(as.matrix(valY),c(1,2),as.numeric))
  lm_pcs <- NULL
  model <- NULL
  k <- 1
  for (outVar in colnames(Y_A)){
    form <- as.formula(paste0(outVar,'~.'))
    #which(colnames(data_modeling) %in% c(outVar,colnames(X_A_B)))
    # lm_model <- train(form,
    #                   data = data_modeling[,c(colnames(X_A_B),outVar)],
    #                   method = 'lm',
    #                   trControl = ctrl,
    #                   metric = "MAE",
    #                   tuneGrid = expand.grid(intercept = c(TRUE,FALSE)))
    lm_model <- train(form,
                      data = data_modeling[,c(colnames(X_A_B),outVar)],
                      method = 'glmnet',
                      trControl = ctrl,
                      metric = "MAE",
                      tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-2, 1, length = 100)))
    # lm_model <- train(form,
    #                   data = data_modeling[,c(colnames(X_A_B),outVar)],
    #                   method = 'svmLinear3',
    #                   trControl = ctrl)
    if (k == 1){
      predicted <- predict(lm_model,data_modeling[,c(colnames(X_A_B),outVar)])
      if (sd(predicted)==0){
        predicted <- predicted + 1e-8*rnorm(length(predicted))
      }
      # validate
      predicted_val <- predict(lm_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)])
      if (sd(predicted_val)==0){
        predicted_val <- predicted_val + 1e-8*rnorm(length(predicted_val))
      }
    }else{
      if (sd(predict(lm_model,data_modeling[,c(colnames(X_A_B),outVar)]))==0){
        predicted <- cbind(predicted, predict(lm_model,data_modeling[,c(colnames(X_A_B),outVar)])+ 1e-8*rnorm(nrow(data_modeling)))
      }else{
        predicted <- cbind(predicted, predict(lm_model,data_modeling[,c(colnames(X_A_B),outVar)]))
      }
      if (sd(predict(lm_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)]))==0){
        predicted_val <- cbind(predicted_val, predict(lm_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)])+ 1e-8*rnorm(nrow(data_modeling_validation)))
      }else{
        predicted_val <- cbind(predicted_val, predict(lm_model,data_modeling_validation[,c(colnames(val_X_A_B),outVar)]))
      }
    }
    # pcs <- rownames(summary(lm_model)$coefficients)[which(summary(lm_model)$coefficients[,4]<=0.01)]
    # pcs <- pcs[which(pcs!="(Intercept)")]
    pcs <- coef(lm_model$finalModel,s=lm_model$finalModel$lambdaOpt)
    pcs <- pcs@Dimnames[[1]][pcs@i]
    lm_pcs[[outVar]] <- pcs[which(pcs!="(Intercept)")]
    model[[outVar]] <- lm_model 
    message(paste0('Finished output ',outVar))
    k <- k+1
  }
  colnames(predicted) <- colnames(Y_A)
  colnames(predicted_val) <- colnames(valY)
  # evaluate 
  mae <- apply(abs(Y_A - predicted),2,mean)
  tau <- cor(predicted,Y_A,method = 'kendall')
  tau <- diag(tau)
  r <- cor(predicted,Y_A,method = 'pearson')
  r <- diag(r)
  # validation
  valmae <- apply(abs(valY - predicted_val),2,mean)
  valtau <- cor(predicted_val,valY,method = 'kendall')
  valtau <- diag(valtau)
  valr <- cor(predicted_val,valY,method = 'pearson')
  valr <- diag(valr)
  
  # Return PC and translation models
  return(list(sPC_A = PC_A, 
              sPC_B = PC_B, 
              X_A_B = X_A_B,
              model_A_B = model,
              lm_pcs = lm_pcs,
              pcs_kept = pcs_2_keep,
              tain_mae = mae,
              train_r= r,
              train_tau = tau,
              val_mae = valmae,
              val_r= valr,
              val_tau = valtau))
}
