library(tidyverse)

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
