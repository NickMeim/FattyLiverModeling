# Auxiliary functions for plotting or other repetitive tasks. Non-specific to a project
library(ggrepel)


# Function 1: Dummy function for running PCA and making plots

# Samples in columns, features in rows. Metadata has sample info to add aesthetics
# aesthetics MUST correspond to a column in metadata
plot_pca <- function(data, metadata = NULL, color_by = NULL, shape_by = NULL, 
                     scale_data = T, return_data = F, label_features = F){
  # Do PCA and extract variances
  pca <- prcomp(t(data), scale. = scale_data)
  pca_var <- round(100*pca$sdev^2/sum(pca$sdev^2),2)
  # Dummy dataframe with metadata
  df <- cbind(metadata, pca$x)
  # Make plot - scores
  plt <- ggplot(df, aes_string(x = "PC1", y = "PC2", 
                               color = color_by, shape = shape_by)) +
          geom_point(size = 1.5) +
          labs(x = paste0("PC1 (",pca_var[1],"%)"),
               y = paste0("PC2 (",pca_var[2],"%)"))
  
  # Make plot - loadings
  df <- data.frame(pca$rotation, feature = rownames(pca$rotation))
  sd_PC1 <- sd(df$PC1)
  sd_PC2 <- sd(df$PC2)
  df <- df %>% mutate(is.big = (abs(PC1) > 2*sd_PC1) | (abs(PC2) > 2*sd_PC2))
  
  # Override if we want tolabel all loadings
  if (label_features == T){
    plt_load <- ggplot(df, aes(x = PC1, y = PC2, label = feature)) +
                geom_point(size = 1.5) +
                geom_text_repel() 
  }
  else {
    plt_load <- ggplot(df, aes(x = PC1, y = PC2, alpha = is.big, label = feature)) +
      geom_point(size = 1.5) +
      geom_vline(xintercept = 2*sd_PC1, color = "red", linetype = 2) +
      geom_vline(xintercept = -2*sd_PC1, color = "red", linetype = 2) +
      geom_hline(yintercept = 2*sd_PC2, color = "red", linetype = 2) +
      geom_hline(yintercept = -2*sd_PC2, color = "red", linetype = 2) +
      geom_text_repel(data = df %>% filter(is.big == T))
    
  }
  plt_load <- plt_load + 
              labs(x = paste0("PC1 (",pca_var[1],"%)"),
                   y = paste0("PC2 (",pca_var[2],"%)")) +
              theme(legend.position = "none")
                
  
  if (return_data){
    return(list(data = pca, plt_score = plt, plt_load = plt_load))
  }
  else {
    return(list(plt_score = plt, plt_load = plt_load))
  }
}
