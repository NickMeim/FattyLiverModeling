# Functions for getting TF activity and pathway activity on bulk RNA seq data
# with the decoupler package. Using these central functions for consistency

# TO DO: Add functions for plotting

# TO DO: Maybe add function for getting resources to avoid reloading every time

# Libraries
library(decoupleR)
library(dplyr)

# Get resources (TF net and pwy net) to pass as defaults.
# Collectri has way more TF with way more targets than Dorothea, so 
# what I am doing is using the TF from dorothea but augment them with the targets 
# from Collectri
#tf_net_dorothea <- get_dorothea(levels = c("A", "B", "C"))
tf_net_collectri <- get_collectri() #%>% 
                    #filter(source %in% tf_net_dorothea$source) %>%
                    #mutate(confidence = "collectri")

#tf_net <- rbind(tf_net_dorothea, tf_net_collectri) %>% 
#          select(!confidence) %>%
#          unique()


#tf_net <- tf_net_dorothea
tf_net <- tf_net_collectri
# Get gene signatures for the pathways. Top 500 genes by default
pwy_net <- get_progeny(organism = "human", top = 500)

# Function 1: Estimate TF activities. Using the ulm method as done in the vignettes
# by the authors, but included an option to run with VIPER. Using human data 

get_TF_activities <- function(data, net = tf_net, method = "ulm", tf_keep = NULL){
  # Prune to only keep TF of interest
  if (!is.null(tf_keep)){
    net <- net %>% filter(source %in% tf_keep)
  }
  # Get activities - only ulm or VIPER for now
  tf_acts <- switch(method,
                    "ulm" = run_ulm(mat = data, net = net, minsize = 5),
                    "viper" = run_viper(mat = data, net = net, minsize = 5),
                    "wmean" = run_wmean(mat = data, net = net, minsize = 5) %>% filter(statistic == "norm_wmean"),
                    "mlm" = run_mlm(mat = data, net = net, minsize = 5)
                    )

  # Pivot to a long table (matrix)
  tf_mat <- tf_acts %>% 
            pivot_wider(id_cols = source, values_from = score, names_from = condition) %>%
            column_to_rownames("source") %>%
            as.matrix()
  
  # Return both matrices
  return(list(MAT = tf_mat, TABLE = tf_acts))
}

# Function 2: Estimate pathway activities. Using the mlm method as done in the vignettes
# by the authors. Using human data 

get_pwy_activities <- function(data, net = pwy_net, method = "ulm"){
  # Get pathway activity
  pwy_acts <- run_mlm(mat = data, net = net, .source='source', .target='target',
                      .mor='weight', minsize = 5)
  
  # Pivot to a long table (matrix)
  pwy_mat <- pwy_acts %>% 
             pivot_wider(id_cols = source, values_from = score, names_from = condition) %>%
             column_to_rownames("source") %>%
             as.matrix()
  
  # Pivot to store pvalues as matrix
  p_mat <- pwy_acts %>% 
            pivot_wider(id_cols = source, values_from = p_value, names_from = condition) %>%
            column_to_rownames("source") %>%
            as.matrix()
          
  # Return both matrices
  return(list(MAT = pwy_mat, TABLE = pwy_acts, MAT_pval = p_mat))
}





