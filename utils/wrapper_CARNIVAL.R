# Functions for running the CARNIVAL pipeline

# TO DO: Verify complexes

# Libraries
library(OmnipathR)
library(CARNIVAL)
library(igraph)
library(dplyr)

# Define cplex path  must be changed for other users or if installation changes
path_cplex <- "C:/Program Files/IBM/ILOG/CPLEX_Studio1210/cplex/bin/x64_win64/cplex.exe"

# Function 1: Build PPI from Omnipath for human data (default)
# Can also prune network based on user input, be it because those proteins are not
# relevant or they do not signal to a TF of interest
build_PKN_Omnipath <- function(organism = 9606, genes_keep = NULL, tf_target = NULL, keep_complexes = T, save_name = NULL){
  # Load omnipath database of PPI
  ppi_net <- import_omnipath_interactions(organism = organism)
  # If genes_keep is not empty, we prune the ppi to only keep relevant nodes
  if (!is.null(genes_keep)){
    nodes <- c(ppi_net$source_genesymbol, ppi_net$target_genesymbol) %>% unique()
    # If keep complexes we don't filter those out unless the target or source of complex is to be removed. 
    # This needs to be improved based on the members of the complex
    ppi_net_complex <- ppi_net %>% filter((grepl("COMPLEX", source) & (target_genesymbol %in% genes_keep))| 
                                          (grepl("COMPLEX", target) & (source_genesymbol %in% genes_keep)))
    # Now remove interactions where either source of target are to be ignored
    ppi_net <- ppi_net %>% filter(source_genesymbol %in% genes_keep & target_genesymbol %in% genes_keep)
    # Add complexes back in if they should be kept
    if (keep_complexes){
      ppi_net <- rbind(ppi_net, ppi_net_complex)
    }
  }
  
  # If tf_target is not empty, we prune the ppi to only keep nodes that can reach TF of interest
  if (!is.null(tf_target)){
    # Build graph to calculate paths
    OPI_g <- interaction_graph(interactions = ppi_net)
    # Get distances of shortest paths between all nodes and the target TFs
    distance_mat <- distances(OPI_g, to = tf_targets[tf_targets %in% ppi_net$target_genesymbol], mode = "out")
    # Define a node "leverage" which is the number of TF of interest it can reach in the directed PPI graph
    node_leverage <- rowSums(!(is.infinite(distance_mat)))
    # Keep nodes that don't have Inf for at least one TF 
    keep_node <- apply(distance_mat,MARGIN = 1, FUN = function(x) !all(is.infinite(x))) 
    nodes_keep <- rownames(distance_mat)[keep_node]
    ppi_net <- ppi_net %>% filter(source_genesymbol %in% nodes_keep & target_genesymbol %in% nodes_keep)
  }
    
  # Prune network - get rid of disconnected nodes and complexes
  OPI_g <- interaction_graph(interactions = ppi_net)
  # Count "in" and "out" degree
  node_degree <- data.frame(n_in = degree(OPI_g, mode = "in"),
                              n_out = degree(OPI_g, mode = "out"))
  # Eliminate complexes that have no incoming nodes
  node_degree <- node_degree %>% filter(!(grepl("_", rownames(node_degree)) & n_in == 0))

  # changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
  ppi_net$consensus_stimulation[which( ppi_net$consensus_stimulation == 0)] = -1
  ppi_net$consensus_inhibition[which( ppi_net$consensus_inhibition == 1)] = -1
  ppi_net$consensus_inhibition[which( ppi_net$consensus_inhibition == 0)] = 1

  # check consistency on consensus sign and select only those in a SIF format
  sif <- ppi_net[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
    filter(consensus_stimulation==consensus_inhibition) %>%
    unique.data.frame()
  sif$consensus_stimulation <- NULL
  colnames(sif) <- c('source', 'interaction', 'target')
  # Rename complexes
  sif$source <- gsub(":", "_", sif$source)
  sif$target <- gsub(":", "_", sif$target)

# We can further simplify here reduce to skip nodes with just one incoming edge 
# and replace edges with the net MOR


  #save SIF if requested
  if (!is.null(save_name)){
  readr::write_tsv(sif, save_name)
  }
  
  # Return ppi in net and sif format
  return(list(NET = ppi_net, SIF = sif))

} # End of function

# Function 2: Get pathway weights. Takes output from decoupler
# The pvalue of the mlm test is based on a two tail test examining differences from zero
# so it can be transformed to a one tailed test and corrected for the value of the 
# score to get a score similar to the percentiles used in the CARNIVAL paper 
# (verified and values match but this is WAY faster)
# Then, assign weights to key nodes in pathways as assigned in the CARNIVAL paper
# Takes as input the long (non-pivoted) output from decoupler
get_pwy_weights <- function(pwy_acts){
  # Get percentiles
  pwy_acts$percentile <- ifelse(pwy_acts$score < 0, 1 - 0.5*pwy_acts$p_value, 0.5*pwy_acts$p_value)
  # Get pathway score based on percentiles
  pwy_acts$scale_score <- 2*(0.5 - pwy_acts$percentile)
  # Pivot wider
  pwy_acts_mat <- pwy_acts %>%
                  pivot_wider(id_cols = c("source"),
                              names_from = "condition", values_from = "scale_score") %>%
                  as.data.frame() %>%
                  column_to_rownames("source") %>%
                  as.matrix()

  # Change JAK-STAT to JAK_STAT
  rownames(pwy_acts_mat) <- gsub("-", "_", rownames(pwy_acts_mat))
  
  # Representative genes per pathway taken from the SI from the CARNIVAL publication
  rep_genes <- list(Androgen = "AR",
                    EGFR = c("EGFR", "ERBB2"),
                    Estrogen = c("ESR1", "ESR2"),
                    Hypoxia = "HIF1A",
                    JAK_STAT = c("JAK1", "JAK2", "JAK3"),
                    MAPK = c("BRAF", "ARAF", "RAF1" ),
                    NFkB = "NFKB1",
                    PI3K = c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"),
                    TGFb = c("TGFBR1", "TGFBR2", "BMPR1A", "BMPR1B", "BMPR2"),
                    TNFa = c("TNFRSF1A", "TNFRSF1B"),
                    Trail = c("CASP8", "CASP10"),
                    VEGF = c("FLT1", "FLT3", "KDR", "PDGFRA", "PDGFRB"),
                    WNT = "DVL1",
                    p53 = "TP53")
  
  # Convert to a dataframe. Might be a more elegant way
  rep_genes_df <- NULL
  for (p in names(rep_genes)){
    rep_genes_df <- rbind(rep_genes_df,
                          data.frame(path = p, genes = rep_genes[[p]]))
  }
  
  # Map pathway scores to representative gene weights. Again, might be a more elegant way 
  # with mapping, but looping for now
  path_weights <- matrix(0, nrow = nrow(rep_genes_df), ncol = ncol(pwy_acts_mat))
  rownames(path_weights) <- rep_genes_df$genes
  colnames(path_weights) <- colnames(pwy_acts_mat)
  
  for (g in rep_genes_df$genes){
    path_weights[g, ] <- pwy_acts_mat[rep_genes_df$path[rep_genes_df$genes == g],]
  }
  
  return(path_weights)
  
}



