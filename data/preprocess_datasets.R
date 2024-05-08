# Structure datasets to a common format and export them
dir_data_save <- "E:/Jose Luis/Documents/GitHub/FattyLiverModeling/data/"
dir_data <- "Y:/users/jcadavid/Datasets"

# Function to add counts of repeated genes
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

# Load gene mart to convert gene names
gene_mart <- read.table("Y:/users/jcadavid/Datasets/Misc/mart_export_GRCh38.txt", header = T, sep = "\t")
gene_mart <- gene_mart[,c(1,5)] %>% unique()
colnames(gene_mart) <- c("X", "gene")

translate_genes <- function(data, gene_mart){
  # Merge gene names and translate to gene symbol
  data <- merge(gene_mart, data, by = "X")
  # Some genes are not translated so we keep their ENSG
  data$gene[data$gene == ""] <- data$X[data$gene == ""]
  # Remove ENSG
  data <- data[,-1]
  # Some genes are duplicated, so add counts - CHECK IF THIS IS CORRECT
  data <- process_duplicated_rows(data, "gene")
  return(data)
}

# TO DO: Wrap up in simpler functions if we decide on how to store the datasets

## Hoang et al (clinical) - Human data
# Load count matrix from Hoang et al retrieved from ARCH4S
data <- read.table(paste0(dir_data, "/Liver/Datasets_ARCH4S/GSE13090_expression_matrix.tsv"), header = T, sep = "\t")
# Some genes are duplicated, so add counts - CHECK IF THIS IS CORRECT
# Order of genes will be shuffled a bit
data <- process_duplicated_rows(data, "X")
# Load metadata
metadata <- read.table(paste0(dir_data, "/Liver/GSE130970_Hoang_etal/METADATA.txt"), header = T, sep = ",")
# Remove outliers based on PCA plots
keep <- !grepl("36311|36349|36347",metadata$Run)
data <- data[,keep]
metadata <- metadata[keep,]
save(data, metadata, file = paste0(dir_data_save,"GSE13090_Hoang_dataset.RData"))

## Govaere et al (clinical) - Human data
# Load count matrix from Govaere et al retrieved from ARCH4S version as well
data <- read.table(paste0(dir_data, "/Liver/Datasets_ARCH4S/GSE135251_expression_matrix.tsv"), header = T, sep = "\t")
# Some genes are duplicated, so add counts - CHECK IF THIS IS CORRECT
# Order of genes will be shuffled a bit
data <- process_duplicated_rows(data, "X")
# Load metadata
metadata <- read.table(paste0(dir_data, "/Liver/GSE135251_Govaere_etal/METADATA.txt"), header = T, sep = ",")
# The dataset from ARCH4S has fewer samples so we'll exclude them from the metadata for now
metadata <- metadata %>% filter(GEO_Accession..exp. %in% colnames(data))
rownames(metadata) <- metadata$GEO_Accession..exp.
metadata <- metadata[colnames(data),]
rownames(metadata) <- NULL
save(data, metadata, file = paste0(dir_data_save,"GSE135251_Govaere_dataset.RData"))

## Nilsson et al (macrophage perturbation as a dummy ctrl ) 
data <- read.table(paste0(dir_data, "/Misc/GSE202515_Nilsson_etal/GSE202515_counts.tsv"), header = T, sep = "\t")
rownames(data) <- data$X
data <- data[,-1]
# Load metadata
metadata <- read.table(paste0(dir_data, "/Misc/GSE202515_Nilsson_etal/GSE202515_metadata.tsv"), header = T, sep = "\t")
exp_factors <- c("ligand", "concentration", "lps.stimulation", "lps.concentration", "donor", "time")
save(data, metadata, exp_factors, file = paste0(dir_data_save,"GSE202515_Nilsson_dataset.RData"))

## Wang et al liver chip
# Load realigned data with Kallisto (ARCH4S pipeline) to conserve more features that might have been filtered
# out during the original processing of the raw data
data <- read.table("Y:/users/jcadavid/Datasets/Liver/GSE166256_Wang_etal_liverchip/GSE166256_kallisto_counts.tsv", header = T, sep = "\t")
rownames(data) <- data$Gene
data <- data[,-1]
# Extract metadata
metadata <- read.table("Y:/users/jcadavid/Datasets/Liver/GSE166256_Wang_etal_liverchip/METADATA_Wangetal2020.txt", header = T, sep = "\t")
exp_factors <- c("scaffold", "media", "condition")
save(data, metadata, exp_factors, file = paste0(dir_data_save,"GSE166256_Wang_dataset.RData"))

## Kostrewski et al liver chip
# Data from GEO as raw data is not available in ARCH4S
data <- read.table(paste0(dir_data, "/Liver/GSE168285_Kostrzewski_etal/GSE168285_raw_counts.txt"), header = T, sep = "\t")
colnames(data)[1] <- "X"
data <- translate_genes(data, gene_mart)
# Load metadata
metadata <- read.table(paste0(dir_data, "/Liver/GSE168285_Kostrzewski_etal/GSE168285_sample_meta_data.txt"), header = T, sep = "\t")
# Experimental cues to be considered as part of a design matrix - these depend on each dataset
# for this dataset I am treating the combination name as a "cocktail" variable, but could break up into 
# individual cues - see how we correct for batches later on
exp_factors <- c("NPC", "Background", "TGF", "LPS",  "Cholesterol", "Fructose")
# Load design matrix to parse combination of cues
exp_matrix <- read.table(paste0(dir_data, "/Liver/GSE168285_Kostrzewski_etal/design_matrix_MPS.txt"), header = T, sep = "\t")
metadata <- cbind(metadata, exp_matrix)
save(data, metadata, exp_factors, file = paste0(dir_data_save,"GSE168285_Kostrzewski_dataset.RData"))

## Feaver et al Hemoshear device
data <- read.table(paste0(dir_data,"/Liver/Datasets_ARCH4S/GSE89063_expression_matrix.tsv"), header = T, sep = "\t")
data <- process_duplicated_rows(data, "X")
# Load metadata
metadata <- read.table(paste0(dir_data,"/Liver/Datasets_ARCH4S/sample_data_Feaver.txt"), header = T, sep = "\t")
exp_factors <- "media"
save(data, metadata, exp_factors, file = paste0(dir_data_save,"GSE89063_Feaver_dataset.RData"))

