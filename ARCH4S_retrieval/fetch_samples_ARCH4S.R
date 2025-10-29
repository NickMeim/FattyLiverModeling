# Download RNAseq datasets from ARCH4S
library("rhdf5")    
# Set directory to export datasets
export_dir <- "Y:/users/jcadavid/Datasets/Liver/Datasets_ARCH4S"
# Define URL to fetch data (can be changed with newer updates from ARCH4S)
url = "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"
# Read samples and genes in archive
samples = h5read(url, "meta/samples/geo_accession", s3 = T)
genes = h5read(url, "meta/genes/symbol", s3 = T)

# Define function to extract expression from a given list of samples
fetch_expression <- function(samples, collection_name, samples_archive = samples, genes_archive = genes,  dir = export_dir){
  # Name expression file
  extracted_expression_file = paste0(collection_name,"_expression_matrix.tsv")
  # Identify columns to be extracted
  sample_locations = which(samples_archive %in% samples)
  # Extract gene expression from compressed data
  expression = t(h5read(url, "data/expression", index=list(sample_locations, 1:length(genes_archive)), s3 = T))
  rownames(expression) = genes_archive
  colnames(expression) = samples_archive[sample_locations]
  # Print file
  write.table(expression, file = extracted_expression_file, sep="\t", quote=FALSE, col.names=NA)
  print(paste0("Expression file was created at ", export_dir, "/", extracted_expression_file))
}

# Load samples to download
source("Y:/users/jcadavid/Datasets/Liver/Datasets_ARCH4S/samples_to_download.R")

# Loop and extract - only if samples haven't been downloaded before
for (ii in 1:length(list_accession)){
  collection_name <- names(list_accession)[ii]
  if(!file.exists(paste0(export_dir, "/", collection_name, "_expression_matrix.tsv"))){
    # Print stuff
    print(paste0("Downloading expression file for ", collection_name ," (",length(list_accession[[ii]]$samples), " samples)"))
    # Get data - append the "GSM" prefix
    fetch_expression(paste0("GSM",list_accession[[ii]]$samples), collection_name, samples_archive = samples, genes_archive = genes, dir = export_dir)
  }
}



