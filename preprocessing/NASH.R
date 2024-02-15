# R script to download selected samples
# Copy code and run on a local machine to initiate download

library(rhdf5)    # can be installed using Bioconductor
library(tidyverse)
library(preprocessCore)
library(sva)
library(ggplot2)

destination_file = "../data/human_gene_v2.2.h5"
extracted_expression_file = "../data/ARCHS4/NASH_expression_matrix.tsv"
extracted_processed_expression_file = "../preprocessing/ARCHS4_processed/NASH_expression_matrix_preprocessed.tsv"
meta_file <- "../data/ARCHS4/NASH_meta_data.tsv"
url = "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    download.file(url, destination_file, quiet = FALSE, mode = 'wb')
}

# Selected samples to be extracted
samp = c("GSM3946330","GSM3946331","GSM5474289","GSM5474290","GSM5474291","GSM5474292","GSM5474293","GSM5474294","GSM5474295","GSM5474296","GSM5474297","GSM5474298","GSM5474299","GSM5474300","GSM5474301","GSM5474302","GSM5474303","GSM5474304","GSM5474305","GSM5474306","GSM5474313","GSM5474314","GSM5474315","GSM5474316","GSM5474317","GSM5474318","")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/symbol")

# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# extract gene expression from compressed data
expression = t(h5read(destination_file, "data/expression", index=list(sample_locations, 1:length(genes))))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE, col.names=NA)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

# Load metadata for these samples
meta <- h5read(destination_file, "meta")
meta <- data.frame(meta[['samples']])
meta <- meta %>% filter(geo_accession %in% samples[sample_locations])
write.table(meta, file=meta_file, sep="\t", quote=FALSE, col.names=NA)
print(paste0("Expression file was created at ", getwd(), meta_file))

### Normalize and batch correct and save pre-processed results ###
# normalize samples and correct for differences in gene count distribution
expression = log2(expression+1)
expression = normalize.quantiles(expression)
rownames(expression) = genes
colnames(expression) = samples[sample_locations]
print(paste0('Meta data and expression matrix samples are in the same order:',all(colnames(expression)==meta$sample)))
# correct batch effects in gene expression
batchid = match(meta$series_id, unique(meta$series_id))
correctedExpression <- ComBat(dat=expression, batch=batchid, par.prior=TRUE, prior.plots=FALSE)

pca <- prcomp(t(correctedExpression),center = T,scale.=F)
df_pca <- pca$x[,1:2]
df_pca <- left_join(as.data.frame(df_pca) %>% rownames_to_column('sample'),meta %>% select(sample,series_id))
ggplot(df_pca,aes(x=PC1,y=PC2,color=series_id)) + geom_point() +
  theme_minimal(base_size=24,base_family='Arial')+
  theme(text = element_text(size=24,family = 'Arial'),
        legend.position = 'none')
# Print pre-processed file
write.table(correctedExpression, file=extracted_processed_expression_file, sep="\t", quote=FALSE, col.names=NA)
print(paste0("Expression file was created at ", getwd(), "/", extracted_processed_expression_file))