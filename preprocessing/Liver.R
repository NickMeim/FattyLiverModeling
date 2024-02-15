# R script to download selected samples
# Copy code and run on a local machine to initiate download

library(rhdf5)    # can be installed using Bioconductor
library(tidyverse)
library(preprocessCore)
library(sva)
library(ggplot2)

destination_file = "../data/human_gene_v2.2.h5"
extracted_expression_file = "../data/ARCHS4/FattyLiver_expression_matrix.tsv"
extracted_expression_fileRDS = "../data/ARCHS4/FattyLiver_expression_matrix.rds"
extracted_processed_expression_file = "../preprocessing/ARCHS4_processed/FattyLiver_expression_matrix_preprocessed.tsv"
meta_file <- "../data/ARCHS4/FattyLiver_meta_data.tsv"
url = "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    download.file(url, destination_file, quiet = FALSE, mode = 'wb')
}

# Relevent NAFL-NASH GEO Accession numbers
geo_nafl_nash <- read.delim2('../data/ARCHS4/NAFL_NASH_meta_data.tsv',row.names = 1)
geo_nafl_nash <- geo_nafl_nash %>% select(series_id) %>% unique()
mannual_geos <- c('GSE162694','GSE130970','GSE126848')
geos <- unique(c(geo_nafl_nash$series_id,mannual_geos))
exclude_geos <- c('GSE180882')
geos <- geos[which(!(geos %in% exclude_geos))]  

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
H5close()
genes = h5read(destination_file, "meta/genes/symbol")
genes <- unique(genes)
H5close()

# Load metadata for these samples
meta <- h5read(destination_file, "meta")
H5close()
meta <- data.frame(meta[['samples']])
meta <- meta %>% filter(series_id %in% geos)
meta <- meta %>% mutate(keep = ifelse(series_id=='GSE134422' & grepl('AlcoholicCirrhosis',title),'no','yes')) %>% filter(keep=='yes') %>% select(-keep)
meta <- meta %>% mutate(keep = ifelse(series_id=='GSE126848' & grepl('Obese',title),'no','yes')) %>% filter(keep=='yes') %>% select(-keep)
write.table(meta, file=meta_file, sep="\t", quote=FALSE, col.names=NA)
print(paste0("Expression file was created at ", getwd(), meta_file))

# Identify columns to be extracted
sample_locations = which(samples %in% meta$sample)


# extract gene expression from compressed data
expression = t(h5read(destination_file, "data/expression", index=list(sample_locations, 1:length(genes))))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE, col.names=NA)
saveRDS(expression,extracted_expression_fileRDS)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_fileRDS))

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


