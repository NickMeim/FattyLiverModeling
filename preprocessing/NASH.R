# R script to download selected samples
# Copy code and run on a local machine to initiate download

library("rhdf5")    # can be installed using Bioconductor

destination_file = "human_gene_v2.2.h5"
extracted_expression_file = "NASH_expression_matrix.tsv"
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

