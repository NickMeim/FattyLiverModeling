# R script to download selected samples
# Copy code and run on a local machine to initiate download

library(rhdf5)    # can be installed using Bioconductor
library(tidyverse)
library(preprocessCore)
library(sva)
library(ggplot2)

destination_file = "../data/human_gene_v2.2.h5"
extracted_expression_file = "../data/ARCHS4/NAFL_NASH_expression_matrix.tsv"
extracted_processed_expression_file = "../preprocessing/ARCHS4_processed/NAFL_NASH_expression_matrix_preprocessed.tsv"
meta_file <- "../data/ARCHS4/NAFL_NASH_meta_data.tsv"
url = "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    download.file(url, destination_file, quiet = FALSE, mode = 'wb')
}

# Selected samples to be extracted
samp = c("GSM3998167","GSM3998168","GSM3998169","GSM3998170","GSM3998171","GSM3998172","GSM3998173","GSM3998174","GSM3998175","GSM3998176","GSM3998177","GSM3998178","GSM3998179","GSM3998180","GSM3998181","GSM3998182","GSM3998183","GSM3998184","GSM3998185","GSM3998186","GSM3998187","GSM3998188","GSM3998189","GSM3998190","GSM3998191","GSM3998192","GSM3998193","GSM3998194","GSM3998195","GSM3998196","GSM3998197",
"GSM3998198","GSM3998199","GSM3998200","GSM3998201","GSM3998202","GSM3998203","GSM3998204","GSM3998205","GSM3998206","GSM3998207","GSM3998208","GSM3998209","GSM3998210","GSM3998211","GSM3998212","GSM3998213","GSM3998214","GSM3998215","GSM3998219","GSM3998225","GSM3998226","GSM3998227","GSM3998228","GSM3998229","GSM3998230","GSM3998231","GSM3998232","GSM3998233","GSM3998234","GSM3998235",
"GSM3998236","GSM3998237","GSM3998238","GSM3998239","GSM3998240","GSM3998241","GSM3998242","GSM3998243","GSM3998244","GSM3998245","GSM3998246","GSM3998247","GSM3998248","GSM3998249","GSM3998250","GSM3998251","GSM3998252","GSM3998253","GSM3998254","GSM3998255","GSM3998256","GSM3998257","GSM3998258","GSM3998259","GSM3998260","GSM3998261","GSM3998262","GSM3998263","GSM3998264","GSM3998265",
"GSM3998266","GSM3998267","GSM3998268","GSM3998269","GSM3998270","GSM3998271","GSM3998272","GSM3998273","GSM3998274","GSM3998275","GSM3998276","GSM3998277","GSM3998278","GSM3998279","GSM3998280","GSM3998281","GSM3998282","GSM3998283","GSM3998284","GSM3998285","GSM3998286","GSM3998287","GSM3998288","GSM3998289","GSM3998290","GSM3998291","GSM3998292","GSM3998293","GSM3998294","GSM3998295",
"GSM3998296","GSM3998297","GSM3998298","GSM3998299","GSM3998300","GSM3998301","GSM3998302","GSM3998303","GSM3998304","GSM3998305","GSM3998306","GSM3998307","GSM3998308","GSM3998309","GSM3998310","GSM3998311","GSM3998312","GSM3998313","GSM3998314","GSM3998315","GSM3998316","GSM3998317","GSM3998319","GSM3998320","GSM3998321","GSM3998322","GSM3998323","GSM3998324","GSM3998325","GSM3998326",
"GSM3998327","GSM3998328","GSM3998329","GSM3998330","GSM3998331","GSM3998332","GSM3998333","GSM3998334","GSM3998335","GSM3998336","GSM3998337","GSM3998338","GSM3998339","GSM3998340","GSM3998343","GSM3998344","GSM3998345","GSM3998346","GSM3998347","GSM3998348","GSM3998349","GSM3998350","GSM3998351","GSM3998352","GSM3998353","GSM3998354","GSM3998355","GSM3998356","GSM3998357","GSM3998358",
"GSM3998359","GSM3998361","GSM3998362","GSM3998363","GSM3998364","GSM3998365","GSM3998366","GSM3998367","GSM3998368","GSM3998369","GSM3998370","GSM3998371","GSM3998372","GSM3998373","GSM3998374","GSM3998375","GSM3998376","GSM3998377","GSM3998378","GSM3998379","GSM3998380","GSM3998381","GSM3998382",
"GSM3946330","GSM3946331","GSM5474289","GSM5474290","GSM5474291","GSM5474292","GSM5474293","GSM5474294","GSM5474295","GSM5474296","GSM5474297","GSM5474298","GSM5474299","GSM5474300","GSM5474301","GSM5474302","GSM5474303","GSM5474304","GSM5474305","GSM5474306","GSM5474313","GSM5474314","GSM5474315","GSM5474316","GSM5474317","GSM5474318","")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/symbol")
genes <- unique(genes)

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
H5close()
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
        legend.position = 'top')
# Print pre-processed file
write.table(correctedExpression, file=extracted_processed_expression_file, sep="\t", quote=FALSE, col.names=NA)
print(paste0("Expression file was created at ", getwd(), "/", extracted_processed_expression_file))
hist(correctedExpression)
