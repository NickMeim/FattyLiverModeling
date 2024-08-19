library(tidyverse)
library(LIVIVTRA)
library(AnnotationDbi)
library(org.Hs.eg.db)
# Function for loading datasets and appending them to a list
load_datasets <- function(dataset_names, dir_data){
  # Look at files in data directory
  files <- dir(dir_data)
  # Trim to datasets
  files <- files[grepl("dataset.RData", files)]
  # Loop for dataset names given and load one by one, appending to a list.
  # dataset names don't have to match perfectly
  data_list <- list(NULL)
  if (length(files) > 0){ # Loop if any match was found
    for (name in dataset_names){
      if (any(grepl(name, files, ignore.case = TRUE))){ # Avoid loading bad matches
        load(file = paste0(dir_data, files[grep(name, files, ignore.case = TRUE)]))
        print(paste0("Loading file: ", paste0(dir_data, files[grep(name, files, ignore.case = TRUE)])))
        data_list[[name]]$counts <- data
        data_list[[name]]$metadata <- metadata
        data_list[[name]]$genes <- rownames(data)
        if (exists("exp_factors")){
          data_list[[name]]$exp_factors <- exp_factors
        }
      }
    }
    # Remove redundant first element
    data_list[[1]] <- NULL
  }

  return(data_list)

} # End function
### Load the all the data to be used-----------------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../../data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
if (ref_dataset=='Hoang'){
  sex_inferred <- data_list[[ref_dataset]]$metadata$sex
}else if (ref_dataset=='Pantano') {
  sex_inferred <- data_list[[ref_dataset]]$metadata$Sex
}else{
  sex_inferred <- apply(as.matrix(Xh[,c('RPS4Y1')]),2,sign)
  sex_inferred <- 1*(sex_inferred>0)
  sex_inferred <- ifelse(sex_inferred==1,'male','female')
}
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

## Run main function with all the data---------------
livivtra_results <- livivtra_run(Xh,Yh,Xm,Wm)
W_translatable <- livivtra_results$W_translatable
W_opt <- livivtra_results$W_opt
W_tot <- livivtra_results$W_tot

### Visualize gene loadings-------------------------
## valta mesa sta functions ayta
library(ggrepel)
library(ggpubr)
p1 <- plot_gene_loadings(loadings = W_opt,
                             selection='V1',
                             y_lab = 'weight in extra LV1',
                             top=20)
p2 <- plot_gene_loadings(loadings = W_opt,
                         selection='V2',
                         y_lab = 'weight in extra LV2',
                         top=20)

### Visualize pathway activity-------------------------
pathway_interpretation_results <- pathway_activity_interpretation(W_opt,Wm)

### Visualize Hallmark genesets enrichment-------------------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(W_opt), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(W_opt[inds,])
rownames(meas) <- entrez_ids
hallmarks_interpretation_results <- HallmarksFastenrichment(signature_ids = colnames(W_opt),
                                                          gene_ids = rownames(meas),
                                                          measurements = meas,
                                                          order_columns = F, # REMOVE THAT FROM THE FUNCTION
                                                          n_permutations=10000)

p1 <- (ggplot(df_msig %>% filter(LV=='V1') %>% arrange(NES) %>%
                filter(padj<=0.1),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV extra 1')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))

## Run cross validation---------------
livivtra_CV_results <- livivtra_run_CV(Xh,Yh,Xm,Wm,num_folds = 10)
