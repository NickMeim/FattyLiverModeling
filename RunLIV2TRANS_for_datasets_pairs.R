### FIRST YOU NEED TO INSTALL LIV2TRANS WITH: 
### remotes::install_github("NickMeim/FattyLiverModeling", subdir = "src/R/LIV2Trans")
### Before that make sure you install the dependencies as described in :
### https://github.com/NickMeim/FattyLiverModeling/tree/main/src/R

library(tidyverse)
library(LIV2Trans)
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
in_vitro_datasets <- c("Kostrzewski","Wang", "Feaver")
in_vivo_datasets <- c("Govaere","Hoang")
dataset_names <- c(in_vivo_datasets,in_vitro_datasets)
# Load
data_list <- load_datasets(dataset_names, dir_data = 'data/')
# Manually load also the other clinical datasets I have from ARCHS4
geo <- 'GSE162694' # only this one has multiple NAS and fibrosis scores
meta_data <- read.delim('data/ARCHS4_Pantano_and_more/FattyLiver_meta_data.tsv',row.names = 1)
meta_data <- meta_data %>% filter(series_id==geo)
old_cols <- colnames(meta_data)
meta_data <- meta_data %>% separate_rows(characteristics_ch1, sep = ",") %>%
  separate(characteristics_ch1, into = c("key", "value"), sep = ":") %>%
  mutate(key = str_trim(key), value = str_trim(value)) %>%
  spread(key, value)
new_cols <- colnames(meta_data)
new_cols <- new_cols[which(!(new_cols %in% old_cols))]
meta_data <- meta_data %>%
  mutate_at(vars(new_cols), ~ ifelse(grepl("\\d", .), as.numeric(.), .))
meta_data <- meta_data %>% filter(!(is.na(`nas score`) & is.na(`fibrosis stage`))) %>%
  mutate(`nas score`=ifelse(`nas score`=='NA',NA,`nas score`)) %>% 
  mutate(`fibrosis stage`=ifelse(`fibrosis stage`=='normal liver histology',0,`fibrosis stage`))
meta_data <- meta_data %>% filter(!is.na(`nas score`)) %>% filter(!is.na(`fibrosis stage`)) %>%
  filter(`nas score`!='NA') %>% filter(`fibrosis stage`!='NA')
expression_matrix <- readRDS('data/ARCHS4_Pantano_and_more/FattyLiver_expression_matrix.rds')
expression_matrix <- expression_matrix[,meta_data$sample]
data_list[['Pantano']] <- list(counts = expression_matrix,
                               metadata = meta_data,
                               genes = rownames(expression_matrix))
in_vivo_datasets <- c(in_vivo_datasets,c("Pantano"))
dataset_names <- c(in_vivo_datasets,in_vitro_datasets)
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list

### Run LIV2TRANS for all combos------------
pathway_results <- data.frame()
hallmark_results <- data.frame()
for (ref_dataset in in_vivo_datasets){
  for (target_dataset in in_vitro_datasets){
    print(paste0('Begun: In vivo : ',ref_dataset,' , In vitro : ',target_dataset))
    ## Define matrices of interest
    if (ref_dataset=='Govaere'){
      Yh <- as.matrix(data_list$Govaere$metadata  %>% dplyr::select(nas_score,Fibrosis_stage))
    }else if (ref_dataset=='Hoang'){
      Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% dplyr::select(nafld_activity_score,Fibrosis_stage))
    }else if (ref_dataset=='Pantano'){
      Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% dplyr::select(`nas score`,`fibrosis stage`))
      Yh <- apply(Yh,c(1,2),as.numeric)
    }
    colnames(Yh) <- c('NAS','fibrosis')
    Xh <- data_list[[ref_dataset]]$data_center %>% t()
    Xm <- data_list[[target_dataset]]$data_center %>% t()
    # Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
    Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()
    
    ## Run main function with all the data
    liv2trans_results <- liv2trans_run(Xh,Yh,Xm,Wm)
    saveRDS(liv2trans_results,paste0('results/liv2trans_results_',ref_dataset,'_',target_dataset,'.rds'))
    W_opt <- liv2trans_results$W_opt
    
    ### Visualize pathway activity
    pathway_interpretation_results <- pathway_activity_interpretation(W_opt,Wm,plotting=FALSE)
    pathway_results <- rbind(pathway_results,
                             pathway_interpretation_results[[2]] %>%
                               mutate(combo = paste0('In vivo:',
                                                     ref_dataset,
                                                     ' , In vitro:',
                                                     target_dataset)))
    
    
    ### Visualize Hallmark genesets enrichment
    entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(W_opt), column = "ENTREZID", keytype = "SYMBOL")
    entrez_ids <- unname(entrez_ids)
    inds <- which(!is.na(entrez_ids))
    entrez_ids <- entrez_ids[inds]
    meas <- as.matrix(W_opt[inds,])
    rownames(meas) <- entrez_ids
    hallmarks_interpretation_results <- HallmarksFastenrichment(signature_ids = colnames(W_opt),
                                                                gene_ids = rownames(meas),
                                                                measurements = meas,
                                                                axis_lab = 'LV extra',
                                                                n_permutations=10000,
                                                                plotting=FALSE)
    hallmark_results <- rbind(hallmark_results,
                              hallmarks_interpretation_results$enrich_results %>%
                               mutate(combo = paste0('In vivo:',
                                                     ref_dataset,
                                                     ' , In vitro:',
                                                     target_dataset)))
    print(paste0('Finished: In vivo : ',ref_dataset,' , In vitro : ',target_dataset))
  }
}

# saveRDS(pathway_results,'results/alldata_combos_pathway_results.rds')
# saveRDS(hallmark_results,'results/alldata_combos_hallmark_results.rds')

### Visualize the results ---------------------
pathway_results <- readRDS('results/alldata_combos_pathway_results.rds')
hallmark_results <- readRDS('results/alldata_combos_hallmark_results.rds')
i <- 1
lim <- 10
pathway_odrer <- c("JAK-STAT",
                   "Hypoxia",
                   "EGFR",
                   "WNT",
                   "MAPK",
                   "NFkB",
                   "Androgen",
                   "VEGF",
                   "TNFa",
                   "Trail",
                   "PI3K",
                   "Estrogen",
                   "TGFb",
                   "p53") 
pathway_results$Pathway <- factor(pathway_results$Pathway,levels = pathway_odrer)
p1 <- ggplot(pathway_results %>% dplyr::select(c('activity'='score'),Pathway,p_value,condition,combo) %>%
               filter (condition==paste0('V',i)),
             aes(x=activity,y=Pathway,fill=activity)) + 
  geom_bar(stat = 'identity') +
        scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                             midpoint = 0,limits = c(-lim,lim))+
        scale_x_continuous(n.breaks = 8,limits = c(-lim,lim))+
        ggtitle(paste0('LV extra ',i))+
        ylab('Pathway')+
        geom_text(aes(label = ifelse(p_value <= 0.0001, "****",
                                     ifelse(p_value <= 0.001,"***",
                                            ifelse(p_value<=0.01,"**",
                                                   ifelse(p_value<=0.05,'*',
                                                          ifelse(p_value<=0.1,'\u2219',
                                                                 'ns'))))),
                      x = ifelse(activity < 0, activity - 0.2, activity + 0.2)),
                  size = 4,
                  color = 'black',
                  angle=90) +
  facet_wrap(~combo)+
        theme_minimal(base_size = 20,base_family = 'Arial')+
        theme(text = element_text(size = 20,family = 'Arial'),
              legend.position = 'right',
              plot.title = element_text(hjust = 0.5))
print(p1)
ggsave(filename = "figures/Supplementary_figure_other_invivo_invitro_pairs.png", 
       # device = cairo_pdf,
       plot = p1, 
       units = "in", 
       width = 16, 
       height = 12,
       dpi = 600)
