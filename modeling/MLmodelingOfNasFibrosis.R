library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(caret)
source("../utils/plotting_functions.R")
source("functions_translation.R")
source("CrossValidationUtilFunctions.R")
memory.limit(size = 32000)

### Load datasets and learned extra basis---------------------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = '../data/')
# Manually load also the other clinical datasets I have from ARCHS4
geo <- 'GSE162694' # only this one has multiple NAS and fibrosis scores
meta_data <- read.delim('../data/ARCHS4/FattyLiver_meta_data.tsv',row.names = 1)
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
expression_matrix <- readRDS('../data/ARCHS4/FattyLiver_expression_matrix.rds')
expression_matrix <- expression_matrix[,meta_data$sample]
data_list[['Pantano']] <- list(counts = expression_matrix,
                               metadata = meta_data,
                               genes = rownames(expression_matrix))


# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list$Govaere$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

# Load previously found extra basis
Wm_opt <- readRDS(paste0('../results/Wm_',target_dataset,'_extra.rds'))
Wm_tot <- readRDS(paste0('../results/Wm_',target_dataset,'_total.rds'))

### Train ML models different from PLSR------------------------
models <- c('knn','lasso','ridge','elasticnet','lm',
            'svmRadial','svmPoly','svmLinear',
            'gaussprRadial','gaussprPoly','gaussprLinear',
            'xgbTree','rf','neuralnet')
mdl_results <- GeneralMLCompletePipeLineCrossVal(model = 'lm',
                                                 cv_location='../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/',
                                                 external_datasets=c("Hoang","Pantano"),
                                                 data_list_all=data_list,
                                                 Wm=Wm,
                                                 Wtot=Wm_tot,
                                                 pheno = 'NAS', # can be 'NAS' or 'fibrosis' (or both for ANN)
                                                 num_folds=10,
                                                 shuffling = NULL # can be only one of NULL,'X','Y'
                                                 )
