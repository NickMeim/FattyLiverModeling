library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ggsignif)
library(patchwork)
library(caret)
library(ropls)
source("../utils/plotting_functions.R")
source("functions_translation.R")
source("CrossValidationUtilFunctions.R")

### Load in-vitro and in-vivo datasets-----------------------------------------
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
  
### Load found extra basis from the Govaere and Kostrzewski datasets and check------------------------------
### if performance can be retrieved when backprojecting for other human datasets
external_clinical_datasets <- c("Hoang","Pantano")
df_results <- data.frame()
for (dataset in external_clinical_datasets){
  message(paste0('Begun ',dataset ,' dataset'))
  X <- data_list[[dataset]]$data_center %>% t()
  if (dataset == 'Hoang'){
    Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(nafld_activity_score,Fibrosis_stage))

  }else if (dataset == 'Pantano'){
    Y <- as.matrix(data_list[[dataset]]$metadata  %>% select(`nas score`,`fibrosis stage`))
    Y <- apply(Y,c(1,2),as.numeric)
  }
  colnames(Y) <- c('NAS','fibrosis')
  ### Run a small 5-fold cross validation procedure for 
  data_splits <- createMultiFolds(y = Y[,2], k = 5, times = 1)
  train_r <- NULL
  val_r <- NULL
  train_r_backproj <- NULL
  val_r_backproj <- NULL
  train_r_extra <- NULL
  val_r_extra <- NULL
  j <- 1
  for (ind in data_splits){
    y_val <- Y[-ind,]
    y_train <- Y[ind,]
    x_val <- X[-ind,]
    x_train <- X[ind,]
    plsr_model <- opls(x = x_train,
                       y = y_train,
                       predI = 8,
                       crossvalI = 1,
                       scaleC = "center",
                       fig.pdfC = "none",
                       info.txtC = "none")
    y_train_hat <- predict(plsr_model,x_train)
    y_val_hat <- predict(plsr_model,x_val)
    train_r[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r[j] <- mean(diag(cor(y_val_hat,y_val)))
    
    # Get Wh of PLSR
    Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
    rownames(Wh) <- colnames(x_train)
    colnames(Wh) <- colnames(plsr_model@weightMN)
    for (ii in 1:nrow(plsr_model@weightMN)){
      Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
    }
    # Get regression coefficients
    Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
    # Define projection matrices to make more readable
    Th_train <- x_train %*% Wh
    Thm_train <- x_train %*% Wm %*% t(Wm) %*% Wh
    Th_val<- x_val %*% Wh
    Thm_val <- x_val %*% Wm %*% t(Wm) %*% Wh
    
    y_train_hat <- cbind(1, x_train %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    y_val_hat <- cbind(1, x_val %*% Wm_tot %*% t(Wm_tot) %*% Wh) %*% rbind(apply(y_train,2,mean),Bh)
    
    train_r_extra[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r_extra[j] <- mean(diag(cor(y_val_hat,y_val)))
    
    ### simple backprojection 
    y_val_hat <- cbind(1, Thm_val)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    y_train_hat <- cbind(1, Thm_train)  %*% rbind(apply(y_train,2,mean),t(plsr_model@weightMN) %*% plsr_model@coefficientMN)
    train_r_backproj[j] <- mean(diag(cor(y_train_hat,y_train)))
    val_r_backproj[j] <- mean(diag(cor(y_val_hat,y_val)))
    
    print(paste0('Finished ',j,' fold out of ',length(data_splits)))
    j <- j+1
  }
  tmp_train <- data.frame(PLSR = train_r,PC = train_r_backproj,'extra basis'=train_r_extra) %>% 
    mutate(dataset = dataset) %>% mutate(set = 'train')
  tmp_val <- data.frame(PLSR = val_r,PC =val_r_backproj,'extra basis'=val_r_extra) %>% 
    mutate(dataset = dataset) %>% mutate(set = 'test')
  df_results <- rbind(df_results,tmp_train,tmp_val)
}
df_results <- df_results %>% gather('model','r',-dataset,-set)
df_results <- df_results %>% mutate(model = ifelse(model=='extra.basis','extra basis',
                                                   ifelse(model!='PLSR','in-vitro PCs','original PLSR')))
df_results$model <- factor(df_results$model,levels = c('original PLSR','extra basis','in-vitro PCs'))


### Visualize 
model_comparisons <- list(
  c('original PLSR', 'extra basis'),
  c('extra basis', 'in-vitro PCs')
)
p <- ggboxplot(df_results,x='set',y='r',color = 'model',add='jitter') +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,NA))+
  xlab('')+
  theme(text = element_text(size = 20,family = 'Arial'),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.minor.y = element_line(linewidth = 0.5))+
  facet_wrap(~dataset)+
  stat_compare_means(comparisons = list(c('original PLSR','extra basis'),
                                        c('extra basis','in-vitro PCs')),
                     method = 'wilcox')
print(p)
ggsave('../results/inter_clinical_dataset_performance.png',
       plot=p,
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)



















