library(tidyverse)
# library(dorothea)
library(ggrepel)
library(ggpubr)

### Load relevent TF perturbations---------------------------
tf_responses <- readRDS("../data/ChemPert/Transcriptional_responses.rds")
colnames(tf_responses)[1] <- "Response_ID"
metadata <- read.csv("../data/ChemPert/Information_for_transcriptional_responses.csv")
metadata_human <- metadata %>% filter(Species == "Human" & grepl("Hepatocy", Cell_Source))
tf_responses_hep <- tf_responses %>% filter(Response_ID %in% metadata_human$Response_ID)
resp_net <- NULL
for (ii in 1:nrow(tf_responses_hep)){
  if (tf_responses_hep$`Up-regulation number`[ii]>0){
    dummy <- data.frame(source = tf_responses_hep$Response_ID[ii],
                        target = strsplit(tf_responses_hep$`Up-regulation`[ii], split = "; ", fixed = T)
                        %>% unlist(),
                        mor = 1)
  }
  if (tf_responses_hep$`Down-regulation number`[ii]>0){
    dummy <- rbind(dummy,
                   data.frame(source = tf_responses_hep$Response_ID[ii],
                              target = strsplit(tf_responses_hep$`Down-regulation`[ii], split = "; ", fixed = T)
                              %>% unlist(),
                              mor = -1))
  }
  resp_net <- rbind(resp_net, dummy)
}
resp_net <- resp_net %>% mutate(Response_ID = source)
resp_net <- merge(resp_net, metadata_human, by = "Response_ID")

### Load PCA loadings--------------------------------------------
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_A <- data$data_A
data_B <- data$data_B
Y_A <- data$Y_A
## normalize
# Intersect datasets - keep common genes
genes_common <- intersect(rownames(data_A), rownames(data_B))
# Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
data_A <- apply(data_A[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_B <- apply(data_B[genes_common,], MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
# Remove features absent from at least 10% in each samples
keep_gene <- (rowSums(data_A >= 1) >= 0.1*ncol(data_A)) &
  (rowSums(data_B >= 1) >= 0.1*ncol(data_B))
# Run PCA on each dataset. Not scaling PCA
# Log and center dataset A and run PCA
X_A <- log2(1 + data_A[keep_gene,])
X_A <- t(X_A - rowMeans(X_A))
# Log and center dataset B and run PCA
X_B <- log2(1 + data_B[keep_gene,])
X_B <- t(X_B - rowMeans(X_B))
### Perform PCA and get loadings
PCA_alldata <- prcomp(X_B, scale. = F, center = T)
gene_loadings <- PCA_alldata$rotation[,c('PC1','PC2','PC8','PC12')]

### Load optimal gene loading and infer TF activities-----------------------------------
Woptim <- readRDS('../results/Wm_opt.rds')
# regulon <- decoupleR::get_dorothea(levels = c("A","B"))
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
minNrOfGenes  <-  5
TF_activities = decoupleR::run_viper(Woptim, dorotheaData,minsize = minNrOfGenes,verbose = TRUE)
# TF_activities <- as.matrix(TF_activities[,1])
# colnames(TF_activities) <- 'optimal direction'
# TF_activities = run_viper(gene_loadings, dorotheaData, options =  settings)

### Get perturbation activity-------------------------
pert_activity <- decoupleR::run_viper(TF_activities %>% select(-statistic,-condition,-p_value) %>% column_to_rownames('source'), 
                                      resp_net %>% unique(),
                                      minsize = 1,
                                      verbose = TRUE)
# pert_activity <- as.matrix(pert_activity[,1])
# colnames(pert_activity) <- 'optimal direction'
# pert_activity <- as.data.frame(pert_activity) %>% rownames_to_column('Response_ID')
colnames(pert_activity)[2] <- 'Response_ID'
pert_activity <- left_join(pert_activity,metadata_human)
pert_activity <- pert_activity %>% select(Chemical_Compound,c('activity'='score'),Duration,Concentration,c('p-value'='p_value'))
pert_activity <- pert_activity %>% mutate(perturbation = ifelse(!is.na(Concentration),
                                                                paste0(toupper(Chemical_Compound),':',Duration,':',Concentration),
                                                                paste0(toupper(Chemical_Compound),':',Duration)))
pert_activity <- pert_activity %>% group_by(perturbation) %>% mutate(min_pvalue=min(`p-value`)) %>%
  ungroup() %>% filter(`p-value`==min_pvalue) %>% select(-min_pvalue) %>% unique()
pert_activity$significant <- ''
pert_activity <- pert_activity %>% mutate(significant = ifelse(`p-value`<=0.05,
                                                               perturbation,
                                                               significant))
if (nrow(pert_activity %>% filter(`p-value`<=0.05))<20){
  top_20s <- pert_activity$perturbation[order(-abs(pert_activity$activity))[1:20]]
  pert_activity <- pert_activity %>% mutate(significant = ifelse((perturbation %in% top_20s) | (`p-value`<=0.05),
                                                                 perturbation,
                                                                 significant))
}
# top_20s <- pert_activity$perturbation[order(-abs(pert_activity$activity))[1:20]]
# pert_activity <- pert_activity %>% mutate(significant = ifelse((perturbation %in% top_20s) & (`p-value`<=0.05),
#                                                                perturbation,
#                                                                significant))
# pert_activity$significant[order(-abs(pert_activity$activity))[1:20]] <- pert_activity$perturbation[order(-abs(pert_activity$activity))[1:20]]
pert_activity <- pert_activity[order(pert_activity$activity),]
pert_activity$perturbation <- factor(pert_activity$perturbation,levels = unique(pert_activity$perturbation))
ggplot(pert_activity,aes(x=perturbation,y=activity,fill = `p-value`)) + geom_point(shape=21,size=2) +
  geom_text_repel(aes(label=significant,color=`p-value`),size=5,max.overlaps=60,box.padding = 0.7)+
  xlab('perturbations') + ylab('inferred activity')+
  scale_fill_gradient(trans='log10',low = "red",high = "white",
                       limits = c(min(pert_activity$`p-value`),1),
                      breaks = c(0.01,0.05,0.1,0.5,1)) +
  scale_color_gradient(trans='log10',low = "red",high = "#f96464",
                      limits = c(min(pert_activity$`p-value`),0.05)) +
  guides(color = "none")+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(text = element_text(family = 'Arial',size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'right')
ggsave('../results/pc_loadings_scores_analysis/optimal_direction_perturbation_activity.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)
