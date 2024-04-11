# load in packages
library(tidyverse)
library(DESeq2)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(pls)
library(pheatmap)
library(patchwork)
library(limma)
library(edgeR)
# library(dorothea)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(hgu133a.db)
library(rstatix)
library(fgsea)
library(topGO)
library(GO.db)
# library(progeny)
library(OmnipathR)
library(EGSEAdata)
library(nichenetr)

#### Load pre-processed data----------------------
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_A <- data$data_A
data_B <- data$data_B
Y_A <- data$Y_A

### read optimal direction
Woptim <- readRDS('../results/Wm_opt.rds')
colnames(Woptim) <- 'w_optimal'

## normalize
## I do not think I need any more common genes to interogate
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
# Transform to log2(cpm + 1) and keep features that are present in at least 10% of each dataset
# data_B <- apply(data_B, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
# Remove features absent from at least 10% in each samples
# keep_gene <- rowSums(data_B >= 1) >= 0.1*ncol(data_B)
# Run PCA on each dataset. Not scaling PCA
# Log and center dataset B and run PCA
# keep_gene <- rownames(Woptim)[which(rownames(Woptim) %in% rownames(data_B))]
# X_B <- log2(1 + data_B[keep_gene,])
# X_B <- t(X_B - rowMeans(X_B))

Woptim <- as.data.frame(Woptim) %>% rownames_to_column('gene')
Woptim$significant <- ''
Woptim$significant[order(-abs(Woptim$w_optimal))[1:20]] <- Woptim$gene[order(-abs(Woptim$w_optimal))[1:20]]
Woptim <- Woptim[order(Woptim$w_optimal),]
Woptim$gene <- factor(Woptim$gene,levels = Woptim$gene)

### Perform PCA and get loadings----------------
PCA_alldata <- prcomp(X_B, scale. = F, center = T)
loadings <- PCA_alldata$rotation
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,PC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$PC12))[1:10]] <- loadings$gene[order(-abs(loadings$PC12))[1:10]]
loadings <- loadings[order(loadings$PC12),]
loadings$gene <- factor(loadings$gene,levels = loadings$gene)
### Make barplot to look at top genes----------------
ggplot(Woptim,aes(x=gene,y=w_optimal,color = significant)) + geom_point() +
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
  xlab('genes') + ylab('optimal direction loadings')+
  scale_x_discrete(expand = c(0.1, 0.1))+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(text = element_text(family = 'Arial',size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none')
ggsave('../results/pc_loadings_scores_analysis/gene_woptimal_loadings.png',
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

### Infer statistical significance for gene loadings of PC12 and PC8----------------------------------
gene_loadings <- PCA_alldata$rotation
gene_loadings <- as.data.frame(gene_loadings) %>% rownames_to_column('gene')
null_loadings <- gene_loadings %>% gather('PC','null_act',-gene) %>% select(-PC)
gene_loadings$significant <- ''
gene_loadings_PC8 <- left_join(gene_loadings %>% select(gene,PC8,significant),
                               null_loadings) %>%
  group_by(gene) %>% mutate(p.value = ifelse(PC8>=0,
                                             sum(null_act>=PC8)/(ncol(gene_loadings)-2),
                                             sum(null_act<=PC8)/(ncol(gene_loadings)-2))) %>%
  ungroup()
adjustement_tmp <- gene_loadings_PC8 %>% select(gene,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'fdr')
adjustement_tmp <- adjustement_tmp %>% select(gene,p.adj)
gene_loadings_PC8 <- left_join(gene_loadings_PC8,adjustement_tmp)%>% select(-null_act) %>% unique()
gene_loadings_PC12 <- left_join(gene_loadings %>% select(gene,PC12,significant),
                                null_loadings) %>%
  group_by(gene) %>% mutate(p.value = ifelse(PC12>=0,
                                             sum(null_act>=PC12)/(ncol(gene_loadings)-2),
                                             sum(null_act<=PC12)/(ncol(gene_loadings)-2))) %>%
  ungroup()
adjustement_tmp <- gene_loadings_PC12 %>% select(gene,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'fdr')
adjustement_tmp <- adjustement_tmp %>% select(gene,p.adj)
gene_loadings_PC12 <- left_join(gene_loadings_PC12,adjustement_tmp) %>% select(-null_act) %>% unique()
gene_loadings_PC8$significant[order(-abs(gene_loadings_PC8$PC8))[1:10]] <- gene_loadings_PC8$gene[order(-abs(gene_loadings_PC8$PC8))[1:10]]
gene_loadings_PC8 <- gene_loadings_PC8[order(gene_loadings_PC8$PC8),]
gene_loadings_PC8$gene <- factor(gene_loadings_PC8$gene,levels = gene_loadings_PC8$gene)
gene_loadings_PC12$significant[order(-abs(gene_loadings_PC12$PC12))[1:10]] <- gene_loadings_PC12$gene[order(-abs(gene_loadings_PC12$PC12))[1:10]]
gene_loadings_PC12 <- gene_loadings_PC12[order(gene_loadings_PC12$PC12),]
gene_loadings_PC12$gene <- factor(gene_loadings_PC12$gene,levels = gene_loadings_PC12$gene)
# change significant
# gene_loadings_PC12 <- gene_loadings_PC12 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                        ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
# gene_loadings_PC8 <- gene_loadings_PC8 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                      ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
gene_loadings_PC12 <- gene_loadings_PC12 %>% mutate(statistical=ifelse(p.value<=0.01,'p-value<=0.01',
                                                                       ifelse(p.value<=0.05,'p-value<=0.05','p-value>0.05')))
gene_loadings_PC8 <- gene_loadings_PC8 %>% mutate(statistical=ifelse(p.value<=0.01,'p-value<=0.01',
                                                                     ifelse(p.value<=0.05,'p-value<=0.05','p-value>0.05')))

# Combine into one data frame
gene_loadings_plot_frame <- rbind(gene_loadings_PC8 %>% dplyr::rename(value = PC8) %>% mutate(PC='PC8'),
                                  gene_loadings_PC12 %>% dplyr::rename(value = PC12) %>% mutate(PC='PC12'))
gene_loadings_plot_frame <- gene_loadings_plot_frame %>% mutate(PC=paste0('used ',PC,' gene loadings'))

p <- ggplot(gene_loadings_plot_frame %>% dplyr::rename(`statistical threshold`=statistical),
            aes(x=value,y=-log10(p.value),color = `statistical threshold`)) + geom_point() +
  scale_color_manual(values = c('red','#fa8e8e','black'))+
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40,point.padding = 0.5)+
  xlab('gene loading') + ylab(expression(-log[10]('p-value'))) +
  geom_hline(yintercept=-log10(0.01), linetype="dashed",color = "#525252", size=1) +
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'top')+
  facet_wrap(~PC,scales = 'free_x')+
  guides(color = guide_legend(
    override.aes = list(
      linetype = NA,
      size = 3
    )
  ))
print(p)
ggsave('../results/pc_loadings_scores_analysis/gene_loadings_volcano.png',
       plot = p,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

p2 <- (ggplot(gene_loadings_PC8,aes(x=gene,y=PC8,color = statistical)) + geom_point() +
         # scale_color_gradient(high = 'red',low='white')+
         scale_color_manual(values = c('red','#fa8e8e','black'))+
         geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
         xlab('TFs') + ylab('PC8 loadings')+
         scale_x_discrete(expand = c(0.1, 0.1))+
         theme_pubr(base_family = 'Arial',base_size = 20)+
         theme(text = element_text(family = 'Arial',size=20),
               axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               legend.position = 'none')) +
  (ggplot(gene_loadings_PC12,aes(x=gene,y=PC12,color = statistical)) + geom_point() +
     scale_color_manual(values =c('red','#fa8e8e','black'))+
     geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
     xlab('TFs') + ylab('PC12 loadings')+
     scale_x_discrete(expand = c(0.1, 0.1))+
     theme_pubr(base_family = 'Arial',base_size = 20)+
     theme(text = element_text(family = 'Arial',size=20),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = 'none'))
print(p2)
ggsave('../results/pc_loadings_scores_analysis/gene_loadings_points_withstats.png',
       plot = p2,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)
### Infer Pathway activity with progeny----------------------------------------------------------------------------
gene_loadings <- PCA_alldata$rotation[,paste0('PC',c(1,2,8,12))]
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
pathways_progeny_pcs <- decoupleR::run_mlm(mat=gene_loadings, net=net_prog, .source='source', .target='target',
                                           .mor='weight', minsize = 1)
# pathways_progeny_pcs <- progeny(gene_loadings,
#                                 perm = 10000,
#                                 z_scores = T,
#                                 top = 500)
# tt <- pathways_progeny_pcs[c('PC12','PC8'),]
Woptim <- readRDS('../results/Wm_opt.rds')
pathways_progeny <- decoupleR::run_mlm(mat=Woptim, net=net_prog, .source='source', .target='target',
                                           .mor='weight', minsize = 1)
# pathways_progeny <- progeny(Woptim,
#                             perm = 10000,
#                             z_scores = T,
#                             top = 500)
# rownames(pathways_progeny) <- 'optimal direction'
ggplot(pathways_progeny %>% select(c('activity'='score'),c('Pathway'='source'),p_value) %>% 
         mutate(Pathway=factor(Pathway,levels=pathways_progeny$source[order(pathways_progeny$score)])),
       aes(x=activity,y=Pathway,fill=activity)) + geom_bar(stat = 'identity',color='black') +
  scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
  geom_text(aes(label = ifelse(p_value <= 0.0001, "****",
                               ifelse(p_value <= 0.001,"***",
                                      ifelse(p_value<=0.01,"**",
                                             ifelse(p_value<=0.05,'*',
                                                    ifelse(p_value<=0.1,'\u2219',
                                                           'ns'))))),
                x = ifelse(activity < 0, activity - 0.2, activity + 0.2)),
            size = 6,
            color = 'black',
            angle=90) +                                   
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size = 20,family = 'Arial'))
ggsave('../results/pc_loadings_scores_analysis/progenies_only_optimal_loadings_barplot.png',
       width = 12,
       height = 9,
       dpi = 600)
pathways_progeny_all <- rbind(pathways_progeny,pathways_progeny_pcs)

### use progeny weights
# progeny_weights <- progeny::getModel()
progeny_weights <- net_prog %>% select(-p_value) %>% spread(target,weight) %>% replace(is.na(.), 0) %>% column_to_rownames('source')
progeny_weights <- t(progeny_weights)
common_genes_wopt_prog <- Reduce(intersect,list(rownames(progeny_weights), rownames(Woptim)))
progeny_weights_subset <- progeny_weights[common_genes_wopt_prog,]
Woptim_subset <- as.matrix(Woptim[common_genes_wopt_prog,])
woptim_pathways <- t(progeny_weights_subset) %*% Woptim_subset
colnames(woptim_pathways) <- 'optimal direction'
common_genes_pc_prog <- Reduce(intersect,list(rownames(progeny_weights), rownames(gene_loadings)))
progeny_weights_subset2 <- progeny_weights[common_genes_pc_prog,]
# gene_loadings_subset <- as.matrix(gene_loadings[common_genes_pc_prog,])
gene_loadings_subset <- as.matrix(PCA_alldata$rotation[common_genes_pc_prog,])
gene_loadings_pathways <- t(progeny_weights_subset2) %*% gene_loadings_subset
pathways_progeny_weights_all <- cbind(woptim_pathways,gene_loadings_pathways)
# colnames(pathways_progeny_all) <- rownames(pathways_progeny_weights_all)
# pheatmap(pathways_progeny_all)
df_pathways_progeny_weights_all <- as.data.frame(pathways_progeny_weights_all) %>% rownames_to_column('Pathway') %>%
  gather('direction','loading',-Pathway)
weight_distribution <- df_pathways_progeny_weights_all$loading
df_pathways_progeny_weights_all <- df_pathways_progeny_weights_all %>% group_by(direction,Pathway) %>%
  mutate(p_value_naive = sum(abs(weight_distribution)>=abs(loading))/length(weight_distribution)) %>% ungroup()

### Create Null distribution for Woptim and PCs by shuffling genes when performing the matrix multiplication
iters <- 10000
df_shuffled_optim <- data.frame()
df_shuffled_pcs <- data.frame()
for (i in 1:iters){
  ### For optimal vector
  shuffled <- as.matrix(Woptim_subset[sample.int(nrow(Woptim_subset)),])
  rownames(shuffled) <- rownames(Woptim_subset)
  woptim_pathways_shuffled <- t(progeny_weights_subset) %*% shuffled
  colnames(woptim_pathways_shuffled) <- 'score'
  df_shuffled_optim <- rbind(df_shuffled_optim,
                             as.data.frame(woptim_pathways_shuffled) %>% rownames_to_column('Pathway') %>% mutate(trial=i) %>%
                               mutate(direction = 'optimal direction') %>% select(Pathway,direction,score,trial))
  ### For PC loadings
  shuffled <- as.matrix(gene_loadings_subset[sample.int(nrow(gene_loadings_subset)),])
  rownames(shuffled) <- rownames(gene_loadings_subset)
  pc_pathways_shuffled <- t(progeny_weights_subset2) %*% shuffled
  df_shuffled_pcs <- rbind(df_shuffled_pcs,
                           as.data.frame(pc_pathways_shuffled)  %>% 
                             rownames_to_column('Pathway') %>% mutate(trial=i) %>% 
                             gather('direction','score',-Pathway,-trial) %>% select(Pathway,direction,score,trial))
  if (i %% 1000 == 0){
    print(paste0('Finished permutation ',i))
  }
}
df_all_nulls <- rbind(df_shuffled_optim,df_shuffled_pcs)
# saveRDS(df_all_nulls,'../results/pc_loadings_scores_analysis/random_progeny_loadings.rds')
df_all_nulls <- readRDS('../results/pc_loadings_scores_analysis/random_progeny_loadings.rds')

df_pathways_progeny_weights_all <- left_join(df_pathways_progeny_weights_all,df_all_nulls)
df_pathways_progeny_weights_all <- df_pathways_progeny_weights_all %>% group_by(direction,Pathway) %>%
  mutate(bool=abs(score)>=abs(loading)) %>% mutate(p_value_null = sum(bool)) %>%
  ungroup() %>% select(-bool) %>% mutate(p_value_null=p_value_null/iters)
df_pathways_progeny_weights_all  <- df_pathways_progeny_weights_all %>% 
  mutate(p.adj = p_value_null*length(unique(df_pathways_progeny_weights_all$Pathway)))
df_pathways_progeny_weights_all <- df_pathways_progeny_weights_all %>% group_by(Pathway) %>%
  mutate(sig1 = sum(p.adj<=0.05)) %>% ungroup() %>%
  group_by(direction) %>% mutate(sig2 = sum(p.adj<=0.05)) %>% ungroup()
keep_paths <- unique(df_pathways_progeny_weights_all$Pathway[which(df_pathways_progeny_weights_all$sig1!=0)])
keep_pcs <- unique(df_pathways_progeny_weights_all$direction[which(df_pathways_progeny_weights_all$sig2!=0)])
# pheatmap(t(pathways_progeny_weights_all[,c('optimal direction','PC12','PC8')]))
significant_pathways_progeny_weights_all <- pathways_progeny_weights_all[keep_paths,keep_pcs]
selected_row_names <- c("PC1", "PC8","PC12","optimal direction")
custom_labels <- rownames(t(pathways_progeny_weights_all)[c('PC1','PC11','PC12','PC13','PC35','PC8','optimal direction'),])
custom_labels[!custom_labels %in% selected_row_names] <- "" 
l <- distinct(df_pathways_progeny_weights_all %>% select(Pathway,direction,p.adj) %>%
  filter(Pathway %in% keep_paths) %>% filter(direction %in% keep_pcs)) %>% 
  mutate(p.adj = ifelse(p.adj<=0.0001,'****',
                        ifelse(p.adj<=0.001,'***',
                               ifelse(p.adj<=0.01,'**',
                               ifelse(p.adj<=0.05,'*',
                               ''))))) %>% spread(direction,p.adj) %>% column_to_rownames('Pathway')
png('../results/pc_loadings_scores_analysis/significant_pathWeights_subset_heatmap.png',width = 6,height = 6,units = 'in',res = 600)
pheatmap(t(significant_pathways_progeny_weights_all)[c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC11','PC12','PC13','PC35','PC8','optimal direction'),],
         display_numbers = t(l)[c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC11','PC12','PC13','PC35','PC8','optimal direction'),]) #,
         # labels_row = custom_labels)
dev.off()

### barplot weights
df_pathways_progeny_weights_all_visualize <- distinct(df_pathways_progeny_weights_all %>% select(Pathway,direction,loading,p.adj))
ggplot(df_pathways_progeny_weights_all_visualize %>% filter(direction=='optimal direction') %>%
         mutate(Pathway=factor(Pathway,
                               levels=df_pathways_progeny_weights_all_visualize[
                                 df_pathways_progeny_weights_all_visualize$direction=='optimal direction',]$Pathway[
                                 order(df_pathways_progeny_weights_all_visualize[
                                   df_pathways_progeny_weights_all_visualize$direction=='optimal direction',]$loading)])),
       aes(x=loading,y=Pathway,fill=loading)) + geom_bar(stat = 'identity',color='black') +
  scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
  geom_text(aes(label = ifelse(p.adj <= 0.0001, "****",
                               ifelse(p.adj <= 0.001,"***",
                                      ifelse(p.adj<=0.01,"**",
                                             ifelse(p.adj<=0.05,'*',
                                                    ifelse(p.adj<=0.1,'\u2219',
                                                           'ns'))))),
                x = ifelse(loading < 0, loading - 0.2, loading + 0.2)),
            size = 6,
            color = 'black',
            angle=90) +                                   
  theme_pubr(base_size = 20,base_family = 'Arial')+
  theme(text = element_text(size = 20,family = 'Arial'))
ggsave('../results/pc_loadings_scores_analysis/progeniesLoadings_optimal_barplot.png',
       width = 12,
       height = 9,
       dpi = 600)

#### Run CARNIVAL--------------------------
library(CARNIVAL)
gene_loadings <- PCA_alldata$rotation[,paste0('PC',c(1,2,3,8,12))]
minNrOfGenes  <-  5
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
settings = list(verbose = TRUE, minsize = minNrOfGenes)
## read optimal W
Woptim <- readRDS('../results/Wm_opt.rds')
TF_activities = run_viper(cbind(Woptim,Woptim), dorotheaData, options =  settings)
TF_activities <- as.matrix(TF_activities[,1])
tfs <- rownames(TF_activities)
## load omnipath network
interactions <- import_omnipath_interactions()
interactions <- interactions  %>% mutate(kegg=grepl(pattern="KEGG",x=sources))%>% 
  filter(n_resources>=3 | kegg==T)  %>% 
  filter(n_references>=3 | kegg==T)  %>% 
  dplyr::select(-kegg) %>%
  dplyr::select(c('source'='source_genesymbol'),
                c('target'='target_genesymbol'),
                is_inhibition,is_stimulation) %>% unique()

interactions <- interactions %>% filter(!(is_inhibition==0 & is_stimulation==0)) %>% unique()
interactions <- interactions %>% mutate(interaction=ifelse(is_stimulation!=0,1,-1)) %>%
  dplyr::select(source,interaction,target) %>% unique()
### Run carnival
log_con <- file("../results/pc_loadings_scores_analysis/log_cavinval.txt", open="a")
for (j in 1:ncol(TF_activities)){
  
  cat(paste0('Iteration ',j,'/',ncol(TF_activities)), file = log_con, sep="\n")
  tf_activities <- TF_activities[which(rownames(TF_activities) %in% c(interactions$source,
                                                                      interactions$target)),
                                 1]
  tf_activities <- tf_activities[order(-abs(tf_activities))[1:20]]
  # Run carnival
  # YOU MUST FIND WHERE CPLEX OR CUROBI IS INSTALLED IN YOUR OWN COMPUTER
  # GurPath <- 'C:/gurobi1100/win64/bin/gurobi_cl.exe'
  # CplexPath <- 'C:/Program Files/IBM/ILOG/CPLEX_Studio_Community2211/cplex/bin/x64_win64/cplex.exe'
  CbcPath <- 'C:/Users/nmeim/Documents/CBC/bin/cbc.exe'
  carnivalOptions <- list("solver"="cbc",
                          "betaWeight"=0.2,
                          "solverPath"=CbcPath,
                          "timelimit"=1800,
                          "poolrelGap"=1e-4,
                          "lpFilename"="",
                          "outputFolder"="",
                          "cleanTmpFiles"=TRUE,
                          "keepLPFiles"=TRUE,
                          "poolCap"=100)
  # carnivalOptions <- defaultCplexCarnivalOptions()
  # carnivalOptions$timelimit <- 1200
  # carnivalOptions$solverPath <- CplexPath
  # Output dir
  Result_dir <- paste0("../results/pc_loadings_scores_analysis/CARNIVAL/Woptim")
  dir.create(Result_dir, showWarnings = FALSE)
  carnivalOptions$outputFolder <- Result_dir
  
  inverseCarnivalResults <- runInverseCarnival( measurements = tf_activities, 
                                                priorKnowledgeNetwork = interactions, 
                                                carnivalOptions = carnivalOptions)
  
  # Save interaction networks
  nets <- inverseCarnivalResults$sifAll
  nodes <- inverseCarnivalResults$attributesAll
  for (i in 1:length(nets)){
    t <- nets[[i]]
    t <- as.data.frame(t)
    t$Node1 <- as.character(t$Node1)
    t$Node2 <- as.character(t$Node2)
    t$Sign <- as.character(t$Sign)
    write_tsv(t,paste0(Result_dir,'/','interactions_1_model',i,'.tsv'))
    t <- nodes[[i]]
    t <- as.data.frame(t)
    t$Nodes <- as.character(t$Nodes)
    t$Activity <- as.numeric(t$Activity)
    write_delim(t,paste0(Result_dir,'/','nodesActivity_1_model',i,'.txt'),delim = '\t')
  }
  t <- as.data.frame(inverseCarnivalResults$weightedSIF)
  t$Node1 <- as.character(t$Node1)
  t$Node2 <- as.character(t$Node2)
  t$Sign <- as.character(t$Sign)
  t$Weight <- as.numeric(t$Weight)
  write_delim(t,paste0(Result_dir,'/','weightedModel_1.txt'),delim = '\t')
  t <- as.data.frame(inverseCarnivalResults[["nodesAttributes"]])
  write_delim(t,paste0(Result_dir,'/','nodesAttributes_1.txt'),delim = '\t')
  message(paste0('Finished ',colnames(TF_activities)[j]))
}
close(log_con)

### Perform Fisher test between KEGG pathways and nodes in the net
Result_dir <- paste0("../results/pc_loadings_scores_analysis/CARNIVAL/PC12")
nodes <- read.delim(paste0(Result_dir,'/','nodesAttributes_1.txt'))
nodes <- nodes %>% filter(AvgAct!=0)
net <- read.delim(paste0(Result_dir,'/','weightedModel_1.txt'))
net <- net %>% filter(Node1 %in% nodes$Node) %>% filter(Node2 %in% nodes$Node) %>% unique()
egsea.data(species = "human",returnInfo = TRUE)
pathways <- kegg.pathways$human
pathways <- pathways$kg.sets
entrez_ids <- mapIds(org.Hs.eg.db, keys = nodes$Node, column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
all_nodes <- unique(interactions$source,interactions$target)
all_entrez_ids <- mapIds(org.Hs.eg.db, keys = all_nodes, column = "ENTREZID", keytype = "SYMBOL")
all_entrez_ids <- unname(all_entrez_ids)
inds <- which(!is.na(all_entrez_ids))
all_entrez_ids <- all_entrez_ids[inds]

## Perform Fisher Exact Test for each KEGG pathway
fisherExactTest <- function(Set,selected_nodes,all_possible_nodes){
  a  <- sum(selected_nodes %in% Set)
  b <- length(selected_nodes) - a
  c <- length(Set) - a
  d <- length(all_possible_nodes) - a -b-c
  contigencyMat <- data.frame('selected'=c(a,b),
                              'not_selected'= c(c,d))
  rownames(contigencyMat) <- c('in_set','not_in_set')
  return(fisher.test(contigencyMat)$p.value)
}
# universe_nodes <- unique(c(unique(unname(unlist(pathways))),all_entrez_ids))
fisher_results <- lapply(pathways,fisherExactTest,entrez_ids,all_entrez_ids)
fisher_results <- as.data.frame(unlist(fisher_results))
colnames(fisher_results) <- 'p.value'
fisher_results <- fisher_results %>% rstatix::adjust_pvalue(method = 'bonferroni')
fisher_results <- fisher_results %>% rownames_to_column('pathway')
my_data_frame <- bind_rows(lapply(names(pathways), function(name) {
  data.frame(names = name, values = list(pathways[[name]]), stringsAsFactors = FALSE)
})) %>% gather('key','node',-names) %>% filter(!is.na(key))
my_data_frame <- distinct(my_data_frame %>% select(c('pathway'='names'),node))
fisher_results <- left_join(fisher_results,my_data_frame)
fisher_results <- fisher_results %>% filter(p.value.adj<0.05)
fisher_results <- fisher_results %>% filter(!is.na(node))
fishered_nodes <- mapIds(org.Hs.eg.db, keys = fisher_results$node, column = "SYMBOL", keytype = "ENTREZID")
fisher_results$node <- fishered_nodes
fisher_results <- fisher_results %>% filter(node %in% nodes$Node)
print(unique(fisher_results$pathway))
fisher_results <- fisher_results %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))
fisher_results <- fisher_results %>% group_by(pathway) %>% mutate(num_nodes = n_distinct(node)) %>% ungroup()
fisher_results <- fisher_results %>% filter(!grepl('cancer',pathway))
fisher_results <- fisher_results %>% filter(!grepl('infection',pathway))
fisher_results <- fisher_results %>% filter(grepl('signaling',pathway) | grepl('pathway',pathway))
fisher_results <- fisher_results %>% filter(num_nodes>3) %>% select(pathway,node) %>% unique()
colnames(fisher_results) <- c('pathway','Node')
nodes <- left_join(nodes,fisher_results)
write_delim(nodes,paste0(Result_dir,'/','edited_node_attributes.txt'),delim = '\t')


### Infer TF activity with Dorothea from loadings of PC12 and PC8--------------------------------------------------
gene_loadings <- PCA_alldata$rotation
# gene_loadings <- gene_loadings[,c('PC8','PC12')]
minNrOfGenes  <-  5
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]

settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities_loadings = run_viper(gene_loadings, dorotheaData, options =  settings)
Woptim <- readRDS('../results/Wm_opt.rds')
TF_activities = run_viper(cbind(Woptim,Woptim), dorotheaData, options =  settings)
TF_activities <- as.matrix(TF_activities[rownames(TF_activities_loadings),1])
colnames(TF_activities) <- 'optimal direction'
TF_activities <- cbind(TF_activities,TF_activities_loadings)
# pheatmap(TF_activities[,paste0('PC',c(1,2,3,4,5,8,11,12,13,35))])
pheatmap(TF_activities)
genes <- rownames(gene_loadings)
# Build distribution of TF activities across all PCs and
# check what TF activities you can anyway observe in the data
iters <- 10000
# null_activity <- data.frame()
# settings = list(verbose = FALSE, minsize = minNrOfGenes)
# for (i in 1:iters){
#   shuffled_genes <-  gene_loadings[sample(nrow(gene_loadings)),]
#   rownames(shuffled_genes) <- genes
#   tmp_tfs = run_viper(shuffled_genes, dorotheaData, options =  settings)
#   tmp_tfs <- as.data.frame(tmp_tfs) %>% rownames_to_column('TF')
#   null_activity <- rbind(null_activity,tmp_tfs)
#   if (i %% 1000 ==0){
#     # saveRDS(null_activity,'../results/pc_loadings_scores_analysis/dorothea_shuffled_gene_loads.rds')
#     message(paste0('Finished iteration ',i))
#   }
# }
# # saveRDS(null_activity,'../results/pc_loadings_scores_analysis/dorothea_shuffled_gene_loads.rds')
# null_activity <- readRDS('../results/pc_loadings_scores_analysis/dorothea_shuffled_gene_loads.rds')

### Null activity for Woptim
null_activity_optimal <- data.frame()
settings = list(verbose = FALSE, minsize = minNrOfGenes)
shuffled_genes <- matrix(nrow=nrow(Woptim),ncol=iters)
for (i in 1:iters){
  shuffled_genes[,i] <- Woptim[sample(nrow(Woptim)),]
}
rownames(shuffled_genes) <- rownames(Woptim)
random_tfs = run_viper(shuffled_genes, dorotheaData, options =  settings)
random_tfs <- as.data.frame(random_tfs) %>% rownames_to_column('TF')
null_activity_optimal <- random_tfs %>% gather('trial','null_activity',-TF)
# saveRDS(null_activity_optimal,'../results/pc_loadings_scores_analysis/dorothea_shuffled_gene_loads_optimal.rds')
null_activity_optimal <- readRDS('../results/pc_loadings_scores_analysis/dorothea_shuffled_gene_loads_optimal.rds')

#plot TF activities that explain a gene being expressed toward the direction of the loading
TF_activities <- as.data.frame(TF_activities) %>% rownames_to_column('TF')
# Build distribution of TF activities across all PCs and 
# check what TF activities you can anyway observe in the data
null_activity <- TF_activities %>% gather('PC','null_act',-TF) %>% select(-PC)
TF_activities$significant <- ''

### Get p.values for Woptimal with shuffling
TF_activities_nulled <- left_join(TF_activities,null_activity_optimal,by='TF')
TF_activities_optim <- TF_activities_nulled %>% select(TF,`optimal direction`,significant,c('null_act'='null_activity')) %>%
  group_by(TF) %>% mutate(p.value = ifelse(`optimal direction`>=0,
                                           sum(null_act>=`optimal direction`)/iters,
                                           sum(null_act<=`optimal direction`)/iters)) %>%
  ungroup()
adjustement_tmp <- TF_activities_optim %>% select(TF,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(TF,p.adj)
TF_activities_optim <- left_join(TF_activities_optim,adjustement_tmp)%>% select(-null_act) %>% unique()
TF_activities_optim$significant[order(-abs(TF_activities_optim$`optimal direction`))[1:20]] <- TF_activities_optim$TF[order(-abs(TF_activities_optim$`optimal direction`))[1:20]]
TF_activities_optim <- TF_activities_optim[order(TF_activities_optim$`optimal direction`),]
TF_activities_optim$TF <- factor(TF_activities_optim$TF,levels = TF_activities_optim$TF)
# #### This only for shuffled
# TF_activities_nulled <- left_join(TF_activities,null_activity,by='TF')
# TF_activities_PC8 <- TF_activities_nulled %>% select(TF,c('PC8'='PC8.x'),significant,c('null_act'='PC8.y')) %>% 
#   group_by(TF) %>% mutate(p.value = ifelse(PC8>=0,
#                                            sum(null_act>=PC8)/iters,
#                                            sum(null_act<=PC8)/iters)) %>%
#   ungroup()
# TF_activities_PC12 <- TF_activities_nulled %>% select(TF,c('PC12'='PC12.x'),significant,c('null_act'='PC12.y')) %>% 
#   group_by(TF) %>% mutate(p.value = ifelse(PC12>=0,
#                                            sum(null_act>=PC12)/iters,
#                                            sum(null_act<=PC12)/iters)) %>%
#   ungroup()
TF_activities_PC8 <- left_join(TF_activities %>% select(TF,PC8,significant),
                               null_activity) %>%
  group_by(TF) %>% mutate(p.value = ifelse(PC8>=0,
                                           sum(null_act>=PC8)/(ncol(TF_activities)-2),
                                           sum(null_act<=PC8)/(ncol(TF_activities)-2))) %>%
  ungroup()
adjustement_tmp <- TF_activities_PC8 %>% select(TF,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(TF,p.adj)
TF_activities_PC8 <- left_join(TF_activities_PC8,adjustement_tmp)%>% select(-null_act) %>% unique()
TF_activities_PC12 <- left_join(TF_activities %>% select(TF,PC12,significant),
                                null_activity) %>%
  group_by(TF) %>% mutate(p.value = ifelse(PC12>=0,
                                           sum(null_act>=PC12)/(ncol(TF_activities)-2),
                                           sum(null_act<=PC12)/(ncol(TF_activities)-2))) %>%
  ungroup()
adjustement_tmp <- TF_activities_PC12 %>% select(TF,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(TF,p.adj)
TF_activities_PC12 <- left_join(TF_activities_PC12,adjustement_tmp) %>% select(-null_act) %>% unique()
TF_activities_PC8$significant[order(-abs(TF_activities_PC8$PC8))[1:20]] <- TF_activities_PC8$TF[order(-abs(TF_activities_PC8$PC8))[1:20]]
TF_activities_PC8 <- TF_activities_PC8[order(TF_activities_PC8$PC8),]
TF_activities_PC8$TF <- factor(TF_activities_PC8$TF,levels = TF_activities_PC8$TF)
TF_activities_PC12$significant[order(-abs(TF_activities_PC12$PC12))[1:20]] <- TF_activities_PC12$TF[order(-abs(TF_activities_PC12$PC12))[1:20]]
TF_activities_PC12 <- TF_activities_PC12[order(TF_activities_PC12$PC12),]
TF_activities_PC12$TF <- factor(TF_activities_PC12$TF,levels = TF_activities_PC12$TF)

# change significant
TF_activities_PC12 <- TF_activities_PC12 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                       ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
TF_activities_PC8 <- TF_activities_PC8 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                     ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
TF_activities_optim <- TF_activities_optim %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                     ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))

# Combine into one data frame
TF_activities_plot_frame <- rbind(TF_activities_PC8 %>% dplyr::rename(value = PC8) %>% mutate(PC='PC8'),
                                  TF_activities_PC12 %>% dplyr::rename(value = PC12) %>% mutate(PC='PC12'),
                                  TF_activities_optim %>% dplyr::rename(value = `optimal direction`) %>% mutate(PC='optimal direction'))
TF_activities_plot_frame <- TF_activities_plot_frame %>% mutate(PC=ifelse(PC=='optimal direction','used optimal gene loadings',paste0('used ',PC,' gene loadings')))
### Make barplot to look at top TFs
# p <- (ggplot(TF_activities_PC8,aes(x=TF,y=PC8,color = statistical)) + geom_point() +
#   # scale_color_gradient(high = 'red',low='white')+
#     scale_color_manual(values = c('#fa8e8e','black'))+
#   geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
#   xlab('TFs') + ylab('activity from PC8 loadings')+
#   scale_x_discrete(expand = c(0.1, 0.1))+
#   theme_pubr(base_family = 'Arial',base_size = 20)+
#   theme(text = element_text(family = 'Arial',size=20),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = 'none')) +
#   (ggplot(TF_activities_PC12,aes(x=TF,y=PC12,color = statistical)) + geom_point() +
#      scale_color_manual(values =c('red','#fa8e8e','black'))+
#      # scale_color_gradient(high = 'red',low='white')+
#      geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
#      xlab('TFs') + ylab('activity from PC12 loadings')+
#      scale_x_discrete(expand = c(0.1, 0.1))+
#      theme_pubr(base_family = 'Arial',base_size = 20)+
#      theme(text = element_text(family = 'Arial',size=20),
#            axis.ticks.x = element_blank(),
#            axis.text.x = element_blank(),
#            legend.position = 'none'))
p <- ggplot(TF_activities_plot_frame %>% dplyr::rename(`statistical threshold`=statistical) %>%
              mutate(p.adj = ifelse(p.adj<1e-16,1e-16,p.adj)),
            aes(x=value,y=-log10(p.adj),color = `statistical threshold`)) + geom_point() +
  scale_color_manual(values = c('red','#fa8e8e','black'))+
  ggbreak::scale_y_break(breaks = c(2,15),space = 0.5,ticklabels = c(15,16))+
  geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
  xlab('activity') + ylab(expression(-log[10]('adjusted p-value'))) +
  # ylab('adjusted p-value') +
  # scale_x_discrete(expand = c(0.1, 0.1))+
  geom_hline(yintercept=-log10(0.1), linetype="dashed",color = "#525252", size=1) +
  # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'top')+
  facet_wrap(~PC,scales = 'free_x')+
  guides(color = guide_legend(
    override.aes = list(
      linetype = NA,
      size = 3
    )
  ))
print(p)
ggsave('../results/pc_loadings_scores_analysis/TFs_from_loadings_optimal_volcano.png',
       plot = p,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)
significant_tfs <- TF_activities_plot_frame %>% filter(p.adj<0.1)
significant_tfs <- unique(as.character(significant_tfs$TF))
TF_activities_loadings = run_viper(gene_loadings, dorotheaData, options =  settings)
TF_activities = run_viper(cbind(Woptim,Woptim), dorotheaData, options =  settings)
TF_activities <- as.matrix(TF_activities[rownames(TF_activities_loadings),1])
colnames(TF_activities) <- 'optimal direction'
TF_activities <- cbind(TF_activities,TF_activities_loadings)
png('../results/pc_loadings_scores_analysis/significant_TFs_from_loadings_optinal_heatmap.png',width = 9,height = 9,units = 'in',res = 600)
pheatmap(TF_activities[which(rownames(TF_activities) %in% significant_tfs),
                       c('optimal direction',paste0('PC',c(1,2,8,12)))])
dev.off()

### Infer protein activity-----------------------------------------
interactions <- import_omnipath_interactions()
interactions <- interactions %>% mutate(kegg=grepl(pattern="KEGG",x=sources)) %>%
  mutate(signoir=grepl(pattern="SIGNOR",x=sources)) %>% filter(kegg==T | signoir==T)
interactions <- interactions  %>% filter(n_resources>1)  %>% #filter(n_references>=2)%>% 
  dplyr::select(c('source'='source_genesymbol'),
                c('target'='target_genesymbol'),
                is_inhibition,is_stimulation) %>% unique()
interactions <- interactions %>% mutate(interaction=ifelse(is_stimulation==1,
                                                           ifelse(is_inhibition!=1,1,0),
                                                           ifelse(is_inhibition!=0,-1,0))) %>%
  dplyr::select(source,interaction,target) %>% unique()
# interactions <- interactions %>% filter(interaction!=0)
minNrOfGenes  <-  5
settings = list(verbose = TRUE, minsize = minNrOfGenes)
# protein_activity = run_viper(PCA_alldata$rotation, 
#                              rbind(interactions %>% mutate(confidence = 'A') %>% 
#                                      select(c('tf'='source'),confidence,target,c('mor'='interaction')) %>%
#                                      filter(!(paste0(tf,' , ',target) %in% paste0(dorotheaData$tf,' , ',dorotheaData$target))),
#                                    dorotheaData) %>% unique(), 
#                              options =  settings)
protein_activity = run_viper(PCA_alldata$rotation, 
                             interactions %>% mutate(confidence = 'A') %>% 
                               select(c('tf'='source'),confidence,target,c('mor'='interaction')), 
                             options =  settings)
pheatmap(protein_activity[,paste0('PC',c(8,11,12,13,35))])
protein_activity <- as.data.frame(protein_activity) %>% rownames_to_column('prot')
# Build distribution of prot activities across all PCs and 
# check what prot activities you can anyway observe in the data
null_activity <- protein_activity %>% gather('PC','null_act',-prot) %>% select(-PC)
protein_activity$significant <- ''
protein_activity_PC8 <- left_join(protein_activity %>% select(prot,PC8,significant),
                                  null_activity) %>%
  group_by(prot) %>% mutate(p.value = ifelse(PC8>=0,
                                             sum(null_act>=PC8)/(ncol(protein_activity)-2),
                                             sum(null_act<=PC8)/(ncol(protein_activity)-2))) %>%
  ungroup()
adjustement_tmp <- protein_activity_PC8 %>% select(prot,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(prot,p.adj)
protein_activity_PC8 <- left_join(protein_activity_PC8,adjustement_tmp)%>% select(-null_act) %>% unique()
protein_activity_PC12 <- left_join(protein_activity %>% select(prot,PC12,significant),
                                   null_activity) %>%
  group_by(prot) %>% mutate(p.value = ifelse(PC12>=0,
                                             sum(null_act>=PC12)/(ncol(protein_activity)-2),
                                             sum(null_act<=PC12)/(ncol(protein_activity)-2))) %>%
  ungroup()
adjustement_tmp <- protein_activity_PC12 %>% select(prot,p.value) %>% unique()
adjustement_tmp$p.adj <- p.adjust(adjustement_tmp$p.value,method = 'BH')
adjustement_tmp <- adjustement_tmp %>% select(prot,p.adj)
protein_activity_PC12 <- left_join(protein_activity_PC12,adjustement_tmp) %>% select(-null_act) %>% unique()
protein_activity_PC8$significant[order(-abs(protein_activity_PC8$PC8))[1:20]] <- protein_activity_PC8$prot[order(-abs(protein_activity_PC8$PC8))[1:20]]
protein_activity_PC8 <- protein_activity_PC8[order(protein_activity_PC8$PC8),]
protein_activity_PC8$prot <- factor(protein_activity_PC8$prot,levels = protein_activity_PC8$prot)
protein_activity_PC12$significant[order(-abs(protein_activity_PC12$PC12))[1:20]] <- protein_activity_PC12$prot[order(-abs(protein_activity_PC12$PC12))[1:20]]
protein_activity_PC12 <- protein_activity_PC12[order(protein_activity_PC12$PC12),]
protein_activity_PC12$prot <- factor(protein_activity_PC12$prot,levels = protein_activity_PC12$prot)
# change significant
protein_activity_PC12 <- protein_activity_PC12 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                             ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))
protein_activity_PC8 <- protein_activity_PC8 %>% mutate(statistical=ifelse(p.adj<=0.05,'p.adj<=0.05',
                                                                           ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1')))

# protein_activity_PC12 <- protein_activity_PC12 %>% mutate(statistical=ifelse(p.adj<=0.01,'p.adj<=0.01',
#                                                                              ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                                     ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1'))))
# protein_activity_PC8 <- protein_activity_PC8 %>% mutate(statistical=ifelse(p.adj<=0.01,'p.adj<=0.01',
#                                                                            ifelse(p.adj<=0.05,'p.adj<=0.05',
#                                                                                   ifelse(p.adj<=0.1,'p.adj<=0.1','p.adj>0.1'))))

# Combine into one data frame
protein_activity_plot_frame <- rbind(protein_activity_PC8 %>% dplyr::rename(value = PC8) %>% mutate(PC='PC8'),
                                     protein_activity_PC12 %>% dplyr::rename(value = PC12) %>% mutate(PC='PC12'))
protein_activity_plot_frame <- protein_activity_plot_frame %>% mutate(PC=paste0('used ',PC,' gene loadings'))
### Make barplot to look at top prots
# p <- (ggplot(protein_activity_PC8,aes(x=prot,y=PC8,color = statistical)) + geom_point() +
#   # scale_color_gradient(high = 'red',low='white')+
#     scale_color_manual(values = c('#fa8e8e','black'))+
#   geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
#   xlab('prots') + ylab('activity from PC8 loadings')+
#   scale_x_discrete(expand = c(0.1, 0.1))+
#   theme_pubr(base_family = 'Arial',base_size = 20)+
#   theme(text = element_text(family = 'Arial',size=20),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = 'none')) +
#   (ggplot(protein_activity_PC12,aes(x=prot,y=PC12,color = statistical)) + geom_point() +
#      scale_color_manual(values =c('red','#fa8e8e','black'))+
#      # scale_color_gradient(high = 'red',low='white')+
#      geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
#      xlab('prots') + ylab('activity from PC12 loadings')+
#      scale_x_discrete(expand = c(0.1, 0.1))+
#      theme_pubr(base_family = 'Arial',base_size = 20)+
#      theme(text = element_text(family = 'Arial',size=20),
#            axis.ticks.x = element_blank(),
#            axis.text.x = element_blank(),
#            legend.position = 'none'))
p <- ggplot(protein_activity_plot_frame %>% dplyr::rename(`statistical threshold`=statistical),
            aes(x=value,y=-log10(p.adj),color = `statistical threshold`)) + geom_point() +
  scale_color_manual(values = c('black'))+
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40,point.padding = 0.5)+
  xlab('activity') + ylab(expression(-log[10]('adjusted p-value'))) +
  geom_hline(yintercept=-log10(0.1), linetype="dashed",color = "#525252", size=1) +
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        legend.position = 'top')+
  facet_wrap(~PC,scales = 'free_x')+
  guides(color = guide_legend(
    override.aes = list(
      linetype = NA,
      size = 3
    )
  ))
print(p)
ggsave('../results/pc_loadings_scores_analysis/protein_activity_from_loadings_volcano.png',
       plot = p,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

### Perform KEGG pathway analysis (GSEA)------------
# rownames(loadings) <- NULL
PCA_alldata <- prcomp(X_B, scale. = F, center = T)
loadings <- PCA_alldata$rotation
Woptim <- readRDS('../results/Wm_opt.rds')
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,PC12,PC8) %>% column_to_rownames('gene')
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(loadings), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- loadings[inds,]
# meas <- as.matrix(loadings[inds,1])
# rownames(meas) <- NULL
rownames(meas) <- entrez_ids
# names(meas) <- entrez_ids

## run gsea
egsea.data(species = "human",returnInfo = TRUE)
kegg_list <-  kegg.pathways$human$kg.sets
# keggs <- fgseaSimple(pathways = kegg_list,
#                      stats=meas,
#                      minSize=10,
#                      maxSize=500,
#                      nperm = 10000)
keggs <- fastenrichment(colnames(meas),
                        rownames(meas),
                        meas,
                        enrichment_space = 'kegg',
                        n_permutations = 10000)
kegg_nes <- as.data.frame(keggs$NES$`NES KEGG`) %>% rownames_to_column('pathway') %>% gather('PC','NES',-pathway)
kegg_pval <- as.data.frame(keggs$Pval$`Pval KEGG`) %>% rownames_to_column('pathway') %>% gather('PC','padj',-pathway)
df_keggs <- left_join(kegg_nes,kegg_pval)
# df_keggs <- keggs %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))
df_keggs <- df_keggs %>% mutate(pathway=strsplit(pathway,"_"))
df_keggs <- df_keggs %>% unnest(pathway) %>% filter(!(pathway %in% c("KEGG","FL1000")))
df_keggs <- df_keggs %>% mutate(pathway=as.character(pathway))
df_keggs <- df_keggs %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))

df_keggs <- df_keggs %>% mutate(sig=ifelse(abs(NES)>1.5 & padj<0.05,'yes','no'))
df_keggs <- df_keggs%>% mutate(label = ifelse(sig=='yes',pathway,NA))
p1 <- ggplot(df_keggs %>% filter(PC=='PC12') %>% arrange(NES),aes(x=NES,y=-log10(padj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Pathways enriched in PC12 loadings')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)
# df_keggs <- df_keggs[order(df_keggs$NES),]
# df_keggs$pathway <- factor(df_keggs$pathway,levels = df_keggs$pathway)
# p2 <- ggplot(df_keggs %>% arrange(NES) %>% filter(padj<0.05),aes(x=NES,y=pathway,fill=padj))+ 
#   geom_bar(stat = 'identity',color='black',size=1.5) +
#   scale_fill_gradient(low = "red",high = "white",limits = c(min(df_keggs$padj),0.05)) +
#   xlab('Normalized Enrichment Score') +
#   ggtitle('Pathways enriched in PC12 loadings')+
#   theme_pubr(base_family = 'Arial',base_size = 20)+
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.key.size = unit(1.5, "lines"),
#         legend.position = 'right',
#         legend.justification = "center")
# print(p2)
# ggsave('../results/pc_loadings_scores_analysis/kegg_on_pc12.png',
#        plot=p2,
#        width=16,
#        height=9,
#        units = 'in',
#        dpi = 600)

### Repeat for Wopt
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Woptim), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Woptim[inds,])
rownames(meas) <- entrez_ids
keggs <- fastenrichment(colnames(meas),
                        rownames(meas),
                        meas,
                        enrichment_space = 'kegg',
                        n_permutations = 10000,
                        order_columns=F)
kegg_nes <- as.data.frame(keggs$NES$`NES KEGG`) %>% rownames_to_column('pathway')
colnames(kegg_nes)[2] <- 'NES'
kegg_pval <- as.data.frame(keggs$Pval$`Pval KEGG`) %>% rownames_to_column('pathway')
colnames(kegg_pval)[2] <- 'padj'
df_keggs <- left_join(kegg_nes,kegg_pval)
df_keggs <- df_keggs %>% mutate(pathway=strsplit(pathway,"_"))
df_keggs <- df_keggs %>% unnest(pathway) %>% filter(!(pathway %in% c("KEGG","FL1000")))
df_keggs <- df_keggs %>% mutate(pathway=as.character(pathway))
df_keggs <- df_keggs %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))
df_keggs <- df_keggs %>% mutate(sig=ifelse(abs(NES)>1.5 & padj<0.05,'yes','no'))
df_keggs <- df_keggs%>% mutate(label = ifelse(sig=='yes',pathway,NA))
p1 <- ggplot(df_keggs %>% arrange(NES),aes(x=NES,y=-log10(padj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Pathways enriched in optimal loadings')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

### Perform GO Terms analysis (GSEA)------------
loadings <- PCA_alldata$rotation
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,PC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$PC12))[1:10]] <- loadings$gene[order(-abs(loadings$PC12))[1:10]]
loadings <- loadings[order(loadings$PC12),]
loadings$gene <- factor(loadings$gene,levels = loadings$gene)
meas <- loadings[,2]
# meas <- as.matrix(loadings[inds,1])
# rownames(meas) <- NULL
# rownames(meas) <- entrez_ids
names(meas) <- loadings$gene

## run gsea
genes <- factor(x = rep(1,length(meas)),levels = c(0,1))
names(genes) <- names(meas)

GOobject <- new("topGOdata",ontology = "BP", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                ID = "symbol", nodeSize = 10)
term.genes <- genesInTerm(GOobject, GOobject@graph@nodes)
gos <- fgseaSimple(pathways = term.genes,
                   stats=meas,
                   minSize=10,
                   maxSize=500,
                   nperm = 10000)
colnames(gos)[1] <- 'GO'
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
df_gos <- left_join(gos,go_annotations %>% filter(ontology=='BP') %>% dplyr::select(GO,`GO Terms`))
df_gos <- df_gos %>% mutate(sig=ifelse(abs(NES)>2 & padj<0.05,'yes','no'))
df_gos <- df_gos%>% mutate(label = ifelse(sig=='yes',`GO Terms`,NA))
p1 <- ggplot(df_gos %>% arrange(NES),aes(x=NES,y=-log10(padj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('GO Terms enriched in PC12 loadings')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

df_gos <- df_gos[order(df_gos$NES),]
df_gos$`GO Terms` <- factor(df_gos$`GO Terms`,levels = df_gos$`GO Terms`)
p2 <- ggplot(df_gos %>% arrange(NES) %>% filter(padj<0.05 & abs(NES)>2),aes(x=NES,y=`GO Terms`,fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(low = "red",high = "white",limits = c(min(df_gos$padj),0.05)) +
  xlab('Normalized Enrichment Score') +
  ggtitle('GO Terms enriched in PC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/pc_loadings_scores_analysis/go_on_pc12.png',
       plot=p2,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)

### Do GO Enrichment for Woptim
gos <- fastenrichment(colnames(Woptim),
                      rownames(Woptim),
                      Woptim,
                      enrichment_space = 'go_bp',
                      gene_id_type="symbol",
                      n_permutations = 10000,
                      order_columns=F)
go_nes <- as.data.frame(gos$NES$`NES GO BP`) %>% rownames_to_column('GO')
colnames(go_nes)[2] <- 'NES'
go_pval <- as.data.frame(gos$Pval$`Pval GO BP`) %>% rownames_to_column('GO')
colnames(go_pval)[2] <- 'padj'
df_gos <- left_join(go_nes,go_pval)
df_gos <- df_gos %>% mutate(GO=strsplit(GO,"_"))
df_gos <- df_gos %>% unnest(GO) %>% filter(!(GO %in% c("GO","FL1000","BP")))
df_gos <- df_gos %>% mutate(GO=as.character(GO))
df_gos <- left_join(df_gos,go_annotations)

### GSEA with Hallmarks MSIGDB------------------------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Woptim), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Woptim[inds,])
rownames(meas) <- entrez_ids
msig <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes <- as.data.frame(msig$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark')
colnames(msig_nes)[2] <- 'NES'
msig_pval <- as.data.frame(msig$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')
colnames(msig_pval)[2] <- 'padj'
df_msig <- left_join(msig_nes,msig_pval)
df_msig <- df_msig %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))

df_msig <- df_msig %>% mutate(Hallmark=strsplit(Hallmark,"_"))
df_msig <- df_msig %>% unnest(Hallmark) %>% filter(!(Hallmark %in% c("HALLMARK","FL1000","MSIG","H")))
df_msig <- df_msig %>% mutate(Hallmark=as.character(Hallmark))
df_msig <- df_msig %>% mutate(Hallmark=gsub('_'," ",Hallmark))
df_msig$Hallmark <- factor(df_msig$Hallmark,levels = df_msig$Hallmark[order(df_msig$NES)])
p1 <- ggplot(df_msig %>% arrange(NES),aes(x=NES,y=Hallmark,fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
  xlab('Normalized Enrichment Score') +
  ggtitle('Hallmarks enriched in optimal loadings')+
  theme_pubr(base_family = 'Arial',base_size = 15)+
  theme(text = element_text(family = 'Arial',size=15),
        axis.text.y = element_text(size=13),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p1)
ggsave('../results/pc_loadings_scores_analysis/hallmark_on_optimal_loadings.png',
       plot=p1,
       width=12,
       height=9,
       units = 'in',
       dpi = 600)
### NicheNet analysis-----------------------------------------
Woptim <- readRDS('../results/Wm_opt.rds')
colnames(Woptim) <- 'w_optimal'
gene_loadings <- PCA_alldata$rotation[,paste0('PC',c(1,2,8,12))]

weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
ligands_all = lr_network %>% pull(from) %>% unique()
print(paste0('Is IFNa in available ligands? : ',any(grepl('IFNA',ligands_all))))
print(paste0('Is IFNg in available ligands? : ',any(grepl('IFNG',ligands_all))))
selected_ligands <- ligands_all[grepl('IFNA',ligands_all) | grepl('IFNG',ligands_all)]
targets_all = names(Woptim[,1][order(-abs(Woptim[,1]))[1:40]])
targets_all <- targets_all[which(targets_all %in% gr_network$to)]

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

DiagrammeR::render_graph(graph_min_max, layout = "tree")