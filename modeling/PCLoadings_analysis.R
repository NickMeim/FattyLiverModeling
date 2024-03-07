# load in packages
library(tidyverse)
library(DESeq2)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(patchwork)
library(pls)
library(pheatmap)
library(patchwork)
library(limma)
library(edgeR)
library(dorothea)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(hgu133a.db)
library(rstatix)
library(fgsea)
library(topGO)
library(GO.db)
library(progeny)
library(OmnipathR)
library(EGSEAdata)

#### Load pre-processed data----------------------
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

### Perform PCA and get loadings----------------
PCA_alldata <- prcomp(X_B, scale. = F, center = T)
loadings <- PCA_alldata$rotation
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,PC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$PC12))[1:10]] <- loadings$gene[order(-abs(loadings$PC12))[1:10]]
loadings <- loadings[order(loadings$PC12),]
loadings$gene <- factor(loadings$gene,levels = loadings$gene)
### Make barplot to look at top genes----------------
ggplot(loadings,aes(x=gene,y=PC12,color = significant)) + geom_point() +
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
  xlab('genes') + ylab('PC12 loadings')+
  scale_x_discrete(expand = c(0.1, 0.1))+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(text = element_text(family = 'Arial',size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none')
ggsave('../results/pc_loadings_scores_analysis/gene_pc12_loadings.png',
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
gene_loadings <- PCA_alldata$rotation[,paste0('PC',c(1,2,3,8,12))]
pathways_progeny <- progeny(gene_loadings,
                            perm = 10000,
                            z_scores = T,
                            top = 500)
pheatmap(pathways_progeny)
png('../results/pc_loadings_scores_analysis/significant_progenies_from_loadings_heatmap.png',width = 6,height = 6,units = 'in',res = 600)
pheatmap(pathways_progeny)
dev.off()

#### Run CARNIVAL--------------------------
library(CARNIVAL)
gene_loadings <- PCA_alldata$rotation[,paste0('PC',c(1,2,3,8,12))]
minNrOfGenes  <-  5
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities = run_viper(gene_loadings, dorotheaData, options =  settings)
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
                                 colnames(TF_activities)[j]]
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
  Result_dir <- paste0("../results/pc_loadings_scores_analysis/",colnames(TF_activities)[j])
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
fisher_results <- fisher_results %>% filter(num_nodes>=10) %>% select(pathway,node) %>% unique()
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
TF_activities = run_viper(gene_loadings, dorotheaData, options =  settings)
pheatmap(TF_activities[,paste0('PC',c(1,2,3,4,5,8,11,12,13,35))])
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

#plot TF activities that explain a gene being expressed toward the direction of the loading
TF_activities <- as.data.frame(TF_activities) %>% rownames_to_column('TF')
# Build distribution of TF activities across all PCs and 
# check what TF activities you can anyway observe in the data
null_activity <- TF_activities %>% gather('PC','null_act',-TF) %>% select(-PC)
TF_activities$significant <- ''
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

# Combine into one data frame
TF_activities_plot_frame <- rbind(TF_activities_PC8 %>% dplyr::rename(value = PC8) %>% mutate(PC='PC8'),
                                  TF_activities_PC12 %>% dplyr::rename(value = PC12) %>% mutate(PC='PC12'))
TF_activities_plot_frame <- TF_activities_plot_frame %>% mutate(PC=paste0('used ',PC,' gene loadings'))
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
p <- ggplot(TF_activities_plot_frame %>% dplyr::rename(`statistical threshold`=statistical),
            aes(x=value,y=-log10(p.adj),color = `statistical threshold`)) + geom_point() +
  scale_color_manual(values = c('red','#fa8e8e','black'))+
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40,point.padding = 0.5)+
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
ggsave('../results/pc_loadings_scores_analysis/TFs_from_loadings_volcano.png',
       plot = p,
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)
significant_tfs <- TF_activities_plot_frame %>% filter(p.adj<0.1)
significant_tfs <- unique(as.character(significant_tfs$TF))
TF_activities = run_viper(gene_loadings, dorotheaData, options =  settings)

png('../results/pc_loadings_scores_analysis/significant_TFs_from_loadings_heatmap.png',width = 6,height = 6,units = 'in',res = 600)
pheatmap(TF_activities[which(rownames(TF_activities) %in% significant_tfs),paste0('PC',c(1,2,3,8,12))])
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
loadings <- loadings %>% select(gene,PC12) %>% column_to_rownames('gene')
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(loadings), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- loadings[inds,1]
# meas <- as.matrix(loadings[inds,1])
# rownames(meas) <- NULL
# rownames(meas) <- entrez_ids
names(meas) <- entrez_ids

## run gsea
egsea.data(species = "human",returnInfo = TRUE)
kegg_list <-  kegg.pathways$human$kg.sets
keggs <- fgseaSimple(pathways = kegg_list,
                     stats=meas,
                     minSize=10,
                     maxSize=500,
                     nperm = 10000)
df_keggs <- keggs %>% mutate(pathway=substr(pathway, 9, nchar(pathway)))
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
  ggtitle('Pathways enriched in PC12 loadings')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

df_keggs <- df_keggs[order(df_keggs$NES),]
df_keggs$pathway <- factor(df_keggs$pathway,levels = df_keggs$pathway)
p2 <- ggplot(df_keggs %>% arrange(NES) %>% filter(padj<0.05),aes(x=NES,y=pathway,fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(low = "red",high = "white",limits = c(min(df_keggs$padj),0.05)) +
  xlab('Normalized Enrichment Score') +
  ggtitle('Pathways enriched in PC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/pc_loadings_scores_analysis/kegg_on_pc12.png',
       plot=p2,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)

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


### Perform TF analysis (GSEA)------------
dorothera_regulon <- read.delim('../data/dorothea.tsv')
dorothera_regulon <- dorothera_regulon %>% filter(confidence %in% c('A','B'))
dorothera_regulon <- dorothera_regulon %>% dplyr::select(tf,target)
dorothera_regulon <- distinct(dorothera_regulon)
tf_list <- NULL
for (ind in unique(dorothera_regulon$tf)){
  tmp <- dorothera_regulon %>% filter(tf==ind)
  tf_list[[ind]] <- tmp$target
}
tfs <- fgseaSimple(pathways = tf_list,
                   stats=meas,
                   minSize=1,
                   maxSize=500,
                   nperm = 10000)
colnames(tfs)[1] <- 'TF'
df_tfs <- tfs %>% mutate(sig=ifelse(abs(NES)>2 & padj<0.05,'yes','no'))
df_tfs <- df_tfs %>% mutate(label = ifelse(sig=='yes',TF,NA))
p1 <- ggplot(df_tfs %>% arrange(NES),aes(x=NES,y=-log10(padj),color=sig))+ 
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

df_tfs <- df_tfs[order(df_tfs$NES),]
df_tfs$TF <- factor(df_tfs$TF,levels = df_tfs$TF)
p2 <- ggplot(df_tfs %>% arrange(NES) %>% filter(padj<0.05),aes(x=NES,y=TF,fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(low = "red",high = "white",limits = c(min(df_tfs$padj),0.05)) +
  xlab('Normalized Enrichment Score') +
  ggtitle('TFs enriched in PC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/pc_loadings_scores_analysis/tfs_on_pc12.png',
       plot=p2,
       width=16,
       height=10,
       units = 'in',
       dpi = 600)
