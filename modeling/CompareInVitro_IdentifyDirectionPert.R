# load in packages
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(pheatmap)
library(patchwork)

# Load metadata
mps_metadata <- read.table('../data/GSE168285/GSE168285_sample_meta_data.txt',header = TRUE, sep = "\t") %>% select(-X)
wang_metadata <- data.table::fread('../data/METADATA_Wangetal2020.txt',header = TRUE, sep = "\t")
feaver_metadata <- read.table('../data/GSE89063_metadata.txt',header = TRUE, sep = ",")# %>% select(-X)

# Load gene expression
data_wang <- read.delim2('../data/GSE166256_kallisto_counts.tsv') %>% column_to_rownames('Gene')
data_wang <- as.matrix(data_wang)
data_wang <- apply(data_wang,MARGIN = c(1,2),as.numeric)
data <- readRDS("../data/preprocessed_NAFLD.rds")
data_mps <- data$data_B
data_feaver <- read.delim2('../data/GSE89063_expression_matrix.tsv') # %>% column_to_rownames('X')
data_feaver <- aggregate(data_feaver[,2:ncol(data_feaver)],by = list(data_feaver$X),median)
gc()
data_feaver <- data_feaver %>% column_to_rownames('Group.1')

### Pre-process MPS
data_mps  <- apply(data_mps, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_mps <- log10(1 + data_mps)
data_mps <- t(scale(t(data_mps),scale = F))
### Pre-process Wang
data_wang  <- apply(data_wang, MARGIN = 2, FUN = function(x) x/sum(x)*1e6) 
data_wang <- log10(1 + data_wang)
data_wang <- t(scale(t(data_wang),scale = F))
## Pre-process feaver
data_feaver  <- apply(data_feaver, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
data_feaver <- log10(1 + data_feaver)
data_feaver <- t(scale(t(data_feaver),scale = F))

### read optimal direction
w_opt <- readRDS('../results/Wm_opt.rds')

### Check how other in-vitro jiggle in the W_opt-------------------------------------
data_mps_proj <- t(data_mps[reduce(list(rownames(data_mps), rownames(w_opt)),intersect),
                            ]) %*%  w_opt[reduce(list(rownames(data_mps), rownames(w_opt)),intersect),]
data_mps_proj <- as.data.frame(data_mps_proj) %>% mutate(dataset = 'MPS')
data_wang_proj <- t(data_wang[reduce(list(rownames(data_wang), rownames(w_opt)),intersect),
]) %*%  w_opt[reduce(list(rownames(data_wang), rownames(w_opt)),intersect),]
data_wang_proj <- as.data.frame(data_wang_proj) %>% mutate(dataset = 'Wang')
data_feaver_proj <- t(data_feaver[reduce(list(rownames(data_feaver), rownames(w_opt)),intersect),
]) %*%  w_opt[reduce(list(rownames(data_feaver), rownames(w_opt)),intersect),]
data_feaver_proj <- as.data.frame(data_feaver_proj) %>% mutate(dataset = 'Feaver')
data_combo <- rbind(data_mps_proj,data_wang_proj,data_feaver_proj)
data_combo <- data_combo %>% mutate(type='in-vitro')

# Project also human data
geos <- c("GSE126848","GSE130970","GSE134422","GSE135251","GSE162694")
expression_matrix <- readRDS('../data/ARCHS4/FattyLiver_expression_matrix.rds')
meta_data <- read.delim('../data/ARCHS4/FattyLiver_meta_data.tsv',row.names = 1)
data_combo_clinical <- data.frame()
for (i in 1:length(geos)){
  geo <- geos[i]
  tmp <- meta_data %>% filter(series_id==geo)
  X <- expression_matrix[,tmp$sample]
  genes <- rownames(X)
  X <- apply(X, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
  X <- log10(X+1)
  X <- t(scale(t(X),scale = F))
  rownames(X) <- genes
  X <- t(X[reduce(list(rownames(X), rownames(w_opt)),intersect),
  ]) %*%  w_opt[reduce(list(rownames(X), rownames(w_opt)),intersect),]
  X <- as.data.frame(X) %>% mutate(dataset = geo) %>% mutate(type='clinical')
  data_combo_clinical <- rbind(data_combo_clinical,X)
  print(paste0('Finished dataset ',i,' out of ',length(geos)))
}

data_combo_all <- rbind(data_combo_clinical,data_combo)
comparisons <- c('Wang','Feaver','GSE126848','GSE130970','GSE134422','GSE135251','GSE162694')
stats_results <- data.frame()
for (comp in comparisons){
  tmp <- data_combo_all %>% filter(dataset %in% c("MPS",comp))
  res <- var.test(V1~dataset,tmp,alternative='two.sided')
  stats_results <- rbind(stats_results,
                         data.frame(group1 = 'MPS',
                                    group2 = comp,
                                    p_value = res$p.value,
                                    variance = var(tmp$V1[which(tmp$dataset==comp)])))
}
stats_results$p.adj = p.adjust(stats_results$p_value,method = 'bonferroni')
stats_results <- rbind(stats_results,
                       data.frame(group1='MPS',
                                  group2 = 'MPS',
                                  p_value = NA,
                                  variance=var(data_combo_all$V1[which(data_combo_all$dataset=='MPS')]),
                                  p.adj = NA))
stats_results <- stats_results %>% 
  mutate(p.sig = ifelse(p.adj<=1e-04,"****",
                        ifelse(p.adj<=1e-03,"***",
                               ifelse(p.adj<=1e-02,"**",
                               ifelse(p.adj<=0.05,"*","ns")))))
stats_results <- stats_results %>% mutate(label = ifelse(is.na(p_value),
                                                         paste0('\n',"variance=",round(variance,3)),
                                                         paste0(p.sig,'\n',"variance=",round(variance,3))))
data_combo_all$dataset <- factor(data_combo_all$dataset,
                                    levels = c("MPS","Feaver",'Wang',geos))
data_combo_all$dataset <- relevel(data_combo_all$dataset, ref = "MPS")
data_combo_all <- data_combo_all[order(data_combo_all$type),]
data_combo_all$type <- factor(data_combo_all$type,levels=c('clinical',"in-vitro"))
ggboxplot(data_combo_all,x='dataset',y='V1',add = 'jitter',
          color = 'type') +
  ylab("score")+
  geom_hline(yintercept = 0,linetype='dashed',lwd=1,color='black')+
  stat_pvalue_manual(stats_results,
                     label = "label",
                     remove.bracket = T,
                     y.position = 2,
                     size=6.5)+
  theme(text = element_text(family = 'Arial',size=20),
        legend.title = element_blank())
ggsave('../results/pc_loadings_scores_analysis/variation_in_optimal_direction.png',
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)

### Infer pathway activity--------------------
keep_gene <- rowSums(data_mps > 0) >= 0.1*ncol(data_mps)
data_mps_aggr <- data_mps[keep_gene,]
data_mps_aggr <- t(data_mps_aggr)
grouped_data <- mps_metadata %>% select(sampleName,treatment) %>% filter(sampleName %in% rownames(data_mps_aggr)) %>% unique()
print(all(grouped_data$sampleName == rownames(data_mps_aggr)))
data_mps_aggr <- aggregate(data_mps_aggr,by=list(grouped_data$treatment),FUN=mean)
data_mps_aggr <- t(data_mps_aggr %>% column_to_rownames('Group.1'))
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
pathways_mps <- decoupleR::run_mlm(mat=data_mps_aggr, net=net_prog, .source='source', .target='target',
                                           .mor='weight', minsize = 1)
jakstat_mps <- pathways_mps %>% filter(source=='JAK-STAT')
jakstat_mps$condition <- factor(jakstat_mps$condition,levels = jakstat_mps$condition[order(jakstat_mps$score)]) 
jakstat_mps <- jakstat_mps %>% filter(p_value<0.05)
# jakstat_mps <- jakstat_mps %>% mutate(Fructose = grepl('Fructose',condition))
ggplot(jakstat_mps,aes(x=score,y=condition,fill=score,color=score)) + 
  geom_bar(stat = 'identity') +
  scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
  scale_color_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(text = element_text(family = 'Arial',size = 20),
        axis.text.y = element_blank())
# annotation <- mps_metadata %>% select(Combination_name,Background,NPC)
# annotation$Combination_name <- as.factor(annotation$Combination_name)
# annotation$NPC <- as.factor(annotation$NPC)
# annotation$Background <- as.factor(annotation$Background)
# mat <- pathways_mps %>% select(condition,source,score) %>% spread(condition,score) %>% column_to_rownames('source')
# ann_colors = list(
#   Background = c(fat="orange", lean="#50C878"),
#   NPC=c(high="red", low="blue"))
# rownames(annotation) <- mps_metadata$sampleName
# breaksList = seq(-20, 20, by = 10)
# pheatmap(mat[,mps_metadata$sampleName],
#          annotation = annotation,annotation_colors = ann_colors,
#          show_colnames = F)

pathways_wang <- decoupleR::run_mlm(mat=data_wang, net=net_prog, .source='source', .target='target',
                                   .mor='weight', minsize = 1)
pathways_feaver <- decoupleR::run_mlm(mat=data_feaver, net=net_prog, .source='source', .target='target',
                                   .mor='weight', minsize = 1)



