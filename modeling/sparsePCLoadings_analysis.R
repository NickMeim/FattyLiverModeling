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
sPC_B <- readRDS('../results/TransCompR_sparsePCA/optimal_sPC_B.rds')
sPC_A <- readRDS('../results/TransCompR_sparsePCA/optimal_sPC_A.rds')
loadings <- sPC_B$loadings
colnames(loadings) <- paste0('sPC',seq(1:ncol(loadings)))
rownames(loadings) <- colnames(X_B)
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,sPC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$sPC12))[1:10]] <- loadings$gene[order(-abs(loadings$sPC12))[1:10]]
loadings <- loadings[order(loadings$sPC12),]
loadings$gene <- factor(loadings$gene,levels = loadings$gene)
### Make barplot to look at top genes-----------
ggplot(loadings,aes(x=gene,y=sPC12,color = significant)) + geom_point() +
  geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
  xlab('genes') + ylab('sPC12 loadings')+
  scale_x_discrete(expand = c(0.1, 0.1))+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(text = element_text(family = 'Arial',size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none')
ggsave('../results/TransCompR_sparsePCA/gene_sPC12_loadings.png',
       width = 14,
       height = 8,
       units = 'in',
       dpi = 600)

### Fisher exact test
fisherExactTest <- function(geneSet,Top,measGenes){
  not_in_geneset = sum(!(measGenes %in% geneSet ))
  contigencyMat <- data.frame('dex'=c(sum(geneSet %in% Top),length(geneSet)-sum(geneSet %in% Top)),
                              'not_dex'= c(sum(!(Top %in% geneSet)),not_in_geneset-sum(!(Top %in% geneSet))))
  rownames(contigencyMat) <- c('in_set','not_in_set')
  return(fisher.test(contigencyMat)$p.value)
}
egsea.data(species = "human",returnInfo = TRUE)
kegg_list <- kegg.pathways$human$kg.sets
top_genes <- loadings %>% filter(sPC12!=0)
top_genes <- top_genes %>% select(gene,sPC12) %>% column_to_rownames('gene')
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(top_genes), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
top_genes <- rownames(top_genes)
all_genes <- loadings %>% select(gene,sPC12) %>% column_to_rownames('gene')
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(all_genes), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
all_genes <- rownames(all_genes)
fisher_results <- unlist(lapply(kegg_list,fisherExactTest,top_genes,all_genes))
fisher_results_adj <- p.adjust(fisher_results,'BH')
inds <- which(fisher_results_adj<0.05)
pathways <- names(fisher_results_adj)[inds]
df_pathways <- data.frame(pathways)
df_pathways <- df_pathways %>% mutate(pathways=substr(pathways, 9, nchar(pathways)))

### Perform KEGG pathway analysis (GSEA)------------
rownames(loadings) <- NULL
loadings <- loadings %>% select(gene,sPC12) %>% column_to_rownames('gene')
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
  ggtitle('Pathways enriched in sPC12 loadings')+
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
  ggtitle('Pathways enriched in sPC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/TransCompR_sparsePCA/kegg_on_sPC12.png',
       plot=p2,
       width=16,
       height=9,
       units = 'in',
       dpi = 600)

### Perform GO Terms analysis (GSEA)------------
loadings <- sPC_B$loadings
colnames(loadings) <- paste0('sPC',seq(1:ncol(loadings)))
rownames(loadings) <- colnames(X_B)
loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,sPC12)
loadings$significant <- ''
loadings$significant[order(-abs(loadings$sPC12))[1:10]] <- loadings$gene[order(-abs(loadings$sPC12))[1:10]]
loadings <- loadings[order(loadings$sPC12),]
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
df_gos <- df_gos %>% mutate(sig=ifelse(abs(NES)>1.5 & padj<0.05,'yes','no'))
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
  ggtitle('GO Terms enriched in sPC12 loadings')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

df_gos <- df_gos[order(df_gos$NES),]
df_gos$`GO Terms` <- factor(df_gos$`GO Terms`,levels = df_gos$`GO Terms`)
p2 <- ggplot(df_gos %>% arrange(NES) %>% filter(padj<0.05 & abs(NES)>1.5),aes(x=NES,y=`GO Terms`,fill=padj))+ 
  geom_bar(stat = 'identity',color='black',size=1.5) +
  scale_fill_gradient(low = "red",high = "white",limits = c(min(df_gos$padj),0.05)) +
  xlab('Normalized Enrichment Score') +
  ggtitle('GO Terms enriched in sPC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/TransCompR_sparsePCA/go_on_sPC12.png',
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
df_tfs <- tfs %>% mutate(sig=ifelse(abs(NES)>1.4 & padj<0.05,'yes','no'))
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
  ggtitle('GO Terms enriched in sPC12 loadings')+
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
  ggtitle('TFs enriched in sPC12 loadings')+
  theme_pubr(base_family = 'Arial',base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.position = 'right',
        legend.justification = "center")
print(p2)
ggsave('../results/TransCompR_sparsePCA/tfs_on_sPC12.png',
       plot=p2,
       width=16,
       height=10,
       units = 'in',
       dpi = 600)
