library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(OmnipathR)
library(EGSEAdata)
library(igraph)
library(org.Hs.eg.db) 
library(tidyverse)

### Load data and each time filter to keep statistically significant pathways-------------------------------
geos <- c("GSE126848","GSE130970","GSE134422","GSE135251","GSE162694")

### OR LOGIC Internally
df_all <- data.frame()
for (geo in geos){
  df_gsea <- readRDS(paste0('../results/bioinformatics_results/',geo,'FinalGSEADF.rds'))
  df_gsea <- df_gsea %>% filter(padj<0.05) %>% filter(abs(NES)>0.5)
  df_gsea <- df_gsea %>% select(pathway) %>% mutate(dataset=geo) %>% unique()
  df_all <- rbind(df_all,df_gsea)
}
# saveRDS(df_all,'../results/bioinformatics_results/all_significant_keggs.rds')
print(paste0('All unique pathways appearing across datasets are ',length(unique(df_all$pathway))))

### Keep pathways appearing in at least 4 out of 5 since this is also in how many datasets the NAFLD pathway appears
df_filt <- df_all %>% group_by(pathway) %>% mutate(counts = n_distinct(dataset)) %>% ungroup() %>% 
  filter(counts>=4) %>% #filter(counts>=0.5*length(geos)) %>%
  select(-dataset) %>% unique()
print(paste0('Filtered pathways appearing across datasets are ',length(unique(df_filt$pathway))))
# saveRDS(df_filt,'../results/bioinformatics_results/all_filtered_significant_keggs.rds')

### Get core signaling network to use-------------------------------
egsea.data(species = "human",returnInfo = TRUE)
pathways <- kegg.pathways$human
pathways <- pathways$kg.sets
pathways <- pathways[df_filt$pathway]
genes <- c()
for (set in pathways){
  genes <- c(genes,set)
}
genes <- unique(genes)
all_gene_symbols <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENTREZID")
all_gene_symbols <- unname(all_gene_symbols)
inds <- which(!is.na(all_gene_symbols))
all_gene_symbols <- all_gene_symbols[inds]

### Get OmniPath and keep only the smallest subgraph containing all these nodes
interactions <- import_omnipath_interactions()
interactions <- interactions  %>% mutate(kegg=grepl(pattern="KEGG",x=sources))%>% 
  filter(n_resources>=3 | kegg==T)  %>%
  filter(n_references>=3 | kegg==T)  %>%
  dplyr::select(-kegg) %>%
  dplyr::select(c('source'='source_genesymbol'),
                c('target'='target_genesymbol'),
                is_inhibition,is_stimulation) %>% unique()
all_gene_symbols <- all_gene_symbols[which(all_gene_symbols %in% c(interactions$source,interactions$target))]

### Find the subgraph
g = graph.data.frame(interactions[,c('source','target')])
needNodes <- character()
## loop through all nodes and calculate the path to each other node
for(i in 1:(length(all_gene_symbols)-1)){
  paths <- shortest_paths(g, from=all_gene_symbols[i], 
                          to=all_gene_symbols[(i+1):length(all_gene_symbols)],
                          mode="all")
  needNodes <- unique(c(needNodes, unlist(lapply(paths$vpath, names))))
}
## subset the graph
subGr <- induced_subgraph(g, vids=needNodes)
saveRDS(subGr,'../results/core_signaling_net_igraph_strict.rds')

interactions2keep <- igraph::as_data_frame(subGr)
interactions2keep <- distinct(interactions2keep)

### Filter graph even more
# interactions_filt <- import_omnipath_interactions() %>% 
#   filter(source_genesymbol %in% interactions2keep$from) %>% filter(target_genesymbol %in% interactions2keep$to) %>% 
#   mutate(kegg=grepl(pattern="KEGG",x=sources))%>% 
#   filter(n_resources>=3 | kegg==T)  %>%
#   filter(n_references>=3 | kegg==T)  %>%
#   dplyr::select(-kegg) %>% 
#   dplyr::select(c('source'='source_genesymbol'),
#                 c('target'='target_genesymbol'),
#                 is_inhibition,is_stimulation) %>% unique()
interactions_filt <- interactions %>% 
  filter(source %in% interactions2keep$from) %>% filter(target %in% interactions2keep$to) %>% 
  unique()

## Now get rid of disconnected parts
g <- graph.data.frame(interactions_filt[,c('source','target')])
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
# ids
vert_ids <- V(g)[components$membership == biggest_cluster_id]
# final subgraph
final_subgraph <- induced_subgraph(g, vert_ids)
final_subgraph_df <- distinct(igraph::as_data_frame(final_subgraph))

final_interactions_filt <- interactions_filt %>% filter(source %in% final_subgraph_df$from) %>% 
  filter(target %in% final_subgraph_df$to)
saveRDS(final_interactions_filt,'../results/filtered_core_signaling_net_strict.rds')
