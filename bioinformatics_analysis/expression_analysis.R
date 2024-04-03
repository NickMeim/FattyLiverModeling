library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(factoextra)
library(org.Hs.eg.db) 
source('enrichment_calculations.R')
library(RColorBrewer)
library("gplots")

### Load desired data set-----------------------
geo <- "GSE130970"# one of c("GSE126848","GSE130970","GSE134422","GSE135251","GSE162694")
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
expression_matrix <- readRDS('../data/ARCHS4/FattyLiver_expression_matrix.rds')
expression_matrix <- expression_matrix[,meta_data$sample]

### Analyze covariates in PCA space------------------------------
keep_gene <- (rowSums(expression_matrix >= 1) >= 0.1*ncol(expression_matrix))
expression_matrix_normed <- expression_matrix[keep_gene,]
expression_matrix_normed <- apply(expression_matrix_normed, MARGIN = 2, FUN = function(x) x/sum(x)*1e6)
# keep_gene <- (rowSums(expression_matrix_normed >= 1) >= 0.1*ncol(expression_matrix))
# expression_matrix_normed <- expression_matrix_normed[keep_gene,]
expression_matrix_normed <- log2(1 + expression_matrix_normed)
expression_matrix_normed <- t(expression_matrix_normed - rowMeans(expression_matrix_normed))
pca_transformed <- prcomp(expression_matrix_normed,scale. = F,center = T)
### check PCs that explain at least 1% of the variance
#fviz_eig(pca_transformed)
perc_var <- 100*pca_transformed$sdev^2/sum(pca_transformed$sdev^2)
inds <- which(perc_var>=1)
df_pcs <- pca_transformed$x[,inds]
df_pcs <- as.data.frame(df_pcs) %>% rownames_to_column('sample')
df_pcs <- left_join(df_pcs,meta_data %>% dplyr::select(sample,all_of(new_cols)))
df_pcs <- df_pcs %>% dplyr::select(-tissue)
colnames(df_pcs)[(inds[length(inds)]+2):ncol(df_pcs)] <- c('age','ballooning','fibrosis','lobular_inflammation','nas','sex','steatosis')
anova_results <- df_pcs  %>% mutate(sex = ifelse(sex=='M',1,0)) %>% select(-sample,-ballooning,-lobular_inflammation,-steatosis) %>%
  gather('PC','value',-age,-fibrosis,-sex,-nas) %>% group_by(PC) %>%
  rstatix::anova_test(value~age+sex+nas+fibrosis) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  ungroup()
# anova_results <- df_pcs  %>% mutate(gender = ifelse(gender=='Male',1,0)) %>% 
#   mutate(disease = ifelse(disease=='healthy',0,
#                           ifelse(disease=='NAFLD',1,2))) %>% 
#   dplyr::select(-sample) %>%
#   gather('PC','value',-gender,-disease) %>% group_by(PC) %>%
#   rstatix::anova_test(value~gender+disease) %>% 
#   rstatix::adjust_pvalue(method = 'BH') %>%
#   ungroup()
anova_results <- as.data.frame(anova_results)
anova_results <- anova_results %>%
  mutate(significance = ifelse(p.adj<0.001,'***',
                               ifelse(p.adj<0.01,'**',
                                      ifelse(p.adj<0.05,'*',
                                             ifelse(p.adj<0.1,'.','')))))
p_F <- ggplot(anova_results,aes(x=PC,y=Effect,fill=`F`)) + 
  geom_tile(color='black') +
  geom_text(aes(label=significance),size=7)+
  scale_fill_gradient(low = 'white',high='red')+
  labs(fill = "F-statistic") + 
  ggtitle('ANOVA results')+
  theme_pubr() + 
  theme(legend.position = 'right',
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab('') + xlab('')
print(p_F)

# ggboxplot( df_pcs %>% 
#              gather('PC','value',-age,-fibrosis,-sex,-nas,-sample,-ballooning,-lobular_inflammation,-steatosis),
#            x='sex',y='value',
#           cor.coef = T,cor.coef.size = 4) + 
#   facet_wrap(~PC) +
#   stat_compare_means()
# ggscatter(df_pcs %>% gather('PC','value',-age,-fibrosis,-sex,-nas,-sample,-ballooning,-lobular_inflammation,-steatosis),
#            x='value',y='age',
#            cor.coef = T,cor.coef.size = 4) + 
#   facet_wrap(~PC)+
#   geom_smooth(se=TRUE,color='blue',method = 'lm') 

### Perform differential expression analysis with DESeq2--------------------------------------------------------------------------
keep_gene <- (rowSums(expression_matrix >= 1) >= 0.1*ncol(expression_matrix))
expression_matrix <- expression_matrix[keep_gene,]
colnames(meta_data)[(ncol(meta_data)-7):ncol(meta_data)] <- c('age','ballooning','fibrosis','lobular_inflammation','nas','sex','steatosis','tissue')
meta_data_normed <- meta_data
meta_data_normed$age <- (meta_data_normed$age-mean(meta_data_normed$age))/sd(meta_data_normed$age)
# Compare NAS 3 vs NAS 0

### I can use DESeq2 with multiple NAS stages and then thake the comparisons for stage X VS healthy 
### Or I can just fit a model only between each pair
meta_data_filt <- meta_data_normed %>% mutate(nas=ifelse(nas==0,'healthy',paste0('NAFL',nas)))
meta_data_filt$fibrosis <- as.factor(meta_data_filt$fibrosis)
dds <- DESeqDataSetFromMatrix(countData = expression_matrix[,meta_data_filt$sample], colData = meta_data_filt, design = ~ disease+gender)
dds <- DESeq(dds,betaPrior = FALSE)
comparisons <- list(c("nas","NAFL1","healthy"),c("nas","NAFL2","healthy"),c("nas","NAFL3","healthy"),
                  c("nas","NAFL4","healthy"),c("nas","NAFL5","healthy"),c("nas","NAFL6","healthy"))
# comparisons <- list(c("disease","NAFLD","healthy"),c("disease","NASH","healthy"))
all_results <- NULL
all_plots <- NULL
for (contrast in comparisons){
  res <- results(dds, contrast = contrast)
  res <- as.data.frame(res)
  res <- res %>% rownames_to_column('gene')
  res <- res %>% mutate(significance = ifelse(padj<0.05,ifelse(log2FoldChange>1,'upregulated',ifelse(log2FoldChange<(-1),'downregulated',NA)),NA))
  res <- res %>% mutate(label = ifelse(!is.na(significance),gene,NA))
  p <- ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=significance)) + geom_point() +
    xlab(expression(log[2]('Fold Change'))) + ylab(expression(-log[10]('adjusted p-value'))) +
    scale_color_manual(values=c("#EC2326","#008c41"))+
    geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "#525252", size=0.5) +
    geom_vline(xintercept=-1, linetype="dashed",color = "#EC2326", size=1) +
    geom_vline(xintercept=1, linetype="dashed",color = "#008c41", size=1)+ 
    # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
    geom_label_repel(aes(label=label),size=6,max.overlaps = 15)+
    ggtitle(paste0(contrast[2],' vs Healthy'))+
    theme_minimal(base_size=18,base_family = 'Arial')+
    theme(text = element_text(size=18,family = 'Arial'),
          plot.title = element_text(hjust=0.5))
  all_results[[paste0(contrast[2],' vs Healthy')]] <- res
  all_plots[[paste0(contrast[2],' vs Healthy')]] <- p
  print(paste0('Finished ',contrast[2],' vs Healthy'))
  
}
saveRDS(all_results,paste0('../results/bioinformatics_results/',geo,'_Deseq2_results.rds'))
saveRDS(all_plots,paste0('../results/bioinformatics_results/',geo,'_volcano_plots.rds'))

volcano_plot <- ggarrange(plotlist=all_plots,ncol=2,nrow=3,common.legend = TRUE,legend = 'bottom')
# print(volcano_plot)
ggsave(paste0('../results/bioinformatics_results/',geo,'_volcano_plots.png'),
       bg = 'white',
       plot=volcano_plot,
       width = 12,
       height = 18,
       units = 'in',
       dpi = 600)
### Run GSEA using logFC across all samples---------------------
#### Ftiaxto me to gProfiler i kati tetoio gt ayto einai thlipsis
df_res <- all_results$`NAFL1 vs Healthy`
df_res <- df_res %>% dplyr::select(gene,stat)
colnames(df_res)[2] <- 'NAFL1 vs Healthy'
for (i in 2:length(all_results)){
  res <- all_results[[i]]
  res <- res %>% dplyr::select(gene,stat)
  colnames(res)[2] <- names(all_results)[i]
  df_res <- left_join(df_res,res)
}

entrez_ids <- mapIds(org.Hs.eg.db, keys = df_res$gene, column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- df_res[inds,]
meas <- meas %>% dplyr::select(-gene)
rownames(meas) <- NULL
rownames(meas) <- entrez_ids
keggs <- fastenrichment(signature_ids = colnames(meas),
                        gene_ids = entrez_ids,
                        measurements=meas,
                        enrichment_space = 'kegg',
                        n_permutations = 10000)

keggs_nes <- as.data.frame(as.matrix(keggs[["NES"]]$`NES KEGG`)) %>% rownames_to_column('KEGG pathway')
keggs_nes <- keggs_nes %>% gather('background','NES',-`KEGG pathway`)
keggs_pval <- as.data.frame(as.matrix(keggs[["Pval"]]$`Pval KEGG`)) %>% rownames_to_column('KEGG pathway')
keggs_pval <- keggs_pval %>% gather('background','p.adj',-`KEGG pathway`)
df_keggs <- left_join(keggs_nes,keggs_pval)
df_keggs <- df_keggs %>% mutate(`KEGG pathway`=strsplit(`KEGG pathway`,"_"))
df_keggs <- df_keggs %>% unnest(`KEGG pathway`) %>% filter(!(`KEGG pathway` %in% c("KEGG","FL1000")))
df_keggs <- df_keggs %>% mutate(`KEGG pathway`=as.character(`KEGG pathway`))
genesets <- lapply(kegg.pathways[["human"]]$kg.sets,length)
names_genesets <- names(genesets)
genesets <- t(data.frame(genesets))
rownames(genesets) <- NULL
colnames(genesets) <- 'size'
genesets <- as.data.frame(genesets)
genesets$`KEGG pathway` <- names_genesets
df_keggs <- left_join(df_keggs,genesets)
df_keggs <- df_keggs %>% mutate(`KEGG pathway`=substr(`KEGG pathway`, 9, nchar(`KEGG pathway`)))
df_keggs <- df_keggs %>% mutate(sig=ifelse(abs(NES)>=1.5 & p.adj<=0.01,'yes','no'))
df_keggs <- df_keggs%>% mutate(label = ifelse(sig=='yes',`KEGG pathway`,NA))
df_keggs <- df_keggs %>% group_by(`KEGG pathway`) %>% mutate(mean_nes=mean(NES)) %>% mutate(sig_counts=sum(sig=='yes')) %>% ungroup()
df_keggs <- df_keggs[order(-df_keggs$mean_nes),]

ggplot(df_keggs %>% filter(sig=='yes'),aes(x=background,y=`KEGG pathway`)) + 
  geom_point( aes(fill = NES, size = size), shape=21) + xlab('Comparison') + 
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0) +
  theme_minimal(base_size = 20,base_family = 'Arial')
df_keggs <- df_keggs %>% mutate(significance = ifelse(p.adj<0.01,ifelse(NES>=2.5,'upregulated',ifelse(NES<(-2),'downregulated',NA)),NA))
show_paths <- c(' Non-alcoholic fatty liver disease (NAFLD)','Alcoholism',' Fructose and mannose metabolism',' Toll-like receptor signaling pathway',
                ' TNF signaling pathway',' TGF-beta signaling pathway','Apoptosis',
                ' Insulin signaling pathway',' Insulin secretion')
tmp <- df_keggs %>% filter(!is.na(significance)) %>% filter(background=='NAFL1 vs Healthy') %>% filter(NES<0)
show_paths <- c(show_paths,unique(tmp$`KEGG pathway`))
show_paths <- unique(show_paths)
df_keggs <- df_keggs %>% mutate(significance = ifelse(p.adj<0.05,ifelse(NES>1.5,'upregulated',ifelse(NES<(-1.5),'downregulated',NA)),NA))
# df_keggs <- df_keggs %>% mutate(label = ifelse(!is.na(significance),`KEGG pathway`,NA))
df_keggs <- df_keggs %>% mutate(label = ifelse(((!is.na(significance)) & (`KEGG pathway` %in% show_paths)),
                                               `KEGG pathway`,NA))

for (comp in unique(df_keggs$background)){
  color_values <- c("#EC2326","#008c41")
  if (comp=='NAFL3 vs Healthy'){
    color_values <- c("#008c41")
  }
  kegg_plot <- ggplot(df_keggs %>% filter(background==comp),aes(x=NES,y=-log10(p.adj),color=significance)) + geom_point() +
    xlab('Normalized Enrichment Score') + ylab(expression(-log[10]('adjusted p-value'))) +
    scale_color_manual(values=color_values)+
    geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "#525252", size=0.5) +
    geom_vline(xintercept=-1.5, linetype="dashed",color = "#EC2326", size=1) +
    geom_vline(xintercept=1.5, linetype="dashed",color = "#008c41", size=1)+ 
    # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
    geom_label_repel(aes(label=label),size=6,max.overlaps = 15,min.segment.length = 0.1,max.iter = 50000)+
    ggtitle(comp)+
    theme_minimal(base_size=24,base_family = 'Arial')+
    theme(text = element_text(size=24,family = 'Arial'),
          plot.title = element_text(hjust=0.5))
  ggsave(paste0('../results/bioinformatics_results/',geo,'_',substr(comp,1,5),'_kegg_pathways.png'),
         bg = 'white',
         plot=kegg_plot,
         width = 12,
         height = 12,
         units = 'in',
         dpi = 600)
}


#### GO-Term enrichment-----------------------------------
meas <- df_res
meas <- meas %>% column_to_rownames('gene')
gos <- fastenrichment(signature_ids = colnames(meas),
                        gene_ids = rownames(meas),
                        measurements=meas,
                        enrichment_space = 'go_bp',
                        gene_id_type = 'symbol',
                        n_permutations = 10000)
gos_nes <- as.data.frame(as.matrix(gos[["NES"]]$`NES GO BP`)) %>% rownames_to_column('GO')
gos_nes <- gos_nes %>% gather('background','NES',-`GO`)
gos_pval <- as.data.frame(as.matrix(gos[["Pval"]]$`Pval GO BP`)) %>% rownames_to_column('GO')
gos_pval <- gos_pval %>% gather('background','p.adj',-`GO`)
df_gos <- left_join(gos_nes,gos_pval)
df_gos <- df_gos %>% mutate(`GO`=strsplit(`GO`,"_"))
df_gos <- df_gos %>% unnest(`GO`) %>% filter(!(`GO` %in% c("GO","BP","FL1000")))
df_gos <- df_gos %>% mutate(`GO`=as.character(`GO`))
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
df_gos <- left_join(df_gos,go_annotations %>% filter(ontology=='BP') %>% dplyr::select(GO,`GO Terms`)) 

df_gos <- df_gos %>% mutate(significance = ifelse(p.adj<0.01,ifelse(NES>=2.5,'upregulated',ifelse(NES<(-2),'downregulated',NA)),NA))
df_gos <- df_gos %>% mutate(label = ifelse(!is.na(significance),
                                           `GO Terms`,NA))
df_gos <- df_gos %>% mutate(significance = ifelse(p.adj<0.05,ifelse(NES>1.5,'upregulated',ifelse(NES<(-1.5),'downregulated',NA)),NA))

for (comp in unique(df_gos$background)){
  color_values <- c("#EC2326","#008c41")
  # if (comp=='NAFL3 vs Healthy'){
  #   color_values <- c("#008c41")
  # }
  go_plot <- ggplot(df_gos %>% filter(background==comp),aes(x=NES,y=-log10(p.adj),color=significance)) + geom_point() +
    xlab('Normalized Enrichment Score') + ylab(expression(-log[10]('adjusted p-value'))) +
    scale_color_manual(values=color_values)+
    geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "#525252", size=0.5) +
    geom_vline(xintercept=-1.5, linetype="dashed",color = "#EC2326", size=1) +
    geom_vline(xintercept=1.5, linetype="dashed",color = "#008c41", size=1)+ 
    # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
    geom_label_repel(aes(label=label),size=6,max.iter = 50000,max.overlaps = 40)+
    ggtitle(comp)+
    theme_minimal(base_size=24,base_family = 'Arial')+
    theme(text = element_text(size=24,family = 'Arial'),
          plot.title = element_text(hjust=0.5))
  ggsave(paste0('../results/bioinformatics_results/',geo,'_',substr(comp,1,5),'_goterms.png'),
         bg = 'white',
         plot=go_plot,
         width = 12,
         height = 12,
         units = 'in',
         dpi = 600)
}
