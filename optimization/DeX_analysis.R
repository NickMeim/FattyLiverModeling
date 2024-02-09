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
# library(decoupleR)
# source('utils-omnipath.R')

# Load new data----------------------------------------------------------
data_lean <- data.table::fread('../results/optimized_mps/pertubed_control_lean.csv')
data_lean <- data_lean %>% column_to_rownames('V1')
data_fatty <- data.table::fread('../results/optimized_mps/pertubed_control_fatty.csv')
data_fatty <- data_fatty %>% column_to_rownames('V1')
conditions <- c(rep('control',12),rep('perturbed',12))
design <- model.matrix(~0+conditions)
colnames(design) <- c('control','perturbed')

# # load dx
# dx_lean <-  data.table::fread('../results/optimized_mps/dx_lean.csv')
# dx_lean$V1 <- 'lean'
# colnames(dx_lean)[1] <- 'background' 
# dx_fatty <-  data.table::fread('../results/optimized_mps/dx_fatty.csv')
# dx_fatty$V1 <- 'fatty'
# colnames(dx_fatty)[1] <- 'background'

# use limma with lean transformed data-------------------------------------
fit_lean <- lmFit(data_lean,design,method = 'ls')
# plotSA(fit_lean, main="Final model: Mean-variance trend")
contr_lean <- makeContrasts(contrasts = 'perturbed-control',levels = colnames(coef(fit_lean)))
contr_fit_lean <- contrasts.fit(fit_lean, contr_lean)
fit_ebayes_lean <- eBayes(contr_fit_lean)
# plotSA(fit_ebayes_lean, main="Final model: Mean-variance trend")
top.table_lean <- topTable(fit_ebayes_lean,number = nrow(data_lean))
# Make adjustements to get the desired plots
top.table_lean <- top.table_lean %>% mutate(significance = ifelse(adj.P.Val<0.05,ifelse(logFC>1.5,'upregulated',ifelse(logFC<(-1.5),'downregulated',NA)),NA))
top.table_lean <- top.table_lean %>% rownames_to_column('gene')
top.table_lean <- top.table_lean %>% mutate(label = ifelse(!is.na(significance),gene,NA))
ggplot(top.table_lean,aes(x=logFC,y=-log10(adj.P.Val),color=significance)) + geom_point() +
  xlab(expression(log[2]('Fold Change'))) + ylab(expression(-log[10]('adjusted p-value'))) +
  scale_color_manual(values=c("#EC2326","#008c41"))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "#525252", size=0.5) +
  geom_vline(xintercept=-1.5, linetype="dashed",color = "#EC2326", size=1) +
  geom_vline(xintercept=1.5, linetype="dashed",color = "#008c41", size=1)+ 
  # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
  geom_label_repel(aes(label=label),size=6)+
  ggtitle('Proposed perturbations for lean-derived samples')+
  theme(text = element_text(size=18,family = 'Arial'),
        plot.title = element_text(hjust=0.5))
ggsave('../results/optimized_mps/lean_volcano.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)

# use limma with fatty transformed data------------------------------------
fit_fatty <- lmFit(data_fatty,design,method = 'ls')
# plotSA(fit_fatty, main="Final model: Mean-variance trend")
contr_fatty <- makeContrasts(contrasts = 'perturbed-control',levels = colnames(coef(fit_fatty)))
contr_fit_fatty <- contrasts.fit(fit_fatty, contr_fatty)
fit_ebayes_fatty <- eBayes(contr_fit_fatty)
# plotSA(fit_ebayes_fatty, main="Final model: Mean-variance trend")
top.table_fatty <- topTable(fit_ebayes_fatty,number = nrow(data_fatty))
# Make adjustements to get the desired plots
top.table_fatty <- top.table_fatty %>% mutate(significance = ifelse(adj.P.Val<0.05,ifelse(logFC>1.5,'upregulated',ifelse(logFC<(-1.5),'downregulated',NA)),NA))
top.table_fatty <- top.table_fatty %>% rownames_to_column('gene')
top.table_fatty <- top.table_fatty %>% mutate(label = ifelse(!is.na(significance),gene,NA))
ggplot(top.table_fatty,aes(x=logFC,y=-log10(adj.P.Val),color=significance)) + geom_point() +
  xlab(expression(log[2]('Fold Change'))) + ylab(expression(-log[10]('adjusted p-value'))) +
  scale_color_manual(values=c("#EC2326","#008c41"))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "#525252", size=0.5) +
  geom_vline(xintercept=-1.5, linetype="dashed",color = "#EC2326", size=1) +
  geom_vline(xintercept=1.5, linetype="dashed",color = "#008c41", size=1)+ 
  # annotate("text",x=-2,y=1.5,label="adjusted p-value=0.05",size=6)+
  geom_label_repel(aes(label=label),size=6)+
  ggtitle('Proposed perturbations for fatty-derived samples')+
  theme(text = element_text(size=18,family = 'Arial'),
        plot.title = element_text(hjust=0.5))
ggsave('../results/optimized_mps/fatty_volcano.png',
       width = 12,
       height = 9,
       units = 'in',
       dpi = 600)

### Perform dorothea analysis using logFCs OR directly dXs ---------------
dex_data <- rbind(top.table_lean %>% mutate(background = 'lean') %>% select(gene,logFC,background),
                  top.table_fatty %>% mutate(background = 'fatty') %>% select(gene,logFC,background))
dex_data <- dex_data %>% spread('gene','logFC')
# dex_data <- rbind(dx_lean,dx_fatty)
dex_data <- dex_data %>% column_to_rownames('background')
dex_data <- t(dex_data)
hist(as.matrix(dex_data))

### get high-confidence dorothea regulon
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]

settings = list(verbose = TRUE, minsize = 5)
TF_activities = run_viper(dex_data, dorotheaData, options =  settings)

### visualize differential activity
TF_activities <- as.data.frame(TF_activities) %>% rownames_to_column('TF') %>% gather('background','activity',-TF) 
TF_activities <- TF_activities %>% mutate(sig = ifelse(abs(activity)>=2,'yes','no'))
TF_activities <- TF_activities %>% mutate(label = ifelse(sig=='yes',TF,NA))
# p1 <- ggbarplot(TF_activities %>% filter(background=='lean') %>% arrange(activity),x='TF',y='activity',fill='sig')+ 
#   scale_fill_manual(values = c('#b30024','#00b374'))+
#   ggtitle('Lean')+
#   theme(text=element_text(family = 'Arial',size=18),
#         legend.position = 'none',
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(size=12,angle = 90))
p1 <- ggplot(TF_activities %>% filter(background=='lean') %>% arrange(activity),aes(x=activity,y=reorder(TF,abs(activity)),color=sig))+ 
  geom_point()+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Lean')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
        panel.grid.major.y = element_blank())
print(p1)

# p2 <- ggbarplot(TF_activities %>% filter(background=='fatty') %>% arrange(activity),x='TF',y='activity',fill='sig')+ 
#   scale_fill_manual(values = c('#b30024','#00b374'))+
#   ggtitle('Fatty')+
#   theme(text=element_text(family = 'Arial',size=18),
#         legend.position = 'none',
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(size=12,angle = 90))
p2 <- ggplot(TF_activities %>% filter(background=='fatty') %>% arrange(activity),aes(x=activity,y=reorder(TF,abs(activity)),color=sig))+ 
  geom_point()+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Fatty')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
        panel.grid.major.y = element_blank())
print(p2)

### do it with GSEA
tfs <- fastenrichment(signature_ids = colnames(dex_data),
                      gene_ids = rownames(dex_data),
                      measurements=dex_data,
                      enrichment_space = 'tf_dorothea',
                      tf_path = '../Artificial-Signaling-Network/DrugsANNSignaling/data/dorothea.tsv',
                      n_permutations = 10000)
tfs_nes <- as.data.frame(as.matrix(tfs[["NES"]]$`NES TF`)) %>% rownames_to_column('TF')
tfs_nes <- tfs_nes %>% gather('background','NES',-TF)
tfs_pval <- as.data.frame(as.matrix(tfs[["Pval"]]$`Pval TF`)) %>% rownames_to_column('TF')
tfs_pval <- tfs_pval %>% gather('background','p.adj',-TF)
df_tfs <- left_join(tfs_nes,tfs_pval)
df_tfs <- df_tfs %>% mutate(`TF`=strsplit(`TF`,"_"))
df_tfs <- df_tfs %>% unnest(`TF`) %>% filter(!(`TF` %in% c("TF","DOROTHEA","FL1000")))
df_tfs <- df_tfs %>% mutate(sig=ifelse(abs(NES)>1.5 & p.adj<0.05,'yes','no'))
df_tfs <- df_tfs%>% mutate(label = ifelse(sig=='yes',TF,NA))
p1 <- ggplot(df_tfs %>% filter(background=='lean') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Lean')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

p2 <-  ggplot(df_tfs %>% filter(background=='fatty') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Fatty')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p2)

p <- p1/p2
print(p)

ggsave('../results/optimized_mps/tfs_enrichment_volcano.png',
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

### Perform GSEA for KEGG pathways ---------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(dex_data), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- dex_data[inds,]
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
df_keggs <- df_keggs %>% mutate(`KEGG pathway`=substr(`KEGG pathway`, 9, nchar(`KEGG pathway`)))
df_keggs <- df_keggs %>% mutate(sig=ifelse(abs(NES)>1.5 & p.adj<0.05,'yes','no'))
df_keggs <- df_keggs%>% mutate(label = ifelse(sig=='yes',`KEGG pathway`,NA))
p1 <- ggplot(df_keggs %>% filter(background=='lean') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Lean')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

p2 <-  ggplot(df_keggs %>% filter(background=='fatty') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Fatty')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p2)

p <- p1/p2
print(p)

ggsave('../results/optimized_mps/kegg_enrichment_volcano.png',
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)


#### GSEA for GO Terms of Biological Processes--------------------------------
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
go_bps <- fastenrichment(signature_ids = colnames(dex_data),
                        gene_ids = rownames(dex_data),
                        measurements=dex_data,
                        enrichment_space = 'go_bp',
                        gene_id_type = 'symbol',
                        n_permutations = 10000)
gos_bp_nes <- as.data.frame(as.matrix(go_bps[["NES"]]$`NES GO BP`)) %>% rownames_to_column('GO')
gos_bp_nes <- gos_bp_nes %>% gather('background','NES',-`GO`)
gos_bp_pval <- as.data.frame(as.matrix(go_bps[["Pval"]]$`Pval GO BP`)) %>% rownames_to_column('GO')
gos_bp_pval <- gos_bp_pval %>% gather('background','p.adj',-`GO`)
df_gos_bp <- left_join(gos_bp_nes,gos_bp_pval)
df_gos_bp <- df_gos_bp %>% mutate(`GO`=strsplit(`GO`,"_"))
df_gos_bp <- df_gos_bp %>% unnest(`GO`) %>% filter(!(`GO` %in% c("GO","BP","FL1000")))
df_gos_bp <- left_join(df_gos_bp,go_annotations %>% filter(ontology=='BP') %>% dplyr::select(GO,`GO Terms`)) 
df_gos_bp <- df_gos_bp %>% mutate(sig=ifelse(abs(NES)>=2 & p.adj<0.05,'yes','no'))
df_gos_bp <- df_gos_bp%>% mutate(label = ifelse(sig=='yes',`GO Terms`,NA))
p1 <- ggplot(df_gos_bp %>% filter(background=='lean') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Lean')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

p2 <-  ggplot(df_gos_bp %>% filter(background=='fatty') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Fatty')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p2)

p <- p1/p2
print(p)

ggsave('../results/optimized_mps/gos_bp_enrichment_volcano.png',
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

#### GSEA for GO Terms of Molecular Functions--------------------------------
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
go_mfs <- fastenrichment(signature_ids = colnames(dex_data),
                         gene_ids = rownames(dex_data),
                         measurements=dex_data,
                         enrichment_space = 'go_mf',
                         gene_id_type = 'symbol',
                         n_permutations = 10000)
gos_mf_nes <- as.data.frame(as.matrix(go_mfs[["NES"]]$`NES GO MF`)) %>% rownames_to_column('GO')
gos_mf_nes <- gos_mf_nes %>% gather('background','NES',-`GO`)
gos_mf_pval <- as.data.frame(as.matrix(go_mfs[["Pval"]]$`Pval GO MF`)) %>% rownames_to_column('GO')
gos_mf_pval <- gos_mf_pval %>% gather('background','p.adj',-`GO`)
df_gos_mf <- left_join(gos_mf_nes,gos_mf_pval)
df_gos_mf <- df_gos_mf %>% mutate(`GO`=strsplit(`GO`,"_"))
df_gos_mf <- df_gos_mf %>% unnest(`GO`) %>% filter(!(`GO` %in% c("GO","MF","FL1000")))
df_gos_mf <- left_join(df_gos_mf,go_annotations %>% filter(ontology=='MF') %>% dplyr::select(GO,`GO Terms`)) 
df_gos_mf <- df_gos_mf %>% mutate(sig=ifelse(abs(NES)>=2 & p.adj<0.05,'yes','no'))
df_gos_mf <- df_gos_mf%>% mutate(label = ifelse(sig=='yes',`GO Terms`,NA))
p1 <- ggplot(df_gos_mf %>% filter(background=='lean') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Lean')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

p2 <-  ggplot(df_gos_mf %>% filter(background=='fatty') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Fatty')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p2)

p <- p1/p2
print(p)

ggsave('../results/optimized_mps/gos_mf_enrichment_volcano.png',
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)

#### GSEA for GO Terms of Cellular Components--------------------------------
go_annotations <- data.frame(GOs = Term(GOTERM),
                             'GO Terms' = GOID(GOTERM),
                             definition = Definition(GOTERM),
                             ontology = Ontology(GOTERM))
colnames(go_annotations) <- c('GO Terms','GO','definition','ontology')
go_ccs <- fastenrichment(signature_ids = colnames(dex_data),
                         gene_ids = rownames(dex_data),
                         measurements=dex_data,
                         enrichment_space = 'go_cc',
                         gene_id_type = 'symbol',
                         n_permutations = 10000)
gos_cc_nes <- as.data.frame(as.matrix(go_ccs[["NES"]]$`NES GO CC`)) %>% rownames_to_column('GO')
gos_cc_nes <- gos_cc_nes %>% gather('background','NES',-`GO`)
gos_cc_pval <- as.data.frame(as.matrix(go_ccs[["Pval"]]$`Pval GO CC`)) %>% rownames_to_column('GO')
gos_cc_pval <- gos_cc_pval %>% gather('background','p.adj',-`GO`)
df_gos_cc <- left_join(gos_cc_nes,gos_cc_pval)
df_gos_cc <- df_gos_cc %>% mutate(`GO`=strsplit(`GO`,"_"))
df_gos_cc <- df_gos_cc %>% unnest(`GO`) %>% filter(!(`GO` %in% c("GO","CC","FL1000")))
df_gos_cc <- left_join(df_gos_cc,go_annotations %>% filter(ontology=='CC') %>% dplyr::select(GO,`GO Terms`)) 
df_gos_cc <- df_gos_cc %>% mutate(sig=ifelse(abs(NES)>=2 & p.adj<0.05,'yes','no'))
df_gos_cc <- df_gos_cc%>% mutate(label = ifelse(sig=='yes',`GO Terms`,NA))
p1 <- ggplot(df_gos_cc %>% filter(background=='lean') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Lean')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p1)

p2 <-  ggplot(df_gos_cc %>% filter(background=='fatty') %>% arrange(NES),aes(x=NES,y=-log10(p.adj),color=sig))+ 
  geom_point()+
  ylab(expression(-log[10]('p.adj')))+
  geom_label_repel(aes(label=label),
                   size=4,
                   box.padding = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 40)+
  scale_x_continuous(n.breaks = 10)+
  ggtitle('Fatty')+
  theme(text=element_text(family = 'Arial',size=18),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())
print(p2)

p <- p1/p2
print(p)

ggsave('../results/optimized_mps/gos_cc_enrichment_volcano.png',
       width = 16,
       height = 12,
       units = 'in',
       dpi = 600)
