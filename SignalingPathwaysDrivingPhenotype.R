library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(ggforce)
library(ggradar)
library(ggiraphExtra)
library(grid)
library(ggExtra)
library(matrixStats)
source('modeling/CrossValidationUtilFunctions.R')
source('modeling/functions_translation.R')
source("utils/plotting_functions.R")
source("modeling/vector_space_interpretation.R")

## Load the all the data to be used----------------
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load
data_list <- load_datasets(dataset_names, dir_data = 'data/')
# Run PCA
tmp <- process_datasets(data_list, filter_variance = F)
data_list <- tmp$data_list
plt_list <- tmp$plt_list
# Define matrices of interest
Yh <- as.matrix(data_list[[ref_dataset]]$metadata  %>% select(nas_score,Fibrosis_stage)) #keep both Fibrosis and NAS
colnames(Yh) <- c('NAS','fibrosis')
Xh <- data_list[[ref_dataset]]$data_center %>% t()
# Get matrices of target dataset
Xm <- data_list[[target_dataset]]$data_center %>% t()
# Get Wm as the PC space of the MPS data when averaging tech replicates to capture variance due to experimental factors
Wm <- data_list[[target_dataset]]$Wm_group %>% as.matrix()

### Run PLSR and find extra basis--------------------------
plsr_model <- opls(x = Xh, 
                   y = Yh,
                   predI = 8,
                   crossvalI = 1,
                   scaleC = "center",
                   fig.pdfC = "none",
                   info.txtC = "none")
vip_scores <- getVipVn(plsr_model)
# Select variables with VIP > 1
important_vars <- names(vip_scores[vip_scores > 1])
cat("Number of important variables:", length(important_vars), "\n")
print(important_vars)

# Subset your X matrix to only important variables
Xh_filtered <- Xh[, important_vars]

# Get Wh of PLSR
Wh <- matrix(data = 0, ncol = ncol(plsr_model@weightMN), nrow = ncol(Xh))
rownames(Wh) <- colnames(Xh)
colnames(Wh) <- colnames(plsr_model@weightMN)
for (ii in 1:nrow(plsr_model@weightMN)){
  Wh[rownames(plsr_model@weightMN)[ii], ] <- plsr_model@weightMN[ii,]
}
# Get regression coefficients
Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN
Wh_red <- plsr_model@weightMN
## Find extra latent variables
phi <- Wh %*% Bh
Wm_opt <- analytical_solution_opt(y=Yh,
                                  W_invitro = Wm,
                                  phi = phi)
Wm_tot <- cbind(Wm, Wm_opt)

### Find translatable LV of the in vitro system
### Run evolutionary algorithm
Wm_combo <- get_translatable_LV_2phenotype(Xh, Yh, Wh, Wm,Bh)
Wm_combo <- Wm_combo$Wm_TC

rownames(Wm_combo) <- rownames(Wm)
rownames(Wm_tot) <- rownames(Wm)
rownames(Wm_opt) <- rownames(Wm)
rownames(Wm_combo) <- rownames(Wm)
rownames(Wh) <- rownames(Wm)

### Analyze gene weights (Wh)----------------------------------
# First find the correlation of each LV of PLSR with phenotype
Zh_plsr <- plsr_model@scoreMN
mcor <- cor(cbind(Zh_plsr,Yh))
mcor[upper.tri(mcor,diag = TRUE)] <- 100
mcor <- reshape2::melt(mcor)
mcor <- mcor %>% filter(value != 100)
mcor <- mcor %>% mutate(keep = ifelse(Var1=='fibrosis' | Var2=='fibrosis' | Var1=='NAS' | Var2=='NAS',TRUE,FALSE)) %>%
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(keep = ifelse((Var1=='fibrosis' & Var2=='NAS') | (Var1=='NAS' & Var2=='fibrosis'),FALSE,TRUE))%>%
  filter(keep==TRUE) %>% select(-keep) 
## Also infer pathway activity of in vivo data
net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
pathway_human  <- decoupleR::run_viper(t(Xh), net_prog,minsize = 1,verbose = TRUE) %>% select(-statistic)

invivo_plsr <- as.data.frame(Zh_plsr[,c('p1','p2','p3')])
colnames(invivo_plsr) <- c('LV1','LV2','LV3')
invivo_plsr <- cbind(invivo_plsr,as.data.frame(Yh))
invivo_plsr <- invivo_plsr %>% rownames_to_column('condition')
invivo_plsr <- left_join(invivo_plsr,pathway_human)

# Test enrichment of pathways in progeny using VIP scores
vip_scores <- getVipVn(plsr_model)
# Background = ALL genes in your model
background_genes <- names(vip_scores)
# Hit list = important genes above threshold
threshold <- 1.0   # or try 0.8, 1.5
hits <- names(vip_scores[vip_scores > threshold])
# Convert to a named list of gene sets (one vector per pathway)
progeny_sets <- net_prog %>%
  group_by(source) %>%
  summarise(genes = list(target)) %>%
  tibble::deframe()
run_fisher_ora <- function(hits, background, gene_sets, min_overlap = 3) {
  
  results <- lapply(names(gene_sets), function(pathway) {
    
    pathway_genes <- intersect(gene_sets[[pathway]], background)
    
    # Skip pathways with insufficient overlap
    if (length(intersect(hits, pathway_genes)) < min_overlap) {
      return(NULL)
    }
    
    # Build 2x2 contingency table
    a <- length(intersect(hits, pathway_genes))          # hit & in pathway
    b <- length(setdiff(hits, pathway_genes))            # hit & not in pathway
    c <- length(intersect(setdiff(background, hits),     # not hit & in pathway
                          pathway_genes))
    d <- length(setdiff(background,                      # not hit & not in pathway
                        union(hits, pathway_genes)))
    
    cont_table <- matrix(c(a, c, b, d), nrow = 2,
                         dimnames = list(c("VIP_hit", "VIP_not"),
                                         c("In_pathway", "Not_pathway")))
    
    ft  <- fisher.test(cont_table, alternative = "greater")
    
    data.frame(
      pathway       = pathway,
      n_pathway     = length(pathway_genes),     # pathway size in background
      n_overlap     = a,                          # key overlap
      n_hits        = length(hits),
      odds_ratio    = ft$estimate,
      p_value       = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine and multiple-testing correct
  results_df <- do.call(rbind, Filter(Negate(is.null), results))
  results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
  results_df %>% arrange(p_adj)
}
ora_results <- run_fisher_ora(
  hits       = hits,
  background = background_genes,
  gene_sets  = progeny_sets,
  min_overlap = 3
)

thresholds <- c(0.7,0.8, 1.0, 1.2, 1.5, 2.0,2.5,3,3.5,4)
threshold_results <- lapply(thresholds, function(thr) {
  hits_thr <- names(vip_scores[vip_scores > thr])
  if (length(hits_thr) < 10) return(NULL)   # skip if too few genes
  
  res <- run_fisher_ora(hits_thr, background_genes, progeny_sets)
  res$threshold <- thr
  res$n_hits    <- length(hits_thr)
  res
})
sensitivity_df <- bind_rows(Filter(Negate(is.null), threshold_results))
# Which pathways are consistently significant?
sensitivity_df <- sensitivity_df %>%
  filter(p_adj < 0.05) %>%
  group_by(pathway) %>%
  summarise(n = n()) %>%
  mutate(perc_sig = n/length(thresholds)) %>%
  arrange(-perc_sig) %>%
  filter(perc_sig>0.5)


pPlsr <- ggplot(invivo_plsr %>% filter(source %in% c(sensitivity_df$pathway,'JAK-STAT')),
       aes(x=LV1,y=LV3,fill=score)) +
  geom_point(size = 1.5, shape = 21, color = "black")+
  scale_fill_viridis_c()+
  labs(fill = 'Score')+
  xlab('Human LV1') +
  ylab('Human LV3') +
  facet_wrap(~source) +
  theme_bw()
print(pPlsr)
ggsave(plot=pPlsr,
       filename = 'figures/suppl_VIP_path_enrichment_3.png',
       units = 'cm',
       width = 16,
       height = 12,
       dpi=600)

  

### Analyze gene weights at the TF activity level--------------
dorotheaData = read.table('data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
extra_basis_TF_activity <- TF_activity_interpretation(Wh,
                                                      Wm,
                                                      dorotheaData)
tf_activity <- extra_basis_TF_activity[[2]] %>% select(TF,condition,score,p_value,significant)
print(mcor %>% filter(Var1=='NAS') %>% arrange(-abs(value)))
print(mcor %>% filter(Var1!='NAS') %>% arrange(-abs(value)))
p <- (ggplot(tf_activity %>% select(c('activity'='score'),TF,p_value,condition,significant) %>% 
               filter(condition=='p1'),
             aes(x=as.numeric(reorder(TF,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
        geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
        scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.001,1))+
        scale_y_continuous(n.breaks = 6,limits = c(-10,10))+
        ggtitle('LV extra 1')+
        xlab('Rank')+
        theme_pubr(base_size = 20,base_family = 'Arial')+
        theme(text = element_text(size = 20,family = 'Arial'),
              legend.position = 'right',
              plot.title = element_text(hjust = 0.5))) 
print(p)

### Analyze gene weights at the Pathway activity level--------------
extra_basis_pathway_activity <- pathway_activity_interpretation(Wh,
                                                                Wm)
pa <- (ggplot(extra_basis_pathway_activity %>% mutate(significant = ifelse(p_value<0.05,Pathway,NA)) %>%
               select(c('activity'='score'),Pathway,p_value,condition,significant) %>% 
               filter(condition=='p1'),
             aes(x=as.numeric(reorder(Pathway,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
        geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
        scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.001,1))+
        scale_y_continuous(n.breaks = 6,limits = c(-15,15))+
        ggtitle('LV1')+
        xlab('Rank')+
        theme_pubr(base_size = 20,base_family = 'Arial')+
        theme(text = element_text(size = 20,family = 'Arial'),
              legend.position = 'right',
              plot.title = element_text(hjust = 0.5))) +
  (ggplot(extra_basis_pathway_activity %>% mutate(significant = ifelse(p_value<0.05,Pathway,NA)) %>%
            select(c('activity'='score'),Pathway,p_value,condition,significant) %>% 
            filter(condition=='p2'),
          aes(x=as.numeric(reorder(Pathway,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
     geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
     scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.001,1))+
     scale_y_continuous(n.breaks = 6,limits = c(-15,15))+
     ggtitle('LV2')+
     xlab('Rank')+
     theme_pubr(base_size = 20,base_family = 'Arial')+
     theme(text = element_text(size = 20,family = 'Arial'),
           legend.position = 'right',
           plot.title = element_text(hjust = 0.5)))
pb <- (ggplot(extra_basis_pathway_activity %>% mutate(significant = ifelse(p_value<0.05,Pathway,NA)) %>%
                select(c('activity'='score'),Pathway,p_value,condition,significant) %>% 
                filter(condition=='p3'),
              aes(x=as.numeric(reorder(Pathway,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
         geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
         scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.001,1))+
         scale_y_continuous(n.breaks = 6,limits = c(-15,15))+
         ggtitle('LV3')+
         xlab('Rank')+
         theme_pubr(base_size = 20,base_family = 'Arial')+
         theme(text = element_text(size = 20,family = 'Arial'),
               legend.position = 'right',
               plot.title = element_text(hjust = 0.5))) +
  (ggplot(extra_basis_pathway_activity %>% mutate(significant = ifelse(p_value<0.05,Pathway,NA)) %>%
            select(c('activity'='score'),Pathway,p_value,condition,significant) %>% 
            filter(condition=='p4'),
          aes(x=as.numeric(reorder(Pathway,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
     geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
     scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.001,1))+
     scale_y_continuous(n.breaks = 6,limits = c(-15,15))+
     ggtitle('LV4')+
     xlab('Rank')+
     theme_pubr(base_size = 20,base_family = 'Arial')+
     theme(text = element_text(size = 20,family = 'Arial'),
           legend.position = 'right',
           plot.title = element_text(hjust = 0.5)))
p <- pa / pb
print(p)
ggsave(plot=p,
       filename = 'figures/suppl_Wh_path_enrichment.png',
       units = 'cm',
       height = 24,
       width = 24,
       dpi=600)

### Analyze gene weights with GSEA on MSIG Hallmarks genesets--------------
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(Wh), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(Wh[inds,])
rownames(meas) <- entrez_ids
msig <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes <- as.data.frame(msig$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark')  #%>% gather('PC','NES',-Hallmark)
msig_nes <- msig_nes %>% gather('LV','NES',-Hallmark)
msig_pval <- as.data.frame(msig$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')#%>% gather('PC','padj',-Hallmark)
msig_pval <- msig_pval %>% gather('LV','padj',-Hallmark)
df_msig <- left_join(msig_nes,msig_pval)
df_msig <- df_msig %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))
df_msig <- df_msig %>% mutate(Hallmark=str_replace_all(Hallmark,'_',' '))
df_msig <- df_msig %>% mutate(Hallmark = tolower(Hallmark)) %>% 
  mutate(Hallmark = paste0(toupper(substr(Hallmark, 1, 1)), tolower(substr(Hallmark, 2, nchar(Hallmark)))))



p1 <- (ggplot(df_msig %>% filter(LV=='p1') %>% arrange(NES) %>%
                filter(padj<=0.05),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV1')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p1)

p2 <- (ggplot(df_msig %>% filter(LV=='p2') %>% arrange(NES) %>%
                filter(padj<=0.05),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV2')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p2)

p3 <- (ggplot(df_msig %>% filter(LV=='p3') %>% arrange(NES) %>%
                filter(padj<=0.05),
              aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+ 
         geom_bar(stat = 'identity') +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
         xlab('Normalized Enrichment Score') + ylab('Hallmark')+
         ggtitle('Hallmarks enriched in LV3')+
         theme_minimal(base_family = 'Arial',base_size = 18)+
         theme(text = element_text(family = 'Arial',size=18),
               axis.text.y = element_text(size=18),
               plot.title = element_text(hjust = 0.5),
               legend.key.size = unit(1.5, "lines"),
               legend.position = 'none'))
print(p3)

### Repeat ORA analysis and scatterplots from PROGENY but with Hallmarks------------
library(msigdbr)

# Pull Hallmark gene sets (human, HGNC symbols)
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")

# Convert to named list of gene vectors — same format as progeny_sets before
hallmark_sets <- hallmark_df %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  tibble::deframe()

# Clean names: remove HALLMARK_ prefix for readability
names(hallmark_sets) <- gsub("HALLMARK_", "", names(hallmark_sets))
names(hallmark_sets) <- str_replace_all(names(hallmark_sets), "_", " ")
names(hallmark_sets) <- paste0(
  toupper(substr(names(hallmark_sets), 1, 1)),
  tolower(substr(names(hallmark_sets), 2, nchar(names(hallmark_sets))))
)

# Check
length(hallmark_sets)       # should be 50
hallmark_sets[["Tnfa signaling via nfkb"]][1:5]   # example

vip_scores       <- getVipVn(plsr_model)
background_genes <- names(vip_scores)
# thresholds       <- c(0.8, 1.0, 1.2, 1.5, 2.0)
thresholds <- c(0.7,0.8, 1.0, 1.2, 1.5, 2.0,2.5,3,3.5,4)

threshold_results_hallmark <- lapply(thresholds, function(thr) {
  hits <- names(vip_scores[vip_scores > thr])
  
  if (length(hits) < 10) return(NULL)
  
  res <- run_fisher_ora(
    hits       = hits,
    background = background_genes,
    gene_sets  = hallmark_sets,
    min_overlap = 3
  )
  res$threshold <- thr
  res$n_hits    <- length(hits)
  res
})

sensitivity_hallmark <- bind_rows(Filter(Negate(is.null), threshold_results_hallmark))

# Which pathways are consistently significant?
sensitivity_hallmark <- sensitivity_hallmark %>%
  filter(p_adj < 0.05) %>%
  group_by(pathway) %>%
  summarise(n = n()) %>%
  mutate(perc_sig = n/length(thresholds)) %>%
  arrange(-perc_sig) %>%
  filter(perc_sig>0.5)

## Also infer Hallmark score of in vivo data
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(t(Xh)), column = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- unname(entrez_ids)
inds <- which(!is.na(entrez_ids))
entrez_ids <- entrez_ids[inds]
meas <- as.matrix(t(Xh)[inds,])
rownames(meas) <- entrez_ids
msig_human <- fastenrichment(colnames(meas),
                       entrez_ids,
                       meas,
                       enrichment_space = 'msig_db_h',
                       n_permutations = 10000,
                       order_columns=F)
msig_nes_human <- as.data.frame(msig_human$NES$`NES MSIG Hallmark`) %>% rownames_to_column('Hallmark')
msig_nes_human <- msig_nes_human %>% gather('condition','NES',-Hallmark)
msig_pval_human <- as.data.frame(msig_human$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column('Hallmark')
msig_pval_human <- msig_pval_human %>% gather('condition','padj',-Hallmark)
df_msig_human <- left_join(msig_nes_human,msig_pval_human)
df_msig_human <- df_msig_human %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))
df_msig_human <- df_msig_human %>% mutate(Hallmark=str_replace_all(Hallmark,'_',' '))
df_msig_human <- df_msig_human %>% mutate(Hallmark = tolower(Hallmark)) %>% 
  mutate(Hallmark = paste0(toupper(substr(Hallmark, 1, 1)), tolower(substr(Hallmark, 2, nchar(Hallmark)))))
invivo_plsr <- as.data.frame(Zh_plsr[,c('p1','p2','p3')])
colnames(invivo_plsr) <- c('LV1','LV2','LV3')
invivo_plsr <- cbind(invivo_plsr,as.data.frame(Yh))
invivo_plsr <- invivo_plsr %>% rownames_to_column('condition')
invivo_plsr <- left_join(invivo_plsr,df_msig_human)


pPlsr2 <- ggplot(invivo_plsr %>% filter(Hallmark %in% sensitivity_hallmark$pathway),
                aes(x=LV1,y=LV3,fill=NES)) +
  geom_point(size = 1.5, shape = 21, color = "black")+
  scale_fill_viridis_c()+
  labs(fill = 'Score')+
  xlab('Human LV1') +
  ylab('Human LV3') +
  facet_wrap(~Hallmark) +
  theme_bw()
print(pPlsr2)
ggsave(plot=pPlsr2,
       filename = 'figures/suppl_VIP_Hallmark_enrichment_3.png',
       units = 'cm',
       width = 16,
       height = 12,
       dpi=600)
