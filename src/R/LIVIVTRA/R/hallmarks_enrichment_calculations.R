HallmarksFastenrichment <- function(signature_ids,
                           gene_ids,
                           measurements,
                           axis_lab,
                           pval_adjustment=T,
                           n_permutations=1000){

  # Load packages
  library(tidyverse)
  library(fgsea)
  library(gage)
  library(EGSEAdata)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggfortify)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(ggrepel)
  library(patchwork)

  sig_ids <- as.character(signature_ids)

  if (length(sig_ids)<2){
    measurements <- as.matrix(measurements)
  }
  colnames(measurements) <- paste0('V',seq(1,ncol(measurements)))
  # load pathway data
  egsea.data(species = "human",returnInfo = TRUE)

  all_genesets <- list()
  # Build the list of pathways, genesets to run gsea
  print("building enrichment space")

  msig_overview <- t(as.data.frame(msigdb[[1]]))
  for (i in 2:(length(msigdb)-1)) {
    msig_overview <- rbind(msig_overview,t(as.data.frame(msigdb[[i]])))
  }
  msig_overview <- as.data.frame(msig_overview)
  rownames(msig_overview) <- msig_overview$STANDARD_NAME

  hallmark <- list()
  h <- 1
  for (i in 1:(length(msigdb)-1)) {
    a <- msigdb[[i]]
    if (a["CATEGORY_CODE"] == "H") {
      nm <- a["STANDARD_NAME"]
      hallmark <- c(hallmark,str_split(a["MEMBERS_EZID"],pattern = ","))
      names(hallmark)[h] <- nm
      h <- h+1
    }
  }
  names(hallmark) <- paste0("FL1000_MSIG_H_",names(hallmark))
  all_genesets <- c(all_genesets,hallmark)

  print("running fgsea for enrichment space")
  genesets_list <-apply(measurements,MARGIN = 2,fgsea,pathways = all_genesets,
                  minSize=10,
                  maxSize=500,
                  nperm = n_permutations)
  print("fgsea finished, preparing outputs")

  # Prepare output

  NES <- genesets_list[[1]]$NES
  if (pval_adjustment==T){
    pval <- genesets_list[[1]]$padj
  }else{
    pval <- genesets_list[[1]]$pval
  }
  if (length(genesets_list)>1){
    for (i in 2:length(genesets_list)) {

      NES <- cbind(NES,genesets_list[[i]]$NES)
      if (pval_adjustment==T){
        pval <- cbind(pval,genesets_list[[i]]$padj)
      }else{
        pval <- cbind(pval,genesets_list[[i]]$pval)
      }
    }
  }else{
    NES <- as.matrix(NES)
    pval <- as.matrix(pval)
  }

  colnames(NES) <- names(genesets_list)
  rownames(NES) <- genesets_list[[1]]$pathway
  colnames(pval) <- names(genesets_list)
  rownames(pval) <- genesets_list[[1]]$pathway

  NES_h <- NES[sapply(rownames(NES),FUN=grepl,pattern="FL1000_MSIG_H_"),]
  pval_h <- pval[sapply(rownames(pval),FUN=grepl,pattern="FL1000_MSIG_H_"),]

  msig_nes <- as.data.frame(NES_h) %>% rownames_to_column('Hallmark')
  msig_nes <- msig_nes %>% gather('LV','NES',-Hallmark)
  msig_pval <- as.data.frame(pval_h) %>% rownames_to_column('Hallmark')
  msig_pval <- msig_pval %>% gather('LV','padj',-Hallmark)
  df_msig <- left_join(msig_nes,msig_pval)
  df_msig <- df_msig %>% mutate(Hallmark=substr(Hallmark, nchar('FL1000_MSIG_H_HALLMARK_')+1, nchar(Hallmark)))
  df_msig <- df_msig %>% mutate(Hallmark=str_replace_all(Hallmark,'_',' '))
  df_msig <- df_msig %>% mutate(Hallmark = tolower(Hallmark)) %>%
    mutate(Hallmark = paste0(toupper(substr(Hallmark, 1, 1)), tolower(substr(Hallmark, 2, nchar(Hallmark)))))

  plot_list <- NULL
  for (i in 1:ncol(measurements)){
    p <- (ggplot(df_msig %>% filter(LV==paste0('V',i)) %>% arrange(NES) %>%
                    filter(padj<=0.1),
                  aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+
             geom_bar(stat = 'identity') +
             scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0)+
             xlab('Normalized Enrichment Score') + ylab('Hallmark')+
             ggtitle(paste0('Hallmarks enriched in ',axis_label ,' ',i))+
             theme_minimal(base_family = 'Arial',base_size = 18)+
             theme(text = element_text(family = 'Arial',size=18),
                   axis.text.y = element_text(size=18),
                   plot.title = element_text(hjust = 0.5),
                   legend.key.size = unit(1.5, "lines"),
                   legend.position = 'none'))
    print(p)
    plot_list[[i]] <- p
  }
  result <- list(enrich_results = df_msig, plots = plot_list)
  return(result)
}
