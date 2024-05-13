# load in packages
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(ggbreak) 
library(patchwork)
library(ggradar)
library(pls)
library(pheatmap)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(hgu133a.db)
library(rstatix)
library(fgsea)
library(topGO)
library(GO.db)
library(OmnipathR)
library(EGSEAdata)
source("../utils/plotting_functions.R")
source("functions_translation_jose.R")
source("CrossValidationUtilFunctions.R")
source('enrichment_calculations.R')

plot_gene_loadings <- function(loadings,selection,y_lab,plotting=TRUE){
  # n_lvs <- ncol(loadings)
  loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,selection)
  loadings$significant <- ''
  loadings$significant[order(-abs(loadings[,selection]))[1:10]] <- loadings$gene[order(-abs(loadings[,selection]))[1:10]]
  loadings <- loadings[order(loadings[,selection]),]
  loadings$gene <- factor(loadings$gene,levels = loadings$gene)
  colnames(loadings)[2] <- 'LV'
  p <- ggplot(loadings,aes(x=gene,y=LV,color = significant)) + geom_point() +
    geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
    xlab('genes') + ylab(y_lab)+
    scale_x_discrete(expand = c(0.1, 0.1))+
    theme_pubr(base_family = 'Arial',base_size = 20)+
    theme(text = element_text(family = 'Arial',size=20),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  if (plotting==TRUE){
    print(p)
  }
  return(p)
}

pathway_activity_interpretation <- function(W,W_PCspace,plotting=TRUE){
  net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
  extra_basis_paths <- decoupleR::run_viper(W, net_prog,minsize = 1,verbose = TRUE) %>% select(-statistic)
  PC_space_paths <- decoupleR::run_viper(W_PCspace, net_prog,minsize = 1,verbose = TRUE)
  PC_space_paths <- PC_space_paths %>% select(source,c('pc_activity'='score'))
  PC_space_paths <- PC_space_paths$pc_activity
  extra_basis_paths <- extra_basis_paths %>% select(-p_value)
  colnames(extra_basis_paths)[1] <- 'Pathway'
  
  extra_basis_paths <- extra_basis_paths %>% group_by(condition,Pathway) %>%
    mutate(p_value=sum(abs(PC_space_paths)>=abs(score))/length(PC_space_paths))
  extra_basis_paths  <- extra_basis_paths %>%
    mutate(p.adj = p_value*length(unique(extra_basis_paths$Pathway))) %>%
    mutate(p.adj = ifelse(p.adj>1,1,p.adj))
  
  p <- (ggplot(extra_basis_paths %>% select(c('activity'='score'),Pathway,p_value,condition) %>% 
            filter (condition=='V1'),
          aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity',color='black') +
      scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
      scale_x_continuous(n.breaks = 8,limits = c(-8,8))+
      ggtitle('Extra LV1')+
      ylab('Pathway')+
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
      theme(text = element_text(size = 20,family = 'Arial'),
            legend.position = 'right',
            plot.title = element_text(hjust = 0.5))) +
    (ggplot(extra_basis_paths %>% select(c('activity'='score'),Pathway,p_value,condition) %>% 
              filter (condition=='V2'),
            aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity',color='black') +
       ggtitle('Extra LV2')+
       ylab('Pathway')+
       scale_fill_gradient2(low='blue',high = 'red',mid = 'white',midpoint = 0)+
       scale_x_continuous(n.breaks = 8,limits = c(-8,8))+
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
       theme_pubr(base_size = 16,base_family = 'Arial')+
       theme(text = element_text(size = 16,family = 'Arial'),
             legend.position = 'right',
             plot.title = element_text(hjust = 0.5),
             axis.title.y = element_blank())) 
  if (plotting==TRUE){
    print(p)
  }
  return(list(figure=p,extra_basis_paths))
}
