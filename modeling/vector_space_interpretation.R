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
dir = getwd()
if (rev.Vector(strsplit(dir,split = '/')[[1]])[1] == "FattyLiverModeling"){
  source("utils/plotting_functions.R")
  source("modeling/functions_translation.R")
  source("modeling/CrossValidationUtilFunctions.R")
  source('modeling/enrichment_calculations.R')
}else{
  source("../utils/plotting_functions.R")
  source("functions_translation.R")
  source("CrossValidationUtilFunctions.R")
  source('enrichment_calculations.R')
}

plot_gene_loadings <- function(loadings,selection,y_lab,plotting=TRUE){
  # n_lvs <- ncol(loadings)
  loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% select(gene,selection)
  loadings$significant <- ''
  loadings$significant[order(-abs(loadings[,selection]))[1:10]] <- loadings$gene[order(-abs(loadings[,selection]))[1:10]]
  loadings <- loadings[order(loadings[,selection]),]
  loadings$gene <- factor(loadings$gene,levels = loadings$gene)
  colnames(loadings)[2] <- 'LV'
  p <- ggplot(loadings,aes(x=as.numeric(gene),y=LV,color = significant)) + geom_point() +
    geom_text_repel(aes(label=significant),size=6,max.overlaps=40)+
    xlab('Rank') + ylab(y_lab)+
    theme_pubr(base_family = 'Arial',base_size = 20)+
    theme(text = element_text(family = 'Arial',size=20),
          legend.position = 'none')
  if (plotting==TRUE){
    print(p)
  }
  return(p)
}

pathway_activity_interpretation <- function(W,W_PCspace,plotting=TRUE,lim=8){
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
          aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity') +
      scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                           midpoint = 0,limits = c(-lim,lim))+
      scale_x_continuous(n.breaks = 8,limits = c(-lim,lim))+
      ggtitle('LV extra 1')+
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
        theme_minimal(base_size = 24,base_family = 'Arial')+
      theme(text = element_text(size = 24,family = 'Arial'),
            legend.position = 'right',
            plot.title = element_text(hjust = 0.5))) +
    (ggplot(extra_basis_paths %>% select(c('activity'='score'),Pathway,p_value,condition) %>% 
              filter (condition=='V2'),
            aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity') +
       ggtitle('LV extra 2')+
       ylab('Pathway')+
       scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                            midpoint = 0,limits = c(-lim,lim))+
       scale_x_continuous(n.breaks = 8,limits = c(-lim,lim))+
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
       theme_minimal(base_size = 24,base_family = 'Arial')+
       theme(text = element_text(size = 24,family = 'Arial'),
             legend.position = 'right',
             plot.title = element_text(hjust = 0.5),
             axis.title.y = element_blank())) 
  if (plotting==TRUE){
    print(p)
  }
  return(list(figure=p,extra_basis_paths))
}

TF_activity_interpretation <- function(W,W_PCspace,regulon,plotting=TRUE,lim=6){
  extra_basis_tfs <- decoupleR::run_viper(W, regulon,minsize = 5,verbose = TRUE) %>% select(-statistic)
  PC_space_tfs <- decoupleR::run_viper(W_PCspace, regulon,minsize = 5,verbose = TRUE)
  PC_space_tfs <- PC_space_tfs %>% select(source,c('pc_activity'='score'))
  PC_space_tfs <- PC_space_tfs$pc_activity
  extra_basis_tfs <- extra_basis_tfs %>% select(-p_value)
  colnames(extra_basis_tfs)[1] <- 'TF'
  
  extra_basis_tfs <- extra_basis_tfs %>% group_by(condition,TF) %>%
    mutate(p_value=sum(abs(PC_space_tfs)>=abs(score))/length(PC_space_tfs))
  extra_basis_tfs  <- extra_basis_tfs %>%
    mutate(p.adj = p_value*length(unique(extra_basis_tfs$TF))) %>%
    mutate(p.adj = ifelse(p.adj>1,1,p.adj))
  extra_basis_tfs$significant <- NA  
  extra_basis_tfs <- extra_basis_tfs %>% mutate(significant = ifelse(p_value<0.1,
                                                                     TF,
                                                                     significant))
  
  p <- (ggplot(extra_basis_tfs %>% select(c('activity'='score'),TF,p_value,condition,significant) %>% 
                 filter (condition=='V1'),
               aes(x=as.numeric(reorder(TF,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
          geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
          scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.01,1))+
          scale_y_continuous(n.breaks = 6,limits = c(-lim,lim))+
          ggtitle('LV extra 1')+
          xlab('Rank')+
          theme_pubr(base_size = 20,base_family = 'Arial')+
          theme(text = element_text(size = 20,family = 'Arial'),
                legend.position = 'right',
                plot.title = element_text(hjust = 0.5))) +
    (ggplot(extra_basis_tfs %>% select(c('activity'='score'),TF,p_value,condition,significant) %>% 
              filter (condition=='V2'),
            aes(x=as.numeric(reorder(TF,activity)),y=activity,fill=p_value)) + geom_point(shape=21,size=2) +
       geom_text_repel(aes(label=significant),size=5,max.overlaps=60,box.padding = 0.7)+
       scale_fill_gradient(low='red',high = 'white',trans = 'log',breaks = c(0.01,0.05,0.1,0.5),limits = c(0.01,1))+
       scale_y_continuous(n.breaks = 6,limits = c(-lim,lim))+
       ggtitle('LV extra 2')+
       xlab('Rank')+
       theme_pubr(base_size = 20,base_family = 'Arial')+
       theme(text = element_text(size = 20,family = 'Arial'),
             legend.position = 'right',
             plot.title = element_text(hjust = 0.5),
             axis.title.y = element_blank())) 
  if (plotting==TRUE){
    print(p)
  }
  return(list(figure=p,extra_basis_tfs))
}

perturnation_activity_inference <- function(W,metadata_human,regulon,perturbation_regulon,plotting=TRUE){
  extra_basis_tfs <- decoupleR::run_viper(W, regulon,minsize = 5,verbose = FALSE) %>% select(-statistic)
  extra_basis_tfs <- extra_basis_tfs %>% select(-p_value)
  colnames(extra_basis_tfs)[1] <- 'TF'
  extra_basis_tfs <- extra_basis_tfs %>% spread('condition','score') %>% column_to_rownames('TF')
  pert_activity <- decoupleR::run_viper(extra_basis_tfs, 
                                        perturbation_regulon %>% unique(),
                                        minsize = 1,
                                        verbose = TRUE)
  colnames(pert_activity)[2] <- 'Response_ID'
  
  pert_activity <- left_join(pert_activity,metadata_human)
  pert_activity <- pert_activity %>% select(Chemical_Compound,c('activity'='score'),Duration,Concentration,c('p-value'='p_value'),condition)
  pert_activity <- pert_activity %>% mutate(perturbation = ifelse(!is.na(Concentration),
                                                                  paste0(toupper(Chemical_Compound),':',Duration,':',Concentration),
                                                                  paste0(toupper(Chemical_Compound),':',Duration)))
  pert_activity <- pert_activity %>% group_by(perturbation,condition) %>% mutate(min_pvalue=min(`p-value`)) %>%
    ungroup() %>% filter(`p-value`==min_pvalue) %>% select(-min_pvalue) %>% unique()
  pert_activity$significant <- ''
  pert_activity <- pert_activity %>% mutate(significant = ifelse(`p-value`<=0.05 & abs(activity)>1.5,
                                                                 perturbation,
                                                                 significant))
  pert_activity <- pert_activity %>% mutate(significant = ifelse(Chemical_Compound=='Ifna',
                                                                 perturbation,
                                                                 significant))
  
  p <- (ggplot(pert_activity %>% filter(condition=='V1'),
               aes(x=as.numeric(reorder(perturbation,activity)),y=activity)) +  # fill = `p-value`
          geom_point(size=2,color = '#CD5C5C') + #shape=21
          # geom_text_repel(aes(label=significant,color=`p-value`),size=8,max.overlaps=60,box.padding = 1.25)+
          geom_text_repel(aes(label=significant),size=8,max.overlaps=60,box.padding = 1.5)+
          xlab('Rank') + ylab('Inferred perturbation activity')+
          # scale_fill_gradient(trans='log10',low = "red",high = "white",
          #                     limits = c(min(pert_activity$`p-value`),1),
          #                     breaks = c(0.01,0.05,0.1,0.5,1)) +
          # scale_color_gradient(trans='log10',low = "red",high = "#f96464",
          #                      limits = c(min(pert_activity$`p-value`),0.05)) +
          # guides(color = "none")+
          ggtitle('LV extra 1')+
          theme_pubr(base_family = 'Arial',base_size = 24)+
          theme(text = element_text(family = 'Arial',size=24),
                panel.grid.major = element_line(),
                # axis.ticks.x = element_blank(),
                # axis.text.x = element_blank(),
                # legend.position = 'right',
                # legend.text = element_text(size=14),
                plot.title = element_text(hjust=0.5))) + 
    (ggplot(pert_activity %>% filter(condition=='V2'),
            aes(x=as.numeric(reorder(perturbation,activity)),y=activity)) + 
       geom_point(size=2,color = '#CD5C5C')  +
       geom_text_repel(aes(label=significant),size=8,max.overlaps=60,box.padding = 1.25)+
       xlab('Rank') + ylab('Inferred perturbation activity')+
       # scale_fill_gradient(trans='log10',low = "red",high = "white",
       #                     limits = c(min(pert_activity$`p-value`),1),
       #                     breaks = c(0.01,0.05,0.1,0.5,1)) +
       # scale_color_gradient(trans='log10',low = "red",high = "#f96464",
       #                      limits = c(min(pert_activity$`p-value`),0.05)) +
       # guides(color = "none")+
       ggtitle('LV extra 2')+
       theme_pubr(base_family = 'Arial',base_size = 24)+
       theme(text = element_text(family = 'Arial',size=24),
             panel.grid.major =  element_line(),
             # axis.ticks.x = element_blank(),
             # axis.title.y = element_blank(),
             # axis.text.x = element_blank(),
             # legend.position = 'right',
             # legend.text = element_text(size=14),
             plot.title = element_text(hjust=0.5)))
  if (plotting==TRUE){
    print(p)
  }
  return(list(figure=p,pert_activity))
}
