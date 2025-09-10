# load in packages
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(patchwork)

plot_gene_loadings <- function(loadings,selection,y_lab,plotting=TRUE,top=10){
  # n_lvs <- ncol(loadings)
  library(tidyverse)
  library(ggfortify)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(ggrepel)
  library(patchwork)
  loadings <- as.data.frame(loadings) %>% rownames_to_column('gene') %>% dplyr::select(gene,selection)
  loadings$significant <- ''
  loadings$significant[order(-abs(loadings[,selection]))[1:top]] <- loadings$gene[order(-abs(loadings[,selection]))[1:top]]
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
  library(tidyverse)
  library(ggfortify)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(ggrepel)
  library(patchwork)
  colnames(W) <- paste0('V',seq(1,ncol(W)))
  net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
  extra_basis_paths <- decoupleR::run_viper(W, net_prog,minsize = 1,verbose = TRUE) %>% dplyr::select(-statistic)
  PC_space_paths <- decoupleR::run_viper(W_PCspace, net_prog,minsize = 1,verbose = TRUE)
  PC_space_paths <- PC_space_paths %>% dplyr::select(source,c('pc_activity'='score'))
  PC_space_paths <- PC_space_paths$pc_activity
  extra_basis_paths <- extra_basis_paths %>% dplyr::select(-p_value)
  colnames(extra_basis_paths)[1] <- 'Pathway'

  extra_basis_paths <- extra_basis_paths %>% group_by(condition,Pathway) %>%
    mutate(p_value=sum(abs(PC_space_paths)>=abs(score))/length(PC_space_paths))

  plots <- NULL
  for (i in 1:ncol(W)){
    p <- (ggplot(extra_basis_paths %>% dplyr::select(c('activity'='score'),Pathway,p_value,condition) %>%
              filter (condition==paste0('V',i)),
            aes(x=activity,y=reorder(Pathway,activity),fill=activity)) + geom_bar(stat = 'identity') +
       scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                            midpoint = 0,limits = c(-lim,lim))+
       scale_x_continuous(n.breaks = 8,limits = c(-lim,lim))+
       ggtitle(paste0('LV extra ',i))+
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
             plot.title = element_text(hjust = 0.5)))
    if (plotting==TRUE){
      print(p)
    }
    plots[[i]] <- p
  }
  return(list(plots=plots,extra_basis_paths))
}
