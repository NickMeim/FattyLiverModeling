### Functions for making elements of figure 2 - PLSR model
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets

root_dir <- "E:/Jose Luis/Documents/GitHub/FattyLiverModeling"
setwd(root_dir)
source("./utils/plotting_functions.R")

#################################################################################
### Panel - human phenotypes correlate
plt_pheno_cor <- ggplot(data.frame(Yh), aes(x = NAS, y = fibrosis)) +
                  geom_point(size = 2*size_dot, alpha = 0.2, color = "steelblue") +
                  labs(x = "MAFLD score (MAS)", y = "Fibrosis stage")
plt_pheno_cor <- add_theme(plt_pheno_cor)
plt_pheno_cor <- ggMarginal(plt_pheno_cor, type = "histogram", fill = "steelblue")

#################################################################################
### Panel - cross-validation of PLSR model
performance_PLSR <- readRDS("./results/performance_df_human_plsr.rds")
performance_PLSR <- performance_PLSR %>% mutate(type = ifelse(type=='model',set,type))
performance_PLSR <- performance_PLSR %>% filter(metric=='r') %>% select(-metric)
performance_PLSR$type <- factor(performance_PLSR$type ,levels=c('train','test','shuffle Y','shuffle X','random X'))
performance_PLSR <- performance_PLSR %>% mutate(phenotype = ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))

plt_PLSR_CV <- ggplot(performance_PLSR, aes(x = type, y = value, color = type)) +
                geom_boxplot(outlier.alpha = 0, size = size_line*0.7, fill = "white", show.legend = F) +
                geom_jitter(width = 0.1, size = size_dot, show.legend = F) +
                stat_compare_means(comparisons = list(c('train','test'),
                                                      c('test','shuffle Y'),
                                                      c('test','shuffle X'),
                                                      c('test','random X')),
                                   method = 'wilcox',#label = 'p.signif',
                                   tip.length = 0.025,
                                   label.y = c(1.1,0.8, 1.1, 0.95),
                                   size = size_annotation*0.7)+
                facet_wrap(~phenotype, nrow = 2) +
                labs(x = NULL, y = "Pearson correlation") +
                scale_y_continuous(breaks = seq(-0.75,1,0.25),limits = c(NA,1.3)) +
                scale_color_brewer(palette = "Dark2")
plt_PLSR_CV <- add_theme(plt_PLSR_CV)

#################################################################################
### Panel - PLSR performance in training data
plt_PLSR_training <- rbind(data.frame(Measured = Yh[,1], Phenotype = "MAS"), data.frame(Measured = Yh[,2], Phenotype = "Fibrosis stage")) %>%
                      mutate(Predicted = c(predict(plsr_model)[,1], predict(plsr_model)[,2])) %>%
                    ggplot(aes(x = Measured, y = Predicted)) +
                      geom_jitter(size = size_dot*0.8, color = "steelblue", width = 0.1, alpha = 0.8) +
                      geom_abline(slope = 1, intercept = 0, color = "black", linewidth = size_line, linetype = 2) +
                      stat_cor(color = "black", size = size_annotation*0.7) +
                      facet_wrap(~Phenotype, scales = "free") +
                      labs(x = "Measured", y = "Predicted")
plt_PLSR_training <- add_theme(plt_PLSR_training)
  
#################################################################################
### Panel - PLSR scores of human data
plt_PLSR_human <- ggplot(invivo_plsr, aes(x=V1,y=V2,color=normed_score)) +
                  geom_point(size = size_dot)+
                  scale_color_viridis_c()+
                  labs(color = 'Score')+
                  xlab(paste0('Human LV1 (',lv_vars[1],'%)')) + ylab(paste0('Human LV2 (',lv_vars[2],'%)')) +
                  facet_wrap(~phenotype)

plt_PLSR_human <- add_theme(plt_PLSR_human)

### Save panels as figures
ggsave(filename = "./Figures/figure2/plt_pheno_cor.pdf", plot = plt_pheno_cor, units = "cm", width = 5, height = 5)
ggsave(filename = "./Figures/figure2/plt_PLSR_CV.pdf", plot = plt_PLSR_CV, units = "cm", width = 8, height = 8)
ggsave(filename = "./Figures/figure2/plt_PLSR_training.pdf", plot = plt_PLSR_training, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure2/plt_PLSR_human.pdf", plot = plt_PLSR_human, units = "cm", width = 8, height = 5)
