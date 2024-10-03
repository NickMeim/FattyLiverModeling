### Functions for making elements of figure 2 - PLSR model
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets
library(tidyverse)
library(matrixStats)
library(ggExtra)
# root_dir <- "C:/Users/nmeim/Documents/LiverModeling/FattyLiverModeling"
# setwd(root_dir)
source("../utils/plotting_functions.R")
target_dataset <- "Kostrzewski"
ref_dataset <-  "Govaere"
processed_data_list <- readRDS(paste0('../results/processed_data_list_',
                               tolower(ref_dataset),'_',
                               tolower(target_dataset),'.rds'))
Xh <- processed_data_list$Xh
Yh <- processed_data_list$Yh
sex_inferred <- processed_data_list$sex_inferred
Xm <- processed_data_list$Xm
Wm <- processed_data_list$Wm
plsr_model <- readRDS( paste0('../results/PLSR_model_',tolower(ref_dataset),'.rds'))
## Get LV projections
Zh_plsr <- plsr_model@scoreMN
mcor <- cor(cbind(Zh_plsr,Yh))
mcor[upper.tri(mcor,diag = TRUE)] <- 100
mcor <- reshape2::melt(mcor)
mcor <- mcor %>% filter(value != 100)
mcor <- mcor %>% mutate(keep = ifelse(Var1=='fibrosis' | Var2=='fibrosis' | Var1=='NAS' | Var2=='NAS',TRUE,FALSE)) %>%
  filter(keep==TRUE) %>% select(-keep) %>%
  mutate(keep = ifelse((Var1=='fibrosis' & Var2=='NAS') | (Var1=='NAS' & Var2=='fibrosis'),FALSE,TRUE))%>%
  filter(keep==TRUE) %>% select(-keep)
mcor <- mcor %>% group_by(Var2) %>% mutate(avg_abs_cor = mean(abs(value))) %>% select(c('LV'='Var2'),avg_abs_cor) %>%
  unique()
lvs <- mcor$LV[order(-mcor$avg_abs_cor)]
lvs <- as.character(lvs[1:2])
lv_vars <- round(100 * colVars(Zh_plsr)/sum(colVars(Xh)),2)
lv1_var <- round(100 * var(Zh_plsr[,c(lvs[1])])/sum(apply(Xh,2,var)),2)
lv2_var <- round(100 * var(Zh_plsr[,c(lvs[2])])/sum(apply(Xh,2,var)),2)
invivo_plsr <- as.data.frame(Zh_plsr[,c(lvs[1],lvs[2])])
invivo_plsr <- cbind(invivo_plsr,as.data.frame(Yh))
invivo_plsr <- invivo_plsr %>% gather('phenotype','Score',-all_of(lvs))
colnames(invivo_plsr)[1:2] <- c('V1','V2')
invivo_plsr <- invivo_plsr %>%
  mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype)) %>%
  group_by(phenotype) %>%
  mutate(normed_score=Score/max(Score))
#################################################################################
### Panel - human phenotypes correlate
plt_pheno_cor <- ggplot(data.frame(Yh), aes(x = NAS, y = fibrosis)) +
                  geom_point(size = 2*size_dot, shape = 21, fill = "steelblue", stroke = 0.1*size_stroke, color = "black", alpha = 0.2) +
                  labs(x = "MAFLD score (MAS)", y = "Fibrosis stage")
plt_pheno_cor <- add_theme(plt_pheno_cor)
plt_pheno_cor <- ggMarginal(plt_pheno_cor, type = "histogram", fill = "steelblue")

#################################################################################
### Panel - cross-validation of PLSR model
performance_PLSR <- readRDS("../results/performance_df_human_plsr.rds")
performance_PLSR <- performance_PLSR %>% mutate(type = ifelse(type=='model',set,type))
performance_PLSR <- performance_PLSR %>% filter(metric=='r') %>% select(-metric)
performance_PLSR$type <- factor(performance_PLSR$type ,levels=c('train','test','shuffle Y','shuffle X','random X'))
performance_PLSR <- performance_PLSR %>% mutate(phenotype = ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype))

plt_PLSR_CV <- ggplot(performance_PLSR, aes(x = type, y = value, color = type, fill = type)) +
                geom_boxplot(outlier.alpha = 0, size = size_line*0.7, fill = "white", show.legend = F) +
                geom_jitter(width = 0.1, size = size_dot, shape = 21, color = "black", stroke = size_stroke, show.legend = F) +
                stat_compare_means(comparisons = list(c('train','test'),
                                                      c('test','shuffle Y'),
                                                      c('test','shuffle X'),
                                                      c('test','random X')),
                                   method = 'wilcox',#label = 'p.signif',
                                   tip.length = 0.025,
                                   label.y = c(1.1,0.8, 1.1, 0.95),
                                   size = size_annotation*0.5)+
                facet_wrap(~phenotype, nrow = 2) +
                labs(x = NULL, y = "Pearson correlation") +
                scale_y_continuous(breaks = seq(-0.75,1,0.25),limits = c(NA,1.3)) +
                scale_color_brewer(palette = "Dark2") +
                scale_fill_brewer(palette = "Dark2") +
                scale_x_discrete(labels = c("train" = "train\ndata", "test" = "test\ndata", "shuffle Y" = "shuffle\nY",
                                           "shuffle X" = "shuffle\nX", "random X" = "random\nX"))

plt_PLSR_CV <- add_theme(plt_PLSR_CV)

#################################################################################
### Panel - PLSR performance in training data
plt_PLSR_training <- rbind(data.frame(Measured = Yh[,1], Phenotype = "MAS"), data.frame(Measured = Yh[,2], Phenotype = "Fibrosis stage")) %>%
                      mutate(Predicted = c(predict(plsr_model)[,1], predict(plsr_model)[,2])) %>%
                    ggplot(aes(x = Measured, y = Predicted)) +
                      geom_jitter(size = size_dot*0.8, shape = 21, fill = "steelblue", color = "black", stroke = size_stroke, width = 0.1) +
                      geom_abline(slope = 1, intercept = 0, color = "black", linewidth = size_line, linetype = 2) +
                      stat_cor(color = "black", size = size_annotation*0.7) +
                      facet_wrap(~Phenotype, scales = "free") +
                      labs(x = "Measured", y = "Predicted")
plt_PLSR_training <- add_theme(plt_PLSR_training)

#################################################################################
### Panel - PLSR scores of human data
plt_PLSR_human <- ggplot(invivo_plsr, aes(x=V1,y=V2,fill=normed_score)) +
                  geom_point(size = size_dot, shape = 21, color = "black", stroke = size_stroke)+
                  scale_fill_viridis_c()+
                  labs(fill = 'Score')+
                  xlab(paste0('Human LV',substr(lvs[1],2,2),' (',lv1_var,'%)')) +
  ylab(paste0('Human LV',substr(lvs[2],2,2),' (',lv2_var,'%)')) +
                  facet_wrap(~phenotype)

plt_PLSR_human <- add_theme(plt_PLSR_human)
print(plt_PLSR_human)
### Save panels as figures
ggsave(filename = "figure2/plt_pheno_cor.pdf", plot = plt_pheno_cor, units = "cm", width = 6, height = 5)
ggsave(filename = "figure2/plt_PLSR_CV.pdf", plot = plt_PLSR_CV, units = "cm", width = 6, height = 6)
ggsave(filename = "figure2/plt_PLSR_training.pdf", plot = plt_PLSR_training, units = "cm", width = 6.5, height = 5)
ggsave(filename = "figure2/plt_PLSR_human.pdf", plot = plt_PLSR_human, units = "cm", width = 7.5, height = 5)
# Save short versions of PLSR for figure 3
ggsave(filename = "figure3/plt_PLSR_training_short.pdf", plot = plt_PLSR_training, units = "cm", width = 6.5, height = 4)
ggsave(filename = "figure3/plt_PLSR_human_short.pdf", plot = plt_PLSR_human, units = "cm", width = 7.5, height = 4)
