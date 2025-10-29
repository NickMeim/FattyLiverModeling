### Functions for making elements of figure 4 - Additional latent variables
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets
library(tidyverse)
library(matrixStats)
library(ggExtra)
root_dir <- "C:/Users/nmeim/Documents/LiverModeling/FattyLiverModeling"
setwd(root_dir)
source("utils/plotting_functions.R")
target_dataset <- "Kostrzewski"
ref_dataset <-  "Govaere"
dataset_names <- c("Govaere", "Kostrzewski", "Wang", "Feaver")
processed_data_list <- readRDS(paste0('results/processed_data_list_',
                                      tolower(ref_dataset),'_',
                                      tolower(target_dataset),'.rds'))
Xh <- processed_data_list$Xh
Yh <- processed_data_list$Yh
sex_inferred <- processed_data_list$sex_inferred
Xm <- processed_data_list$Xm
#################################################################################
### Panel - Human data on LV extras
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
Zh <- as.data.frame(Xh %*% Wm_opt)
Zh$NAS <- Yh[,1]
Zh$fibrosis <- Yh[,2]
Zh <- Zh %>% rownames_to_column('sample')
combined_plot_nas <- scatter_box_plot_V2(Zh %>% mutate(pheno=NAS),'MAS',
                                      font_size = size_text, point_size = size_dot, point_stroke = size_stroke)
combined_plot_fibrosis <- scatter_box_plot_V2(Zh %>% mutate(pheno=fibrosis),'Fibrosis stage',
                                           font_size = size_text, point_size = size_dot, point_stroke = size_stroke)

combined_plot_nas <- add_theme(combined_plot_nas)
combined_plot_fibrosis <- add_theme(combined_plot_fibrosis) # Fix issue with PDF export. Maybe use cowplot for panels

#################################################################################
### Panel - Cross validation of extended performance with LV extra
performance_all_plot <- readRDS('results/performanceall_plot_spearman.rds')

plt_p_test <- performance_all_plot %>% filter(set=='test' & approach!='Wopt') %>%
  ggplot(aes(x = approach, y = rho, color = approach, fill = approach)) +
  geom_boxplot(outlier.alpha = 0, size = size_line*0.7, fill = "white", show.legend = F) +
  geom_jitter(width = 0.1, size = size_dot, shape = 21, color = "black", stroke = size_stroke, show.legend = F) +
  stat_compare_means(comparisons = list(c('human genes','optimized MPS'),
                                        c('human genes','truncated'),
                                        c('optimized MPS','truncated'),
                                        c('optimized MPS','shuffle X')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.025,
                     label.y = c(1.1,0.9,0.9,0.9),
                     size = size_annotation*0.75)+
  facet_wrap(~phenotype, nrow = 2) +
  labs(x = NULL, y = "Spearman`s rank correlation") +
  scale_y_continuous(breaks = seq(-0.75,1,0.25),limits = c(NA,1.3)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("human genes" = "Human genes", "truncated" = "Truncated data", "optimized MPS" = "Optimized MPS",
                              "shuffle X" = "shuffle X"))

plt_p_test <- add_theme(plt_p_test) +
                theme(legend.position = 'none')

#################################################################################
### Panel - Improvement in performance with each extra LV
all_performance_res <- readRDS('results/all_performance_res_spearman.rds')
mu_plsr_train <- mean(all_performance_res$rho[all_performance_res$set=='train' & all_performance_res$model=='PLSR'])
sd_plsr_train<- sd(all_performance_res$rho[all_performance_res$set=='train'& all_performance_res$model=='PLSR'])
mu_plsr_test<- mean(all_performance_res$rho[all_performance_res$set=='test'& all_performance_res$model=='PLSR'])
sd_plsr_test<- sd(all_performance_res$rho[all_performance_res$set=='test'& all_performance_res$model=='PLSR'])
test_r_shuffled<- readRDS('results/test_r_shuffled_cumulative_lvs_spearman.rds')
num_folds <- length(unique(all_performance_res$fold))
plt_required_extra_basis <-
  ggplot(all_performance_res %>% filter(model!='PLSR') %>% select(model,set,mu,std),
         aes(x=model,y=mu,color = set , group =set))+
  geom_point(size=size_dot)+
  geom_line(lwd = size_line)+
  geom_errorbar(aes(ymax = mu + std/sqrt(num_folds),ymin = mu - std/sqrt(num_folds)),
                width = 0.05,size=0.75)+
  ylim(NA,1) +
  ylab('Spearman`s rank correlation')+
  ### train performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mu_plsr_train - sd_plsr_train/sqrt(num_folds),
                  ymax = mu_plsr_train + sd_plsr_train/sqrt(num_folds)),
              fill = "#01B8BB", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mu_plsr_train,
           yend = mu_plsr_train,
           color = "#01B8BB", linewidth = size_line,linetype='dashed') +
  annotate('text',x=1,
           y=mu_plsr_train + 0.03,
           label="human PLSR train performance",
           hjust = 0 ,
           size=0.75*size_annotation)+
  ### test performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mu_plsr_test - sd_plsr_test/sqrt(num_folds),
                  ymax = mu_plsr_test + sd_plsr_test/sqrt(num_folds)),
              fill = "#E0766D", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mu_plsr_test,
           yend = mu_plsr_test,
           color = "#E0766D", linewidth = size_line,linetype='dashed') +
  annotate('text',x=1,
           y=mu_plsr_test + 0.03,
           label="human PLSR test performance",
           hjust = 0 ,
           size=0.75*size_annotation)+
  ### random performance shaded area
  geom_ribbon(inherit.aes = FALSE,
              xmin=1,xmax=3,
              aes(x=seq(1,3,length.out=nrow(all_performance_res %>% filter(model!='PLSR'))),
                  ymin = mean(test_r_shuffled) - sd(test_r_shuffled)/sqrt(num_folds),
                  ymax = mean(test_r_shuffled) + sd(test_r_shuffled)/sqrt(num_folds)),
              fill = "#F564E3", alpha = 0.25) +  # Shaded area
  annotate("segment", x = 1, xend = 3, y = mean(test_r_shuffled),
           yend = mean(test_r_shuffled),
           color = "#F564E3", linewidth = size_line,linetype='dashed') +
  annotate('text',x=1,
           y=mean(test_r_shuffled) + 0.03,
           label="shuffled model performance",
           hjust = 0 ,
           size=0.75*size_annotation)+
  geom_hline(yintercept = 0,linetype = 'solid',color='black',lwd=0.75) +
  scale_x_discrete(labels=c('translatable\n PCs',
                            'PCs +\n extra LV1',
                            'PCs +\n extra LV1 and LV2'))



plt_required_extra_basis <- add_theme(plt_required_extra_basis) +
                            theme(axis.title.x = element_blank(),
                                  axis.line.x = element_blank(),
                                  axis.line.y = element_line(linewidth = 0.75),
                                  panel.grid.major = element_line(),
                                  legend.position = "top")

#################################################################################
### Panel - LV extras are somewhat generalizable
results_external_data <- readRDS("results/external_clinical_differences_results_of_extra_vector_spearman.rds") %>%
                          pivot_longer(cols = 1:3, names_to = "model", values_to = "Spearman")

results_external_data$fold <- sapply(results_external_data$fold_ids, FUN = function(x){strsplit(x, ".", fixed = T) %>% unlist()})[1,]

results_external_data <- results_external_data %>%
                          group_by(dataset, set, model, fold) %>%
                          summarize(Spearman = mean(Spearman)) %>%
                          mutate(dataset_GEO = ifelse(dataset == "Hoang", "GSE13090", "GSE162694"),
                                 set = factor(set, levels = c("train", "test"))) # Hoang et al or Pantano et al

# Plot
plt_external_data <- results_external_data %>%
  ggplot(aes(x = model, y = Spearman, color = model, fill = model)) +
  geom_boxplot(outlier.alpha = 0, size = size_line*0.7, fill = "white", show.legend = F) +
  geom_jitter(width = 0.1, size = size_dot, shape = 21, color = "black", stroke = size_stroke, show.legend = F) +
  stat_compare_means(comparisons = list(c('extra.basis','PC'),
                                        c('extra.basis','PLSR'),
                                        c('PC','PLSR')),
                     method = 'wilcox',#label = 'p.signif',
                     tip.length = 0.025,
                     #label.y = c(0.9,1.1, 0.9),
                     size = size_annotation*0.75)+
  facet_grid(rows = vars(set), cols = vars(dataset_GEO)) +
  labs(x = NULL, y = "Spearman`s rank correlation") +
  scale_y_continuous(breaks = seq(0.25,1,0.25),limits = c(NA,1.2)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(limits = c("PLSR", "PC", "extra.basis"),
                   labels = c("extra.basis" = "Optimized\nMPS", "PC" = "Truncated\ndata", "PLSR" = "Human\ngenes"))

plt_external_data <- add_theme(plt_external_data)



#################################################################################
### Panel - Pathway activity of TCs

### Save panels as figures
ggsave(filename = "figures/figure4/plt_combined_NAS.pdf", plot = combined_plot_nas, units = "cm", width = 7, height = 5.5)
ggsave(filename = "figures/figure4/plt_combined_Fib.pdf", plot = combined_plot_fibrosis, units = "cm", width = 7, height = 5.5)
ggsave(filename = "figures/figure4/plt_required_extra_basis.pdf", plot = plt_required_extra_basis, units = "cm", width = 6, height = 6)
ggsave(filename = "figures/figure4/plt_p_test.pdf", plot = plt_p_test, units = "cm", width = 8, height = 9)
ggsave(filename = "figures/figure4/plt_external_data.pdf", plot = plt_external_data, units = "cm", width = 9, height = 6)


