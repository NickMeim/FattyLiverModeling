### Functions for making elements of figure 6 - General optimization of variance
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets

root_dir <- "E:/Jose Luis/Documents/GitHub/FattyLiverModeling"
setwd(root_dir)
source("./utils/plotting_functions.R")

#################################################################################
### Panel - Human data on LV extras

#################################################################################
### Panel - Cross validation of extended performance with LV extra

#################################################################################
### Panel - Improvement in performance with each extra LV

#################################################################################
### Panel - LV extras are somewhat generalizable



#################################################################################
### Panel - PLSR labelled by sex
plt_PLSR_sex <- ggplot(cbind(as.data.frame(Zh_plsr[,c(lvs[1],lvs[2])]), 
                             as.data.frame(sex_inferred) %>% select(c('sex'='V1'))),
                       aes(x=p1,y=p2,color=sex)) +
                geom_point(size = size_dot)+
                scale_color_viridis_d()+
                stat_ellipse(size = size_line, show.legend = F, linetype = 2) +
                labs(color = 'Sex')+
                xlab(paste0('Human LV1 (',lv_vars[1],'%)')) + ylab(paste0('Human LV2 (',lv_vars[2],'%)'))
plt_PLSR_sex <- add_theme(plt_PLSR_sex)



### Save panels as figures
ggsave(filename = "./Figures/figure3/plt_PCA_MPS.pdf", plot = plt_PCA_MPS, units = "cm", width = 5, height = 5)
ggsave(filename = "./Figures/figure6/plt_PLSR_sex.pdf", plot = plt_PLSR_sex, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure3/plt_PLSR_training_back.pdf", plot = plt_PLSR_training_back, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure3/plt_TC_prediction.pdf", plot = plt_TC_prediction, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure3/plt_TC_MPS.pdf", plot = plt_TC_MPS, units = "cm", width = 5, height = 6)


