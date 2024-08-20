### Functions for making elements of figure 3 - MPS and backprojection
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets

root_dir <- "E:/Jose Luis/Documents/GitHub/FattyLiverModeling"
setwd(root_dir)
source("./utils/plotting_functions.R")
target_dataset <- "Kostrzewski"

#################################################################################
### Load results
Wm_tot <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_TC <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_combo.rds'))
Wh <- readRDS(file = paste0('results/Wh_',tolower(target_dataset),'.rds'))
plsr_model <- readRDS(file = paste0('results/PLSR_model_',tolower(target_dataset),'.rds'))
Bh <- t(plsr_model@weightMN) %*% plsr_model@coefficientMN

#################################################################################
### Panel - MPS PCA plot
# Project MPS data onto reduced PCA space defined based on averaged experimental treatments
Zm <- Xm %*% Wm
per_var <- round(100 * colVars(Zm)/sum(colVars(Xm)), 2)
plt_PCA_MPS <- ggplot(data.frame(Zm), aes(x = PC1,y = PC2)) +
                geom_point(fill = '#FD0714',size = size_dot, shape = 21, stroke = size_stroke, color = "black")+
                xlab(paste0('MPS PC1 (',per_var[1],'%)')) + ylab(paste0('MPS PC2 (',per_var[2],'%)')) 
plt_PCA_MPS <- add_theme(plt_PCA_MPS)

#################################################################################
### Panel - PLSR scores of backprojected human data
# Extract the necessary matrices from the opls object
P <- plsr_model@loadingMN  # Loadings matrix (P)
# plsr_scores <- plsr_model@scoreMN    # Scores matrix (T)
W <- plsr_model@weightMN   # Weights matrix (W)
# C <- plsr_model@loadingYMN # Y-loadings matrix (C)
# E <- plsr_model@residualsMN # Residuals matrix (E)
# U <- plsr_model@orthoScoreMN # Orthogonal scores matrix (U)
# Manually calculate the scores (T)
# T = X * W * (P' * W)^-1
Zh_plsr_backprj <- Xh %*% Wm %*% t(Wm)  %*% Wh %*% solve(t(P) %*% W)
invivo_plsr_backprj <- as.data.frame(Zh_plsr_backprj[,1:2])
invivo_plsr_backprj <- cbind(invivo_plsr_backprj,as.data.frame(Yh))
invivo_plsr_backprj <- invivo_plsr_backprj %>% gather('phenotype','Score',-all_of(lvs))
colnames(invivo_plsr_backprj)[1:2] <- c('V1','V2')
invivo_plsr_backprj <- invivo_plsr_backprj %>%
                        mutate(phenotype=ifelse(phenotype=='fibrosis','Fibrosis stage',phenotype)) %>%
                        group_by(phenotype) %>% 
                        mutate(normed_score=Score/max(Score))

plt_PLSR_backproject <- ggplot(invivo_plsr_backprj,aes(x=V1,y=V2,fill=normed_score)) +
                          geom_point(size = size_dot, shape = 21, stroke = size_stroke, color = "black")+
                          scale_fill_viridis_c()+
                          labs(color = 'Score')+
                          xlab(paste0('Human LV1 (',lv_vars[1],'%)')) + ylab(paste0('Human LV2 (',lv_vars[2],'%)')) +
                          facet_wrap(~phenotype)
plt_PLSR_backproject <- add_theme(plt_PLSR_backproject)

#################################################################################
### Panel - PLSR performance in backprojected data
Yh_back <- Xh %*% Wm %*% t(Wm) %*% Wh %*% Bh 
plt_PLSR_training_back <- rbind(data.frame(Measured = Yh[,1], Phenotype = "MAS"), data.frame(Measured = Yh[,2], Phenotype = "Fibrosis stage")) %>%
                          mutate(Predicted = c(Yh_back[,1] + mean(Yh[,1]), Yh_back[,2] + mean(Yh[,2]))) %>%
                          ggplot(aes(x = Measured, y = Predicted)) +
                          geom_jitter(size = size_dot, shape = 21, stroke = size_stroke, fill = "steelblue", color = "black", width = 0.1) +
                          geom_segment(data = data.frame(Phenotype = c("Fibrosis stage", "MAS"),
                                                         low = c(0,0),
                                                         high = c(4,8)), 
                                       aes(x = low, xend = high, y = low, yend = high),
                                       color = "black", linewidth = size_line, linetype = 2) +
                          stat_cor(color = "black", size = size_annotation*0.7) +
                          facet_wrap(~Phenotype, scales = "free") +
                          labs(x = "Measured", y = "Predicted")
plt_PLSR_training_back <- add_theme(plt_PLSR_training_back)

#################################################################################
### Panel - Predictions of TCs vs all PCs
Yh_back_TC <- Xh %*% Wm_TC %*% t(Wm_TC) %*% Wh %*% Bh 
plt_TC_prediction <- rbind(data.frame(Measured = Yh_back[,1] + mean(Yh[,1]), Phenotype = "MAS"), 
                           data.frame(Measured = Yh_back[,2] + mean(Yh[,2]), Phenotype = "Fibrosis stage")) %>%
                          mutate(Predicted = c(Yh_back_TC[,1] + mean(Yh[,1]), Yh_back_TC[,2] + mean(Yh[,2]))) %>%
                          ggplot(aes(x = Measured, y = Predicted)) +
                          geom_jitter(size = size_dot, shape = 21, stroke = size_stroke, fill = "steelblue", color = "black", width = 0.1) +
                          geom_abline(intercept = 0, slope = 1, linewidth = 1, color = "black", linetype = 2) +
                          stat_cor(color = "black", size = size_annotation*0.7) +
                          facet_wrap(~Phenotype, scales = "free") +
                          labs(x = "Prediction with all PCs", y = "Prediction with TCs only")
plt_TC_prediction <- add_theme(plt_TC_prediction)


#################################################################################
### Panel - Projection of MPS onto TCs colored by TGFb
plt_TC_MPS <- data.frame(x = Xm %*% Wm_TC[,1], y = Xm %*% Wm_TC[,2], TGFb = data_list$Kostrzewski$metadata$TGF) %>%
              ggplot(aes(x = x, y= y, color = TGFb, fill = TGFb)) +
                geom_point(size = size_dot, shape = 21, stroke = size_stroke, color = "black") +
                scale_fill_manual(values = c('#66C2A5','#FC8D62')) +
                scale_color_manual(values = c('#66C2A5','#FC8D62')) +
                stat_ellipse(linewidth = size_line, show.legend = F, linetype = 2) +
                labs(x = 'MPS TC1', y = 'MPS TC2')  

plt_TC_MPS <- add_theme(plt_TC_MPS) + theme(legend.position = "top")

#################################################################################
### Panel - Pathway activity of TCs
translatable_components_progenies <- pathway_activity_interpretation(Wm_TC, Wm)

plt_pwy_TC1 <- plot_pwy_activity(translatable_components_progenies %>% filter(condition == "LV_data1"), plt_lim = 16, show_fill_legend = T)
plt_pwy_TC1 <- add_theme(plt_pwy_TC1)

plt_pwy_TC2 <- plot_pwy_activity(translatable_components_progenies %>% filter(condition == "LV_data2"), plt_lim = 16, show_fill_legend = T)
plt_pwy_TC2 <- add_theme(plt_pwy_TC2)

### Save panels as figures
ggsave(filename = "./Figures/figure3/plt_PCA_MPS.pdf", plot = plt_PCA_MPS, units = "cm", width = 5, height = 5)
ggsave(filename = "./Figures/figure3/plt_PLSR_backproject.pdf", plot = plt_PLSR_backproject, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure3/plt_PLSR_training_back2.pdf", plot = plt_PLSR_training_back, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure3/plt_TC_prediction.pdf", plot = plt_TC_prediction, units = "cm", width = 8, height = 5)
ggsave(filename = "./Figures/figure3/plt_TC_MPS2.pdf", plot = plt_TC_MPS, units = "cm", width = 5, height = 6)
ggsave(filename = "./Figures/figure3/plt_pwy_TC1.pdf", plot = plt_pwy_TC1, units = "cm", width = 7, height = 6)
ggsave(filename = "./Figures/figure3/plt_pwy_TC2.pdf", plot = plt_pwy_TC2, units = "cm", width = 7, height = 6)


