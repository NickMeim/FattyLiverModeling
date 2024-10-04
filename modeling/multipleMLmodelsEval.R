library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(patchwork)

### Load results-----------------
tested_models <- list.files('../results/MLresults/',include.dirs = TRUE)
files <- list.files('../results/MLresults/',recursive = TRUE)
files <- files[grepl('.csv',files)]
cv_files <- files[grepl('_cv',files)]
clinical_files <- files[grepl('_external',files)]
df_res_cv <- data.frame() 
df_res_external <- data.frame() 
for (i in 1:length(clinical_files)){
  df_res_cv <- rbind(df_res_cv,
                     data.table::fread(paste0('../results/MLresults/',cv_files[i])) %>% select(-V1))
  df_res_external <- rbind(df_res_external,
                           data.table::fread(paste0('../results/MLresults/',clinical_files[i])) %>% select(-V1))
}
df_res_cv <- df_res_cv %>% gather('input','r',-model,-set,-fold)
num_folds <- length(unique(df_res_cv$fold))

### Visualize results---------------------------------
df_res <- rbind(df_res_cv,df_res_external %>% mutate(set=dataset) %>% select(all_of(colnames(df_res_cv))))
df_plot <- df_res %>% mutate(input = ifelse(input=='back-projected','truncated',input)) %>%
  group_by(model,set,input) %>% mutate(mu = mean(r)) %>% mutate(se = sd(r)/sqrt(num_folds)) %>% ungroup() %>%
  mutate(input = factor(input,levels=c('optimized MPS','truncated','human genes')))
df_plot$model <- factor(df_plot$model,
                           levels = c('PLSR','neuralNet','rf','xgboost','svmRBF','svmPoly','svmLinear','knn','elasticNet','lasso','ridge'))
df_plot$model <-factor(df_plot$model, levels = rev(levels(df_plot$model)))
df_plot$set <- factor(df_plot$set,levels = c('train','test','Hoang','Pantano'))
p1 <- ggplot(df_plot,
       aes(x=model,y=mu,color=input,fill=input)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymax = mu + se, ymin = mu - se),color='black', position = position_dodge(width = 0.8), width = 0.25, size = 1) +
  geom_hline(yintercept = 0,size=1,color='black',linetype='dashed')+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(NA,1))+
  ylab('correlation')+
  theme_pubr(base_size = 24,base_family = 'Arial') +
  theme(text = element_text(size = 24,family = 'Arial'),
        # axis.line.y = element_blank(),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size = 24,family = 'Arial')) +
  facet_wrap(~set)+
  coord_flip()
print(p1)

ggsave('../figures/multiple_ml_models_eval.png',
       plot = p1,
       width = 12,
       height = 9,
       units = 'in',
       dpi=600)
