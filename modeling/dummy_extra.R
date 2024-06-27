# Plot translatable components
Wm_combo <- readRDS("../results/Wm_kostrzewski_combo.rds")

Wh <- read.table("../results/Wh_govaere.csv", header = T, sep = ",", row.names = 1)

plot(predict(plsr_model, Xh %*% Wm %*% t(Wm))[,1], predict(plsr_model, Xh %*% Wm_combo %*% t(Wm_combo))[,1])

# Project human data on these components
plt_plsr_trans_human <- rbind(data.frame(x = Xh %*% Wm_combo[,1], y = Xh %*% Wm_combo[,2], phenotype = "Fibrosis stage", Score = Yh[,2]),
                             data.frame(x = Xh %*% Wm_combo[,1], y = Xh %*% Wm_combo[,2], phenotype = "NAS", Score = Yh[,1])) %>%
                        ggplot(aes(x = x, y = y, color = Score)) +
                        geom_point(size = size_dot) +
                        labs(x = "MPS TC1", y = "MPS TC2") +
                        scale_color_viridis_c() +
                        facet_wrap(~factor(phenotype), nrow = 1)

plt_plsr_trans_human <- add_theme(plt_plsr_trans_human)

ggsave('../results/plot_files/plt_plsr_trans_human.png', 
       plot = plt_plsr_trans_human,
       device = "png",
       height = 6,
       width=10.5,
       units = 'cm')

# Project MPS on trans components
plt_plsr_trans_MPS <- data.frame(x = Xm %*% Wm_combo[,1], y = Xm %*% Wm_combo[,2]) %>%
                        ggplot(aes(x = x, y = y)) +
                        geom_point(size = size_dot, color = "indianred") +
                        labs(x = "MPS TC1", y = "MPS TC2")

ggsave('../results/plot_files/plt_plsr_trans_MPS.png', 
       plot = add_theme(plt_plsr_trans_MPS),
       device = "png",
       height = 6,
       width=6,
       units = 'cm')

# Color by TGFb
plt_plsr_trans_MPS2 <- data.frame(x = Xm %*% Wm_combo[,1], y = Xm %*% Wm_combo[,2], data_list$Kostrzewski$metadata) %>%
                        ggplot(aes(x = x, y = y, color = TGF)) +
                        geom_point(size = size_dot) +
                        labs(x = "MPS TC1", y = "MPS TC2") +
                        scale_color_brewer(palette = color_palette) 

ggsave('../results/plot_files/plt_plsr_trans_MPS2.png', 
       plot = add_theme(plt_plsr_trans_MPS2),
       device = "png",
       height = 6,
       width=7,
       units = 'cm')

# Load extra LVs
Wm_extra <- readRDS("../results/Wm_kostrzewski_extra.rds")
Wm_tot <- cbind(Wm, Wm_extra)

all_scatter_plot_extra <- left_join(data.frame(Yh) %>% 
                                mutate(id = seq(1,nrow(Yh))) %>% 
                                gather('phenotype','true',-id),
                              data.frame(predict(plsr_model, Xh %*% Wm_tot %*% t(Wm_tot))) %>% 
                                mutate(id = seq(1,nrow(Yh))) %>% 
                                gather('phenotype','prediction',-id)) %>%
                                select(-id)

all_cor_results_extra <- all_scatter_plot_extra %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))


plt_plsr_scatter_extra <- all_scatter_plot_extra %>% 
                              ggplot(aes(x = true,y=prediction)) +
                              geom_jitter(width = 0.05, color = "steelblue", size = size_dot) + 
                              geom_abline(slope=1,intercept = 0, linetype = 2)+
                              geom_text(data = all_cor_results_extra, 
                                        aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
                                        hjust = 0, vjust =  1.5, size = size_annotation) +
                              facet_wrap(~phenotype,scales = 'free', labeller = labeller(phenotype = c(fibrosis = "Fibrosis stage", NAS = "NAS"))) +
                              labs(x = "Measured", y = "Prediction")



ggsave('../results/plot_files/plt_plsr_scatter_extra.png', plot = add_theme(plt_plsr_scatter_extra),
       device = "png",
       height = 6,
       width=9,
       units = 'cm')


# Project human data onto these extra variables
plt_PCA_human_extra <- rbind(data.frame(x = Xh %*% Wm_extra[,1], y = Xh %*% Wm_extra[,2], phenotype = "Fibrosis stage", Score = Yh[,2]),
                               data.frame(x = Xh %*% Wm_extra[,1], y = Xh %*% Wm_extra[,2], phenotype = "NAS", Score = Yh[,1])) %>%
  ggplot(aes(x = x, y = y, color = Score)) +
  geom_point(size = size_dot) +
  labs(x = "LV extra 1", y = "LV extra 2") +
  scale_color_viridis_c() +
  facet_wrap(~factor(phenotype), nrow = 1)



ggsave('../results/plot_files/plt_PCA_human_extra.png', plot = add_theme(plt_PCA_human_extra),
       device = "png",
       height = 6,
       width=9,
       units = 'cm')

# CV of LV opt
p_test3 <- ggboxplot(performance_all_plot %>% filter(set=='test') %>% filter(approach %in% c("PLSR", "backprojected", "analytical Wopt", "shuffle X")),
                     x='approach',y='r',add='jitter') +
  scale_y_continuous(breaks = seq(-0.5,1,0.1))+
  labs(x = NULL, y = "Mean Pearson correlation") +
  stat_compare_means(comparisons = list(c('PLSR','backprojected'),
                                        c('analytical Wopt','backprojected'),
                                        # c('PLSR','Wopt'),
                                        c('PLSR','analytical Wopt')),
                     method = 'wilcox',
                     tip.length = 0.01,
                     label.y = c(0.85, 0.85, 1),
                     size = size_annotation) +
  ylim(c(-0.6, 1.2)) +
  scale_x_discrete(name = NULL, labels = c("PLSR" = "PLSR \n human", 
                                           "backprojected" = "MPS \n backprojected", 
                                           "analytical Wopt" = "MPS \n extended", 
                                           "shuffle X" = "shuffle X"))

ggsave('../results/plot_files/plt_CV_plsr_test_all3.png', plot = add_theme(p_test3),
       device = "png",
       height = 6,
       width=8,
       units = 'cm')





##### Make radial plot from scrath
df_radial <- data.frame(theta = pi*seq(0,0.5,0.02), cor_NAS = 0, cor_Fib = 0)
df_pwy_radial <- NULL

for (ii in 1:nrow(df_radial)){
  theta <- df_radial$theta[ii]
  Wm_theta <- Wm_extra[,1]*cos(theta) + Wm_extra[,2]*sin(theta)
  Wm_tot_theta <- cbind(Wm, Wm_theta)
  Y_pred_theta <- predict(plsr_model, Xh %*% Wm_tot_theta %*% t(Wm_tot_theta))
  df_radial$cor_NAS[ii] <- cor(Y_pred_theta[,1], Yh[,1])
  df_radial$cor_Fib[ii] <- cor(Y_pred_theta[,2], Yh[,2])
  
  pwy_acts <- run_mlm(mat = Wm_theta, net = pwy_net, .source='source', .target='target',
                      .mor='weight', minsize = 5)
  
  df_pwy_radial <- rbind(df_pwy_radial, 
                         data.frame(theta = theta, source = pwy_acts$source, score = pwy_acts$score))
  
}
# Transform to cartesian
df_radial <- df_radial %>% 
              pivot_longer(cols = 2:3, names_to = "phenotype", values_to = "cor") %>%
              mutate(x = cor*cos(theta), y = cor*sin(theta)) %>%
              mutate(phenotype = ifelse(phenotype == "cor_NAS", "NAS", "Fibrosis stage"))



plt_cor_radial <-   ggplot(df_radial, aes(x = theta, y = cor, color = phenotype, group = phenotype)) + 
                    geom_line(linewidth = 1) +
                    coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F, direction = 1) +
                    labs(y = "Pearson correlation", x = NULL) +
                    theme(text = element_text(size = size_text),
                          axis.text.y = element_text(color = "black"),
                          axis.text.x = element_text(color = "black"),
                          legend.key.size = unit(size_legend, "cm")) +
                    scale_color_brewer(palette = "Dark2")


ggsave('../results/plot_files/plt_cor_radial.pdf', plot = plt_cor_radial,
       device = "pdf",
       height = 6,
       width=8,
       units = 'cm')

  
plt_cor_radial2 <-   ggplot(df_radial, aes(x = theta/(0.5*pi), y = cor, color = phenotype, group = phenotype)) + 
  geom_line(linewidth = 1) +
  theme(text = element_text(size = size_text),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.key.size = unit(size_legend, "cm")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Mixing proportion", y = "Pearson correlation", color = "Phenotype")

ggsave('../results/plot_files/plt_cor_radial2.png', plot = add_theme(plt_cor_radial2),
       device = "png",
       height = 6,
       width=8,
       units = 'cm')



plt_cor_radial_pwy <-  df_pwy_radial %>% filter(source %in% c("JAK-STAT", "p53")) %>%
                        ggplot(aes(x = theta/(0.5*pi), y = score, color = source)) + 
                        geom_line(linewidth = 1) +
                        theme(text = element_text(size = size_text),
                              axis.text.y = element_text(color = "black"),
                              axis.text.x = element_text(color = "black"),
                              legend.key.size = unit(size_legend, "cm")) +
                        labs(x = "Mixing proportion", y = "Inferred pathway activity", color = "Pathway") +
                        scale_color_brewer(palette = "Dark2")

ggsave('../results/plot_files/plt_cor_radial_pwy.png', plot = add_theme(plt_cor_radial_pwy),
       device = "png",
       height = 6,
       width=8,
       units = 'cm')


# Retrain regression coefficients
mod_NAS <- lm(Yh[,1] ~ Xh %*% Wm %*% t(Wm) %*% Wh)



# Project human data on translatable components
mod1 <- lm(Yh[,1] ~ Xh %*% Wm_combo)
mod2 <- lm(Yh[,2] ~ Xh %*% Wm_combo)

all_scatter_project_TC <-  rbind(data.frame(Ytrue = Yh[,2], Ypred = predict(mod2), Phenotype = "Fibrosis stage"),
                                 data.frame(Ytrue = Yh[,1], Ypred = predict(mod1), Phenotype = "NAS")) %>%
                            mutate(Phenotype = factor(Phenotype, levels = c("Fibrosis stage", "NAS")))

all_scatter_project_TC_res <- all_scatter_project_TC %>%
  group_by(Phenotype) %>%
  summarise(cor_test = list(cor.test(Ytrue, Ypred))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))




plt_cor_combo <- all_scatter_project_TC %>% 
  ggplot(aes(x = Ytrue,y=Ypred)) +
  geom_jitter(width = 0.05, color = "steelblue", size = size_dot) + 
  geom_abline(slope=1,intercept = 0, linetype = 2)+
  geom_text(data = all_scatter_project_TC_res, 
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
            hjust = 0, vjust =  1.5, size = size_annotation) +
  facet_wrap(~Phenotype,scales = 'free') +
  labs(x = "Measured", y = "Prediction") +
  ggh4x::facetted_pos_scales(y = list(scale_y_continuous(limits = c(0, 4)),
                                      scale_y_continuous(limits = c(0, 8))))

ggsave('../results/plot_files/plt_cor_combo.png', plot = add_theme(plt_cor_combo),
       device = "png",
       height = 6,
       width=9,
       units = 'cm')
