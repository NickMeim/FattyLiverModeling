theta1 <- 0.5*pi*seq(0,1,0.05)
theta2 <- theta1

cor_theta_NAS <- matrix(data = 0, nrow = length(theta1), ncol = length(theta1))
cor_theta_Fib <- cor_theta_NAS

cor_theta_long <- NULL

for (ii in 1:length(theta1)){
  for (jj in 1:length(theta2)){
    V1 <- sin(theta1[ii])*Wm_opt[,1] + cos(theta1[ii])*Wm_opt[,2]
    V2 <- sin(theta2[jj])*Wm_opt[,1] + cos(theta2[jj])*Wm_opt[,2]
    W_dummy <- cbind(Wm, V1,V2)
    Y_theta <- Xh %*% W_dummy %*% t(W_dummy) %*% Wh %*% Bh
    
    cor_theta_NAS[ii,jj] <- cor(Y_theta[,1], Yh[,1])
    cor_theta_Fib[ii,jj] <- cor(Y_theta[,2], Yh[,2])
    
    cor_theta_long <- rbind(cor_theta_long,
                            data.frame(theta1 = theta1[ii],
                                       theta2 = theta2[jj],
                                       NAS = cor_theta_NAS[ii,jj],
                                       Fibrosis = cor_theta_Fib[ii,jj]))
    
    
    
  }
}


corrplot::corrplot(cor_theta_NAS, is.corr = F)


cor_theta_long %>%
  ggplot(aes(x = theta1/(0.5*pi), y = theta2/(0.5*pi), fill = NAS)) +
    geom_point(color = "white", shape = 22, size = 10) +
    scale_fill_viridis_c()


pwy_net <- decoupleR::get_progeny()
colnames(Wm_opt) <- c("LV1", "LV2")
colnames(pwy_net)[3] <- "mor"
pwy_Wextra <- decoupleR::run_mlm(Wm_opt, pwy_net)


plt_pwy_LV <- pwy_Wextra %>% filter(condition == "LV1") %>%
              ggplot(aes(x = score, y = reorder(source, score), fill = score)) +
                geom_col(show.legend = F) +
                labs(x = "Pathway score", y = NULL) +
                scale_fill_gradient2(low='darkblue',high = 'indianred',
                                     mid = 'whitesmoke',midpoint = 0)
  
  
ggsave('./results/plot_files/plt_pwy_LV.png', plot = add_theme(plt_pwy_LV),
       device = "png",
       height = 10,
       width=6,
       units = 'cm')

plt_msig_LV <- df_msig %>% 
                filter(LV == "LV1") %>% 
                top_n(12, abs(NES)) %>%
                ggplot(aes(x = NES, y = reorder(Hallmark, -NES), fill = NES)) + 
                  geom_col(show.legend = F) +
                labs(x = "Normalized enrichment score", y = NULL) +
                scale_fill_gradient2(low='darkblue',high = 'indianred',
                                     mid = 'whitesmoke',midpoint = 0)

ggsave('./results/plot_files/plt_msig_LV.png', plot = add_theme(plt_msig_LV),
       device = "png",
       height = 10,
       width=10,
       units = 'cm')

extra_basis_inferred_perts <- extra_basis_inferred_perts %>% 
                              group_by("condition") %>% 
                              mutate(rank = rank(activity))

plt_extra_basis_perts <- extra_basis_inferred_perts %>% 
  filter(condition == "LV1") %>% 
  ggplot(aes(x = rank, y = activity)) + 
  geom_point(size = size_dot, color = "indianred") + 
  geom_text_repel(data = extra_basis_inferred_perts %>% 
               top_n(8, abs(activity)), 
             aes(label = perturbation),
             size = size_annotation) +
  labs(x = "Rank", y = "Inferred perturbation activity")

ggsave('./results/plot_files/plt_extra_basis_perts.png', 
       plot = add_theme(plt_extra_basis_perts),
       device = "png",
       height = 10,
       width=7,
       units = 'cm')




phi <- Wh %*% Bh[,c(2,1)]
Wm_opt2 <- analytical_solution_opt(y=Yh[,c(2,1)],
                                  W_invitro = Wm,
                                  phi = phi)

Wm_tot2 <- cbind(Wm, Wm_opt2)



cor_theta_dummy <- theta1 
for (ii in 1:length(theta1)){
  V1 <- cos(theta1[ii])*Wm_opt[,1] + sin(theta1[ii])*Wm_opt2[,1]
  Wdummy <- cbind(Wm, V1)
  cor_theta_dummy[ii] <- cor(Yh[,2], Xh %*% Wdummy %*% t(Wdummy) %*% Wh %*% Bh[,2])
}



Wm_opt_dummy <- analytical_solution_opt(y=Yh[,2],
                                   W_invitro = cbind(Wm,V1),
                                   phi = Wh %*% Bh[,2])

phi <- Wh %*% Bh[,2]
W_invitro <- cbind(Wm, V1)
alpha <- t(W_invitro) %*% phi
Wopt_dummy <- (phi - W_invitro %*% alpha)/sqrt(sum(phi^2) -sum(alpha^2))

Wm_tot_dummy <- cbind(W_invitro, Wopt_dummy)


# Project-backproject with retrain
Y_pred <- Yh
Y_pred[,1] <- predict(lm(Yh[,1] ~ Xh %*% Wm %*% t(Wm) %*% Wh))
Y_pred[,2] <- predict(lm(Yh[,2] ~ Xh %*% Wm %*% t(Wm) %*% Wh))


all_scatter_plot <- left_join(data.frame(Yh) %>% 
                                mutate(id = seq(1,nrow(Yh))) %>% 
                                gather('phenotype','true',-id),
                              data.frame(Y_pred) %>% 
                                mutate(id = seq(1,nrow(Y_pred))) %>% 
                                gather('phenotype','prediction',-id)) %>%
  select(-id)
all_cor_results <- all_scatter_plot %>%
  group_by(phenotype) %>%
  summarise(cor_test = list(cor.test(true, prediction))) %>%
  mutate(cor_coef = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value))


plt_plsr_scatter_retrain <- all_scatter_plot %>% 
  ggplot(aes(x = true,y=prediction)) +
  geom_jitter(width = 0.05, color = "steelblue", size = size_dot) + 
  geom_abline(slope=1,intercept = 0, linetype = 2)+
  geom_text(data = all_cor_results, 
            aes(x = 0, y = Inf, label = sprintf("r = %.2f, p = %.2g", cor_coef, p_value)),
            hjust = 0, vjust =  1.5, size = size_annotation) +
  facet_wrap(~phenotype,scales = 'free', labeller = labeller(phenotype = c(fibrosis = "Fibrosis stage", NAS = "NAS"))) +
  labs(x = "Measured", y = "Prediction")



ggsave('./results/plot_files/plt_plsr_scatter_retrain.pdf', 
       plot = add_theme(plt_plsr_scatter_retrain),
       device = "pdf",
       height = 6,
       width=9,
       units = 'cm')


