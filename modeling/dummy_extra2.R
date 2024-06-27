# Naive comparison - do correlations based on RNAseq
cor_Xh_Xm <- cor(t(Xh), t(Xm))

cor_Xh_Xm <- cor_Xh_Xm %>% as.data.frame() %>% mutate(human_sample = rownames(.)) %>%
              pivot_longer(cols = 1:179, names_to = "mps_sample", values_to = "cor")


cor_Xh_Xm <- merge(cor_Xh_Xm, data_frame(human_sample = rownames(Xh), NAS = Yh[,1], fibrosis = Yh[,2]), by = "human_sample")


cor_Xh_Xm %>% as_tibble() %>% 
  heatmap(mps_sample, human_sample,cor) %>%
    add_tile(NAS)

# Group by MPS sample and get the median NAS and fibrosis for the top 10 correlations
cor_Xh_Xm_top <- cor_Xh_Xm %>% 
                  group_by(mps_sample) %>% 
                  top_n(10, cor) %>% 
                  summarize(NAS = median(NAS), fibrosis = median(fibrosis)) %>%
                  mutate(X = mps_sample)

cor_Xh_Xm_top <- merge(cor_Xh_Xm_top, data_list$Kostrzewski$metadata, by = "X")


plt_cor_NAS <- cor_Xh_Xm_top %>%
              ggplot(aes(x = NAS, y = reorder(Combination_name, NAS))) +
                geom_boxplot() +
                labs(x = "Top correlated human NAS", y = NULL) +
                xlim(0,8)

plt_cor_fib <- cor_Xh_Xm_top %>%
                ggplot(aes(x = fibrosis, y = reorder(Combination_name, fibrosis))) +
                geom_boxplot() +
                labs(x = "Top correlated human Fibrosis", y = NULL) +
                xlim(0,4)

ggsave('../results/plot_files/plt_cor_NAS.png', 
       plot = add_theme(plt_cor_NAS),
       device = "png",
       height = 6,
       width=9,
       units = 'cm')

ggsave('../results/plot_files/plt_cor_fib.png', 
       plot = add_theme(plt_cor_fib),
       device = "png",
       height = 6,
       width=9,
       units = 'cm')

# Signatures
gNAS <- c("ARPC5", "YWHAH", "ARF4", "TNFRSF12A", "ADHFE1",
          "USP33", "CD52", "ACVR2B", "ING5", "ASB3", "IFI30",
          "ERVW-1", "YWHAZ", "ERBB3", "KPNA2", "COQ10B", "MAGI1",
          "MAPRE1", "ABCA6")

gFib <- c("TIMP1", "MYLI2B", "LUM", "ZNF395", "AKAP9",
          "ACTR2", "LGALS3", "MAPRE1", "FRK", "ANKRD28",
          "IGFBP7", "YWHAZ", "USP33", "CD59", "TAX1BP3",
          "FAM221A", "ADHFE1", "TNFRSF12A")

sig_regulon <- rbind(
  data.frame(source = "NAS", target = gNAS, mor = c(1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1)),
  data.frame(source = "Fibrosis stage", target = gFib, mor = c(1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,1,-1,-1,-1)))

sig_viper <- decoupleR::run_viper(t(Xh), net = sig_regulon)
colnames(sig_viper)[3] <- "GEO_Accession..exp."

sig_viper <- merge(sig_viper, data_list$Govaere$metadata[, c("GEO_Accession..exp.", "Fibrosis_stage", "nas_score")], by = "GEO_Accession..exp.")

sig_viper <- sig_viper %>% 
              pivot_longer(cols = c("Fibrosis_stage", "nas_score"), 
                           names_to = "Phenotype", 
                           values_to = "val") %>%
              mutate(Phenotype = factor(ifelse(Phenotype == "nas_score", "NAS", "Fibrosis stage")))

plt_sig_viper <- sig_viper %>%
                  ggplot(aes(x = val, y = score)) +
                  geom_jitter(width = 0.1, color = "steelblue", size = size_dot) +
                  facet_wrap(~Phenotype, scales = "free") +
                  geom_smooth(method = "lm", se = T, color = "indianred", linewidth = size_line, linetype = 2) +
                  labs(x = "Histological score", y = "Gene signature enrichment score")


ggsave('../results/plot_files/plt_sig_viper.png', 
       plot = add_theme(plt_sig_viper),
       device = "png",
       height = 6,
       width=9,
       units = 'cm')





sig_viper_MPS <- decoupleR::run_viper(t(Xm), net = sig_regulon)
colnames(sig_viper_MPS)[3] <- "X"
sig_viper_MPS <- merge(sig_viper_MPS, data_list$Kostrzewski$metadata, by = "X")


plt_sig_viper_MPS <- sig_viper_MPS %>%
                      ggplot(aes(x = score, y = reorder(Combination_name, score))) +
                        geom_boxplot() +
                        facet_wrap(~source) +
                        labs(x = "Gene signature enrichment score", y = NULL)




ggsave('../results/plot_files/plt_sig_viper_MPS.png', 
       plot = add_theme(plt_sig_viper_MPS),
       device = "png",
       height = 6,
       width=12,
       units = 'cm')
