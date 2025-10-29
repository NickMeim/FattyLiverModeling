### Functions for making elements of figure 5 - Interpretation of LV extra
### This function will load results from the pipeline to generate plots.
### These will be specific for the chosen datasets
library(tidyverse)
library(matrixStats)
library(ggplot2)
library(ggExtra)
root_dir <- "C:/Users/nmeim/Documents/LiverModeling/FattyLiverModeling"
setwd(root_dir)
source("utils/plotting_functions.R")
source("modeling/functions_translation.R")
source("modeling/vector_space_interpretation.R")
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
### Panel - Radial plots

# Load results and process
df_radial <- readRDS(paste0("results/df_correlation_radial_", tolower(target_dataset), "_speamaned.rds"))

df_mu_radial <- df_radial %>% group_by(theta) %>% mutate(mu = mean(corr)) %>% ungroup() %>%
  select(-corr,-phenotype) %>% unique() %>%
  mutate(phenotype = 'average') %>% select(theta,phenotype,c('corr'='mu')) %>%
  mutate(x = corr*cos(theta), y = corr*sin(theta))
df_radial <- df_radial %>%
  mutate(x = corr*cos(theta), y = corr*sin(theta)) %>%
  mutate(phenotype = ifelse(phenotype == "fibrosis", "Fibrosis stage",phenotype)) %>%
  mutate(phenotype = ifelse(phenotype == "NAS", "MAS",phenotype))

# Plot
plt_cor_radial <-   ggplot(df_radial %>% filter(theta<=90),
                           aes(x = theta, y = corr, color = phenotype, group = phenotype)) +
  geom_line(linewidth = 1.2*size_line) +
  geom_line(data = df_mu_radial %>% filter(theta<=90),
            aes(x = theta, y = corr),
            color='black',lwd= 1.2*size_line,linetype='dashed',
            inherit.aes = FALSE)+
  scale_x_continuous(breaks = seq(0,90,15))+
  geom_hline(yintercept = 0.4,color='black',lwd = size_line)+
  geom_vline(xintercept = 90,color='black',lwd = size_line)+
  coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F) +
  labs(y = "Spearman`s rank correlation", x = 'LV extra 2') +
  ggtitle('LV extra 1') +
  scale_color_brewer(palette = "Dark2")

# Add theme
plt_cor_radial <- add_theme(plt_cor_radial)

# Adjust theme
plt_cor_radial <- plt_cor_radial +
  theme(text = element_text(size = size_text),
        plot.title = element_text(vjust = -12,hjust = 0.1,size= size_text,face = 'bold'),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(hjust = 0.66),
        axis.text.x = element_blank(),
        axis.title.x = element_text(vjust = 8,hjust = 1.05,size= size_text,face='bold'),
        axis.line.x = element_line(linewidth = size_line),
        axis.line.y = element_line(linewidth = size_line),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -1.5, unit = "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = size_line),
        panel.grid.minor = element_line(linewidth = size_line))

#################################################################################
### Panel - Pathway radial plots

df_pwy_plot <- readRDS(paste0("results/df_pwy_radial_", tolower(target_dataset),"_.rds"))

plt_pwy_radial <- ggplot(df_pwy_plot  %>% gather('theta','activity',-Pathway) %>%
                          mutate(theta=as.numeric(theta))%>%
                          filter(Pathway %in% c('JAK-STAT','p53','NFkB')),
                          aes(x = theta, y = activity, color = Pathway,fill=Pathway, group = Pathway)) +
  geom_line(linewidth = size_line, show.legend = F) +
  geom_area(alpha=0.2, show.legend = F)+
  scale_x_continuous(breaks = seq(0,90,15))+
  scale_y_continuous(limits = c(0,8))+
  scale_colour_manual(values = c('#53B400','#00C094','#00BFC4'))+
  scale_fill_manual(values = c('#53B400','#00C094','#00BFC4'))+
  geom_hline(yintercept = 0,color='black',lwd=size_line)+
  geom_vline(xintercept = 90,color='black',lwd=size_line)+
  labs(y = "Absolute pathway activity", x = NULL) +
  coord_radial(start = 0, end = 0.5*pi, inner.radius = 0.4, expand = F) +
  facet_wrap(~Pathway)

plt_pwy_radial <- add_theme(plt_pwy_radial)

plt_pwy_radial <- plt_pwy_radial +
                    theme(axis.text.y = element_text(color = "black"),
                          axis.title.y = element_text(hjust = 0.66),
                          axis.text.x = element_blank(),
                          # axis.title.x = element_text(vjust = 9,hjust = 1.05,size=18,face='bold'),
                          axis.line.x = element_line(linewidth = size_line),
                          axis.line.y = element_line(linewidth = size_line),
                          # legend.key.size = unit(0.3, "cm"),
                          # legend.margin = margin(t = 0, r = 0, b = 0, l = -1.5, unit = "cm"),
                          legend.position = 'none',
                          panel.border = element_blank(),
                          panel.grid.major = element_line(linewidth = size_line),
                          panel.grid.minor = element_line(linewidth = size_line))

#################################################################################
### Panel - Pathway activity on LV1
Wm <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_total.rds'))
Wm_opt <- readRDS(paste0('results/Wm_',tolower(target_dataset),'_extra.rds'))
colnames(Wm_opt) <- c("LV_extra1", "LV_extra2")

translatable_components_progenies <- pathway_activity_interpretation(Wm_opt, Wm)

plt_pwy_LV1 <- plot_pwy_activity(translatable_components_progenies %>% filter(condition == "LV_extra1"),
                                 plt_lim = 9, show_fill_legend = T, offset_annot = 0.8, n.breaks = 6)

plt_pwy_LV1 <- add_theme(plt_pwy_LV1)

plt_pwy_LV2 <- plot_pwy_activity(translatable_components_progenies %>% filter(condition == "LV_extra2"),
                                 plt_lim = 9, show_fill_legend = T, offset_annot = 0.8, n.breaks = 6)

plt_pwy_LV2 <- add_theme(plt_pwy_LV2)

#################################################################################
### Panel - Hallmark enrichment on LV1

df_msig <- readRDS(paste0("results/hallmark_enrichment_", tolower(target_dataset),"_LVopt.rds" ))


plt_hallmarks_LV1 <- ggplot(df_msig %>% filter(LV=='V1') %>% arrange(NES) %>%
                filter(padj<=0.1), aes(x=NES,y=reorder(Hallmark,-NES),fill=NES))+
         geom_bar(stat = 'identity', size = size_col, color = "black", show.legend = F) +
         # scale_fill_gradient(trans='log10',low = "red",high = "white",limits = c(min(df_msig$padj),1)) +
         scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',midpoint = 0) +
         labs(x = 'Normalized Enrichment Score', y = 'Hallmark')

plt_hallmarks_LV1 <- add_theme(plt_hallmarks_LV1)

#################################################################################
### Panel - ChemPert inference
extra_basis_inferred_perts <- readRDS(paste0("results/utilityRDS/extra_basis_", tolower(target_dataset),"_inferred_perts.rds"))

plt_LV1_extra_chempert <- extra_basis_inferred_perts[[2]] %>% filter(condition == "V1") %>%
                              ggplot(aes(x=as.numeric(reorder(perturbation,activity)),y=activity)) +
                              geom_point(size=size_dot,color = '#CD5C5C') +
                              geom_text_repel(aes(label=significant),size=0.75*size_annotation, max.overlaps=80,
                                              box.padding = 0.5, nudge_x = 1, nudge_y = 0.2)+
                              labs(x = 'Rank', y = 'Inferred perturbation activity')

plt_LV1_extra_chempert <- add_theme(plt_LV1_extra_chempert)


#################################################################################
### Panel - Pathway activity of TCs

### Save panels as figures
ggsave(filename = "figures/figure5/plt_cor_radial.pdf", plot = plt_cor_radial, units = "cm", width = 6, height = 6)
ggsave(filename = "figures/figure5/plt_pwy_radial.pdf", plot = plt_pwy_radial, units = "cm", width = 11, height = 7)
ggsave(filename = "figures/figure5/plt_pwy_LV1.pdf", plot = plt_pwy_LV1, units = "cm", width = 6, height = 6)
ggsave(filename = "figures/figure5/plt_hallmarks_LV1.pdf", plot = plt_hallmarks_LV1, units = "cm", width = 7, height = 6)
ggsave(filename = "figures/figure5/plt_LV1_extra_chempert.pdf", plot = plt_LV1_extra_chempert, units = "cm", width = 5, height = 6)
ggsave(filename = "figures/figure5/plt_pwy_LV2.pdf", plot = plt_pwy_LV2, units = "cm", width = 5, height = 6)




