library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
source('CrossValidationUtilFunctions.R')
source('functions_translation.R')
source("../utils/plotting_functions.R")
source("vector_space_interpretation.R")
library(EGSEAdata)
library(org.Hs.eg.db)
library(fgsea)
egsea.data(species = "human",returnInfo = TRUE)

## Initial LIV2TRANS run loading-------------------------------
Wm_tot <- readRDS('../results/Wm_kostrzewski_total.rds')
Wm_opt <- readRDS('../results/Wm_kostrzewski_extra.rds')
Wm_combo <- readRDS('../results/Wm_kostrzewski_combo.rds')
Wm_opt <- as.data.frame(Wm_opt)
colnames(Wm_opt) <- c('extraLV1','extraLV2')

### Identify P53 and IFNa Genesets using Hallmark and Dorothea-------------------

## Load and use Hamllmarks
inds2keep <- lapply(msigdb, function(x){return(x["CATEGORY_CODE"]=="H")}) 
inds2keep <- which(do.call(c,inds2keep))
hallmarks <- data.frame()
for (i in inds2keep) {
  hallmarks <- rbind(hallmarks,t(as.data.frame(msigdb[[i]])))
}
hallmarks <- as.data.frame(hallmarks)
rownames(hallmarks) <- hallmarks$STANDARD_NAME

p53_set <- hallmarks %>% filter(grepl('P53',STANDARD_NAME))
p53_set <- str_split(p53_set$MEMBERS_SYMBOLIZED,pattern = ",")[[1]]
ifna_set <- hallmarks %>% filter(grepl('INTERFERON_ALPHA',STANDARD_NAME))
ifna_set <- str_split(ifna_set$MEMBERS_SYMBOLIZED,pattern = ",")[[1]]


## Load and use dorothea (for IFNa use its downstream maker TFs)
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
colnames(dorotheaData)[1] <- 'source' 
ifna_regulon <- dorotheaData %>% filter(source %in% c('IRF1','IRF9','STAT1','STAT2'))
p53_regulon <- dorotheaData %>% filter(source %in% c('TP53'))

## Get interection of those
p53_geneset <- unique(c(p53_regulon$target,p53_set))
p53_geneset <- p53_geneset[p53_geneset %in% rownames(Wm_combo)]
ifna_geneset <- unique(c(ifna_regulon$target,ifna_set))
ifna_geneset <- ifna_geneset[ifna_geneset %in% rownames(Wm_combo)]

## Check which of those are highly ranked in extra LV1---------------
data_p53 <- Wm_opt %>% rownames_to_column('name')
data_p53 <- data_p53 %>% mutate(geneset = ifelse(name %in% p53_geneset,'P53','other'))
data_p53 <- data_p53 %>% mutate(
  rank1 = rank(-abs(extraLV1), na.last = "keep"),
  rank2 = rank(-abs(extraLV2), na.last = "keep")
)
p1 <- ggplot(data_p53, aes(x =rank1, y = abs(extraLV1),color=geneset)) +
  geom_point(size=0.5,alpha=0.5)+
  # Annotate specific genes
  geom_text_repel(data = data_p53 %>% filter(geneset=='P53') %>% filter(rank1<=260),
                  aes(label = name),
                  size = 5,
                  fontface = "bold",
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "black",
                  segment.size = 1,
                  max.overlaps = 50) +
  labs(title = "Ranking of genes in the P53 geneset or downstream of the TP53 TF",
       x = "Rank",
       y = "Absolute weight in extra LV1") +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(size = 16,face = "bold"))
print(p1)
ggsave("../figures/p53_downstream_gene_candidates.png", p1, width = 10, height = 8, dpi = 600)


data_ifna <- Wm_opt %>% rownames_to_column('name')
data_ifna <- data_ifna %>% mutate(geneset = ifelse(name %in% ifna_geneset,'IFNA','other'))
data_ifna <- data_ifna %>% mutate(
  rank1 = rank(-abs(extraLV1), na.last = "keep"),
  rank2 = rank(-abs(extraLV2), na.last = "keep")
)
p2 <- ggplot(data_ifna, aes(x =rank1, y = abs(extraLV1),color=geneset)) +
  geom_point(size=0.5,alpha=0.5)+
  # Annotate specific genes
  geom_text_repel(data = data_ifna %>% filter(geneset=='IFNA') %>% filter(rank1<=330),
                  aes(label = name),
                  size = 5,
                  fontface = "bold",
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "black",
                  segment.size = 1,
                  max.overlaps = 50) +
  labs(title = "Ranking of genes in the IFNA geneset or downstream of the TF makers of IFNA",
       x = "Rank",
       y = "Absolute weight in extra LV1") +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(size = 16,face = "bold"))
print(p2)
ggsave("../figures/ifna_downstream_gene_candidates.png", p2, width = 10, height = 8, dpi = 600)
