library(tidyverse)
library(ggrepel)
library(ggsignif)
library(ggpubr)
library(patchwork)
library(matrixStats)
library(ropls)
library(LIV2Trans)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(LIV2Trans)
source('enrichment_calculations.R')
fig_dir <- "../figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# Function for loading datasets and appending them to a list
load_datasets <- function(dataset_names, dir_data){
  # Look at files in data directory
  files <- dir(dir_data)
  # Trim to datasets
  files <- files[grepl("dataset.RData", files)]
  # Loop for dataset names given and load one by one, appending to a list.
  # dataset names don't have to match perfectly
  data_list <- list(NULL)
  if (length(files) > 0){ # Loop if any match was found
    for (name in dataset_names){
      if (any(grepl(name, files, ignore.case = TRUE))){ # Avoid loading bad matches
        load(file = paste0(dir_data, files[grep(name, files, ignore.case = TRUE)]))
        print(paste0("Loading file: ", paste0(dir_data, files[grep(name, files, ignore.case = TRUE)])))
        data_list[[name]]$counts <- data
        data_list[[name]]$metadata <- metadata
        data_list[[name]]$genes <- rownames(data)
        if (exists("exp_factors")){
          data_list[[name]]$exp_factors <- exp_factors
        }
      }
    }
    # Remove redundant first element
    data_list[[1]] <- NULL
  }
  
  return(data_list)
  
}
run_phenotype_models <- function(X_h, X_v, Y, meta_v, space_label, fig_dir,
                                 phenos = c("NAS", "fibrosis")) {
  
  stopifnot(!is.null(rownames(X_h)), !is.null(rownames(X_v)))
  
  # Drop human samples missing a phenotype label
  complete <- complete.cases(Y)
  X_h <- X_h[complete, , drop = FALSE]
  Y   <- Y[complete,   , drop = FALSE]
  
  # ---- LOOCV harness ----
  loocv_predict <- function(X, y, train_fn, predict_fn) {
    preds <- numeric(length(y))
    for (i in seq_along(y)) {
      fit <- train_fn(X[-i, , drop = FALSE], y[-i])
      preds[i] <- predict_fn(fit, X[i, , drop = FALSE])
    }
    preds
  }
  
  # ---- Model definitions ----
  train_ord   <- function(X, y) {
    df <- as.data.frame(X); df$y <- ordered(y)
    MASS::polr(y ~ ., data = df, Hess = TRUE)
  }
  predict_ord <- function(fit, X_new) {
    p <- predict(fit, newdata = as.data.frame(X_new), type = "probs")
    if (is.null(dim(p))) p <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
    as.numeric(p %*% as.numeric(colnames(p)))
  }
  train_svm   <- function(X, y) e1071::svm(x = X, y = y, type = "eps-regression", kernel = "radial")
  predict_svm <- function(fit, X_new) as.numeric(predict(fit, X_new))
  train_lm    <- function(X, y) lm(y ~ ., data = data.frame(y = y, as.data.frame(X)))
  predict_lm  <- function(fit, X_new) as.numeric(predict(fit, newdata = as.data.frame(X_new)))
  
  models <- list(
    ordinal = list(train = train_ord, predict = predict_ord),
    svm     = list(train = train_svm, predict = predict_svm),
    linear  = list(train = train_lm,  predict = predict_lm)
  )
  
  # ---- Run all combos ----
  loocv_list  <- list()
  scores_list <- list()
  for (mname in names(models)) {
    for (ph in phenos) {
      y <- Y[, ph]
      pred_loocv <- loocv_predict(X_h, y, models[[mname]]$train, models[[mname]]$predict)
      loocv_list[[paste(mname, ph)]] <- data.frame(
        sample = rownames(X_h), actual = y, predicted = pred_loocv,
        model = mname, phenotype = ph
      )
      fit  <- models[[mname]]$train(X_h, y)
      pred <- models[[mname]]$predict(fit, X_v)
      scores_list[[paste(mname, ph)]] <- data.frame(
        sample = rownames(X_v), score = pred,
        model = mname, phenotype = ph
      )
    }
  }
  df_loocv  <- bind_rows(loocv_list)
  df_scores <- bind_rows(scores_list) %>%
    left_join(meta_v %>% dplyr::rename(sample = sampleName), by = "sample")
  
  # ---- Plot helpers ----
  loocv_plot <- function(model_name) {
    d <- df_loocv %>% filter(model == model_name)
    ann <- d %>% group_by(phenotype) %>%
      summarise(r   = cor(actual, predicted, method = "pearson"),
                rho = cor(actual, predicted, method = "spearman"),
                .groups = "drop") %>%
      mutate(label = sprintf("r = %.2f\n\u03C1 = %.2f", r, rho))
    ggplot(d, aes(actual, predicted)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 3, alpha = 0.7, color = "steelblue") +
      geom_smooth(method = "lm", se = FALSE, color = "indianred", linewidth = 0.8) +
      geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
                hjust = -0.1, vjust = 1.2, size = 4.2, inherit.aes = FALSE) +
      facet_wrap(~ phenotype, scales = "free", ncol = 2) +
      labs(title = sprintf("LOOCV on Govaere %s \u2014 %s model", space_label, model_name),
           x = "Actual", y = "Predicted") +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  scores_plot <- function(model_name) {
    d <- df_scores %>% filter(model == model_name)
    ggplot(d, aes(x = condition, y = score, fill = condition)) +
      geom_jitter(width = 0.12, height = 0, size = 3.5, shape = 21, color = "black") +
      stat_summary(fun = mean, geom = "crossbar",
                   width = 0.5, color = "grey20", linewidth = 0.4) +
      stat_compare_means(comparisons = list(c('AlphaBeta','Beta'),c("T2D", "AlphaBeta")),
                         method = "t.test",
                         method.args = list(
                           alternative = "less",
                           var.equal = FALSE
                         ))+
      facet_wrap(~ phenotype, scales = "free_y", ncol = 2) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = sprintf("RNAval predicted scores (%s) \u2014 %s model", space_label, model_name),
           x = "Condition", y = "Predicted score") +
      theme_bw(base_size = 14) +
      theme(plot.title  = element_text(hjust = 0.5),
            legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1))
  }
  
  for (mname in names(models)) {
    ggsave(file.path(fig_dir, sprintf("loocv_%s_%s.png", space_label, mname)),
           loocv_plot(mname),
           width = 10, height = 5,   dpi = 600, units = "in")
    ggsave(file.path(fig_dir, sprintf("rnaval_scores_%s_%s.png", space_label, mname)),
           scores_plot(mname),
           width = 10, height = 5.5, dpi = 600, units = "in")
  }
  
  invisible(list(loocv = df_loocv, scores = df_scores))
}
### Load LIV2TRANS extra LVs and TCs----------
dataset_names <- c("Govaere", "Kostrzewski","Wang", "Feaver",'Hoang')
ref_dataset <- "Govaere"
target_dataset <- "Kostrzewski"
# Load previously found extra basis
Wm_opt <- readRDS(paste0('../results/Wm_',target_dataset,'_extra.rds'))
Wm_TC <- readRDS(paste0('../results/Wm_',target_dataset,'_combo.rds'))

### Load pre-processed RNAseq results---------
gene_expression <- read.delim('../RNAseq_validation/2C72ND-expression-matrix.tsv')
gene_expression <- gene_expression %>% dplyr::select('gene_id','gene_name','gene_biotype',
                                              "X2C72ND_10_cpm","X2C72ND_11_cpm","X2C72ND_12_cpm",
                                              "X2C72ND_13_cpm", "X2C72ND_14_cpm" , "X2C72ND_15_cpm" ,
                                              "X2C72ND_1_cpm" ,"X2C72ND_2_cpm","X2C72ND_3_cpm",
                                              "X2C72ND_4_cpm","X2C72ND_5_cpm","X2C72ND_6_cpm",
                                              "X2C72ND_7_cpm","X2C72ND_8_cpm","X2C72ND_9_cpm")
colnames(gene_expression) <- c('gene_id','gene_name','gene_biotype',
                               'BetaDMSO1','BetaDMSO2','BetaDMSO3',
                               'BetaM1','BetaM2','BetaM3',
                               'T2D1','T2D2','T2D3',
                               'Beta1','Beta2','Beta3',
                               'AlphaBeta1','AlphaBeta2','AlphaBeta3')
metadata <- data.frame(sampleID=c('BetaDMSO1','BetaDMSO2','BetaDMSO3',
                                  'BetaM1','BetaM2','BetaM3',
                                  'T2D1','T2D2','T2D3',
                                  'Beta1','Beta2','Beta3',
                                  'AlphaBeta1','AlphaBeta2','AlphaBeta3'),
                       condition = c('BetaDMSO','BetaDMSO','BetaDMSO',
                                     'BetaM','BetaM','BetaM',
                                     'T2D','T2D','T2D',
                                     'Beta','Beta','Beta',
                                     'AlphaBeta','AlphaBeta','AlphaBeta'))
gene_expression <- gene_expression %>% mutate(gene_name = ifelse(gene_name=="",gene_id,gene_name))
gene_expression <- gene_expression %>% filter(gene_name %in% rownames(Wm_opt))
# inds <- which(apply(gene_expression[,4:18], 1,sum)>0)
# gene_expression <- gene_expression[inds,]
# keep <- log2(matrixStats::rowVars(as.matrix(gene_expression[,4:18]))) > -5
# gene_expression <- gene_expression[keep,]
genes <- list(gene_expression$gene_name)
gene_expression <- gene_expression %>% dplyr::select(-gene_id,-gene_biotype,-gene_name)
gene_expression <- aggregate(gene_expression,by = genes,sum)
Xval <- log2(gene_expression %>% column_to_rownames('Group.1')+1)
Xval <- t(Xval)
Xval <- Xval - colMeans(Xval)

pca_diagnostic <- prcomp(Xval,scale. = FALSE,center = TRUE)
factoextra::fviz_screeplot(pca_diagnostic,ncp=20)

df_pca <- as.data.frame(pca_diagnostic$x[,1:2])
df_pca <- left_join(df_pca %>% rownames_to_column('sampleID'),metadata)
ggplot(df_pca,aes(x=PC1,y=PC2,fill=condition))+
  geom_point(size=2.5,shape=21,stroke =1.2) +
  theme_bw()

### Comparison with W_opt
Wm_opt <- Wm_opt[colnames(Xval),]
Wm_TC <- Wm_TC[colnames(Xval),]
# missing_genes <- setdiff(rownames(Wm_opt), colnames(Xval))
# if (length(missing_genes) > 0) {
#   zeros <- matrix(
#     0,
#     nrow = nrow(Xval),
#     ncol = length(missing_genes),
#     dimnames = list(rownames(Xval), missing_genes)
#   )
#
#   Xval <- cbind(Xval, zeros)
# }
# Xval <- Xval[, rownames(Wm_opt), drop = FALSE]

# Z <- Xval %*% Wm_opt
Z <- Xval %*% cbind(Wm_TC,Wm_opt)[,c(1,3)]
colnames(Z) <- c('TC1','extraLV1')
df_proj <- left_join(as.data.frame(Z) %>% rownames_to_column('sampleID'),metadata)
p_opt <- ggplot(df_proj,
                     aes(x = TC1, y = extraLV1,
                         fill  = condition)) +
  geom_point(size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Valdiation RNAseq projected into Extra LVs",
       x    = "TC1",
       y    = "Extra LV1") +
  theme_bw()
print(p_opt)
ggsave('../figures/RNAval_mpsTC1_mpsExtraLV1_cpmed.png',p_opt,
       width = 7.5, height = 5.5, dpi = 600, units = "in")

Xh <- readRDS("../results/processed_data_list_govaere_kostrzewski.rds")$Xh
Yh <- readRDS("../results/processed_data_list_govaere_kostrzewski.rds")$Yh
# Genes common to both sides AND present in the basis
shared <- Reduce(intersect, list(colnames(Xh), colnames(Xval), rownames(Wm_opt)))
res_opt <- run_phenotype_models(
  X_h         = Xh[,   shared] %*% Wm_opt[shared, ],
  X_v         = Xval[, shared] %*% Wm_opt[shared, ],
  Y           = Yh,
  meta_v      = metadata %>% dplyr::rename(sampleName=sampleID),
  space_label = "extraLV",
  fig_dir     = fig_dir
)

### Load pre-processed RNAseq and perform differential gene expression analysis---------------
gene_expression <- read.delim('../RNAseq_validation/2C72ND-expression-matrix.tsv')
gene_expression <- gene_expression %>% dplyr::select('gene_id','gene_name','gene_biotype',
                                              "X2C72ND_10_count","X2C72ND_11_count","X2C72ND_12_count",
                                              "X2C72ND_13_count", "X2C72ND_14_count" , "X2C72ND_15_count" ,
                                              "X2C72ND_1_count" ,"X2C72ND_2_count","X2C72ND_3_count",
                                              "X2C72ND_4_count","X2C72ND_5_count","X2C72ND_6_count",
                                              "X2C72ND_7_count","X2C72ND_8_count","X2C72ND_9_count")
colnames(gene_expression) <- c('gene_id','gene_name','gene_biotype',
                               'BetaDMSO1','BetaDMSO2','BetaDMSO3',
                               'BetaM1','BetaM2','BetaM3',
                               'T2D1','T2D2','T2D3',
                               'Beta1','Beta2','Beta3',
                               'AlphaBeta1','AlphaBeta2','AlphaBeta3')
metadata <- data.frame(sampleID=c('BetaDMSO1','BetaDMSO2','BetaDMSO3',
                                  'BetaM1','BetaM2','BetaM3',
                                  'T2D1','T2D2','T2D3',
                                  'Beta1','Beta2','Beta3',
                                  'AlphaBeta1','AlphaBeta2','AlphaBeta3'),
                       condition = c('BetaDMSO','BetaDMSO','BetaDMSO',
                                     'BetaM','BetaM','BetaM',
                                     'T2D','T2D','T2D',
                                     'Beta','Beta','Beta',
                                     'AlphaBeta','AlphaBeta','AlphaBeta'))
gene_expression <- gene_expression %>% mutate(gene_name = ifelse(gene_name=="",gene_id,gene_name))
# gene_expression[,4:18] <- apply(gene_expression[,4:18],c(1,2),round)
exp_factors <- unique(metadata$condition)

# Build count matrix (genes x samples), with gene_id as rownames
# Sum raw counts across gene_ids that share a gene_name
counts_by_name <- gene_expression %>%
  dplyr::select(-gene_id, -gene_biotype) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum), .groups = "drop")
RNAval_dataset <- list(
  counts = counts_by_name %>% column_to_rownames('gene_name'),
  metadata = metadata %>% 
    dplyr::select(sampleName = sampleID, condition),
  genes = counts_by_name$gene_name,
  exp_factors = exp_factors
)
# save(RNAval_dataset, file = "../RNAseq_validation/RNAval_dataset.RData")

counts_by_name[,2:16] <- apply(counts_by_name[,2:16],c(1,2),round)

# Build a gene_name-level annotation table.
# If multiple gene_ids share a name, collapse their ids/biotypes into one string
gene_annot <- gene_expression %>%
  group_by(gene_name) %>%
  summarise(
    gene_ids      = paste(sort(unique(gene_id)),       collapse = ";"),
    gene_biotypes = paste(sort(unique(gene_biotype)),  collapse = ";"),
    n_ids         = n_distinct(gene_id),
    .groups       = "drop"
  )

# Sanity check — how many gene_names had >1 gene_id?
table(gene_annot$n_ids)

# Build count matrix with gene_name as rownames
count_matrix <- as.matrix(counts_by_name[, -1])
rownames(count_matrix) <- counts_by_name$gene_name
storage.mode(count_matrix) <- "integer"

# --- DESeq2 setup ---
rownames(metadata) <- metadata$sampleID
metadata <- metadata[colnames(count_matrix), ]

metadata$condition <- factor(metadata$condition,
                             levels = c("T2D", "BetaDMSO", "BetaM", "Beta", "AlphaBeta"))

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData   = metadata,
                              design    = ~ condition)

# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep, ]

dds <- DESeq(dds)

# Results helper, now keyed by gene_name
get_results <- function(dds, condition_name, annot, ref = "T2D", alpha = 0.05) {
  res <- results(dds,
                 contrast = c("condition", condition_name, ref),
                 alpha    = alpha)
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_name") %>%
    left_join(annot, by = "gene_name") %>%
    arrange(padj)
}
# Run all four contrasts vs T2D
res_BetaDMSO_vs_T2D  <- get_results(dds, "BetaDMSO",  gene_annot)  %>%
  filter(gene_name %in% rownames(Wm_opt)) %>% filter(!is.na(stat))
res_BetaM_vs_T2D     <- get_results(dds, "BetaM",     gene_annot) %>%
  filter(gene_name %in% rownames(Wm_opt)) %>% filter(!is.na(stat))
res_Beta_vs_T2D      <- get_results(dds, "Beta",      gene_annot) %>%
  filter(gene_name %in% rownames(Wm_opt)) %>% filter(!is.na(stat))
res_AlphaBeta_vs_T2D <- get_results(dds, "AlphaBeta", gene_annot) %>%
  filter(gene_name %in% rownames(Wm_opt)) %>% filter(!is.na(stat))


#### Re-run contrasts and directly infer pathway activity and geneset enrichment comparing conditions------
## Contrast IFNA+TGFB vs T2d, Contrast TGFB vs T2d, and Contrast IFNA+TGFB vs TGFB
AlphaBeta_vs_T2D  <- get_results(dds, "AlphaBeta", gene_annot)                %>% filter(!is.na(stat))
Beta_vs_T2D       <- get_results(dds, "Beta",      gene_annot)                %>% filter(!is.na(stat))
BetaM_vs_T2D <- get_results(dds, "BetaM",      gene_annot)                %>% filter(!is.na(stat))
AlphaBeta_vs_Beta <- get_results(dds, "AlphaBeta", gene_annot, ref = "Beta")  %>% filter(!is.na(stat))
BetaM_vs_Beta <- get_results(dds, "BetaM",      gene_annot, ref = "Beta")                %>% filter(!is.na(stat))
contrasts_list <- list(
  `IFNA+TGFB vs T2D`  = AlphaBeta_vs_T2D,
  `TGFB vs T2D`       = Beta_vs_T2D,
  `BetaM vs T2D`      = BetaM_vs_T2D,
  `IFNA+TGFB vs TGFB` = AlphaBeta_vs_Beta,
  `BetaM vs TGFB`     = BetaM_vs_Beta
)
common_genes <- Reduce(intersect, lapply(contrasts_list, `[[`, "gene_name"))
stat_mat <- sapply(contrasts_list, function(df) {
  df$stat[match(common_genes, df$gene_name)]
})
rownames(stat_mat) <- common_genes

net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
pwy_act <- decoupleR::run_viper(stat_mat, net_prog,minsize = 1,verbose = TRUE) %>% dplyr::select(-statistic)

# progeny returns a (contrasts × pathways) matrix
df_pwy <- pwy_act %>%
  dplyr::rename(contrast = condition,
                Pathway = source)

pathway_order <- c("JAK-STAT","Hypoxia","EGFR","WNT","MAPK","NFkB","Androgen",
                   "VEGF","TNFa","Trail","PI3K","Estrogen","TGFb","p53")
lim_pwy <- max(abs(df_pwy$score), na.rm = TRUE) * 1.15

p_pwy_contrasts <- ggplot(
  df_pwy %>% mutate(Pathway  = factor(Pathway,  levels = pathway_order),
                    contrast = factor(contrast, levels = names(contrasts_list))),
  aes(x = score, y = Pathway, fill = score)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = 0,
                       limits = c(-lim_pwy, lim_pwy)) +
  scale_x_continuous(n.breaks = 6, limits = c(-lim_pwy, lim_pwy)) +
  facet_wrap(~ contrast, ncol = 3) +
  labs(title = "PROGENy pathway activity by contrast",
       x = "Pathway activity (z-score)", y = "Pathway") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# Hallmark enrichment per contrast
ids  <- mapIds(org.Hs.eg.db, keys = rownames(stat_mat),
               column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
keep <- which(!is.na(ids))
meas <- stat_mat[keep, , drop = FALSE]
rownames(meas) <- ids[keep]

hm <- fastenrichment(colnames(meas), rownames(meas), meas,
                     enrichment_space = "msig_db_h",
                     n_permutations   = 10000,
                     order_columns    = FALSE)

nes_mat  <- as.matrix(hm$NES$`NES MSIG Hallmark`);   colnames(nes_mat)  <- colnames(meas)
pval_mat <- as.matrix(hm$Pval$`Pval MSIG Hallmark`); colnames(pval_mat) <- colnames(meas)

df_hm <- left_join(
  as.data.frame(nes_mat)  %>% rownames_to_column("Hallmark") %>%
    pivot_longer(-Hallmark, names_to = "contrast", values_to = "NES"),
  as.data.frame(pval_mat) %>% rownames_to_column("Hallmark") %>%
    pivot_longer(-Hallmark, names_to = "contrast", values_to = "padj"),
  by = c("Hallmark", "contrast")) %>%
  mutate(Hallmark = sub("^FL1000_MSIG_H_HALLMARK_", "", Hallmark),
         Hallmark = str_replace_all(Hallmark, "_", " "),
         Hallmark = paste0(toupper(substr(tolower(Hallmark), 1, 1)),
                           substring(tolower(Hallmark), 2)))

# per-facet ordering, with top-15 by |NES| fallback if nothing passes padj <= 0.05
sig_label <- function(p) {
  ifelse(is.na(p),    "",
         ifelse(p <= 0.0001, "****",
                ifelse(p <= 0.001,  "***",
                       ifelse(p <= 0.01,   "**",
                              ifelse(p <= 0.05,   "*", "")))))
}
df_hm_plot <- df_hm %>%
  group_by(contrast) %>%
  group_modify(~{
    sig <- .x %>% filter(padj <= 0.05)
    if (nrow(sig) == 0) .x %>% slice_max(abs(NES), n = 15) else sig
  }) %>% ungroup() %>%
  mutate(contrast = factor(contrast, levels = names(contrasts_list))) %>%
  arrange(contrast, NES) %>%
  mutate(HM_ord = paste(Hallmark, contrast, sep = "___"),
         HM_ord = factor(HM_ord, levels = unique(HM_ord)))
offset <- max(abs(df_hm_plot$NES), na.rm = TRUE) * 1.15 * 0.04

p_hm_contrasts <- ggplot(df_hm_plot, aes(x = NES, y = HM_ord, fill = NES)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = 0) +
  geom_text(aes(label = sig_label(padj),
                x = ifelse(NES < 0, NES - offset, NES + offset)),
            size = 5, color = "black") +
  facet_wrap(~ contrast, ncol = 3, scales = "free_y") +
  labs(title = "Hallmark enrichment by contrast",
       x = "Normalized Enrichment Score", y = "Hallmark") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave(file.path(fig_dir, "rnaval_pathway_contrasts.png"),  p_pwy_contrasts,
       width = 14, height = 10, dpi = 600, units = "in")
ggsave(file.path(fig_dir, "rnaval_hallmark_contrasts.png"), p_hm_contrasts,
       width = 16, height = 14, dpi = 600, units = "in")

contrasts_list_volcano <- lapply(contrasts_list, function(df) {
  df %>% dplyr::select(gene_name, log2FoldChange, padj)
})

df_volcano <- bind_rows(
  lapply(names(contrasts_list_volcano), function(nm) {
    contrasts_list_volcano[[nm]] %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      mutate(contrast = nm)
  })
)

# Direction classification
lfc_thresh  <- 1
padj_thresh <- 0.05
df_volcano <- df_volcano %>%
  mutate(direction = case_when(
    padj <= padj_thresh & log2FoldChange >  lfc_thresh ~ "Up",
    padj <= padj_thresh & log2FoldChange < -lfc_thresh ~ "Down",
    TRUE                                                ~ "ns"
  ),
  contrast = factor(contrast, levels = names(contrasts_list_volcano)))

# Top 10 most significant genes per contrast for labeling
df_top_labels <- df_volcano %>%
  filter(direction != "ns") %>%
  group_by(contrast) %>%
  slice_min(padj, n = 10, with_ties = FALSE) %>%
  ungroup()

p_volcano <- ggplot(df_volcano, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = direction), alpha = 0.55, size = 1.2) +
  scale_color_manual(values = c(Up = "indianred", Down = "steelblue", ns = "grey75")) +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype  = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = df_top_labels,
                           aes(label = gene_name),
                           size = 3, max.overlaps = 20,
                           box.padding = 0.3, segment.size = 0.3) +
  facet_wrap(~ contrast, ncol = 3, scales = "free") +
  labs(title = "Volcano plots by contrast",
       x     = expression(Log[2]~"fold change"),
       y     = expression(-Log[10]~"adjusted p-value"),
       color = "Direction") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(hjust = 0.5),
        legend.position  = "bottom")

ggsave(file.path(fig_dir, "rnaval_volcano_contrasts.png"), p_volcano,
       width = 16, height = 10, dpi = 600, units = "in")


### project on MPS TCs and extraLVs
all_genes_intesection <- intersect(intersect(intersect(res_BetaDMSO_vs_T2D$gene_name,
                                   res_BetaM_vs_T2D$gene_name),
                                   res_Beta_vs_T2D$gene_name),
                                   res_AlphaBeta_vs_T2D$gene_name)
res_BetaDMSO_vs_T2D  <- res_BetaDMSO_vs_T2D  %>%
  filter(gene_name %in% all_genes_intesection) %>% column_to_rownames('gene_name')
res_BetaM_vs_T2D     <- res_BetaM_vs_T2D%>%
  filter(gene_name %in%  all_genes_intesection) %>% column_to_rownames('gene_name')
res_Beta_vs_T2D      <- res_Beta_vs_T2D %>%
  filter(gene_name %in% all_genes_intesection) %>% column_to_rownames('gene_name')
res_AlphaBeta_vs_T2D <- res_AlphaBeta_vs_T2D %>%
  filter(gene_name %in% all_genes_intesection) %>% column_to_rownames('gene_name')
res_BetaDMSO_vs_T2D <- res_BetaDMSO_vs_T2D[all_genes_intesection,]
res_BetaM_vs_T2D <- res_BetaM_vs_T2D[all_genes_intesection,]
res_Beta_vs_T2D <- res_Beta_vs_T2D[all_genes_intesection,]
res_AlphaBeta_vs_T2D <- res_AlphaBeta_vs_T2D[all_genes_intesection,]
res_all <- data.frame(BetaDMSO=res_BetaDMSO_vs_T2D$stat,
                      BetaM = res_BetaM_vs_T2D$stat,
                      Beta=res_Beta_vs_T2D$stat,
                      AlphaBeta=res_AlphaBeta_vs_T2D$stat)
rownames(res_all) <- all_genes_intesection
res_all <- as.matrix(res_all)
res_all <- res_all - rowMeans(res_all)
## Compare with Wm_opt
W <- Wm_opt[all_genes_intesection,]
colnames(W) <- c('extraLV1','extraLV2')
sim <- lsa::cosine(cbind(res_all,W))
sim[!upper.tri(sim)] <- -100
sim <- as.data.frame(sim) %>% rownames_to_column("condition1") %>% gather(key='condition2',value='cosine_sim',-condition1)
sim <- sim %>% filter(cosine_sim>(-100))

Z <- t(res_all) %*% W
# colnames(Z) <- c('TC1','extraLV1')
df_proj <- as.data.frame(Z) %>% rownames_to_column('condition')
p_opt <- ggplot(df_proj,
                aes(x = extraLV1, y = extraLV2,
                    fill  = condition)) +
  geom_point(size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Valdiation RNAseq projected into Extra LVs using Wald statistic",
       x    = "Extra LV1",
       y    = "Extra LV2") +
  theme_bw()
print(p_opt)
ggsave('../figures/RNAval_in_MPS_extraLVs.png',p_opt,
       width = 7.5, height = 5.5, dpi = 600, units = "in")


Wtc <- Wm_TC[all_genes_intesection,]
Z <- t(res_all) %*% Wtc
colnames(Z) <- c('TC1','TC2')
df_proj_tc <- as.data.frame(Z) %>% rownames_to_column('condition')
p_tc <- ggplot(df_proj_tc,
                aes(x = TC1, y = TC2,
                    fill  = condition)) +
  geom_point(size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Valdiation RNAseq projected into TCs using Wald statistic",
       x    = "TC1",
       y    = "TC2") +
  theme_bw()
print(p_tc)
ggsave('../figures/RNAval_in_MPS_TCs.png',p_tc,
       width = 7.5, height = 5.5, dpi = 600, units = "in")

## Re-run LIV2TRANS---------------
# RNAval is already a packaged list:
load("../RNAseq_validation/RNAval_dataset.RData")    # -> object `RNAval_dataset`

# Govaere is the load_datasets-style file (separate `data` and `metadata`):
load("../data/GSE135251_Govaere_dataset.RData")      # -> objects `data`, `metadata`
Govaere <- list(counts   = data,
                metadata = metadata,
                genes    = rownames(data))
rm(data, metadata)
load("../data/GSE13090_Hoang_dataset.RData")      # -> objects `data`, `metadata`
Hoang <- list(counts   = data,
                metadata = metadata,
                genes    = rownames(data))
rm(data, metadata)
load("../data/GSE168285_Kostrzewski_dataset.RData")      # -> objects `data`, `metadata`
Kostrzewski <- list(counts   = data,
                metadata = metadata,
                genes    = rownames(data))
rm(data, metadata)

load("../data/GSE166256_Wang_dataset.RData")      # -> objects `data`, `metadata`
Wang <- list(counts   = data,
                    metadata = metadata,
                    genes    = rownames(data))
rm(data, metadata)

load("../data/GSE89063_Feaver_dataset.RData")      # -> objects `data`, `metadata`
Feaver <- list(counts   = data,
             metadata = metadata,
             genes    = rownames(data))
rm(data, metadata)

# -----------------------------------------------------------------------
# 1. Make both counts matrices numeric with gene-symbol rownames
# -----------------------------------------------------------------------

prep_counts <- function(counts, genes) {
  m <- as.matrix(counts)
  storage.mode(m) <- "numeric"
  rownames(m) <- genes
  m
}
Govaere_counts <- prep_counts(Govaere$counts,         Govaere$genes)
RNAval_counts  <- prep_counts(RNAval_dataset$counts,  RNAval_dataset$genes)

# -----------------------------------------------------------------------
# 2. Common-gene intersection (the only "filtering" process_datasets does
#    when filter_variance = FALSE)
# -----------------------------------------------------------------------

genes_common <- intersect(rownames(Govaere_counts), rownames(RNAval_counts))
genes_common <- intersect(genes_common,Hoang$genes)
genes_common <- intersect(genes_common,Kostrzewski$genes)
genes_common <- intersect(genes_common,Wang$genes)
genes_common <- intersect(genes_common,Feaver$genes)
cat("Common genes:", length(genes_common), "\n")

Govaere_counts <- Govaere_counts[genes_common, ]
RNAval_counts  <- RNAval_counts[genes_common, ]

# -----------------------------------------------------------------------
# 3. Helpers that mirror pca_center() exactly
# -----------------------------------------------------------------------

log2cpm    <- function(m) apply(m, 2, function(x) log2(1 + x / sum(x) * 1e6))
row_center <- function(m) m - rowMeans(m)

# -----------------------------------------------------------------------
# 4. Govaere (in-vivo reference): log2(1+CPM) + per-gene centering -> Xh
#    Phenotype matrix Yh from its metadata
# -----------------------------------------------------------------------

Govaere_lcpm   <- log2cpm(Govaere_counts)
Govaere_center <- row_center(Govaere_lcpm)
Xh <- t(Govaere_center)                 # samples x genes

# Check the column names first -- adjust if your Govaere metadata uses different ones
# colnames(Govaere$metadata)
Yh <- as.matrix(Govaere$metadata %>% dplyr::select(nas_score, Fibrosis_stage))
colnames(Yh) <- c("NAS", "fibrosis")

# -----------------------------------------------------------------------
# 5. RNAval (in-vitro target): same log2(1+CPM) + centering -> Xm
# -----------------------------------------------------------------------

RNAval_lcpm   <- log2cpm(RNAval_counts)
RNAval_center <- row_center(RNAval_lcpm)
Xm <- t(RNAval_center)                  # samples x genes

# -----------------------------------------------------------------------
# 6. Average per condition -> re-center -> PCA -> Wm
# -----------------------------------------------------------------------

meta <- RNAval_dataset$metadata
# align metadata to the columns of the centered matrix
meta <- meta[match(colnames(RNAval_center), meta$sampleName), ]
stopifnot(identical(meta$sampleName, colnames(RNAval_center)))

groups        <- meta$condition
groups_unique <- unique(groups)

# group means
data_grouped <- sapply(groups_unique, function(g) {
  rowMeans(RNAval_center[, groups == g, drop = FALSE])
})
rownames(data_grouped) <- rownames(RNAval_center)
colnames(data_grouped) <- groups_unique

# re-center (no log this time, matches pca_center(..., log_cpm = FALSE, center = TRUE))
data_grouped_c <- row_center(data_grouped)

# PCA on samples x genes
pca_group <- prcomp(t(data_grouped_c), scale. = FALSE)

# Wm = rotation matrix without the final PC (variance zero after centering K groups)
Wm <- pca_group$rotation[, -ncol(pca_group$rotation)]
liv2trans_results <- liv2trans_run(Xh, Yh, Xm, Wm)
Wm_opt <- liv2trans_results$W_opt
Wm_combo <- liv2trans_results$W_translatable

# calculate effective capture of Wm
compute_alpha <- function(W, phi) {
  genes_common <- intersect(rownames(W), rownames(phi))
  if (length(genes_common) == 0) {
    stop("No common genes between W and phi.")
  }
  W_c   <- W[genes_common, , drop = FALSE]
  phi_c <- phi[genes_common, , drop = FALSE]
  alpha <- sapply(1:ncol(phi_c), function(j) {
    phi_j <- phi_c[, j]
    proj  <- as.vector(t(W_c) %*% phi_j)
    proj^2 / sum(phi_j^2)
  })
  if (is.null(dim(alpha))) alpha <- matrix(alpha, ncol = ncol(phi_c))
  rownames(alpha) <- if (is.null(colnames(W))) seq_len(ncol(W)) else colnames(W)
  colnames(alpha) <- colnames(phi_c)
  alpha
}
effective_rank <- function(alpha_vec) {
  s1 <- sum(alpha_vec)
  s2 <- sum(alpha_vec^2)
  if (s2 < .Machine$double.eps) return(NA_real_)
  (s1^2) / s2
}
concentration <- function(alpha_vec) {
  r <- effective_rank(alpha_vec)
  k <- length(alpha_vec)
  if (is.na(r) || k <= 1) return(NA_real_)
  1 - (r / k)
}
# phi_PLS from real PLSR model
model_real <- liv2trans_results$model
Wh_real    <- liv2trans_results$Wh
Bh_real    <- t(model_real@weightMN) %*% model_real@coefficientMN
phi_real   <- Wh_real %*% Bh_real          # p x q
colnames(phi_real) <- colnames(Yh)
k_mps   <- ncol(Wm)
p_genes <- nrow(Wm)
alpha_real    <- compute_alpha(Wm, phi_real)            # k x q
r_eff_real    <- apply(alpha_real, 2, effective_rank)
conc_real     <- apply(alpha_real, 2, concentration)
rho_real_vec  <- colSums(alpha_real)

shared_opt   <- intersect(colnames(Xm), rownames(Wm_opt))
shared_combo <- intersect(colnames(Xm), rownames(Wm_combo))

Zm_opt   <- Xm[, shared_opt]   %*% Wm_opt[shared_opt, ]
Zm_combo <- Xm[, shared_combo] %*% Wm_combo[shared_combo, ]
colnames(Zm_opt)   <- c("extraLV1", "extraLV2")
colnames(Zm_combo) <- c("TC1",      "TC2")

df_rnaval_tc <- as.data.frame(Zm_combo) %>%
  rownames_to_column("sample") %>%
  left_join(meta %>% dplyr::select(c('sample' = 'sampleName'),condition), by = "sample")

p_rnaval_tc <- ggplot(df_rnaval_tc, aes(TC1, TC2, fill = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 4, shape = 21, color = "black") +
  geom_text_repel(aes(label = sample), size = 3, color = "black", max.overlaps = Inf) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "RNAval projected onto Translatable Components",
       x = "TC1", y = "TC2", fill = "Condition") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
print(p_rnaval_tc)
ggsave(file.path(fig_dir, "rnaval_TC_projection.png"), p_rnaval_tc,
       width = 7.5, height = 5.5, dpi = 600, units = "in")

# Pathway-activity barplots (PROGENy) for TCs and Extra LVs
pathway_order <- c("JAK-STAT","Hypoxia","EGFR","WNT","MAPK","NFkB","Androgen",
                   "VEGF","TNFa","Trail","PI3K","Estrogen","TGFb","p53")

pathway_barplot <- function(df, title, axis_label) {
  # df columns: Pathway, score, p_value, condition (V1 / V2)
  # lim <- max(abs(df$score), na.rm = TRUE) * 1.15
  lim <- 25
  df  <- df %>%
    mutate(Pathway   = factor(Pathway, levels = pathway_order),
           condition = factor(condition,
                              levels = c("V1","V2"),
                              labels = paste(axis_label, c(1, 2))))
  ggplot(df, aes(x = score, y = Pathway, fill = score)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred",
                         mid = "whitesmoke", midpoint = 0,
                         limits = c(-lim, lim)) +
    scale_x_continuous(n.breaks = 8, limits = c(-lim, lim)) +
    geom_text(aes(
      label = ifelse(p_value <= 0.0001, "****",
                     ifelse(p_value <= 0.001,  "***",
                            ifelse(p_value <= 0.01,   "**",
                                   ifelse(p_value <= 0.05,   "*",
                                          ifelse(p_value <= 0.1,    "\u2219", "ns"))))),
      x = ifelse(score < 0, score - lim*0.03, score + lim*0.03)),
      size = 4, color = "black", angle = 90) +
    facet_wrap(~ condition, ncol = 2) +
    labs(title = title, x = "Activity", y = "Pathway") +
    theme_minimal(base_size = 16) +
    theme(plot.title  = element_text(hjust = 0.5),
          legend.position = "right")
}

pwy_extra <- pathway_activity_interpretation(Wm_opt,   Wm, plotting = FALSE)[[2]]
pwy_tc    <- pathway_activity_interpretation(Wm_combo, Wm, plotting = FALSE)[[2]]

p_pwy_extra <- pathway_barplot(pwy_extra, "Pathway activity — Extra LVs", "Extra LV")
p_pwy_tc    <- pathway_barplot(pwy_tc,    "Pathway activity — Translatable Components", "TC")

ggsave(file.path(fig_dir, "rnaval_pathway_extraLV.png"), p_pwy_extra,
       width = 11, height = 6.5, dpi = 600, units = "in")
ggsave(file.path(fig_dir, "rnaval_pathway_TC.png"), p_pwy_tc,
       width = 11, height = 6.5, dpi = 600, units = "in")

# Hallmark enrichment barplots for TCs and Extra LVs
hallmark_for_W <- function(W) {
  ids  <- mapIds(org.Hs.eg.db, keys = rownames(W),
                 column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  keep <- which(!is.na(ids))
  meas <- as.matrix(W[keep, , drop = FALSE]); rownames(meas) <- ids[keep]
  
  res <- fastenrichment(colnames(meas), rownames(meas), meas,
                        enrichment_space = "msig_db_h",
                        n_permutations   = 10000,
                        order_columns    = FALSE)
  
  nes  <- as.data.frame(res$NES$`NES MSIG Hallmark`)   %>% rownames_to_column("Hallmark") %>%
    gather("LV", "NES",  -Hallmark)
  pval <- as.data.frame(res$Pval$`Pval MSIG Hallmark`) %>% rownames_to_column("Hallmark") %>%
    gather("LV", "padj", -Hallmark)
  
  left_join(nes, pval, by = c("Hallmark", "LV")) %>%
    mutate(Hallmark = sub("^FL1000_MSIG_H_HALLMARK_", "", Hallmark),
           Hallmark = str_replace_all(Hallmark, "_", " "),
           Hallmark = paste0(toupper(substr(tolower(Hallmark),1,1)),
                             substring(tolower(Hallmark),2)))
}

plot_hallmark <- function(df, title, axis_label, padj_thresh = 0.1, fallback_n = 15) {
  # see what LV values actually exist
  lvs <- sort(unique(df$LV))
  message("LV values found: ", paste(lvs, collapse = ", "))
  
  # keep significant hallmarks per LV; if none, fall back to top |NES|
  df_plot <- df %>%
    group_by(LV) %>%
    group_modify(~{
      sig <- .x %>% filter(padj <= padj_thresh)
      if (nrow(sig) == 0) .x %>% slice_max(abs(NES), n = fallback_n) else sig
    }) %>%
    ungroup()
  
  # readable facet labels: <axis_label> 1, <axis_label> 2, ...
  df_plot <- df_plot %>%
    mutate(LV_label = factor(LV, levels = lvs,
                             labels = paste(axis_label, seq_along(lvs))),
           # unique per-facet y identifier (Hallmark + LV) so reorder works inside each panel
           HM_ord = paste(Hallmark, LV, sep = "___"))
  
  # per-facet ordering of HM_ord by NES, then strip the suffix for the y axis text
  df_plot <- df_plot %>%
    arrange(LV, NES) %>%
    mutate(HM_ord = factor(HM_ord, levels = unique(HM_ord)))
  
  ggplot(df_plot, aes(x = NES, y = HM_ord, fill = NES)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(limits = c(-2.7,2.7))+
    scale_fill_gradient2(low = "darkblue", high = "indianred",
                         mid = "whitesmoke", midpoint = 0) +
    facet_wrap(~ LV_label, ncol = 2, scales = "free_y") +
    labs(title = title, x = "Normalized Enrichment Score", y = "Hallmark") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

df_hm_extra <- hallmark_for_W(Wm_opt)
df_hm_tc    <- hallmark_for_W(Wm_combo)

p_hm_extra <- plot_hallmark(df_hm_extra, "Hallmark enrichment — Extra LVs", "Extra LV",padj_thresh=0.05)
p_hm_tc    <- plot_hallmark(df_hm_tc,    "Hallmark enrichment — Translatable Components", "TC",padj_thresh=0.05)

ggsave(file.path(fig_dir, "rnaval_hallmark_extraLV.png"), p_hm_extra,
       width = 13, height = 8, dpi = 600, units = "in")
ggsave(file.path(fig_dir, "rnaval_hallmark_TC.png"), p_hm_tc,
       width = 13, height = 8, dpi = 600, units = "in")


# Project HUMAN samples (Govaere) onto both spaces, coloured by NAS / fibrosis
shared_h_opt   <- intersect(colnames(Xh), rownames(Wm_opt))
shared_h_combo <- intersect(colnames(Xh), rownames(Wm_combo))

Zh_extra <- Xh[, shared_h_opt]   %*% Wm_opt[shared_h_opt, ]
Zh_tc    <- Xh[, shared_h_combo] %*% Wm_combo[shared_h_combo, ]
colnames(Zh_extra) <- c("extraLV1", "extraLV2")
colnames(Zh_tc)    <- c("TC1",      "TC2")

df_h_extra <- as.data.frame(Zh_extra) %>% rownames_to_column("sample") %>%
  mutate(NAS = Yh[, "NAS"], fibrosis = Yh[, "fibrosis"])
df_h_tc <- as.data.frame(Zh_tc) %>% rownames_to_column("sample") %>%
  mutate(NAS = Yh[, "NAS"], fibrosis = Yh[, "fibrosis"])

plot_human_proj <- function(df, x, y, color_col, title, xlab, ylab) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[color_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = title, x = xlab, y = ylab, color = color_col) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))
}

p_h_extra_fib <- plot_human_proj(df_h_extra, "extraLV1", "extraLV2", "fibrosis",
                                 "Govaere on Extra LVs — fibrosis", "Extra LV1", "Extra LV2")
p_h_extra_nas <- plot_human_proj(df_h_extra, "extraLV1", "extraLV2", "NAS",
                                 "Govaere on Extra LVs — NAS",     "Extra LV1", "Extra LV2")
p_h_tc_fib    <- plot_human_proj(df_h_tc,    "TC1",      "TC2",      "fibrosis",
                                 "Govaere on TCs — fibrosis",       "TC1", "TC2")
p_h_tc_nas    <- plot_human_proj(df_h_tc,    "TC1",      "TC2",      "NAS",
                                 "Govaere on TCs — NAS",            "TC1", "TC2")

ggsave(file.path(fig_dir, "govaere_extraLV_fibrosis.png"), p_h_extra_fib,
       width = 7, height = 5.5, dpi = 600, units = "in")
ggsave(file.path(fig_dir, "govaere_extraLV_NAS.png"),      p_h_extra_nas,
       width = 7, height = 5.5, dpi = 600, units = "in")
ggsave(file.path(fig_dir, "govaere_TC_fibrosis.png"),      p_h_tc_fib,
       width = 7, height = 5.5, dpi = 600, units = "in")
ggsave(file.path(fig_dir, "govaere_TC_NAS.png"),           p_h_tc_nas,
       width = 7, height = 5.5, dpi = 600, units = "in")

message("All figures saved to ", normalizePath(fig_dir))


# ===========================================================================
#  Predict NAS / fibrosis from TC scores
#    - Ordinal regression (MASS::polr)
#    - SVM regression  (e1071::svm)
#    - Linear regression (lm)
#    For each: LOOCV scatter on Govaere, then score RNAval samples by condition.
# ===========================================================================

library(MASS)
library(e1071)

phenos <- c("NAS", "fibrosis")

# Inputs (already in workspace):
#   Zh_tc    : Govaere TC scores  (samples x 2)  with colnames TC1, TC2
#   Zm_combo : RNAval TC scores   (samples x 2)  with colnames TC1, TC2
#   Yh       : Govaere phenotype matrix (samples x 2) with NAS, fibrosis
# Drop human samples missing a phenotype label
complete  <- complete.cases(Yh)
X_h       <- Zh_tc[complete, , drop = FALSE]
Yh_clean  <- Yh[complete,    , drop = FALSE]
X_v       <- Zm_combo

# ---- LOOCV harness ---------------------------------------------------------
loocv_predict <- function(X, y, train_fn, predict_fn) {
  preds <- numeric(length(y))
  for (i in seq_along(y)) {
    fit <- train_fn(X[-i, , drop = FALSE], y[-i])
    preds[i] <- predict_fn(fit, X[i, , drop = FALSE])
  }
  preds
}

# ---- Three model train/predict pairs --------------------------------------
# Ordinal: predict expected value across the ordered classes -> continuous-ish
train_ord <- function(X, y) {
  df <- as.data.frame(X); df$y <- ordered(y)
  polr(y ~ ., data = df, Hess = TRUE)
}
predict_ord <- function(fit, X_new) {
  p <- predict(fit, newdata = as.data.frame(X_new), type = "probs")
  if (is.null(dim(p))) p <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
  as.numeric(p %*% as.numeric(colnames(p)))
}

train_svm   <- function(X, y) svm(x = X, y = y, type = "eps-regression", kernel = "radial")
predict_svm <- function(fit, X_new) as.numeric(predict(fit, X_new))

train_lm    <- function(X, y) lm(y ~ ., data = data.frame(y = y, as.data.frame(X)))
predict_lm  <- function(fit, X_new) as.numeric(predict(fit, newdata = as.data.frame(X_new)))

models <- list(
  ordinal = list(train = train_ord, predict = predict_ord),
  svm     = list(train = train_svm, predict = predict_svm),
  linear  = list(train = train_lm,  predict = predict_lm)
)

# ---- Run all model x phenotype combos ------------------------------------
loocv_list  <- list()
scores_list <- list()

for (mname in names(models)) {
  for (ph in phenos) {
    y <- Yh_clean[, ph]
    
    # LOOCV on Govaere
    pred_loocv <- loocv_predict(X_h, y,
                                models[[mname]]$train,
                                models[[mname]]$predict)
    loocv_list[[paste(mname, ph)]] <- data.frame(
      sample    = rownames(X_h),
      actual    = y,
      predicted = pred_loocv,
      model     = mname,
      phenotype = ph
    )
    
    # Refit on full human data, score RNAval
    fit  <- models[[mname]]$train(X_h, y)
    pred <- models[[mname]]$predict(fit, X_v)
    scores_list[[paste(mname, ph)]] <- data.frame(
      sample    = rownames(X_v),
      score     = pred,
      model     = mname,
      phenotype = ph
    )
  }
}

df_loocv  <- bind_rows(loocv_list)
df_scores <- bind_rows(scores_list) %>%
  left_join(meta %>% dplyr::rename(sample = sampleName), by = "sample")

# ---- LOOCV scatter (one figure per model) --------------------------------
loocv_plot <- function(model_name) {
  d <- df_loocv %>% filter(model == model_name)
  ann <- d %>% group_by(phenotype) %>%
    summarise(r   = cor(actual, predicted, method = "pearson"),
              rho = cor(actual, predicted, method = "spearman"),
              .groups = "drop") %>%
    mutate(label = sprintf("r = %.2f\n\u03C1 = %.2f", r, rho))
  
  ggplot(d, aes(actual, predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "indianred", linewidth = 0.8) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2, size = 4.2, inherit.aes = FALSE) +
    facet_wrap(~ phenotype, scales = "free", ncol = 2) +
    labs(title = paste0("LOOCV on Govaere TCs \u2014 ", model_name, " model"),
         x = "Actual", y = "Predicted") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))
}

# ---- Dot plot of RNAval scores per condition (one figure per model) ------
scores_plot <- function(model_name) {
  d <- df_scores %>% filter(model == model_name)
  ggplot(d, aes(x = condition, y = score, fill = condition)) +
    geom_jitter(width = 0.12, height = 0, size = 3.5, shape = 21, color = "black") +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.5, color = "grey20", linewidth = 0.4) +
    stat_compare_means(comparisons = list(c('AlphaBeta','Beta'),c("T2D", "AlphaBeta")),
                       method = "t.test",
                       method.args = list(
                         alternative = "less",
                         var.equal = FALSE
                       ))+
    facet_wrap(~ phenotype, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste0("RNAval predicted scores \u2014 ", model_name, " model"),
         x = "Condition", y = "Predicted score") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- Save all six figures ------------------------------------------------
for (mname in names(models)) {
  ggsave(file.path(fig_dir, paste0("loocv_TC_",        mname, ".png")),
         loocv_plot(mname),
         width = 10, height = 5,   dpi = 600, units = "in")
  ggsave(file.path(fig_dir, paste0("rnaval_scores_TC_", mname, ".png")),
         scores_plot(mname),
         width = 10, height = 5.5, dpi = 600, units = "in")
}

# Optional: keep the underlying tables around for sanity-checking
saveRDS(df_loocv,  "../results/loocv_TC_predictions.rds")
saveRDS(df_scores, "../results/rnaval_TC_scores.rds")
