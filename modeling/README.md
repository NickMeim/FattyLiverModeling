## Machine learning and statistical analysis algorithms of the project

## Scripts
1. functions_translation_nikos.R: Functions used for translation, using TransCompR-based models, as modified and used by Nikos.
2. functions_translation_jose.R : Original functions used in translation provided by Jose.
3. aux_functions.R : Auxilary functions used in multiple scripts.
4. DeX_analysis.R : Differential gene expression analysis of in-silico perturbed data using Limma.
5. multiPhenoModeling.R : TransCompR-based translation using **standard PCA** for multiple human phenotypes from Hoang et al. 2019, Scientific reports.
6. SparsePCAMultiPheno.R : TransCompR-based translation using **sparse PCA** for multiple human phenotypes from Hoang et al. 2019, Scientific reports.
7. PCLoadings_analysis.R : Derive insights using the gene loading of **PCs** of interest, as identified from multiPhenoModeling.R.
8. sparsePCLoadings_analysis.R : Derive insights using the gene loading of **sparse PCs** of interest, as identified from SparsePCAMultiPheno.R.
9. enrichment_calculations.R : Function for performing Geneset Enrichment Analysis (GSEA) for GO Terms, KEGG pathways or TFs (treated as sets of genes by using the dorothea regulon) using the **fgsea** R library.
