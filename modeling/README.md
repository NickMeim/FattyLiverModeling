## Machine learning and statistical analysis algorithms of the manuscript case study

## Scripts
1. functions_translation.R : Contains functions to load and pre-process the raw data (normalize, scale, center, etc), perform PCA, and identify translatable or extra latent components using an evolutionary algorithm. ** Additionally it contains the function to identify extra latent vectors to span the human space, using an analytical solution.**
2. enrichment_calculations.R : Function for performing Geneset Enrichment Analysis (GSEA) for GO Terms, KEGG pathways or TFs (treated as sets of genes by using the dorothea regulon) using the **fgsea** R library.
3. vector_space_interpretation.R: Contains functions to infer pathways and TFs activities based on a vector of gene "loadings"/"weights", and visualize the results.
4. CrossValidationUtilFunctions.R: Contains function to perform cross-fold validation of various approaches.
5. CrossValidationProcess.R: The script that performs all cross-fold validation and evaluation tests.
7. distance_scores.R: The function used to calculate GSEA-based distance between to list of molecular signatures.
8. InterClinicalValidation.R: Validate the performance in predicting fibrosis and NAS in two external clinical datasets, using the extra basis found with the Govaere clinical dataset and the Kostrzewski MPS liver model.
9. MLBenchmarksOfNasFibrosis.py: Python script to bechmark the performance of multiple other linear and non-linear ML models, as well as the loss of performance when backprojecting and the retrieval of performance when using the identified extra basis.
10. multipleMLmodelsEval.R: R script to load the results of the previous script and visualize them.
11. analyzescRNAseq_Gribben_etal.R: R script for processing single cell data from GSE202379 and map it to translatable and extra latent variables.
