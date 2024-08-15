## Scripts for performing in-silico maximization of the human variance captured by the MPS datase and analyzing the results.
1. maximizeExplainedVariance.ipynb : Jupyter notebook to perform the maximization of the variance of Govaere clinical dataset captured by the MPS data, and save the indetified differenatial gene expression vector.
2. maximizeExplainedVarianceMultiClinical.ipynb : Jupyter notebook to perform the maximization of the variance of Govaere clinical dataset captured by the MPS data, and save the indetified differenatial gene expression vector.
3. DeX_analysis.R : Differential gene expression analysis of in-silico perturbed data using Limma.
4. enrichment_calculations.R : Function to perform gene set enrichment for different types of genesets: KEGG pathways, GO Terms, MSIG genesets, even TFs but by using GSEA.
5. functions_translation.R: The same as the one in the "modeling" folder.
6. vector_space_interpretation.R: The same as the one in the "modeling" folder.
7. prepareData4Optimization.R: Prepare all clinical and in vitro datasets to be ready for use in the optimization algorithm.

## Files used.
1. All the 'X' files contain the normalized and centered gene expression data for the dataset that is denoted by a name (e.g. X_Govaere.csv).
2. All the 'Y' files contain phenotypic data of clinical datasets which are denoted a name (e.g. Y_Govaere.csv).
