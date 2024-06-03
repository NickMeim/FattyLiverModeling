## Scripts for performing in-silico maximization of the human variance captured by the MPS datase and analyzing the results.
1. maximizeExplainedVariance.ipynb : Jupyter notebook to perform the maximization of the variance, perturb in-silico the gene expression and save the results.
2. DeX_analysis.R : Differential gene expression analysis of in-silico perturbed data using Limma.
3. enrichment_calculations.R : Function to perform gene set enrichment for different types of genesets: KEGG pathways, GO Terms, MSIG genesets, even TFs but by using GSEA.
4. functions_translation.R: The same as the one in the "modeling" folder.
5. vector_space_interpretation.R: The same as the one in the "modeling" folder.
6. prepareData4Optimization.R: Prepare all clinical and in vitro datasets to be ready for use in the optimization algorithm.
