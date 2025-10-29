# Systems biology framework for the rational design of operational conditions for in vitro / in vivo translation of tissue models
A computational system biology approach for identifying and optimizing in-vitro models that better recapitulate human patients, with a case study in MAFLD.
Github repository of the study:
> Systems biology framework for rational design of operational conditions for in vitro / in vivo translation of microphysiological systems <br>
> Jose L. Cadavid<sup>1+</sup>, Nikolaos Meimetis<sup>1+</sup>, Tyler Matsuzaki<sup>1</sup>, Erin N. Tevonian<sup>1</sup>, Avlant Nilsson<sup>1,2,3</sup>, Douglas A. Lauffenburger<sup>1*</sup>
> 1) Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139, USA
> 2) Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139, USA
> 3) Center for Gynepathology Research, Massachusetts Institute of Technology, Cambridge, MA 02139, USA
> * *Corresponding author, lauffen@mit.edu
> + +These authors share joint first co-authorship, and they contributed equally to this work.

doi: https://doi.org/10.1101/2025.01.17.633624

This repository is administered by @NickMeim. For questions contact both jcadavid@mit.edu and meimetis@mit.edu

## Plug-and-play scripts
In the `src` folder, there are scripts to run the approach as a package in R or Python. 

## Steps to reproduce the manuscript's analyses and results.
In each folder, the functionality of each script is explained.

To reproduce all the results of the manuscript, you need to run sequentially the following steps, **which will automatically save the results that are required by the scripts in the `figures` folder to make the manuscript's figures** :

0. (**Optional** as everything is already downloaded) First, download the metadata files for each dataset from GEO (and for GSE168285 and GSE166256, which was also their re-aligned data). The dataset from Pantano et al. is also included in the `data` folder and it was fetched manually from the ARCHS4 database. We already include them in the `data` folder, make sure you get them.
1. Get everything that is in the `ARCH4S_retrieval` folder, and run the `fetch_samples_ARCH4S.R` script in the R programming language to fetch all samples from ARCHS4. Specifically, for the GSE166256 dataset, which is not available on ARCHS4, we downloaded the raw fastq files from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166256) using the SRA Run Selector. Then we re-aligned them for consistency by uploading them to ARCHS4, and downloaded the generated count matrix.
2. Download everything in the `data` folder of this repository. Run `preprocess_datasets.R` in the `data` folder, which will generate the final .RData files that we use to load our data for all case studies.
3. This step is optional, only for the supplementary snRNA analysis. Make sure you are now using an environment with tidyverse and Seurat packages installed. In the `scRNA_Gribben_et_al_analysis` folder, run the `analysis_Gribben_etal.Rmd` Rmarkdown. This will process the data (the Seurat object needs to be downloaded from [GSE202379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202379)) and produce the results that are required for Supplementary Figure 11.  
4. Run the `get_MAS_Fibrosis_enrichment_score.R` R script.
5. Run `main_example_run.R`, as it will generate many results required for the figures, and it is the R script which shows line-by-line analysis of the whole approach.
6. Run `RunLIV2TRANS_for_datasets_pairs.R` to apply LIV2TRANS in different in vitro / in vivo dataset pairs and generate Supplementary Figure 10. **Important Note** : First you **need to install LIV2TRANS as an R package** using remotes::install_github("NickMeim/FattyLiverModeling", subdir = "src/R/LIV2Trans"). But first, make sure you have installed all dependencies as described in the `src/R` folder.
7. Go to the `modeling` folder and run the `CrossValidationProcess.R` R script. This will generate train-validation-test splits and also perform the evaluation of LIV2TRANS in various scenarios. It will generate many files that are required to make the figures of the manuscript. It will also plot some figures while it runs for diagnostic and monitoring purposes.
8. In the `modeling` folder, run `InterClinicalValidation.R` to evaluate the performance of the identified extra basis in other clinical datasets. It generates the results required for Supplementary figures and Figure 4.
9. In the `modeling` folder, run the `MLBenchmarksOfNasFibrosis.py` Python script, to benchmark many different machine learning models.
10. In the `modeling` folder, run `multipleMLmodelsEval.R` to make Supplementary Figure 5.
11. In the `optimization` folder, run `prepareData4Optimization.R` script to prepare the data to be used in pytorch for performing the general variance optimization approach, described in the manuscript.
12. In the `optimization` folder, open and run the `maximizeExplainedVariance.ipynb` jupyter notebook. It will run the optimization approach and generate results required for figure 6.
13. In the `optimization` folder, run the `postprocessOptimizedMPS.R` to generate results required for figure 6.
14. Go to the `figures` folder and run `MakeFigure2.R` **to make Figure 2**.
15. In the `figures` folder, run `MakeFigure3.R` **to make Figure 3**.
16. In the `figures` folder, run `MakeFigure4.R` **to make Figure 4**.
17. In the `figures` folder, run `MakeFigure5.R` **to make Figure 5**.
18. In the `figures` folder, run `MakeFigure6.R` **to make Figure 6**.

## Scripts in this folder 
`main_example_run.R`: An R script to run a line-by-line analysis of the whole approach, beginning by using a clinical dataset and an in-vitro dataset, and identifying and interpreting in the end an extra basis for better spanning the human phenotypes.
It produces results relevant to making the study's figures.

`get_MAS_Fibrosis_enrichment_score.R` : Supplementary script for gene signature enrichment using VIPER.

`RunLIV2TRANS_for_datasets_pairs.R` : Supplementary script for running LIV2TRANS for all available combos of in vivo / in vitro datasets.

## Folders' description
1. data : Folder where you should put your raw data. (Be careful! Larger files cannot be pushed into GitHub and can potentially crash it)
2. preprocessing : Pre-processed data along with algorithms for the pre-processing of the raw data. (Be careful! Larger files cannot be pushed into GitHub and can potentially crash it)
3. modeling : Machine learning and statistical analysis algorithms of the project.
4. optimization : Scripts for performing in-silico maximization of the human variance captured by the MPS datase and analyzing the results.
5. results : Storage of the results and figures of the studies.
6. utils: Auxiliary functions used by other scripts in this study.
7. src: Scripts to run the approach as a package in R or Python.
8. figures : Folder containing the scripts to reproduce the figures of the study. You can also store the produced figures there.
9. ARCH4S_retrieval : Folder containing scripts for fethcing data from ARCHS4 and the fecthed data.
10. scRNA_Gribben_et_al_analysis: Folder that contains the scripts and the data required to perform the supplementary snRNA analysis and thus make Supplementary Figure 11.
