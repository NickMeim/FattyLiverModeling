# FattyLiverModeling
A computational system biology approach for identifying and optimizing in-vitro models that better recapitulate human patients, with a case study in MAFLD.

## Plug-and-play scripts
In the `src` folder there are scripts to run the approach as a package in R or Python. 

## Run the manuscript's case study analysis.
`main_example_run.R`: An R script to run a line-by-line analysis of the whole approach, beginning by using a clinical dataset and an in-vitro dataset, and identifying and interpreting in the end an extra basis for better spaning the human phenotypes.
It produces results relevent for making the study's figures. After running you make go to the `figure` folder and run the corresponding scripts to re-produce the figures of the study.

`get_MAS_Fibrosis_enrichment_score.R` : Supplementary script for gene signature enrichment using VIPER.

`RunLIV2TRANS_for_datasets_pairs.R` : Supplementary script for running LIV2TRANS for all available combos of in vivo / in vitro datasets.

## Folders
1. data : Folder where you should put your raw data. (Be careful larger files cannot be pushed into GitHub and can potentially crash it)
2. preprocessing : Pre-processed data along with algorithms for the pre-processing of the raw data. (Be careful larger files cannot be pushed into GitHub and can potentially crash it)
3. modeling : Machine learning and statistical analysis algorithms of the project.
4. optimization : Scripts for performing in-silico maximization of the human variance captured by the MPS datase and analyzing the results.
5. results : Storage of the results and figures of the studies.
6. utils: Auxilary functions used by other scripts in this study.
7. src: Scripts to run the approach as a package in R or Python.
8. figures : Folder containing the scripts to re-produce the figures of the study. You can also store the produced figures there.
