#### Process figure files
#### This script opens up ggplot files created during the analysis pipeline
#### and homogenizes their style and fixes their dimensions for the final figure

root_dir <- "E:/Jose Luis/Documents/GitHub/FattyLiverModeling"
setwd(root_dir)

# Load plotting functions that include custom themes
source("./utils/plotting_functions.R")

