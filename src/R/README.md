# Installation
First install dependencies by starting R and entering:
install.packages("BiocManager")
BiocManager::install(c("cmapR","rhdf5","dorothea"," org.Hs.eg.db ","hgu133a.db","AnnotationDbi","fgsea","topGO","EGSEAdata","GO.db"))
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("caret")
install.packages("patchwork")
install.packages("devtools")
devtools::install_github("thomasp85/patchwork")
BiocManager::install('saezlab/decoupleR')

Then start R and enter:
``` r
remotes::install_github("NickMeim/FattyLiverModeling", subdir = "src/R/LIVIVTRA")
```

# Run

# Example
Run the R_tutorial.R to see an example usage of the package, from training a PLSR model to predict Y from X, to identifying an extra basis and interpreting the results.
You can also open the example in RStudio.