This repository contains the required R code to produce all of the analyses and plots in:  
‘The Impact of Tip Age Distribution on Reconstructing Trait Evolution Using Phylogenetic Comparative Methods’  
William Gearty, Bethany J. Allen, Pedro L. Godoy, and Alfio Alessandro Chiarenza

## Before running the code
Before running this R script, we suggest that you download this folder and set it as your R working directory. Plot files will be saved within the "/figures" folder. Note that many of the steps of the script can be skipped by loading the included simulation data within the "/data" folder.

## Running the code
Run `R/analyses.R` to perform all analyses and produce all plots from the manuscript. Note that the analyses can take several days to run completely.

## Required R packages
The following packages can all be installed from CRAN:

ape  
deeptime  
dplyr  
forcats  
FossilSim  
future  
geiger  
ggh4x  
ggplot2  
mvMORPH  
pbapply  
pcmtools  
phytools  
TreeSim  
tibble  
tidyr  

  
The pcmtools package is currently only available on GitHub and will need to be installed from there to reproduce
these analyses and plots. This can be achieved by running the following commands in your R console (ignore the first line 
if you already have devtools installed).  
```r
install.packages("devtools")  
devtools::install_github("willgearty/pcmtools")
```

  
All other packages can be installed from CRAN. These scripts have been tested using R version 4.4.0 - 
Copyright (C) 2024 The R Foundation for Statistical Computing.
