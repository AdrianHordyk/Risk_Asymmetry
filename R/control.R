
# Load libraries
library(DLMtool); library(dplyr); library(ggplot2)
library(MSEtool); library(flextable); library(forcats);
library(grid); library(gridExtra)



# Directory names
OMdir <- "OMs"
Resultsdir <- "Results"

nsim <- 100

source("R/Functions.R")


# Create MP 
SCA_MP <- make_MP(SCA, HCR_MSY, fix_h=TRUE, fix_tau=TRUE)