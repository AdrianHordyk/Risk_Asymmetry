
loadPackages <- function(pcks=NUL) {
  for (pk in pcks) {
    Chk <- require(pk, character.only=TRUE)
  }
  if (!Chk) {
    install.packages(pk)
    require(pk, character.only=TRUE)
  }
}

pcks <- c('DLMtool', 'MSEtool', 'dplyr', 'ggplot2', 
          'flextable', 'forcats', 'grid', 'gridExtra')

# Load libraries
loadPackages(pcks)

# Directory names
OMdir <- "OMs"
Resultsdir <- "Results"

nsim <- 100

source("R/Functions.R")

# Create MP 
SCA_MP <- make_MP(SCA, HCR_MSY, fix_h=TRUE, fix_tau=TRUE)