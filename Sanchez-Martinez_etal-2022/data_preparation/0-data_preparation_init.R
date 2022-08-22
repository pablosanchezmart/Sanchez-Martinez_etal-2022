#### DATA PREPARATION INITIALITATION -------------------------------------------------------------------------------------- ####

remove(list = ls())

#### PACKAGES -------------------------------------------------------------------------------- ####

# install.packages("V.PhyloMaker", dependencies = T)
library(dplyr)
library(plyr)
library(stringr)
library(humidity)

library(ggpubr)

library(stats)
library(factoextra)

# install.packages("Taxonstand", dependencies = T)
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("wcornwell/taxonlookup")
# install.packages("taxize", dependencies = T)

library(Taxonstand) # step 1: clean up nomenclature
#https://cran.r-project.org/web/packages/Taxonstand/Taxonstand.pdf
#http://viktoriawagner.weebly.com/blog/cleaning-species-names-with-r-i-taxonstand
# install.packages("devtools") # if necessary
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup) # step 2: look it up
#https://github.com/traitecoevo/taxonlookup/blob/master/README.md
#library(phyndr) # step 3: swap species in tree for others with traits
#https://github.com/traitecoevo/phyndr-ms/blob/master/vignette/vignette.md

library(taxize) # OPTION 2: clean nomenclature
#http://viktoriawagner.weebly.com/blog/cleaning-species-names-with-r-ii-taxize

library(ape)
library(sp)
library(rgeos)


library(rgdal)
library(stars)
library(raster)
library(snow)
library(medfate)
library(foreach)
library(parallel)
library(doParallel)

# devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)

library(sf)

library(fasterize)

#### DATASET USED ----------------------------------------------------------------------------------- ####

hammond_xft_data <- T
hammond_xft_data_P88 <- F

if(isTRUE(hammond_xft_data)){
  #### All available data (Hammond) or just previous version of XFT? ####
  processed.data <- "data_processed/hammond_xt_data/"
  output.dir <- "outputs/hammond_xt_data/"
  print("Using Hammond et al. xft data (not published)")
} else{
  if(isTRUE(hammond_xft_data_P88)){
    #### All available data (Hammond) or just previous version of XFT? ####
    processed.data <- "data_processed/hammond_xt_data_P88/"
    output.dir <- "outputs/hammond_xt_data_P88/"
    print("Using Hammond et al. xft data (not published) with P88 values for angiosperms")
  } else{
    processed.data <- "data_processed/"
    output.dir <- "outputs/"
    print("Using Sanchez-Martinez et al. 2020 data")
  }
}
