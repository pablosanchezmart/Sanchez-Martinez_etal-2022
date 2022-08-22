#### PREDICTING FROM PHYLOGENY -------------------------------------------------------------------------------------------------------------------####

# P Sanchez-Marinez

remove(list = ls())
options(scipen = 999)
### Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

# install.packages("itertools", dependencies = T)
# install.packages("rJava", dependencies = T)
# if(!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtreeExtra")
# BiocManager::install("ggtree")
# }


library(dplyr)
library(plyr)
library(tidyverse)
library(caret)
library(stringr)

# library(MCMCglmm)
library(phytools)

# Parallelization
library(parallel)
library(foreach)

# Projections
library(raster)
library(sp)
library(sf)
library(fasterize)
library(SpaDES)

library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(tictoc)

library(ape)
library(missForest)
library(doParallel)

# phylogenetic plot
library(ggtree)
library(egg)
library(ggtreeExtra)
library(ggstance)
library(ggnewscale)
library(ggpubr)

# logit model with r2
library(rms)
# to plot it
library(jtools)

library(emmeans)
# maxent
library(rJava)
options(java.parameters = "-Xmx1g" ) # for a higher memory usage
# only run if the maxent.jar file is available, in the right folder
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
library(dismo)

library(randomForest)
require(caTools)

#discretize
library(arules)

# SEM
library(lavaan)
# library(semPlot)
library(OpenMx)
library(GGally)
library(corrplot)

# Piecewise sem
library(piecewiseSEM)

library(ade4)

# Workspace

print("loading workspace ...")

options(scipen = 999, digits = 3)

nClusters <- 4

hammond_xft_data <- T
hammond_xft_data_P88 <- F


## Directories

print("preparing gneral data and modelling specifications...")

if(isTRUE(hammond_xft_data)){
  #### All available data (Hammond) or just previous version of XFT? ####
  processed.data <- "data_processed/hammond_xt_data/"
  output.dir <- "outputs/hammond_xt_data/"
  rslts.dir <- "results/hammond_xt_data/"
  print("Using Hammond et al. xft data (not published)")
} else{
  if(isTRUE(hammond_xft_data_P88)){
    processed.data <- "data_processed/hammond_xt_data_P88/"
    output.dir <- "outputs/hammond_xt_data_P88/"
    rslts.dir <- "results/hammond_xt_data_P88/"
    print("Using Hammond et al. xft data (not published) with P88 values for angiosperms")
  } else {
  processed.data <- "data_processed/"
  output.dir <- "outputs/"
  rslts.dir <- "results/"
  print("Using Sanchez-Martinez et al. 2020 data")
}
}

# Main projection for results
behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs") # Projection

#### Data ------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021.csv"), header = T)

# load(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021_georeferenced.RData")) # georeferenced, object called te (tree evolution)

length(ht$Species) # 47540

a <- ht %>% filter(!is.na(MinWP_md))
length(a$Species) # 536 (685)
b <- ht %>% filter(!is.na(P50))
length(b$Species) # 863 (1376)

# Genus lvl tree
gen.tr <- read.tree("data/phylogeny/treevol_phylogeny_genus_lvl_2020.tre")

# Species lvl trees
# sp.tr <- read.tree("data/phylogeny/treevol_phylogeny_species_lvl_2020.tre") # From Svenning group
# sp.tr$node.label <- NULL

sp.tr <- read.tree("data/phylogeny/treevol_phylogeny_species_lvl_2021.tre") # v.phylomaker scenario 1
sp.tr$node.label <- NULL

# Random forest dataset

rf.df <- ht
rf.df$group <- as.factor(rf.df$group)
rf.df$leaf.form <- as.factor(rf.df$leaf.form)

# Observed data

obs.rf.df <- rf.df %>% filter(!is.na(P50) | !is.na(MinWP_md))
length(obs.rf.df$genus) # 1079

#### Model parameters -------------------------------------------------------------------------------------------------------------------------------------------------------- ####

CVnIterations <- 100
impnIterations <- 100

cat("Number of clusters used when parallelizing:", nClusters,
    "\nN Iterations for crossVal: ", CVnIterations, 
    "\nIterations for imputation: ", impnIterations, "\n")


#### Models for cross validation --------------------------------------------------------------------------------------------------------------------------------------------- ####

RFmodels <- data.frame(rbind(
  # Climate
  cbind("P50_climate_3", "P50", toString(paste0("PC", 1:3))),
  cbind("MinWP_md_climate_3", "MinWP_md", toString(paste0("PC", 1:3))),
  cbind("bi_climate_3", "MinWP_md, P50", toString(paste0("PC", 1:3))),
  
  cbind("P50_climate_10", "P50", toString(paste0("PC", 1:10))),
  cbind("MinWP_md_climate_10", "MinWP_md", toString(paste0("PC", 1:10))),
  cbind("bi_climate_10", "MinWP_md, P50", toString(paste0("PC", 1:10))),
  
  # Phylogeny
  cbind("P50_phylo_5", "P50", toString(c(paste0("Phylo_axis_", 1:5)))),
  cbind("MinWP_md_phylo_5", "MinWP_md", toString(c(paste0("Phylo_axis_", 1:5)))),
  cbind("bi_phylo_5", "MinWP_md, P50", toString(c(paste0("Phylo_axis_", 1:5)))),
  
  cbind("P50_phylo_20", "P50", toString(c(paste0("Phylo_axis_", 1:20)))),
  cbind("MinWP_md_phylo_20", "MinWP_md", toString(c(paste0("Phylo_axis_", 1:20)))),
  cbind("bi_phylo_20", "MinWP_md, P50", toString(c(paste0("Phylo_axis_", 1:20)))),
  
  cbind("P50_phylo_50", "P50", toString(c(paste0("Phylo_axis_", 1:50)))),
  cbind("MinWP_md_phylo_50", "MinWP_md", toString(c(paste0("Phylo_axis_", 1:50)))),
  cbind("bi_phylo_50", "MinWP_md, P50", toString(c(paste0("Phylo_axis_", 1:50)))),
  
  cbind("P50_phylo_100", "P50", toString(c(paste0("Phylo_axis_", 1:100)))),
  cbind("MinWP_md_phylo_100", "MinWP_md", toString(c(paste0("Phylo_axis_", 1:100)))),
  cbind("bi_phylo_100", "MinWP_md, P50", toString(c(paste0("Phylo_axis_", 1:100)))),
  
  # Combination including group
  
  cbind("P50_climate_3_phylo_5_group", "P50", toString(c("group", paste0("PC", 1:3), paste0("Phylo_axis_", 1:5)))),
  cbind("MinWP_md_climate_3_phylo_5_group", "MinWP_md", toString(c("group", paste0("PC", 1:3), paste0("Phylo_axis_", 1:5)))),
  cbind("bi_climate_3_phylo_5_group", "MinWP_md, P50", toString(c("group", paste0("PC", 1:3), paste0("Phylo_axis_", 1:5)))),
  
  cbind("P50_climate_5_phylo_5_group", "P50", toString(c("group", paste0("PC", 1:5), paste0("Phylo_axis_", 1:5)))),
  cbind("MinWP_md_climate_5_phylo_5_group", "MinWP_md", toString(c("group", paste0("PC", 1:5), paste0("Phylo_axis_", 1:5)))),
  cbind("bi_climate_5_phylo_5_group", "MinWP_md, P50", toString(c("group", paste0("PC", 1:5), paste0("Phylo_axis_", 1:5)))),
  
  # cbind("P50_climate_5_phylo_10_group_lf", "P50", toString(c("group", "leaf.form", paste0("PC", 1:5), paste0("Phylo_axis_", 1:10)))),
  # cbind("MinWP_md_climate_5_phylo_10_group_lf", "MinWP_md", toString(c("group", "leaf.form", paste0("PC", 1:5), paste0("Phylo_axis_", 1:10)))),
  # cbind("bi_climate_5_phylo_10_group_lf", "P50, MinWP_md", toString(c("group", "leaf.form", paste0("PC", 1:5), paste0("Phylo_axis_", 1:10)))),
  
  cbind("P50_climate_10_phylo_10_group", "P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:10)))),
  cbind("MinWP_md_climate_10_phylo_10_group", "MinWP_md", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:10)))),
  cbind("bi_climate_10_phylo_10_group", "MinWP_md, P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:10)))),
  
  cbind("P50_climate_10_phylo_20_group", "P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:20)))),
  cbind("MinWP_md_climate_10_phylo_20_group", "MinWP_md", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:20)))),
  cbind("bi_climate_10_phylo_20_group", "MinWP_md, P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:20)))),
  
  # cbind("P50_climate_10_phylo_20_group_lf", "P50", toString(c("group", "leaf.form", paste0("PC", 1:10), paste0("Phylo_axis_", 1:20)))),
  # cbind("MinWP_md_climate_10_phylo_20_group_lf", "MinWP_md", toString(c("group", "leaf.form", paste0("PC", 1:10), paste0("Phylo_axis_", 1:20)))),
  # cbind("bi_climate_5_phylo_20_3448_group_lf", "P50, MinWP_md", toString(c("group", "leaf.form", paste0("PC", 1:5), paste0("Phylo_axis_", 1:10), paste0("Phylo_axis_", 90:100)))),
  
  cbind("P50_climate_10_phylo_50_group", "P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:50)))),
  cbind("MinWP_md_climate_10_phylo_50_group", "MinWP_md", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:50)))),
  cbind("bi_climate_10_phylo_50_group", "MinWP_md, P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:50)))),
  
  cbind("P50_climate_10_phylo_100_group", "P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:100)))),
  cbind("MinWP_md_climate_10_phylo_100_group", "MinWP_md", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:100)))),
  cbind("bi_climate_10_phylo_100_group", "MinWP_md, P50", toString(c("group", paste0("PC", 1:10), paste0("Phylo_axis_", 1:100))))
  
  # cbind("P50_climate_10_phylo_100_group_lf", "P50", toString(c("group", "leaf.form", paste0("PC", 1:10), paste0("Phylo_axis_", 1:100)))),
  # cbind("MinWP_md_climate_10_phylo_100_group_lf", "MinWP_md", toString(c("group", "leaf.form", paste0("PC", 1:10), paste0("Phylo_axis_", 1:100)))),
  # cbind("bi_climate_10_phylo_100_group_lf", "P50, MinWP_md", toString(c("group", "leaf.form", paste0("PC", 1:10), paste0("Phylo_axis_", 1:100)))),
  
  # cbind("bi_climate_5_phylo_5_group_MinWP_soil_SX", "P50, MinWP_md", toString(c("group", "MinWP_soil_SX", paste0("PC", 1:5), paste0("Phylo_axis_", 1:5))))
), stringsAsFactors = F)
names(RFmodels) <- c("Model", "Response","Predictors")

