##### MAPPING DATA PREPARATION RF ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #####

remove(list = ls())

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")


#### FUNCTIONS ---------------------------------------------------------------------------------------------------------------------- ####

# Extract variable continuous values for each species

extrRastFun <- function(spAlpha = new_all_species_polygons[i, ], st = tile, tile = T){
  e <- extent(spAlpha)
  if(!is.null(intersect(e, st))){                # Ensure that we only work with species present in the tile
    # if(tile == F){    # Crop by the extent is useful when using large maps, if not, gives problems due to extent overlap, so better not to crop when using tiles
    #   r <- crop(st, e)  
    # }
    mask.r <- mask(x = st, mask = spAlpha)
    if(isTRUE(cellStats(mask.r$MinWP_soil_SX, "countNA") == length(mask.r$MinWP_soil_SX))){
      mask.df <- as.data.frame(coordinates(mask.r))
      mask.df[, names(st)] <- NA
    } else {
      mask.df <- as.data.frame(as(mask.r, "SpatialPixelsDataFrame"))
      mask.df <- mask.df[, c("x", "y", names(st))]
    }
    mask.lst <- list(mask.df)
    names(mask.lst) <- as.character(spAlpha$sp)
    return(mask.lst)
  } else {
    print(paste0("tile with no species: ", n))
  }
}

# HSM calculation for species (continuous)
continuousHSMFun <- function(minWP_soil.lst = sp.var.lst, species_traits.df = RFimp.df){
  minWP_soil.lst <- minWP_soil.lst[names(minWP_soil.lst) %in% species_traits.df$sp]
  # HSM calculation per pixel ...
  for(i in 1:length(minWP_soil.lst)){
    print(i)
    sp.name <- names(minWP_soil.lst[i])
    # P50 assignation
    minWP_soil.lst[[i]]$P50 <- species_traits.df[species_traits.df$sp == sp.name, "P50"] 
    minWP_soil.lst[[i]]$HSM  <- minWP_soil.lst[[i]]$MinWP_soil_SX - minWP_soil.lst[[i]]$P50
    minWP_soil.lst[[i]]$Species  <- sp.name
    minWP_soil.lst[[i]][which(is.na(minWP_soil.lst[[i]]$MinWP_soil_SX)), c("Species", "P50")] <- NA  # To not count for species where there is no data
  }
  return(minWP_soil.lst)
}

# HSM calculation for species (mean per species)
meanHSMFun <- function(MinWP_s = MinWP_soil_SX, spAlpha = all_species_traits){
  # HSM calculation per pixel ...
  for(sp in spAlpha$sp){
    # mean MinWP
    sp.pol <- spAlpha[spAlpha$sp == sp, ]
    e <- extent(sp.pol)
    r <- crop(MinWP_s, e)
    mask.r <- raster::mask(x = r, mask = sp.pol)
    varValues <- cellStats(x = mask.r, stat = mean, na.rm = T)
    # P50 assignation
    hsms_values <- varValues - spAlpha[spAlpha$sp == sp, ]$P50
    spAlpha[spAlpha$sp == sp, "HSMs_SaxtonRawls"] <- mean(hsms_values)
  }
  return(spAlpha)
}


#### DATA  -------------------------------------------------------------------------------------------------------------------------- ####

# Load predictions
RFimp.df <- read.csv(paste0(output.dir, "RF_imputations/RFimputation_100_bi_climate_5_phylo_5_group.csv"), header = T)
RFimp.df$sp <- str_replace(RFimp.df$Species, " ", "_")

length(RFimp.df$Species) # 45289 (45475)

summary(RFimp.df$MinWP_md)
RFimp.df[which(RFimp.df$MinWP_md == min(RFimp.df$MinWP_md)), c("Species", "family")]
RFimp.df[which(RFimp.df$MinWP_md == max(RFimp.df$MinWP_md)), c("Species", "family")]

summary(RFimp.df$P50)
RFimp.df[which(RFimp.df$P50 == min(RFimp.df$P50)), c("Species", "family")]
RFimp.df[which(RFimp.df$P50 == max(RFimp.df$P50)), c("Species", "family")]

# Load predictions observed data only
RFimp_obs.df <- read.csv(paste0(output.dir, "RF_imputations/RFimputation_obsData_100_bi_climate_5_phylo_5_group.csv"), header = T)
length(RFimp_obs.df$Species) # 1051 (1469)
RFimp_obs.df$sp <- str_replace(RFimp_obs.df$Species, " ", "_")

# Load predictions using different subset of species each time (80% each time) (uncertainity checking)

RFimp_unc.df <- read.csv("outputs/RF_imputations/RFimputation_uncertainity_80_100_bi_climate_5_phylo_5_group.csv", header =  T)
length(RFimp_unc.df$Species)
RFimp_unc.df$sp <- str_replace(RFimp_unc.df$Species, " ", "_")

# Load predictions using only phylogenetic data (uncertainity checking)
RFimp_phylo.df <- read.csv("outputs/RF_imputations/RFimputation_100_bi_phylo_20.csv", header = T)
summary(RFimp_phylo.df)
length(RFimp_phylo.df$Species)
RFimp_phylo.df$sp <- str_replace(RFimp_phylo.df$Species, " ", "_")

# Load predictions using only climatic data (uncertainity checking)

RFimp_clim.df <- read.csv("outputs/RF_imputations/RFimputation_100_bi_climate_10.csv", header = T)
summary(RFimp_clim.df)
length(RFimp_clim.df$Species)
RFimp_clim.df$sp <- str_replace(RFimp_clim.df$Species, " ", "_")

# Load species distributions
load(file = paste0(processed.data, "species_distributions/all_species_2020_sf.RData"))  # Species alpha-hulls With Behrmann projection, object called new_all_sopecies_polygons
new_all_species_polygons

# MinWP_soil data
MinWP_soil_SX <- raster(paste0("data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif"))
varNames <- names(MinWP_soil_SX)


#### GEOGRAPHICALLY REFERENCED OBSERVED SPECIES MEANS ------------------------------------------------------------------------------- ####

obs_species <- new_all_species_polygons[new_all_species_polygons$sp %in% obs.rf.df$sp, ]
length(obs_species$sp)# 792 (1204) species with observed data + distribution data

obs_species_traits <- merge(obs_species, RFimp_obs.df, by = "sp", all.x = T)

### Calculate HSMs from mean soil minimum water potential ###

obs_species_traits$HSMs <- obs_species_traits$MinWP_soil_SX - obs_species_traits$P50

# Save observed values referenced geographically
save(obs_species_traits, file = paste0(output.dir, "/species_traits/obs_species_traits.RData"))
print(paste0(output.dir, "/species_traits/obs_species_traits.RData"))


#### GEOGRAPHICALLY REFERENCED PREDICTED SPECIES MEANS ------------------------------------------------------------------------------- ####

all_species <- new_all_species_polygons[(new_all_species_polygons$sp %in% RFimp.df$sp), ]
length(all_species$sp) # 44892 (44679) species with observed and predicted data

all_species_traits <- merge(all_species, RFimp.df, by = "sp", all.x = T)
length(all_species_traits$sp)

### Calculate HSMs from mean soil minimum water potential ###

all_species_traits$HSMs <- all_species_traits$MinWP_soil_SX - all_species_traits$P50

# Best model
save(all_species_traits, file = paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))
print(paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))


#### GEOGRAPHICALLY REFERENCED MODEL UNCERTAINITY (USING DIFFERENT SUBSETS OF SPECIES, 80% EACH TIME) ------------------------------------------------ ####

all_species_uncert <- new_all_species_polygons[(new_all_species_polygons$sp %in% RFimp_unc.df$sp), ]
length(all_species_uncert$sp) # 44753 species with observed and predicted data

all_species_traits <- merge(all_species_uncert, RFimp_unc.df, by = "sp", all.x = T)

save(all_species_traits, file = "outputs/species_traits/all_species_traits_uncertainity_RFimp.RData")
print("outputs/species_traits/all_species_traits_uncertainity_RFimp.RData")


#### GEOGRAPHICALLY REFERENCE PREDICTED SPECIES MEANS (ONLY PHYLOGENETIC PREDICTORS) ------------------------------------------------ ####

all_species_phylo <- new_all_species_polygons[(new_all_species_polygons$sp %in% RFimp_phylo.df$sp), ]
length(all_species_phylo$sp) # 44753 species with observed and predicted data

all_species_traits <- merge(all_species_phylo, RFimp_phylo.df, by = "sp", all.x = T)

# Only phylogeny
save(all_species_traits, file = "outputs/species_traits/all_species_traits_phylo20_RFimp.RData")
print("outputs/species_traits/all_species_traits_phylo20_RFimp.RData")


#### GEOGRAPHICALLY REFERENCE PREDICTED SPECIES MEANS (ONLY CLIMATIC PREDICTORS) ---------------------------------------------------- ####

all_species_clim <- new_all_species_polygons[(new_all_species_polygons$sp %in% RFimp_clim.df$sp), ]
length(all_species_clim$sp) # 44753 species with observed and predicted data

all_species_traits <- merge(all_species_clim, RFimp_clim.df, by = "sp", all.x = T)

# Only climate
save(all_species_traits, file = "outputs/species_traits/all_species_traits_climate10_RFimp.RData")
print("outputs/species_traits/all_species_traits_climate10_RFimp.RData")


#### HSM CALCULATION USING CONTINUOUS SOIL MINIMUM WATER POTENTIAL (takes time, do it once, NOT USED FOR NOW) ----------------------------------------------------- ####

if(!file.exists(paste0(output.dir, "/HSMs_min/HSMs_min.tif"))){
  
extent(MinWP_soil_SX) <- extent(new_all_species_polygons)
tiles <- splitRaster(MinWP_soil_SX, nx = 10, ny = 10)

# Ymin per species and per tile

for(n in 1:length(tiles)){
  print(paste0("tile ", n, " out of ", length(tiles)))
  tile <- tiles[[n]]
  # plot(tile)
  # Extract minWP_soil data per species (continuous data)
  cl <- parallel::makeCluster(nCluster)
  registerDoParallel(cl)
  parallel::clusterExport(cl, varlist=c("new_all_species_polygons", "extrRastFun", "extent", "crop", "mask", "tile", "cellStats", "coordinates", "intersect"))
  
  # Extracting ...
  sp.var.lst <- list()
  sp.var.lst <- foreach(i=1:length(new_all_species_polygons$sp), .combine= "c") %dopar% {
    extrRastFun(spAlpha = new_all_species_polygons[i, ], st = tile)
  }
  stopCluster(cl)
  # Save R object for further use
  save(sp.var.lst, file = paste0("outputs/HSMs_tiles/all_species_raster_MinWP_soil_SX_", n, ".RData"))
  print(paste0("outputs/HSMs_tiles/all_species_raster_MinWP_soil_SX_", n, ".RData"))
  remove(sp.var.lst)
  gc()
  removeTmpFiles()
}

### HSM max, min and mean calculation ####

HSMs_mean_tiles <- list()
HSMs_min_tiles <- list()
HSMs_max_tiles <- list()

files <- list.files("outputs/HSMs_tiles", pattern = ".RData", full.names = T)
for(n in 1:length(files)){  # number 23 gives problems (too heavy), so it should be calculated manually (first mean, then min and finally max)
  file <- files[which(parse_number(files) == n)]
  load(file = file)
  ### HSM calculation
  all_species_hsms <- continuousHSMFun(minWP_soil.lst = sp.var.lst, species_traits.df = RFimp.df)
  all_species_hsms <- dplyr::bind_rows(all_species_hsms)
  all_species_hsms <- st_as_sf(all_species_hsms, coords = c("x", "y"), crs = behrmann)
  gc()
  removeTmpFiles()
  
  # Rasterize (mean) sf
  
  HSMs_mean <-  rasterize(all_species_hsms, tiles[[parse_number(file)]], field = "HSM", fun = mean)
  writeRaster(HSMs_mean, paste0("outputs/HSMs_mean/hsms_mean_", n, ".tif"), overwrite = T)
  print(paste0("outputs/HSMs_mean/hsms_mean_", n, ".tif"))
  remove(HSMs_mean)
  gc()
  removeTmpFiles()
  
  # Rasterize (min) sf
  
  HSMs_min <- rasterize(all_species_hsms, tiles[[parse_number(file)]], field = "HSM", fun = min)
  writeRaster(HSMs_min, paste0("outputs/HSMs_min/hsms_min_", n, ".tif"), overwrite = T)
  print(paste0("outputs/HSMs_tiles/hsms_min_", n, ".tif"))
  # HSMs_min_tiles <- c(HSMs_min_tiles, list(HSMs_min))
  remove(HSMs_min)
  gc()
  removeTmpFiles()
  
  # Rasterize (max) sf
  
  HSMs_max <- rasterize(all_species_hsms, tiles[[parse_number(file)]], field = "HSM", fun = max)
  writeRaster(HSMs_max, paste0("outputs/HSMs_max/hsms_max_", n, ".tif"), overwrite = T)
  print(paste0("outputs/HSMs_max/hsms_max_", n, ".tif"))
  # HSMs_max_tiles <- c(HSMs_max_tiles, list(HSMs_max))
  remove(HSMs_max)
  gc()
  removeTmpFiles()
}

# Merge tiles to get the final result

# Mean
mean_files <- list.files(path = "outputs/HSMs_mean/", pattern = "\\d+.tif", full.names = T)
HSMs_mean_tiles <- lapply(mean_files, FUN = raster)
HSMs_mean <- mergeRaster(HSMs_mean_tiles)
writeRaster(HSMs_mean, "outputs/HSMs_mean/HSMs_mean.tif", overwrite = T)

# Max
max_files <- list.files(path = "outputs/HSMs_max/", pattern = "\\d+.tif", full.names = T)
HSMs_max_tiles <- lapply(max_files, FUN = raster)
HSMs_max <- mergeRaster(HSMs_max_tiles)
writeRaster(HSMs_max, "outputs/HSMs_max/HSMs_max.tif", overwrite = T)

# Min
min_files <- list.files(path = "outputs/HSMs_min/", pattern = "\\d+.tif", full.names = T)
HSMs_min_tiles <- lapply(min_files, FUN = raster)
HSMs_min <- mergeRaster(HSMs_min_tiles)
writeRaster(HSMs_min, "outputs/HSMs_min/HSMs_min.tif", overwrite = T)


### Number of species with HSMs data per pixel ####

files <- list.files("outputs/HSMs_tiles", pattern = ".RData", full.names = T)
for(n in 1:length(files)){
  file <- files[which(parse_number(files) == n)]
  load(file = file)
  ### HSM calculation
  all_species_hsms <- continuousHSMFun(minWP_soil.lst = sp.var.lst, species_traits.df = RFimp.df)
  all_species_hsms <- dplyr::bind_rows(all_species_hsms)
  all_species_hsms <- st_as_sf(all_species_hsms, coords = c("x", "y"), crs = behrmann)
  gc()
  removeTmpFiles()
  
  # Rasterize (count) sf
  HSMs_count <-  rasterize(all_species_hsms, tiles[[parse_number(file)]], field = "HSM", fun = "count")
  
  writeRaster(HSMs_count, paste0("outputs/HSMs_count/hsms_count_", n, ".tif"), overwrite = T)
  print(paste0("outputs/HSMs_count/hsms_count_", n, ".tif"))
  
  remove(HSMs_count)
  gc()
  removeTmpFiles()
}

# Merge tiles to get the final result

# total HSMs species count
HSMs_count_files <- list.files(path = "outputs/HSMs_count/", pattern = "\\d+.tif", full.names = T)
HSMs_count_tiles <- lapply(HSMs_count_files, FUN = raster)

beginCluster(4)
HSMs_count <- mergeRaster(HSMs_count_tiles, fun = max)
endCluster()
HSMs_count[HSMs_count == 0] <- NA
writeRaster(HSMs_count, "outputs/HSMs_count/HSMs_count.tif", overwrite = T)
print("outputs/HSMs_count/HSMs_count.tif")


### Number of species with negative HSMs per pixel ####

files <- list.files("outputs/HSMs_tiles", pattern = ".RData", full.names = T)
for(n in 1:length(files)){
  file <- files[which(parse_number(files) == n)]
  load(file = file)
  ### HSM calculation
  all_species_neghsms <- continuousHSMFun(minWP_soil.lst = sp.var.lst, species_traits.df = RFimp.df)
  all_species_neghsms <- dplyr::bind_rows(all_species_neghsms)
  all_species_neghsms[which(all_species_neghsms$HSM > 0), "HSM"] <- NA
  all_species_neghsms <- st_as_sf(all_species_neghsms, coords = c("x", "y"), crs = behrmann)
  gc()
  removeTmpFiles()
  
  # Rasterize (count) sf
  neg_HSMs_count <-  rasterize(all_species_neghsms, tiles[[parse_number(file)]], field = "HSM", fun = "count")
  
  writeRaster(neg_HSMs_count, paste0("outputs/neg_HSMs_count/neg_hsms_count_", n, ".tif"), overwrite = T)
  print(paste0("outputs/neg_HSMs_count/neg_hsms_count_", n, ".tif"))
  
  remove(neg_HSMs_count)
  gc()
  removeTmpFiles()
}

# Merge tiles to get the final result

# neg HSMs count
negHSMs_files <- list.files(path = "outputs/neg_HSMs_count/", pattern = "\\d+.tif", full.names = T)
neg_HSMs_count_tiles <- lapply(negHSMs_files, FUN = raster)

beginCluster(4)
neg_HSMs_count <- mergeRaster(neg_HSMs_count_tiles, fun = max)
endCluster()
neg_HSMs_count[neg_HSMs_count == 0] <- NA
writeRaster(neg_HSMs_count, "outputs/neg_HSMs_count/neg_HSMs_count.tif", overwrite = T)
print("outputs/neg_HSMs_count/neg_HSMs_count.tif")
}

### Mean value per species ####

load(paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))
dim(all_species_traits)

# all_species_traits$HSMs_SaxtonRawls <- meanHSMFun(MinWP_s = MinWP_soil_SX, spAlpha = all_species_traits)
# 
# HSMs_SaxtonRawls <- all_species_traits[, c("sp", "HSMs_SaxtonRawls")]
# save(HSMs_SaxtonRawls, file = paste0(output.dir, "species_traits/all_species_HSMs_RFimp.RData"))
# print(paste0(output.dir, "species_traits/all_species_HSMs_RFimp.RData"))
# 
# load(paste0(output.dir, "species_traits/all_species_HSMs_RFimp.RData"))
# 
# all_species_traits <- merge(all_species_traits, HSMs_SaxtonRawls, by = "sp")


# Order and save
all_species_traits <- all_species_traits[, c("sp", "Species", "P50", "MinWP_md" , "HSM", "MinWP_md_sd", "P50_sd", "HSM_sd", "MinWP_md_cv", "P50_cv", "HSM_cv", "genus", "family" , "order", "group", 
                                             "PC1", "PC2", "PC3", "PC4", "PC5", "Phylo_axis_1", "Phylo_axis_2", "Phylo_axis_3", "Phylo_axis_4", "Phylo_axis_5", "MinWP_soil_SX", "HSMs_SaxtonRawls", "geometry")]

all_species_traits
save(all_species_traits, file = paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))
print(paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))
