#### SOIL MINIMUIM WATER POTENTIAL ------------------------------------------------------------------------------------------------------------------------------------------ ####

print("Soil minimum water potential calculation...")

source("code/data_preparation/0-data_preparation_init.R")


#### FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

# Unify years of the same layer
unifyLayerYearsRastFun <- function(layer = 1){
  r1 <- read_stars(paste0("data/species_data/soil_data/swc/swc", layer, "_1981-2000.grib"), proxy = T)
  r1 <- as(r1, "Raster")
  r2 <- read_stars(paste0("data/species_data/soil_data/swc/swc", layer, "_2000-2019.grib"), proxy = T)
  r2 <- as(r2, "Raster")
  st <- stack(r1, r2)
  remove(list = c("r1", "r2"))
  removeTmpFiles()
  return(st)
}

# Reproject swc data
projSwcFun <- function(target = swc_min, ref){
  target_negative <- raster::crop(target, raster::extent(c(180-0.05, 359.95, -90.05, 90.05)))
  target_positive <- raster::crop(target, raster::extent(c(-0.05, 179.95, -90.05, 90.05)))
  
  # 2. modify the extent
  raster::extent(target_negative) <- raster::extent(c(-180, 0, -90.05, 90.05))
  raster::extent(target_positive) <- raster::extent(c(0, 180, -90.05, 90.05))
  
  # 3. merge
  target_fixed <- raster::merge(target_negative, target_positive)
  
  # 4 reproject
  
  crsToProj <- crs(ref)
  target_reproj <- projectRaster(target_fixed, crs = crsToProj, res = res(target_fixed))
  # target_reproj_agr <- aggregate(target_reproj, fact = 3, fun=min)
  target_reproj <- crop(target_reproj, extent(ref))
  # target_reproj <- projectRaster(from = target_fixed, to = ref)
  
  return(target_reproj)
}

# Soil Hydraulic mean parameters
soilParametersMean <- function(parameters = c("alpha_fit", "n_fit", "mean_theta_s", "mean_theta_r", "var_scaling")){ #  "mean_L", "mean_Ks"
  depths <- c("0cm", "5cm", "15cm", "30cm", "60cm", "100cm") # "200cm"
  for(parameter in parameters){
    st <- stack()
    for(i in 1:length(depths)){
      depth <- depths[i]
      layername <- paste0(parameter, "_", depth)
      r <- raster(paste0("data/species_data/soil_data/soil_hydraulics/Hydraul_Param_SoilGrids_Schaap_", depth,
                         ".nc"), varname = layername)
      st <- stack(st, r)
    }
    par.r <- clusterR(st, overlay, arg = list(fun = mean))
    names(par.r) <- parameter
    writeRaster(par.r, paste0("data/species_data/soil_data/soil_hydraulics/", parameter, ".tif"), overwrite = T)
    print(paste0("data/species_data/soil_data/soil_hydraulics/", parameter, ".tif"))
  }
  f.meanPar <- paste0("data/species_data/soil_data/soil_hydraulics/", list.files("data/species_data/soil_data/soil_hydraulics/",
                                                                                 pattern = ".tif"))
  soilPar.st <- stack(f.meanPar)
  return(soilPar.st)
}

# VG function for rasters
MinWP_soilVGFun <- function(theta_sat, theta_res, theta, n, alpha){
  a <- ( ( ( (theta_sat-theta_res) / (theta-theta_res) )^( n/(n-1) ) ) -1)
  res <- -( (a^(1/n))/alpha)
  # res[res$layer < -40] <- -40 # soil_theta2psiVG seems to saturate at -40
  # res[theta < theta_res] <- -40
  names(res) <- "MinWP_soil"
  res <- stack(res, theta_sat, theta_res, theta, n, alpha)
  return(res)
}

# Curve plot
textReclFun <- function(){
  # Reclassify clay and sand rasters
  x <- raster("data/species_data/soil_data/1km/CLYPPT_M_sl4_1km_ll.tif")
  x <- aggregate(x, fact = 12)
  y <- raster("data/species_data/soil_data/1km/SNDPPT_M_sl4_1km_ll.tif")
  y <-  aggregate(y, fact = 12)
  
  reclass_clay <- matrix(c(0, 18, 1, # coarse
                           18, 35, 2, # medium
                           0, 35, 3, # medium fine
                           35, 60, 4, # fine
                           60, 100, 5), ncol = 3, byrow = T) # very fine
  
  reclass_sand <- matrix(c(65, 100, 1, # coarse
                           15, 65, 2, # medium
                           0, 15, 3), ncol = 3, byrow = T) # medium fine
  
  x.recl <- reclassify(x, reclass_clay)
  y.recl <- reclassify(y, reclass_sand)
  # Texture raster
  
  text.df <- as(stack(x.recl, y.recl), "SpatialPixelsDataFrame")
  text.df$texture <- as.numeric(paste0(text.df$CLYPPT_M_sl4_1km_ll, text.df$SNDPPT_M_sl4_1km_ll))
  text.r <- raster(text.df, layer = "texture")
  
  reclass_text <- matrix(c(01, 1, # Coarse
                           02, 2, # Mdium
                           11, 1, #Coarse
                           12, 2, # Medium
                           13, 3, # Medium-fine
                           21, 2, # Medium
                           22, 2, # medium
                           23, 3, # Medium fine
                           42, 4, # Fine
                           43, 4, # Fine
                           52, 5, # Very fine
                           53, 5 #Very fine
  ), ncol = 2, byrow = T)
  text.recl <- reclassify(text.r, reclass_text)
  writeRaster(text.recl, "data/species_data/soil_data/swc_processed/soil_texture.tif", overwrite = T)
  print("Textures global map saved in data/species_data/soil_data/soil_texture.tif")
  removeTmpFiles() # remove temporal raster files
  return(text.recl)
}

swcCurve <-  function(MinWP_soil, ext = extent(-80, -30, -50, 10)){
  if(!file.exists("data/species_data/soil_data/swc_processed/soil_texture.tif")){
    text.r <- textReclFun() 
  } else {
    text.r <- raster("data/species_data/soil_data/swc_processed/soil_texture.tif")
  }
  text.r <- projectRaster(text.r, MinWP_soil$MinWP_soil_MPa, method = "ngb")
  text.r <- crop(text.r, ext)
  MinWP_soilc <- crop(MinWP_soil, ext)
  MinWP_soilc$texture <- text.r
  MinWP.df <- as.data.frame((as(MinWP_soilc, "SpatialPixelsDataFrame")))
  MinWP.df$texture <- as.factor(MinWP.df$texture)
  p <- ggplot(MinWP.df) + geom_point(aes(x = swc_min, y = MinWP_soil_MPa, color = texture), alpha = 0.5) 
  print(p)
}

# Extract species data (parallelizing)

extrFun <- function(df = sp.df[i, ], spAlpha = all_species[i, ]){
  if(any(is.na(df))){
    e <- extent(spAlpha)
    r <- crop(st, e)
    mask.r <- mask(x = r, mask = spAlpha)
    varValues <- cellStats(x = mask.r, stat = mean, na.rm = T)
    return(varValues)
  }
}

#### DTA AND SPECIFICATIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ####

# treevol data
ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_5_phylo_PCs.csv"), header = T)

# All species polygons
load(file = paste0(processed.data, 'species_distributions/all_species_2020.r'))  # Species alpha-hulls With Behrmann projection, object called new_all_species

# Number of cores to use when paralellizing
nClusters <- 2

#### SWC CALCULATION (ERA from monthly mean) ------------------------------------------------------------------------------------------------------------------------------------------------------- ####
if(!file.exists("data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif")){
### Unify layers (up to 1 m deep) ####

swc1 <- unifyLayerYearsRastFun(layer = 1)
names(swc1) <- paste0(1:length(names(swc1)), "_a")
swc2 <- unifyLayerYearsRastFun(layer = 2)
names(swc2) <- paste0(1:length(names(swc2)), "_b")
swc3 <- unifyLayerYearsRastFun(layer = 3)
names(swc3) <- paste0(1:length(names(swc3)), "_c")
# swc4 <- unifyLayerYearsRastFun(layer = 4)
# names(swc4) <- paste0("month_", 1:length(names(swc4)))

### Calculate mean per pixel (unify layers into one per month) ####

months <- seq(0, 468, 12)
swc_month <- stack()

beginCluster(n = nClusters)
for(i in 1:39){
  start <- months[i] + 1
  finish <- months[i + 1]
  swc_month_layers <- stack(swc1[[start:finish]], swc2[[start:finish]], swc3[[start:finish]])
  swc_month_year <- clusterR(swc_month_layers, stackApply, 
                        args = list(indices = paste0(1:length(names(swc1[[start:finish]]))),
                                    fun = mean), progress='text')
  swc_month <- stack(swc_month, swc_month_year)
  removeTmpFiles()
}
endCluster()

### Annual minimum calculation ####

# Minimum
swc_min_year <- stack()
beginCluster(n = nClusters)
for(i in 1:39){
  start <- months[i] + 1
  finish <- months[i + 1]
  minLayer_annual <- clusterR(swc_month[[start:finish]], overlay, arg = list(fun = min))
  names(minLayer_annual) <- paste0("Y_", as.character(i))
  swc_min_year <- stack(swc_min_year, minLayer_annual)
  remove(minLayer_annual)
  removeTmpFiles()
}
endCluster()

### Absolute minimum per layer ####

beginCluster(n = nClusters)
swc_min <- clusterR(swc_min_year, overlay, arg = list(fun = min))
endCluster()

### Reproject mean from the minimum of each depth layer per pixel ####

beginCluster(n = nClusters)
ref.r <-  raster("data/species_data/soil_data/soil_hydraulics/Hydraul_Param_SoilGrids_Schaap_0cm.nc", varname = "alpha_fit_0cm")
swc_min_proj <- projSwcFun(target = swc_min, ref = ref.r)
endCluster()
plot(swc_min_proj)
writeRaster(swc_min_proj, "data/species_data/soil_data/swc_processed/min_swc_monthly_1980_2019_0_100cm.tif", overwrite = T)
print("data/species_data/soil_data/swc_processed/min_swc_monthly_1980_2019_0_100cm.tif")

save.image(file = "data/RData/minWP_soil.RData")

removeTmpFiles() # remove temporal raster files


### MEDFATE PEDOTRANSFER (Using Saxton and Saxton and Rawls equations) ####

### Soilgrids data on sand, clay and organic matter ####

## Aggregate by cell

## Soil data composition

# Load min swc
swc_min_proj <- raster("data/species_data/soil_data/swc_processed/min_swc_monthly_1980_2019_0_100cm.tif") # minimum of the monthly averaged for the last 40 years 0 - 100 cm averaging layers first

{
  beginCluster(n = nClusters)
  clay.st <- stack(paste0("data/species_data/soil_data/1km/clay/",
                          list.files("data/species_data/soil_data/1km/clay/")))
  
  sand.st <- stack(paste0("data/species_data/soil_data/1km/sand/",
                          list.files("data/species_data/soil_data/1km/sand/")))
  
  soc.st <- stack(paste0("data/species_data/soil_data/1km/soc/",
                         list.files("data/species_data/soil_data/1km/soc/")))
  
  clay <- clusterR(clay.st, overlay, arg = list(fun = mean))
  sand <- clusterR(sand.st, overlay, arg = list(fun = mean))
  soc <- clusterR(soc.st, overlay, arg = list(fun = mean))
  soc <- clusterR(soc, overlay, arg = list(fun = function(x){ x/10}))
  endCluster()
}

soilPart <- stack(clay, sand, soc)
names(soilPart) <- c("clay", "sand", "soc")

save(soilPart, file = "data/species_data/soil_data/1km/processed_sand_clay_soc.RData")
print("data/species_data/soil_data/1km/processed_sand_clay_soc.RData")

# Project and add swc

load("data/species_data/soil_data/1km/processed_sand_clay_soc.RData")

soilPart <- projectRaster(soilPart, swc_min_proj)
soilPart$min_swc <- swc_min_proj

### MinWP by Saxton equations ####

soil.df <- as.data.frame(as(soilPart, "SpatialPixelsDataFrame"))

cl <- makeCluster(nClusters) #not to overload your computer
registerDoParallel(cl)
parallel::clusterExport(cl, varlist=c("soil.df", "soil_theta2psiSX"))

minWP_sx <- foreach(i=1:length(soil.df$x), .combine='c') %dopar% {
  soil_theta2psiSX(sand = soil.df[i, "sand"], 
                   clay = soil.df[i, "clay"],
                   om = soil.df[i, "soc"],
                   theta = soil.df[i, "min_swc"])
}
stopCluster(cl)

soil.df$minWP_sx <- minWP_sx
summary(soil.df)

### Plot ####
minWP_sx <- rasterFromXYZ(soil.df[, c("x", "y", "minWP_sx")], crs = crs(soilPart))

plotMapFun(minWP_sx, "minWP_sx")
log_abs_minWP_sx <- log(-minWP_sx)
plotLogMapFun(log_abs_minWP_sx, "layer")

# Parameters
plotMapFun(soilPart, "min_swc")
plotMapFun(soilPart, "clay")
plotMapFun(soilPart, "sand")
plotMapFun(soilPart, "soc")

### Project to behrmann and save ####

beginCluster(n = nClusters)
minWP_sx_proj <- projectRaster(minWP_sx, crs = behrmann, resol = res(minWP_sx))
endCluster()
writeRaster(minWP_sx_proj, "data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif", overwrite = T)
print("data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif")


}
#### MEAN VALUES PER SPECIES --------------------------------------------------------------------------------------- ####

if(!file.exists("data/species_data/all_species_env_MinWP_soil_SX.csv")){
cl <- parallel::makeCluster(nClusters)
registerDoParallel(cl)
{
  pth = paste0("data/species_data/soil_data/soil_hydraulics/minimum_WP")
  vars = c("MinWP_soil_SX")
  all_species = new_all_species
  # Variables
  vars.df <- data.frame("species" = all_species$Species)
  for(var in vars){
    vars.df[, var] <- NA
    print(paste0("extracting ", var, " values per species..."))
    var.list <- list.files(path = pth, pattern = var, full.names = T)
    print(var.list)
    st <- stack(var.list)
    varNames <- names(st)
    
    # Create empty datafame
    sp.df <- data.frame("species" = all_species$Species)
    sp.df[, varNames] <- rep(NA, length(sp.df$species))
    
    # Extracting ...
    parallel::clusterExport(cl, varlist=c("sp.df", "all_species", "extrFun", "extent", "crop", "mask", "cellStats", "st"))
    
    sp.val <- foreach(i=1:length(sp.df$species), .combine='c') %dopar% {
      extrFun(df = sp.df[i, ], spAlpha = all_species[i, ])
    }
    sp.df$MinWP_soil_SX <- sp.val
    write.csv(sp.df, paste0("data/species_data/all_species_env_", var, ".csv"), row.names = F)
    print(paste0("data/species_data/all_species_env_", var, ".csv"))
    minWP_SX.df <- sp.df
  }
}
stopCluster(cl)
}

# Add to dataset
  
minWP_SX.df <- read.csv("data/species_data/all_species_env_MinWP_soil_SX.csv", header = T)
length(minWP_SX.df$species)
minWP_SX.df$Species <- str_replace_all(minWP_SX.df$species, "_", " ")
minWP_SX.df$species <- NULL

ht$MinWP_soil_SX <- NULL
ht <- merge(ht, minWP_SX.df, by = "Species", all.x = T)

length(ht$Species) # (47197)

write.csv(ht, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_6_mean_MinWP_soil.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_6_mean_MinWP_soil.csv"))
