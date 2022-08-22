#### DATA PREPARATION (BIOCLIM, SOILGRID AND ET)  ------------------------------------------------------------------------------------------------------------------- ####

## P Sanchez-Martinez

print("Species data preparation ...")

#### PACKAGES -------------------------------------------------------------------------------------- ####

# install.packages("rgdal", dependencies = T)

library(raster)
library(stringr)
library(sp)
library(doParallel)
library(snow)
library(rgdal)

#### DATA ----------------------------------------------------------------------------------------- ####

processed.dir <- "data/species_data/env_data/"
load(file='data/species_distributions/all_species.r')  # Species alpha-hulls With Behrmann projection

nClusters <- 12 # set the number of clusters you want to work with when parallelizing

#### FUNCTIONS ------------------------------------------------------------------------------------ #####

# Get WC variables, unzip and delete zip
getWCFun <- function(vars =  c("tmin", "tmax", "tavg", "prec", "srad", "wind", "vapr", "bio"), res = "2.5m"){
  for(var in vars){
    # Download
    print(paste0("Downloading ", var, " ", res, " data..."))
    download.file(paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_", res, "_", var, ".zip"),
                  destfile=paste0("data/species_data/worldclim/wc", res, "/", var, ".zip"), mode="wb")
    # Unzip
    print(paste0("Downloading ", var, " data..."))
    unzip(paste0("data/species_data/worldclim/wc", res, "/", var, ".zip"), exdir = paste0("data/species_data/worldclim/wc", res), overwrite = T)
    # Delete zip file
    file.remove(paste0("data/species_data/worldclim/wc", res, "/", var, ".zip"))
  }
}

# Get SG variables, unzip and delete zip
getSGFun <- function(vars =  c("BDTICM_M_1km_ll.tif", "CECSOL_M_sl4_1km_ll.tif", "ORCDRC_M_sl4_1km_ll.tif", "PHIHOX_M_sl4_1km_ll.tif", "CLYPPT_M_sl4_1km_ll.tif", "SLTPPT_M_sl4_1km_ll.tif", 
                               "SNDPPT_M_sl4_1km_ll.tif", "TAXNWRB_1km_ll.tif", "WWP_M_sl7_1km_ll.tif")){
  for(var in vars){
    print(paste0("Downloading ", var, " data..."))
    download.file(paste0("https://files.isric.org/soilgrids/data/aggregated/1km/", var), destfile = paste0("data/species_data/soil_data/1km/sg_", var), mode = "wb")
    }
}

# Get AI data, unzip and delete zip
getAiFun <- function(dfile = "https://ndownloader.figshare.com/articles/7504448/versions/3", sfile = "data/species_data/aridity_index/AI_30s.zip", 
                     fdir = "data/species_data/aridity_index"){
    # Download
    print(dfile)
    download.file(dfile, destfile = sfile, mode="wb")
    # Unzip
    unzip(sfile, exdir = fdir, overwrite = T)
    unzip("data/species_data/aridity_index/global-ai_et0.zip", exdir = fdir, overwrite = T)
    # Delete zip file
    file.remove(sfile)
  }


# Annualize variables
anuMeanWCFun <- function(vars =  c("tmin", "tmax", "srad", "wind", "vapr"), res = "2.5m"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/worldclim/wc", res, "/wc2.1_", res, "_", var, "_mean_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = mean, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      rm(var.st, anuVar.r) 
    }
  }
}

anuMaxWCFun <- function(vars =  c("wind"), res = "2.5m"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/worldclim/wc", res, "/wc2.1_", res, "_", var, "_max_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = max, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      rm(var.st, anuVar.r) 
    }
  }
}

anuMinWCFun <- function(vars =  c("srad"), res = "2.5m"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/worldclim/wc", res, "/wc2.1_", res, "_", var, "_min_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = min, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      rm(var.st, anuVar.r) 
    }
  }
}

anuSumWCFun <- function(vars =  c("prec"), res = "2.5m"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/worldclim/wc", res, "/wc2.1_", res, "_", var, "_sum_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = sum, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      rm(var.st, anuVar.r) 
    }
  }
}

# Project variables (Behrmann)

projectWCBehFun <- function(res = "2.5m", vars = c("bio", "tmax", "vapr", "annual")){
  for(var in vars){
    pth <- paste0("data/species_data/worldclim/wc", res)
    ptrn <- var
    behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs") # Projection
    var.list <- list.files(path = pth, pattern = ptrn, full.names = T)
    for(var.file in var.list){
      r <- raster(var.file)
      # If file exist
      if(file.exists(paste0(processed.dir, "/", res, "/", names(r), ".tif"))){
        next()
      } else {
        # If alrady reprojected
        if(!is.na(str_extract(projection(r), "cea"))){
          writeRaster(r, paste0(processed.dir, "/", res, "/", names(r), ".tif"), overwrite = T)
          print(paste0("==== >", processed.dir, res, "/", names(r), ".tif"))
        } else {
        # Reproject
        print(paste0("reprojecting...", var.file))
        if(res == "2.5m"){
          resol <- c(4.02, 5.32)
        } else { if(res == "30s"){
          resol <- c(0.008333333333333, 0.008333333333333)
        }
          }
        r <- projectRaster(from = r, crs = behrmann, res = resol)
        writeRaster(r, paste0(processed.dir, "/", res, "/", names(r), ".tif"), overwrite = T)
        print(paste0("==== >", processed.dir, "/",res, "/", names(r), ".tif"))
        removeTmpFiles() # remove temporal raster files
        }
      }
    }
  }
}

projectRastFun <- function(res = "2.5m",pth = "data/species_data/soil_data/1km", ptrn = ".tif", mthd = "bilinear"){
  ref.r <- raster(paste0(processed.dir,res, "/", "wc2.1_2.5m_bio_1.tif"))
  var.list <- list.files(path = pth, pattern = ptrn, full.names = T)
  for(var.file in var.list){
    print(var.file)
    r <- raster(var.file)
    if(file.exists(paste0(processed.dir, "/", res, "/", names(r), ".tif"))){
      next()
    } else {
    if(!is.na(str_extract(projection(r), "cea"))){
      writeRaster(r, paste0(processed.dir, res, "/", names(r), ".tif"), overwrite = T)
      print(paste0("==== >", processed.dir, res, "/", names(r), ".tif"))
    } else {
      r <- projectRaster(from = r, to = ref.r, method = mthd)
      writeRaster(r, paste0(processed.dir, "/", res, "/", names(r), ".tif"), overwrite = T)
      print(paste0("==== >", processed.dir, "/", res, "/", names(r), ".tif"))
      removeTmpFiles() # remove temporal raster files 
    }
    }
  }
}


# Extract environmental variables species mean values
extrSpFun <- function(pth = paste0("data/species_data/env_data/2.5m"), vars = c("ai", "wc", "sg"), sttstc = mean, prev.rslts = F){
  # Variables
  for(var in vars){
    print(paste0("extracting ", var, " values per species..."))
    var.list <- list.files(path = pth, pattern = var, full.names = T)
    print(var.list)
    st <- stack(var.list)
    varNames <- names(st)
    # Load previous results, if exist
    if(isTRUE(prev.rslts)){
      sp.df <- read.csv(paste0("data/species_data/all_species_env_", var, ".csv"), header = T)
    } else {
      sp.df <- data.frame("species" = all_species$SOURCE_SHP)
      sp.df[, varNames] <- rep(NA, length(sp.df$species))
    }
    # Extracting ...
    for(i in 1:length(sp.df$species)){
      if(any(is.na(sp.df[i, varNames]))){
        print(paste0(i, ": ", sp.df$species[i]))
        spAlpha <- all_species[i, ]
        e <- extent(spAlpha)
        r <- crop(st, e)
        # mask.r <- clusterR(x = st, mask, args=list(mask = spAlpha))
        mask.r <- mask(x = r, mask = spAlpha)
        varValues <- cellStats(x = mask.r, stat = sttstc, na.rm = T)
        names(varValues) <- varNames
        for(varName in varNames){
          sp.df[i, varName] <- varValues[varName] 
        }
        write.csv(sp.df, paste0("data/species_data/all_species_env_", var, ".csv"), row.names = F)
        removeTmpFiles(h=1) # remove temporal raster files older than 1 hour
      } else {
        next()
      }
    }
  }
  return(sp.df)
}


#### GET DATA ------------------------------------------------------------------------------------------------------------------------------ ####

# Get worldclim data
getWCFun(vars =  c("tmin", "tmax", "tavg", "prec", "srad", "wind", "vapr", "bio"), res = "2.5m")

# Annualize climate data
anuMeanWCFun(vars = c("tmin", "tmax", "srad", "wind", "vapr"), res = "2.5m")
anuSumWCFun(vars = c("prec"), res = "2.5m")
anuMaxWCFun(vars = c("wind"), res = "2.5m")
anuMinWCFun(vars = c("srad"), res = "2.5m")

beginCluster(n = nClusters)
anuMeanWCFun(vars = c("tmax", "srad", "wind", "vapr"), res = "30s")
anuSumWCFun(vars = c("prec"), res = "30s")
anuMaxWCFun(vars = c("wind"), res = "30s")
anuMinWCFun(vars = c("srad"), res = "30s")
endCluster()

# Get soilgrids data
getSGFun(vars =  c("BDTICM_M_1km_ll.tif", "CECSOL_M_sl4_1km_ll.tif", "ORCDRC_M_sl4_1km_ll.tif", "PHIHOX_M_sl4_1km_ll.tif", "CLYPPT_M_sl4_1km_ll.tif", "SLTPPT_M_sl4_1km_ll.tif", 
                   "SNDPPT_M_sl4_1km_ll.tif", "TAXNWRB_1km_ll.tif", "WWP_M_sl7_1km_ll.tif"))

# Aridity index
getAiFun(dfile = "https://ndownloader.figshare.com/articles/7504448/versions/3", sfile = "data/species_data/aridity_index/AI_30s.zip", 
         fdir = "data/species_data/aridity_index")

#### PROJECT AND SAVE -------------------------------------------------------------------------------------------------------------------- ####

# Reproject WC data (we will use it later as reference)
beginCluster(n = nClusters)
projectWCBehFun(res = "2.5m", vars = c("bio", "tmax", "vapr", "annual"))
endCluster()

beginCluster(n = nClusters)
projectWCBehFun(res = "30s", vars = c("bio", "tmax", "vapr", "annual"))
endCluster()

# Reproject soilgrids
beginCluster(n = nClusters)
projectRastFun(res = "2.5m", pth = "data/species_data/soil_data/1km/", ptrn = ".tif")
endCluster()

beginCluster(n = nClusters)
projectRastFun(res = "30s", pth = "data/species_data/soil_data/1km", ptrn = ".tif")
endCluster()

# Reproject AI

beginCluster(n = nClusters)
projectRastFun(res = "2.5m", pth = "data/species_data/aridity_index/ai_et0/", ptrn = "ai_et0.tif")
endCluster()

beginCluster(n = nClusters)
projectRastFun(res = "30s", pth = "data/species_data/aridity_index/ai_et0/", ptrn = "ai_et0.tif")
endCluster()

r <- raster("data/species_data/env_data/2.5m/ai_et0.tif")
r
#### SPECIES DATA EXTRACTION USING ALPHA HULLS -------------------------------------------------------------------------------------------- ####

# "ai", "bio", "tmax", "vapr", "annual", "sg" # all these need to be run, depending on the computer power run them individually, it might take time.

envSpp.df <- extrSpFun(pth = paste0("data/species_data/env_data/2.5m"), vars = c("annual"), sttstc = mean, prev.rslts = T)
