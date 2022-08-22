##### MAPPING RESULTS ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #####

remove(list = ls())

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

#### DATA AND PARAMETERS ---------------------------------------------------------------- ####

# Figure dimensions in inches

h <- 2.61
w <- 7.87

model.uncertainity.phylo <- F # Using only phylogenetic variables
model.uncertainity.clim <- F # using only climate variables
model.uncertainity.subs <- F # using only a subset of the data to train

# Data regarding uncertainity checking or results (if false and false)

# mortality locations

md <- read.csv("data/mortality/Hammond_etal_database_withspecies_transformed.csv", header = T)

# Observed only
load(paste0(output.dir, "species_traits/obs_species_traits.RData"))

# Observed + Imputed values

if(isFALSE(model.uncertainity.clim) && isFALSE(model.uncertainity.phylo)  && isFALSE(model.uncertainity.subs)){
  load(paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))
  rslts.dir <- paste0(rslts.dir, "figures/maps/")
}

if(isTRUE(model.uncertainity.phylo)){
  load(paste0(output.dir, "species_traits/all_species_traits_phylo20_RFimp.RData"))
  rslts.dir <- paste0(rslts.dir, "figures/maps/only_phylo/")
}

if(isTRUE(model.uncertainity.clim)){
  load(paste0(output.dir, "species_traits/all_species_traits_climate10_RFimp.RData"))
  rslts.dir <- paste0(rslts.dir, "figures/maps/only_clim/")
}

if(isTRUE(model.uncertainity.subs)){
  load("outputs/species_traits/all_species_traits_uncertainity_RFimp.RData") # all_species_traits_uncert
  rslts.dir <- paste0(rslts.dir, "figures/maps/species_uncertainity/")
}

# convert all geometries to MULTIPOLYGON
length(all_species_traits$sp) 
length(obs_species_traits$sp) 

all_species_traits_polygons <- sf::st_cast(all_species_traits, 'MULTIPOLYGON')
obs_species_traits_polygons <- sf::st_cast(obs_species_traits, 'MULTIPOLYGON')

# Load the reference raster
# raster_ref <- raster('data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif') # coarse
raster_ref <- raster("data/species_data/env_data/2.5m/ai_et0.tif") # fine

### Points to zoom in 

zoomPoints <- read.csv(paste0(rslts.dir, "zoomPointsForPlotting.csv"), header = T)
zoomPoints$point_id <- c("a", "b", "c")
coordinates(zoomPoints) <- zoomPoints[, c("x", "y")]
crs(zoomPoints) <- crs(raster_ref)

zoomPoints <- st_as_sf(zoomPoints)
pointsVal <- st_intersection(zoomPoints, all_species_traits_polygons)

### MAPPING ----------------------------------------------------------------------------- ####

### Number of species ####

p50_count_obs <- fasterize(obs_species_traits_polygons, raster_ref, field = 'P50', fun = 'count')
p50_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'P50', fun = 'count')

# Obs count
png(paste0(rslts.dir, "countObsSppMap.png"), height = h, width = w, units = "in", res = 1000)
discPlotMapFun(p50_count_obs)
dev.off()
cat(rslts.dir, "countObsSppMap.png")

# All count
png(paste0(rslts.dir, "countSppMap_rf.png"), height = h, width = w, units = "in", res = 1000)
discPlotMapFun(p50_count)
dev.off()
cat(rslts.dir, "countSppMap_rf.png")

writeRaster(p50_count, paste0(output.dir, "HSM/N_species.tif"), overwrite = T) # To further use
cat(output.dir, "HSM/N_species.tif")


### MinWP_md ####

# Statistics
MinWP_md_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'MinWP_md', fun = 'sum')
MinWP_md_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'MinWP_md', fun = 'count')
MinWP_md_max <- fasterize(all_species_traits_polygons, raster_ref, field = 'MinWP_md', fun = 'max')
MinWP_md_min <- fasterize(all_species_traits_polygons, raster_ref, field = 'MinWP_md', fun = 'min')

# Model performance
MinWP_md_sd_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'MinWP_md_sd', fun = 'sum')

### Statistics calculation

MinWP_md_mean <- overlay(MinWP_md_sum, MinWP_md_count, fun = function(sum, n) {sum/n})

writeRaster(MinWP_md_mean, paste0(output.dir, "Pmin/Pmin_mean.tif"), overwrite = T) # To further use
cat(output.dir, "Pmin/Pmin_mean.tif")

MinWP_md_sd_mean <- overlay(MinWP_md_sd_sum, MinWP_md_count, fun = function(sum, n) {sum/n})

MinWP_md_FRic <- overlay(MinWP_md_max, MinWP_md_min, fun = function(max, min) {max - min})


writeRaster(MinWP_md_FRic, paste0(output.dir, "Pmin/Pmin_FRic.tif"), overwrite = T) # To further use
cat(output.dir, "Pmin/Pmin_FRic.tif")

# beginCluster(n = 2)
# MinWP_md_variance <- rasterize(all_species_traits_polygons, raster_ref, field = 'MinWP_md', fun = 'var', na.rm = T)
# endCluster()

writeRaster(MinWP_md_variance, paste0(output.dir, "diversity_index/MinWP_md_variance.tif"), overwrite = T) # To further use
cat(output.dir, "diversity_index/MinWP_md_variance.tif")

### Plot

# Mean
png(paste0(rslts.dir, "MinWP_md_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
MinWP_md_mean_leg <- plotDataMapFun(mean.r = MinWP_md_mean, max.r = MinWP_md_max, min.r = MinWP_md_min, zoomPnts = T)
dev.off()
cat(rslts.dir, "MinWP_md_mean_rf.png")

pdf(paste0(rslts.dir, "MinWP_md_mean_rf_leg.pdf"),height = 1.1, width = 0.6)
plot(MinWP_md_mean_leg[[1]])
dev.off()
cat(rslts.dir, "MinWP_md_mean_rf_leg.pdf")

# x axis of data plot (not included i the function, because it causes problems with plots alineation)
png(paste0(rslts.dir, "MinWP_md_mean_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = MinWP_md_mean, max.r = MinWP_md_max, min.r = MinWP_md_min, xText = T)
dev.off()
cat(rslts.dir, "MinWP_md_mean_dataplot_xAxisRef.png")

# Histogram zoom points
pointsVal$point_id <- as.factor(pointsVal$point_id)
for(id in levels(pointsVal$point_id)){
  pdf(paste0(rslts.dir, "histograms/MinWP_md_mean_", id, ".pdf"), height = h/1.5, width = h/1.5)
  pointsZoomHistPlot(pointsVal, id = id, var = "MinWP_md")
  dev.off()
  print(paste0(rslts.dir, "histograms/P50_mean"))
}

# Mean sd
png(paste0(rslts.dir, "MinWP_md_sd_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
plotMapFun(MinWP_md_sd_mean)
dev.off()
cat(rslts.dir, "MinWP_md_sd_mean_rf.png")

# Diversity (individual FRIc)

png(paste0(rslts.dir, "MinWP_md_FRic_rf.png"), height = h, width = w, units = "in", res = 1000)
MinWP_md_FRic_leg <- plotDataMapFun(MinWP_md_FRic, var = "layer")
dev.off()
cat(rslts.dir, "MinWP_md_FRic_rf.png")

pdf(paste0(rslts.dir, "MinWP_md_FRic_leg.pdf"), height = 1.1, width = 0.5)
plot(MinWP_md_FRic_leg[[1]])
dev.off()
cat(rslts.dir, "MinWP_md_FRic_leg.pdf")

# X axis text
png(paste0(rslts.dir, "MinWP_md_FRic_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = MinWP_md_FRic, xText = T)
dev.off()
cat(rslts.dir, "MinWP_md_FRic_dataplot_xAxisRef.png")

# MinWP_md_variance

png(paste0(rslts.dir, "MinWP_md_variance_rf.png"), height = h, width = w, units = "in", res = 1000)
MinWP_md_var_leg <- plotDataMapFun(MinWP_md_variance, var = "layer")
dev.off()
cat(rslts.dir, "MinWP_md_variance_rf.png")

pdf(paste0(rslts.dir, "MinWP_md_variance_leg.pdf"), height = 1.4, width = 0.6)
plot(MinWP_md_var_leg[[1]])
dev.off()
cat(rslts.dir, "MinWP_md_variance_leg.pdf")

# X axis text
png(paste0(rslts.dir, "MinWP_md_variance_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = MinWP_md_variance, xText = T)
dev.off()
cat(rslts.dir, "MinWP_md_variance_dataplot_xAxisRef.png")


### P50 ####

p50_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'P50', fun = 'sum')
p50_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'P50', fun = 'count')
p50_max <- fasterize(all_species_traits_polygons, raster_ref, field = 'P50', fun = 'max')
p50_min <- fasterize(all_species_traits_polygons, raster_ref, field = 'P50', fun = 'min')

p50_sd_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'P50_sd', fun = 'sum')

### Statistics calculation

p50_mean <- overlay(p50_sum, p50_count, fun = function(sum, n) {sum/n})

writeRaster(p50_mean, paste0(output.dir, "P50/P50_mean.tif"), overwrite = T) # To further use
cat(output.dir, "P50/P50_mean.tif")

p50_sd_mean <- overlay(p50_sd_sum, p50_count, fun = function(sum, n) {sum/n})

P50_FRic <- overlay(p50_max, p50_min, fun = function(max, min) {max - min})

writeRaster(P50_FRic, paste0(output.dir, "P50/P50_FRic.tif"), overwrite = T) # To further use
cat(output.dir, "P50/P50_FRic.tif")

# beginCluster(n = 2)
# P50_variance <- rasterize(all_species_traits_polygons, raster_ref, field = 'P50', fun = 'var', na.rm = T)
# endCluster()

writeRaster(P50_variance, paste0(output.dir, "diversity_index/P50_sd_variance.tif"), overwrite = T) # To further use
cat(output.dir, "diversity_index/P50_sd_variance.tif")

### Plot

# Mean Map
png(paste0(rslts.dir, "p50_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
p50_mean_leg <- plotDataMapFun(mean.r = p50_mean, max.r = p50_max, min.r = p50_min, zoomPnts = T)
dev.off()
cat(rslts.dir, "p50_mean_rf.png")

pdf(paste0(rslts.dir, "p50_mean_rf_leg.pdf"), height = 1.2, width = 0.6)
plot(p50_mean_leg[[1]])
dev.off()
cat(rslts.dir, "p50_mean_rf_leg.pdf")

pdf(paste0(rslts.dir, "biomes_leg.pdf"))
plot(p50_mean_leg[[2]])
dev.off()
cat(rslts.dir, "biomes_leg.pdf")

# x axis of data plot
png(paste0(rslts.dir, "P50_mean_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = p50_mean, max.r = p50_max, min.r = p50_min, xText = T)
dev.off()
cat(rslts.dir, "P50_mean_dataplot_xAxisRef.png")

# Histogram zoom points
pointsVal$point_id <- as.factor(pointsVal$point_id)
for(id in levels(pointsVal$point_id)){
  pdf(paste0(rslts.dir, "histograms/P50_mean_", id, ".pdf"), height = h/1.5, width = h/1.5)
  pointsZoomHistPlot(pointsVal, id = id, var = "P50")
  dev.off()
  print(paste0(rslts.dir, "histograms/P50_mean"))
}

# sd
png(paste0(rslts.dir, "p50_sd_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
plotMapFun(p50_sd_mean)
dev.off()
cat(rslts.dir, "p50_sd_mean_rf.png")

# Diversity (individual FRIc)

png(paste0(rslts.dir, "P50_FRic_rf.png"), height = h, width = w, units = "in", res = 1000)
P50_FRic_leg <- plotDataMapFun(P50_FRic, var = "layer")
dev.off()
cat(rslts.dir, "P50_FRic_rf.png")

pdf(paste0(rslts.dir, "P50_FRic_leg.pdf"), height = 1.4, width = 0.6)
plot(P50_FRic_leg[[1]])
dev.off()
cat(rslts.dir, "P50_FRic_leg.pdf")

# X axis text
png(paste0(rslts.dir, "P50_FRic_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = P50_FRic, xText = T)
dev.off()
cat(rslts.dir, "P50_FRic_dataplot_xAxisRef.png")

# P50 Diversity (variance)

png(paste0(rslts.dir, "P50_variance_rf.png"), height = h, width = w, units = "in", res = 1000)
P50_var_leg <- plotDataMapFun(P50_variance, var = "layer")
dev.off()
cat(rslts.dir, "P50_variance_rf.png")

pdf(paste0(rslts.dir, "P50_variance_leg.pdf"), height = 1.4, width = 0.6)
plot(P50_var_leg[[1]])
dev.off()
cat(rslts.dir, "P50_variance_leg.pdf")

# X axis text
png(paste0(rslts.dir, "P50_variance_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = P50_variance, xText = T)
dev.off()
cat(rslts.dir, "P50_variance_dataplot_xAxisRef.png")


### HSM ####

HSM_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSM', fun = 'count')
HSM_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSM', fun = 'sum')
HSM_min <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSM', fun = 'min')
HSM_max <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSM', fun = 'max')

writeRaster(HSM_min, paste0(output.dir, "HSM/HSM_min.tif"), overwrite = T) # To further use
cat(output.dir, "HSM/HSM_min.tif")

writeRaster(HSM_max, paste0(output.dir, "HSM/HSM_max.tif"), overwrite = T) # To further use
cat(output.dir, "HSM/HSM_max.tif")

HSM_sd_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSM_sd', fun = 'sum')

## Statistics calculation 

HSM_mean <- overlay(HSM_sum, HSM_count, fun = function(sum, n) {sum/n})

writeRaster(HSM_mean, paste0(output.dir, "HSM/HSM_mean.tif"), overwrite = T) # To furhter use
cat(output.dir, "HSM/HSM_mean.tif")

HSM_FRic <- overlay(HSM_max, HSM_min, fun = function(max, min) {max-min})

writeRaster(HSM_FRic, paste0(output.dir, "HSM/HSM_FRic.tif"), overwrite = T) # To furhter use
cat(output.dir, "HSM/HSM_FRic.tif")


HSM_sd_mean <- overlay(HSM_sd_sum, HSM_count, fun = function(sum, n) {sum/n})

# beginCluster(n = 2)
# HSM_variance <- rasterize(all_species_traits_polygons, raster_ref, field = 'HSM', fun = 'var', na.rm = T)
# endCluster()

writeRaster(HSM_variance, paste0(output.dir, "diversity_index/HSM_variance.tif"), overwrite = T) # To further use
cat(output.dir, "diversity_index/HSM_variance.tif")

HSM_variance <- raster(paste0(output.dir, "diversity_index/HSM_variance.tif"))

writeRaster(HSM_variance, paste0(output.dir, "HSM/HSM_variance.tif"), overwrite = T) # To further use
cat(output.dir, "HSM/HSM_varirance.tif")

### Plot

# Mean
png(paste0(rslts.dir, "HSM_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
HSM_mean_leg <- plotDataMapFun(mean.r = HSM_mean, max.r = HSM_max, min.r = HSM_min, zoomPnts = T)
dev.off()
cat(rslts.dir, "HSM_mean_rf.png")

pdf(paste0(rslts.dir, "HSM_mean_rf_leg.pdf"), height = 1.4, width = 0.6)
plot(HSM_mean_leg[[1]])
dev.off()
cat(rslts.dir, "HSM_mean_rf_leg.pdf")

# x axis of data plot
png(paste0(rslts.dir, "HSM_mean_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = HSM_mean, max.r = HSM_max, min.r = HSM_min, xText = T)
dev.off()
cat(rslts.dir, "HSM_mean_dataplot_xAxisRef.png")

# Histogram zoom points
pointsVal$point_id <- as.factor(pointsVal$point_id)
for(id in levels(pointsVal$point_id)){
  pdf(paste0(rslts.dir, "histograms/HSM_mean_", id, ".pdf"), height = h/1.5, width = h/1.5)
  pointsZoomHistPlot(pointsVal, id = id, var = "HSM")
  dev.off()
  print(paste0(rslts.dir, "histograms/P50_mean"))
}

# Mean sd
png(paste0(rslts.dir, "HSM_sd_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
plotMapFun(HSM_sd_mean)
dev.off()
cat(rslts.dir, "HSM_sd_mean_rf.png")

# Diversity (individual FRIc)

# HSM_FRic <- div.st$HSM_FRic

png(paste0(rslts.dir, "HSM_FRic_rf.png"), height = h, width = w, units = "in", res = 1000)
HSM_FRic_leg <- plotDataMapFun(HSM_FRic, var = "layer")
dev.off()
cat(rslts.dir, "HSM_FRic_rf.png")

pdf(paste0(rslts.dir, "HSM_FRic_rf_leg.pdf"), height = 1.1, width = 0.6)
plot(HSM_FRic_leg[[1]])
dev.off()
cat(rslts.dir, "HSM_FRic_rf_leg.pdf")


# X axis text
png(paste0(rslts.dir, "HSM_FRic_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = HSM_FRic, xText = T)
dev.off()
cat(rslts.dir, "HSM_FRic_dataplot_xAxisRef.png")


# Minimum
png(paste0(rslts.dir, "HSM_min_rf.png"), height = h, width = w, units = "in", res = 1000)
HSM_min_leg <- plotDataMapFun(mean.r = HSM_min, var = "layer")
dev.off()
cat(rslts.dir, "HSM_min_rf.png")

pdf(paste0(rslts.dir, "HSM_min_rf_leg.pdf"), height = 1.4, width = 0.6)
plot(HSM_min_leg[[1]])
dev.off()
cat(rslts.dir, "HSM_min_rf_leg.pdf")

# x axis of data plot
png(paste0(rslts.dir, "HSM_min_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = HSM_min, xText = T)
dev.off()
cat(rslts.dir, "HSM_min_dataplot_xAxisRef.png")


# HSM Diversity (variance)

png(paste0(rslts.dir, "HSM_variance_rf.png"), height = h, width = w, units = "in", res = 1000)
HSM_var_leg <- plotDataMapFun(HSM_variance, var = "layer")
dev.off()
cat(rslts.dir, "HSM_variance_rf.png")

pdf(paste0(rslts.dir, "HSM_variance_leg.pdf"), height = 1.4, width = 0.6)
plot(HSM_var_leg[[1]])
dev.off()
cat(rslts.dir, "HSM_variance_leg.pdf")

# X axis text
png(paste0(rslts.dir, "HSM_variance_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = HSM_variance, xText = T)
dev.off()
cat(rslts.dir, "HSM_variance_dataplot_xAxisRef.png")

### Negative HSM ####


all_species_traits_polygons$neg_HSM_count <- ifelse(all_species_traits_polygons$HSM < 0, 1, 0)

HSM_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSM', fun = 'count')
negHSM_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'neg_HSM_count', fun = 'sum')

writeRaster(negHSM_count, paste0(output.dir, "HSM/negative_HSM_count.tif"), overwrite = T)  # write to further use
cat(output.dir, "HSM/negative_HSM_count.tif")


### Plot

# Neg HSM count

png(paste0(rslts.dir, "negHSM_count.png"), height = h, width = w, units = "in", res = 1000)
negHSM_count_leg <- plotDataMapFun(mean.r = negHSM_count)
dev.off()
cat(rslts.dir, "negHSM_count.png")

pdf(paste0(rslts.dir, "negHSM_count_leg.pdf"),height = 1.1, width = 0.6)
plot(negHSM_count_leg[[1]])
dev.off()
cat(rslts.dir, "negHSM_count_leg.pdf")

# X axis text
png(paste0(rslts.dir, "neg_HSM_count_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = negHSM_count, xText = T)
dev.off()
cat(rslts.dir, "neg_HSM_count_dataplot_xAxisRef.png")


### Minimum water potential (Saxton) ####

minWP_sx <- raster("data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif")

# x axis of data plot
png(paste0(rslts.dir, "minWP_sx_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(minWP_sx, xText = T)
dev.off()
cat(rslts.dir, "minWP_sx_dataplot_xAxisRef.png")


png(paste0(rslts.dir, "minWP_sx.png"), height = h, width = w, units = "in", res = 1000)
minWP_sx_plot_leg <- plotDataMapFun(minWP_sx, var = "MinWP_soil_SX")
dev.off()
cat(rslts.dir, "minWP_sx.png")

pdf(paste0(rslts.dir, "minWP_sx_leg.pdf"),height = 1.4, width = 0.6)
plot(minWP_sx_plot_leg[[1]])
dev.off()
cat(rslts.dir, "minWP_sx_leg.pdf")

# log map
log_minWP_sx <- log(-minWP_sx)
minWP_sx_plot <- plotLogMapFun(log_minWP_sx, var = "layer")
minWP_sx_plot <- minWP_sx_plot + theme(legend.position = "right")
MinWP_soil_sx_leg <- get_legend(minWP_sx_plot)
minWP_sx_plot <- minWP_sx_plot + theme(legend.position = "none")

png(paste0(rslts.dir, "log_abs_minWP_sx.png"), height = h, width = w, units = "in", res = 1000)
cat(minWP_sx_plot)
dev.off()
cat(rslts.dir, "log_abs_minWP_sx.png")

pdf(paste0(rslts.dir, "log_abs_minWP_sx_leg.pdf"),height = 1.4, width = 0.6)
plot(MinWP_soil_sx_leg[[1]])
dev.off()
cat(rslts.dir, "log_abs_minWP_sx_leg.pdf")
### HSMs (mean per species) ####

HSMs_count <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSMs', fun = 'count')
HSMs_sum <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSMs', fun = 'sum')
HSMs_min <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSMs', fun = 'min')
HSMs_max <- fasterize(all_species_traits_polygons, raster_ref, field = 'HSMs', fun = 'max')


HSMs_mean <- overlay(HSMs_sum, HSMs_count, fun = function(sum, n) {sum/n})

writeRaster(HSMs_mean, paste0(output.dir, "HSM/HSMs_species_mean.tif"), overwrite = T) # To further use
cat(output.dir, "HSM/HSMs_species_mean.tif")

# HSMs_mean[HSMs_mean < -20] <- -20
# HSMs_max[HSMs_max < -20] <- -20
# HSMs_min[HSMs_min < -20] <- -20

### Plot

# Mean
png(paste0(rslts.dir, "HSMs_species_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
HSMs_species_mean_leg <- plotDataMapFun(mean.r = HSMs_mean, max.r = HSMs_max, min.r = HSMs_min)
dev.off()
cat(rslts.dir, "HSMs_species_mean_rf.png")

pdf(paste0(rslts.dir, "HSMs_species_mean_rf_leg.pdf"), height = 1.1, width = 0.6)
plot(HSMs_species_mean_leg[[1]])
dev.off()
cat(rslts.dir, "HSMs_species_mean_rf_leg.pdf")

# x axis of data plot
png(paste0(rslts.dir, "HSMs_species_mean_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
plotDataFun(mean.r = HSMs_mean, max.r = HSMs_max, min.r = HSMs_min, xText = T)
dev.off()
cat(rslts.dir, "HSMs_species_mean_dataplot_xAxisRef.png")


### Plot mortality locations ####

png(paste0(rslts.dir, "mortality_vs_bg_HSM_mean.png"), height = h, width = w, units = "in", res = 1000)
plotMapPointsFun(mean.r = HSM_mean, var = "layer", points = md[, c("x", "y")])
dev.off()
cat(rslts.dir, "mortality_vs_bg_HSM_mean.png")


#### NOT USED FOR NOW -->

# ### Continuous HSMs (based on SX soil water potential calculation) ####
# 
# HSMs_mean <- raster(paste0(output.dir,"HSMs_mean/HSMs_mean.tif"))
# HSMs_max <- raster(paste0(output.dir, "HSMs_max/HSMs_max.tif"))
# HSMs_min <- raster(paste0(output.dir, "HSMs_min/HSMs_min.tif"))
# 
# HSMs_mean[HSMs_mean < -20] <- -20
# HSMs_max[HSMs_max < -20] <- -20
# HSMs_min[HSMs_min < -20] <- -20
# 
# png(paste0(rslts.dir, "HSMs_mean_rf.png"), height = h, width = w, units = "in", res = 1000)
# HSMs_mean_leg <- plotDataMapFun(mean.r = HSMs_mean, max.r = HSMs_max, min.r = HSMs_min, var = "HSMs_mean")
# dev.off()
# cat(rslts.dir, "HSMs_mean_rf.png")
# 
# png(paste0(rslts.dir, "HSMs_mean_rf_leg.png"), height = 1.1, width = 0.6, units = "in", res = 1000)
# plot(HSMs_mean_leg[[1]])
# dev.off()
# cat(rslts.dir, "HSMs_mean_rf_leg.png")
# 
# # x axis of data plot
# png(paste0(rslts.dir, "HSMs_continuous_dataplot_xAxisRef.png"), height = h, width = w, units = "in", res = 1000)
# plotDataFun(mean.r = HSMs_mean, max.r = HSMs_max, min.r = HSMs_min, xText = T)
# dev.off()
# cat(rslts.dir, "HSMs_continuous_dataplot_xAxisRef.png")
# 

### Negative HSMs ####
# 
# HSMs_count <- raster(paste0(output.dir, "neg_HSMs_count/HSMs_count.tif"))
# negHSMs_count <- raster(paste0(output.dir, "neg_HSMs_count/neg_HSMs_count.tif"))
# 
# negHSMs_prop <- overlay(negHSMs_count, HSMs_count, fun = function(negn, n) {negn/n})
# plot(negHSMs_prop)
# 
# writeRaster(negHSMs_prop, paste0("neg_HSMs_count/neg_HSMs_prop.tif", overwrite = T))  # write to further use (extension)
# cat(output.dir, "neg_HSMs_count/neg_HSMs_prop.tif")
# 
# ### Plot
# 
# # Neg HSM count
# 
# negHSMs_count_plot <- plotMapFun(negHSMs_count, var = "neg_HSMs_count")
# negHSMs_count_plot <- negHSMs_count_plot + theme(legend.position = "right")
# negHSMs_count_leg <- get_legend(negHSMs_count_plot)
# negHSMs_count_plot <- negHSMs_count_plot + theme(legend.position = "none")
# 
# png(paste0(rslts.dir, "negHSMs_count.png"), height = h, width = w, units = "in", res = 1000)
# cat(negHSMs_count_plot)
# dev.off()
# 
# png(paste0(rslts.dir, "negHSMs_count_leg.png"),height = 1.1, width = 0.6, units = "in", res = 1000)
# plot(negHSMs_count_leg)
# dev.off()
# 
# # Neg HSM prop
# negHSMs_prop_plot <- plotMapFun(negHSMs_prop)
# negHSMs_prop_plot <- negHSMs_prop_plot + theme(legend.position = "right")
# negHSMs_prop_leg <- get_legend(negHSMs_prop_plot)
# negHSMs_prop_plot <- negHSMs_prop_plot + theme(legend.position = "none")
# 
# png(paste0(rslts.dir, "negHSMs_prop.png"), height = h, width = w, units = "in", res = 1000)
# cat(negHSMs_prop_plot)
# dev.off()
# 
# png(paste0(rslts.dir, "negHSMs_prop_leg.png"),height = 1.3, width = 0.6, units = "in", res = 1000)
# plot(negHSMs_prop_leg)
# dev.off()


#### INTERSECTION BETWEEN SPECIES POLYGONS AND SPECIES MORTALITY POINTS (SPECIES BY POINT, SPECIES-WISE DATABASE) ------------------------------------------------- ####

# Intersection between points and polygons (extract each polygon attributes for each point, ref_id would allow to summaryze data by point if needed)
# 
# # mortality points by species
# load(file = "data/mortality/Hammond_etal_database_withspecies_transformed_species_wise_geo.RData") # object called md_species.sf
# 
# 
# ## HSM and mortality congruence at specie level
# 
# # Negative HSM
# neg_HSM_all_species_traits_polygons <- all_species_traits_polygons %>% filter(HSM < 0)
# length(neg_HSM_all_species_traits_polygons$sp) # 8846 species with negative HSM
# length(all_species_traits_polygons$sp)
# 
# length(neg_HSM_all_species_traits_polygons$sp)/ length(all_species_traits_polygons$sp)
# 
# length(unique(md_species.sf$species))
# shared_spp <- all_species_traits_polygons %>% filter(Species %in% md_species.sf$species) %>% dplyr::select(Species)
# length(unique(shared_spp$Species)) # 471 species for which we have trait data and mortality
# 
# 
# shared_mortality_spp <- neg_HSM_all_species_traits_polygons %>% filter(Species %in% md_species.sf$species) %>% dplyr::select(Species)
# length(shared_mortality_spp$Species) # 79 species with negative HSM and mortality
# 
# length(shared_mortality_spp$Species) / length(unique(shared_spp$Species)) #☺17% of species with mortality present negative HSM
# 
# 
# # low HSM (<0.5)
# 
# low_HSM_all_species_traits_polygons <- all_species_traits_polygons %>% filter(HSM < 0.5)
# length(low_HSM_all_species_traits_polygons$sp) # 35445 species with HSM < 0.5
# length(low_HSM_all_species_traits_polygons$sp)/length(all_species_traits_polygons$sp)
#   
# shared_spp <- all_species_traits_polygons %>% filter(Species %in% md_species.sf$species) %>% dplyr::select(Species)
# length(unique(shared_spp$Species)) # 471 species for which we have trait data and mortality
# 
# shared_mortality_spp <- low_HSM_all_species_traits_polygons %>% filter(Species %in% md_species.sf$species) %>% dplyr::select(Species)
# length(shared_mortality_spp$Species) # 304 species with negative HSM and mortality
# 
# length(shared_mortality_spp$Species) / length(unique(shared_spp$Species)) #64.5% %  of species with mortality present HSM < 0.5
# 
# 
# # low HSM (<0.1)
# 
# low_HSM_all_species_traits_polygons <- all_species_traits_polygons %>% filter(HSM < 1)
# length(low_HSM_all_species_traits_polygons$sp) # 44065 species with HSM < 1
# 
# shared_spp <- all_species_traits_polygons %>% filter(Species %in% md_species.sf$species) %>% dplyr::select(Species)
# length(unique(shared_spp$Species)) # 471 species for which we have trait data and mortality
# 
# shared_mortality_spp <- low_HSM_all_species_traits_polygons %>% filter(Species %in% md_species.sf$species) %>% dplyr::select(Species)
# length(shared_mortality_spp$Species) # 413 species with negative HSM and mortality
# 
# length(shared_mortality_spp$Species) / length(unique(shared_spp$Species)) #87.7%  of species with mortality present HSM < 1
# 
# 
# ## HSM and mortality congruence at specie x site level
# 
# # md_species_all_species_traits <- st_intersection(md_species.sf, all_species_traits_polygons) (takes time, run once and save)
# 
# # md_species_all_species_traits$sp <- NULL
# # md_species_all_species_traits <- md_species_all_species_traits %>% dplyr::rename(species_mortality = species, species_range = Species, species_phylogeny = sp.1) %>% dplyr::select(!c(genus.1, family.1, order.1, group.1))
# # 
# # save(md_species_all_species_traits, file = paste0(output.dir, "mortality/spatial_species_wise_mortality_traits.RData")) # object name: md_species_all_species_traits
# # print(output.dir, "mortality/spatial_species_wise_mortality_traits.RData"))
# 
# load(paste0(output.dir, "mortality/spatial_species_wise_mortality_traits.RData"))
# md_species_all_species_traits
# length(md_species_all_species_traits$ref_id)
# 
# # List of mortality species that coicide with species range
# 
# md_species_all_species_traits_overlap <- md_species_all_species_traits %>% filter(species_mortality == species_range)
# length(md_species_all_species_traits_overlap$species_mortality) # 2628
# length(unique(md_species_all_species_traits_overlap$species_mortality)) # 409 species
# 
# 
# 
# plotMapMdPointsFun <- function(mean.r = raster_ref/10000, var = "ai_et0", points = md_species_all_species_traits_overlap){
#   
#   df <- data.frame(as(mean.r, "SpatialPixelsDataFrame"))
#   cols <- brewer.pal(10, "Spectral")
#   
#   p <- ggplot(df) + geom_raster(aes(x=x, y=y, fill=df[, var]))
#   p <- p + scale_fill_gradientn(colours= cols)
#   p <- p + coord_equal()
#   p <- p + theme_light() +  theme(panel.grid.minor.x = element_blank(),
#                                   panel.grid.minor.y = element_blank(),
#                                   axis.title = element_blank(),
#                                   axis.text.x = element_blank(),
#                                   axis.text.y = element_blank(),
#                                   axis.ticks.x = element_blank(),
#                                   axis.ticks.y = element_blank(),
#                                   legend.title = element_blank(),
#                                   legend.position = "left", #legend.position = c(0.21, 0.19)
#                                   plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
# 
#   points$HSM_disc <- discretize(points$HSM, method = "interval", breaks = 5)
#   
#   p <- p + geom_point(data = points, aes(x = x, y = y, color = HSM_disc), size = 1.5)
#   p <- p + scale_color_discrete()
#   print(p)
#   return(p)
# }
# 
