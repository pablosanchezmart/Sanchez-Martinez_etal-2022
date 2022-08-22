#### BIOMES AND FUNCTIONAL TYPES -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

print("Biomes and funcional types reclassification...")

source("code/data_preparation/0-data_preparation_init.R")

#### DATA AND SPECIFICATIONS ------------------------------------------------------------------------------------------------- ####

# Ecoregions (https://ecoregions2017.appspot.com/)
ecoRegions <- st_read("data/species_data/ecoregions/Ecoregions2017.shp")

# Functional types (from copernicus landcover (https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-land-cover?tab=overview))
lc.r <- raster("data/species_data/ecoregions/land_cover_copernicus/land_cover.nc")

# Reference raster (for projection)

ref.r <- raster("data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif")

# Number of clusters to use when paralellizing
nClusters <- 2

#### BIOMES AND ECOREGIONS --------------------------------------------------------------------------------------------------- ####

### ECOREGIONS DATA PREPARATION ####

ecoRegions_proj <- st_transform(ecoRegions, behrmann)
ecoRegions.sf <- st_as_sf(ecoRegions_proj)
ecoRegions.sf$BIOME_NAME <- as.factor(ecoRegions.sf$BIOME_NAME)
ecoRegions.sf$BIOME_NUM <- as.factor(ecoRegions.sf$BIOME_NUM)

biome.r <- fasterize(sf = ecoRegions.sf, raster = ref.r, field = "BIOME_NUM")
biome.r <- projectRaster(from = biome.r, to = ref.r, method = "ngb")
plot(biome.r)

writeRaster(biome.r, "outputs/biomes_ecoregions/biomes.tif", overwrite = T)
print("outputs/biomes_ecoregions/biomes.tif")

## Dataframe

df <- as.data.frame(as(biome.r, "SpatialPixelsDataFrame"))
df$biome <- as.factor(df$layer)

biome.recl <- data.frame("BIOME_NAME" = unique(as.factor(ecoRegions.sf$BIOME_NAME))[-15], 
                         "BIOME_NUM" = unique(as.factor(ecoRegions.sf$BIOME_NUM)))

for(biome_num in levels(df$biome)){
  biome_name <- as.character(biome.recl[biome.recl$BIOME_NUM == biome_num, "BIOME_NAME"])
  df[which(df$biome == biome_num), "biome.name"] <- biome_name
}
df$biome.name <- as.factor(df$biome.name)

ggplot(df) + geom_raster(aes(x = x, y = y, fill = biome.name))


#### BIOME RECLASSIFICATION ####

biome.recl$BIOME_SYNTH <- c("Boreal", "Tropical and subtropical moist", "Mediterranean", "Desert and xeric", "Temperate", 
                            "Boreal", "Temperate", 
                            "Temperate", "Others", "Tropical and subtropical dry",
                            "Others", "Tropical and subtropical dry", "Tropical and subtropical dry", "Tropical and subtropical dry")

for(biome_num in levels(df$biome)){
  biome_name <- as.character(biome.recl[biome.recl$BIOME_NUM == biome_num, "BIOME_NAME"])
  df[which(df$biome == biome_num), "biome.name"] <- biome_name
  
  biome_synth <- as.character(biome.recl[biome.recl$BIOME_NUM == biome_num, "BIOME_SYNTH"])
  df[which(df$biome == biome_num), "biome.synth"] <- biome_synth
  
}

df$biome.name <- as.factor(df$biome.name)
df$biome.synth <- as.factor(df$biome.synth)
summary(df)

write.csv(df, "outputs/biomes_ecoregions/biomes_dataframe.csv", row.names = F)
print("outputs/biomes_ecoregions/biomes_dataframe.csv")


#### LAND COVER COPERNICUS DATA PREPARATION (FOR FUNCTIONAL TYPES) ####

beginCluster(n = nClusters)
lc_proj.r <- projectRaster(lc.r, ref.r, method = "ngb")
endCluster()

writeRaster(lc_proj.r, "data/species_data/ecoregions/land_cover_copernicus/land_cover.tif", overwrite = T)
print("data/species_data/ecoregions/land_cover_copernicus/land_cover.tif")

#### FUNCTIONAL TYPE RECLASSIFICATION ####

lc_proj.r <- raster("data/species_data/ecoregions/land_cover_copernicus/land_cover.tif")

df <- as.data.frame(as(lc_proj.r, "SpatialPixelsDataFrame"))

ft.recl <- data.frame("FT_NAME" = c("Crops and Grassland", "Crops and Grassland", "Crops and Grassland", "Crops and Grassland", "Mosaic", "Mosaic", "Broadleaved evergreen",  "Broadleaved deciduous", "Broadleaved deciduous", "Broadleaved deciduous", "Needleleaved",  "Needleleaved", "Needleleaved",  "Needleleaved", "Needleleaved", "Mixed forest", "Mosaic", "Mosaic", "Shrubland",  "Shrubland",  "Shrubland", "Crops and Grassland", "Others", "Others", "Others",   "Others", "Others", "Others", "Others", "Others", "Others", "Others", "Others","Others",  "Others"),
                      "LC_NUM" = c("10",                         "11",                       "12",              "20",             "30",     "40",            "50",                      "60",                  "61",                   "62",                  "70",             "71",              "72",               "80",       "81",              "90",    "100",   "110",      "120",       "121",      "122",                 "130",     "140",    "150",   "152",      "153",     "160",   "170",     "180",     "190", "200",  "201", "202",      "210",   "220")
)

df$ft <- character(length = length(df$land_cover))
df$land_cover <- as.factor(df$land_cover)

for(lc in levels(df$land_cover)){
  ft_name <- as.character(ft.recl[ft.recl$LC_NUM == lc, "FT_NAME"])
  df[which(df$land_cover == lc), "ft"] <- ft_name
}

df$ft <- as.factor(df$ft)
summary(df$ft)

write.csv(df, "outputs/functional_types/ft_dataframe.csv", row.names = F)
print("outputs/outputs/functional_types/ft_dataframe.csv")
