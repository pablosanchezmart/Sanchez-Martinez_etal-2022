library(XML)
library(rgdal)
library(gdalUtils)


# bb=c(-337500.000,1242500.000,152500.000,527500.000) # Example bounding box (homolosine) for Ghana
sg_url="/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"
igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection

vars <- c("clay", "sand", "soc","cfvo", "bdod")
depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")

var <- "clay"
depth = "0-5cm"


for(var in vars){
  for(depth in depths){
    voi_layer = paste0(var,"_", depth, "_mean") # layer of interest

    gdal_translate(paste0(sg_url,'clay/clay_0-5cm_mean.vrt'),
                   paste0("C:/Users/p.sanchez/OneDrive - CREAF/Doctorado/Papers/Sanchez_etal_2020/predictingFromPhylo/data/species_data/soil_data/sg2020/", voi_layer, ".tif"),
                   s_src=igh,
                   verbose=TRUE)
    
    gdalwarp(paste0("/data/species_data/soil_data/sg2020/", voi_layer, ".tif"),
             paste0("/data/species_data/soil_data/sg2020/", voi_layer, ".tif"), 
             s_src=igh, 
             t_srs="EPSG:4326")
  }
}

