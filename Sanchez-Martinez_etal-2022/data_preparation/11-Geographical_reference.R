#### INCLUDE SPECIES GEOGRAPHICAL RANGE SPECIES ----------------------------------------------------------------------------------------------------------------------------------------- ####

remove(list = ls())


print("Spatially referenced version of the database...")

source("code/data_preparation/0-data_preparation_init.R")

#### FUNCITONS ---------------------------------------------------------------------------------------------------------------------------------------- ####

st_to_sfFun <- function(polygons.df = all_species_traits){
  sf <- sf::st_as_sf(polygons.df[1, ])
  for(n in 2:length(polygons.df$Species)){
    print(n)
    sf[n, ] <- sf::st_as_sf(polygons.df[n, ])
  }
  object.size(sf)
  return(sf)
}

if(!file.exists(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021_georeferenced.RData"))){
#### DATA --------------------------------------------------------------------------------------------------------------------------------------------- ####

# Traits and predictors
ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021.csv"), header = T)

# Load species distributions
load(paste0(processed.data, "species_distributions/all_species_2020.r"))

#### GEOGREAPHICAL DATA PREPARATION DATA ---------------------------------------------------------------------------------------------------------------------------------------------- ####

# Main projection for results
behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs") # Projection

# Removing problematic polygons

n <- which(new_all_species$Species == "Maoutia_australis") # Give problems, thre is a line!, sp number 26624
plot(new_all_species[new_all_species$Species == "Maoutia_australis",])
length(new_all_species)
new_all_species <- new_all_species[-n,]
length(new_all_species)

new_all_species_polygons <- st_to_sfFun(polygons.df = new_all_species)
new_all_species_polygons <- sf::st_make_valid(new_all_species_polygons)
new_all_species_polygons
a <- new_all_species_polygons

length(new_all_species_polygons$sp)

names(new_all_species_polygons)[1] <- "sp"

save(new_all_species_polygons, file = paste0(processed.data, "species_distributions/all_species_2020_sf.RData"))
print(paste0(processed.data, "species_distributions/all_species_2020_sf.RData"))

#### MERGE ---------------------------------------------------------------------------------------------------------------------------------------------- ####

ht$sp <- as.character(ht$sp)
new_all_species_polygons$sp <- as.character(new_all_species_polygons$sp)

te <- merge(new_all_species_polygons, ht, by = "sp", all = T)
length(te$sp) # 47058 (47462)

## Some checkings

a <- te %>% filter(!is.na(MinWP_md))
length(a$Species) # 536
a <- te %>% filter(!is.na(P50))
length(a$Species) # 863

save(te, file = paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021_georeferenced.RData"))
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021_georeferenced.RData"))
}

