#### MORTALITY DATABASE  ----------------------------------------------------------------------------------------------------------------------------------------- ####

print("Mortality database preparation...")

source("code/data_preparation/0-data_preparation_init.R")

#### FUNCITONS ---------------------------------------------------------------------------------------------------------------------------------------- ####

st_to_sfFun <- function(polygons.df = all_species_traits){
  sf <- sf::st_as_sf(polygons.df[1, ])
  for(n in 2:length(polygons.df$SOURCE_SHP)){
    print(n)
    sf[n, ] <- sf::st_as_sf(polygons.df[n, ])
  }
  object.size(sf)
  return(sf)
}

#### DATA ----------------------------------------------------------------------------------------------------------------------------------------------- ####

# Main projection for results
behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs") # Projection

# Mortality data
md <- read.csv("data/mortality/Hammond_etal_database_withspecies.csv", header = T)

# Load species distributions with trait values
load(paste0(output.dir, "species_traits/all_species_traits_RFimp.RData")) # all_species_traits
  
#### MORTALITY DATABASE PREPARATION -------------------------------------------------------------------------------------------------------------------- ####

if(!file.exists("data/mortality/Hammond_etal_database_withspecies_transformed_species_wise_geo.RData")){
  
  
### Geographical database ------------------------------------------------------------------------------------ ####

### Change coordinates projection ###

md$id <- 1:length(md$ref_id)
md.sf <- md
coordinates(md.sf) <- md.sf[, c("long", "lat")]
proj4string(md.sf) <- CRS("+init=EPSG:4326")
summary(md.sf)

md.sf_trans <- spTransform(md.sf, behrmann)

md_coord <- as.data.frame(coordinates(md.sf_trans))
names(md_coord) <- c("x", "y")

md.sf_trans <- st_as_sf(md.sf_trans)
length(md.sf_trans$ref_id)
save(md.sf_trans, file = "data/mortality/Hammond_etal_database_withspecies_transformed_geo.RData")
print("data/mortality/Hammond_etal_database_withspecies_transformed.RData")

md <- cbind(md, md_coord)
md <- md[, c("id","ref_id", "long", "lat", "x", "y", "mortality_year", "species")]

write.csv(md, "data/mortality/Hammond_etal_database_withspecies_transformed.csv", row.names = F)
print("data/mortality/Hammond_etal_database_withspecies_transformed.csv")


### Species-wise database ------------------------------------------------------------------------------------ ####

md <- read.csv("data/mortality/Hammond_etal_database_withspecies_transformed.csv", header = T)

# Species extraction for those studies with any (many sp.)

many_sp_md <- readxl::read_excel("data/mortality/many_sp_ref/mortality_species_extraction_manySpp_PSM_RGV.xlsx")
length(unique(many_sp_md$ref_id)) # we extracted data on species suffering mortality from 9 studies from 23 with many spp.

many_sp_md <- distinct(many_sp_md)
length(many_sp_md$ref_id) # 544 species x site observations
length(unique(many_sp_md$species)) # 422 species

# Remove any duplicate row or observation

length(many_sp_md$ref_id) # 544 new observations at the species level (species x site combinations)
head(many_sp_md)
# coordinate transformation

coordinates(many_sp_md) <- many_sp_md[, c("long", "lat")]
proj4string(many_sp_md) <- CRS("+init=EPSG:4326")

many_sp_md_trans <- spTransform(many_sp_md, behrmann)

many_sp_md_trans <- as.data.frame(many_sp_md_trans)
many_sp_md_trans <- dplyr::rename(many_sp_md_trans, x = long.1, y = lat.1)

# To avoid to include duplicate data coming from the same study, we eliminate all observations for those studies having any "many sp." for which we extracted data

studiesWithSpData <- unique(many_sp_md$ref_id)
length(md$id)
unique(md$ref_id)
md <- md %>% filter(!ref_id %in% studiesWithSpData) %>% dplyr::select(!id) 
length(md$ref_id) # 1273 observations from 1303 comin from studies from which any many sp. observation was recorded

# Species wise (md)

md_species <- data.frame()
for(i in 1:length(md$species)){
  print(i)
  sp <- md$species[i]
  spp <- unlist(str_split(sp, pattern = ", "))
  md_spp <- cbind("species" = spp, md[i, c(-7)])
  md_species <- rbind(md_species, md_spp)
}

# Join new species (from many spp. with species data in their original source)

many_sp_md_trans <- many_sp_md_trans[, names(md_species)]
names(many_sp_md_trans)
length(md_species$species) + length(many_sp_md_trans$ref_id) # 3900 species x site in total (3356 + 544)

md_species_many_sp <- rbind(md_species, many_sp_md_trans)
length(md_species_many_sp$species)

md_species_many_sp <- md_species_many_sp %>% filter(species != "(many spp.)")
md_species_many_sp$species <- as.factor(md_species_many_sp$species)

length(md_species_many_sp$species) # 3809 species x site combinations (once remaining many sp. exculded)
length(unique(levels(md_species_many_sp$species))) # 570 different species

write.csv(md_species_many_sp, "data/mortality/Hammond_etal_database_withspecies_transformed_species_wise.csv", row.names = F)
print("data/mortality/Hammond_etal_database_withspecies_transformed_species_wise.csv")

# Taxon scrubing

md_sp <- str_replace_all(md_species_many_sp$species, "_", " ")
# Scrub_Spp <- TPL(md_sp) # Look for new names in The Plant List
# save(Scrub_Spp, file = "data/mortality/mortality_database_taxonScrub.RData")
# print("data/mortality/mortality_database_taxonScrub.RData")

load("data/mortality/mortality_database_taxonScrub.RData")
sp_noVald <- Scrub_Spp %>% filter(Plant.Name.Index == FALSE)
write.csv(sp_noVald[, c("Taxon", "New.Genus", "New.Species")], "data/mortality/mortality_spNoValid_taxonScrub.csv", row.names = F)
print("data/mortality/mortality_spNoValid_taxonScrub.csv")
length(sp_noVald$Taxon) # most of them are due to genus level data ("spp.")

md_species_many_sp$species <- paste(Scrub_Spp$New.Genus, Scrub_Spp$New.Species, sep=" ") # Add new names to dataframe
md_species_many_sp$sp <- str_replace_all(md_species_many_sp$species, " ", "_")

# How much it changes?

dif <- md_species_many_sp %>% filter(sp != md_sp) %>% dplyr::select(species) # no differences
length(dif$Species)
head(md_species_many_sp)

# Taxon lookup (genus)

taxVars <- c("family", "order", "group")

# Taxonlookup
species <- as.character(unique(md_species_many_sp$species))
taxonomy <- lookup_table(species, by_species = TRUE, family.tax = "plantlist")
taxonomy <- cbind("species" = rownames(taxonomy), data.frame(taxonomy, row.names=NULL))
head(taxonomy)

md_species_many_sp <- merge(md_species_many_sp, taxonomy, by = "species", all.x = T)

# Save

length(md_species_many_sp$species) # 3809 species x site
length(unique(md_species_many_sp$species)) # 554 different species

write.csv(md_species_many_sp, "data/mortality/Hammond_etal_database_withspecies_transformed_species_wise.csv", row.names = F)
print("data/mortality/Hammond_etal_database_withspecies_transformed_species_wise.csv")

# geographically explicit species-wise database
md_species.sf <- md_species_many_sp
coordinates(md_species.sf) <- md_species_many_sp[, c("x", "y")]
proj4string(md_species.sf) <- behrmann
summary(md_species.sf)

md_species.df <- as.data.frame(md_species.sf)
md_species.df <- md_species.df %>% select(species, genus, family, order, group, long, lat)
head(md_species.df)
} else {
  load(file = "data/mortality/Hammond_etal_database_withspecies_transformed_species_wise_geo.RData")
}


#### INTERSECTION BETWEEN SPECIES AND MORTALITY POINTS (SPECIES PRESENT IN EACH MORTAILITY POINT) (not used for now) ------------------------------------------------- ####

# Intersection between points and polygons (extract each polygon attributes for each point, ref_id would allow to summaryze data by point if needed)

md.sf_trans <- st_as_sf(md.sf_trans)
md_all_species_traits <- st_intersection(md.sf_trans, all_species_traits)

md_all_species_traits$sp <- NULL
dim(md_all_species_traits)
md_all_species_traits <- dplyr::rename(md_all_species_traits, species_mortality = species, species_range = Species, species_phylogeny = SOURCE_SHP)

save(md_all_species_traits, file = "outputs/mortality/spatial_mortality_traits.RData") # object name: md_all_species_traits


load("outputs/mortality/spatial_mortality_traits.RData")
md_all_species_traits

print("outputs/mortality/spatial_mortality_traits.RData")

summary(md_all_species_traits)
head(md_all_species_traits)
length(md_all_species_traits$ref_id)

# List of species range per point

md_all_species_traits.df <- as.data.frame(md_all_species_traits)
species_range_byPoints <- as.list(tapply(X = md_all_species_traits.df$species_range, INDEX = md_all_species_traits.df$id, FUN = c))

# Check
# species_range_byPoints$`1`
# unique(md_all_species_traits.df[md_all_species_traits.df$id == 1, "species_range"])

species_range_byPoints.df <- data.frame()
for(i in 1:length(species_range_byPoints)){
  species_range_byPoints[[i]] <- split(species_range_byPoints[[i]], "species_range")
  df <- as.data.frame(species_range_byPoints[[i]]$species_range)
  df$id <- names(species_range_byPoints[i])
  names(df)[1] <- "Species_range"
  species_range_byPoints.df <- cbind(rbind(species_range_byPoints.df, df))
}

length(species_range_byPoints.df$Species_range)

species_range_byPoints_md.df <- join(x = species_range_byPoints.df, y = md, by = "id")
head(species_range_byPoints_md.df)
species_range_byPoints_md.df <- dplyr::rename(species_range_byPoints_md.df, species_mortality = species, point_id = id)
head(species_range_byPoints_md.df)

# Save dataframe of species range present in each point
write.csv(species_range_byPoints_md.df, "outputs/mortality/species_range_list_per_mortality_point.csv", row.names = F)
print("outputs/mortality/species_range_list_per_mortality_point.csv")

# Save array with species present in each point (to handle in R)
save(species_range_byPoints, file = "outputs/mortality/species_range_list_per_mortality_point.RData")
print("outputs/mortality/species_range_list_mortality.RData")

}