#### SPECIES NAMES CHECKING AND TAXON SCRUBBING ---------------------------------------------------------------------------------------------- ####

## P Sanchez-Martinez

print("Taxon scrubbing and taxonomic data filling ...")

source("code/data_preparation/0-data_preparation_init.R")


### FUNCTIONS ------------------------------------------------------------------------------------------------ ####

## Mode function (for categorical variables)
mode.fun <- function(v) {
  v <- as.character(v)
  uniqv <- unique(v)
  mode <- uniqv[which.max(tabulate(match(v, uniqv)))]
  if(is.na(mode)){                                            # To avoid NAs
    uniqv <- uniqv[which(!is.na(uniqv))]
    mode <- uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  return(mode)
}

# Maximum, minimum and mean function to avoid problems with NAs. 
mean.fun <- function(x){
  y <- mean(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}

# Combination function
comb.fun <- function(df, var){
  var.x <- paste(var, ".x", sep = "")
  var.y <- paste(var, ".y", sep = "")
  # choose default value from 'x' in the case that both '.x' and '.y' have values
  df[!is.na(df[, var.x]) & !is.na(df[, var.y]), var] <- df[!is.na(df[, var.x]) & !is.na(df[, var.y]), var.x]
  # replace NA values in 'x' with values from 'y' and vice-versa
  df[is.na(df[, var.x]), var] <- df[is.na(df[, var.x]), var.y]
  df[is.na(df[, var.y]), var] <- df[is.na(df[, var.y]), var.x]
  # Delete ".x" and ".y" columns
  df[, var.x] <- NULL
  df[, var.y] <- NULL
  return(df)
}

#### DATA ---------------------------------------------------------------------------------------------------- ####

ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_1.csv"), header = T)
# tr <- read.tree("data/phylogeny/bigphylo2.tre")   #not used for now, as it is just for plotting, we build it from scractch by using v.phylomaker package
load("data/species_distributions/all_species.r")


### SPECIES NOMENCLATURE ------------------------------------------------------------------------------------- ####

sp <- str_replace_all(ht$Species, "_", " ")
tr_sp <- str_replace_all(tr$tip.label, "_", " ")
dist_sp <- str_replace_all(all_species$SOURCE_SHP, "_", " ")

if(file.exists(paste0(processed.data, "RData/ht_2020_taxon_scrubing.RData"))){
  load(paste0(processed.data, "RData/ht_2020_taxon_scrubing.RData"))
} else {
  htSpp <- TPL(sp) # Look for new names in The Plant List
  save.image(file = paste0(processed.data, "RData/ht_2020_taxon_scrubing.RData"))
  print(paste0(processed.data, "RData/ht_2020_taxon_scrubing.RData"))
}

if(file.exists(paste0(processed.data, "RData/phylogeny_2020_taxon_scrubing.RData"))){
  load(paste0(processed.data, "RData/phylogeny_2020_taxon_scrubing.RData"))
} else {
  trSpp <- TPL(tr_sp) # Look for new names in The Plant List
  save.image(file = paste0(processed.data, "RData/phylogeny_2020_taxon_scrubing.RData"))
  print(paste0(processed.data, "RData/phylogeny_2020_taxon_scrubing.RData"))
}

if(file.exists(paste0(processed.data, "RData/distributions_2020_taxon_scrubing.RData"))){
  load(paste0(processed.data, "RData/distributions_2020_taxon_scrubing.RData"))
} else {
  distSpp <- TPL(dist_sp)
  save.image(file = paste0(processed.data, "RData/distributions_2020_taxon_scrubing.RData"))
  print(paste0(processed.data, "RData/distributions_2020_taxon_scrubing.RData"))
}

ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_1.csv"), header = T)

# New names (corrected)

ht$Species.x <- paste(htSpp$New.Genus, htSpp$New.Species, sep="_") # Add new names to dataframe
ht$Species.x[which(ht$Species.x == "NA_NA")] <- NA
ht$Species.y <- ht$Species
head(ht)
ht$Species <- NULL

ht <- comb.fun(ht, var = "Species")
ht$Species <- str_replace(ht$Species, "_", " ")
ht <- ht[, c("Species", names(ht)[-length(ht)])]

### Manually edit ####

tr_noVald <- trSpp %>% filter(Plant.Name.Index == FALSE)
length(tr_noVald$Taxon)

ht_noVald <- htSpp %>% filter(Plant.Name.Index == FALSE)
length(ht_noVald$Taxon)
write.csv(ht_noVald[, c("Taxon", "Genus", "Species")], paste0(processed.data, "treevol_workflow/species_taxonScrub_nameIndexF.csv"), row.names = F)
print(paste0(processed.data, "RData/species_taxonScrub_nameIndexF.csv"))

dist_noVald <- distSpp %>% filter(Plant.Name.Index == FALSE)
length(dist_noVald$Taxon)
write.csv(dist_noVald[, c("Taxon", "Genus", "Species")], paste0(processed.data, "species_distributions/species_taxonScrub_nameIndexF.csv"), row.names = F)
print(paste0(processed.data,"species_distributions/species_taxonScrub_nameIndexF.csv"))


ht$Species[ht$Species == "Alnus viridis_sinuata"] <-  "Alnus viridis_sinuata"
ht$Species[ht$Species == "Baccaurea ramilfora"] <-  "Baccaurea ramiflora"
ht$Species[ht$Species == "Beureria cumanensis"] <- "Bourreria cumanensis"
ht$Species[ht$Species == "Brachyleana neriifolia"] <- "Brachylaena neriifolia"
ht$Species[ht$Species == "Capparis aristiguetae"] <- "Capparis aristigueta"
ht$Species[ht$Species == "Craibiodendron scleranthum_kwangtungense"] <- "Caryodendron scleranthum_kwangtungense"
ht$Species[ht$Species == "Cyclobalanopsis myrsinaefolia"] <- "Cyclobalanopsis myrsinifolia"
ht$Species[ht$Species == "Dryandra cirisoides"] <- "Dryandra cytisoides"
ht$Species[ht$Species == "Eucalyptus grandis_camaldulensis"] <- "Eucalyptus grandis camaldulensis"
ht$Species[ht$Species == "Eucalyptus grandis_urophylla"] <- "Eucalyptus grandis urophylla"
ht$Species[ht$Species == "Fagus grandiflora"] <- "Fagus grandifolia"
ht$Species[ht$Species == "Homolanthus novoguinensis"] <- "Homalanthus novoguineensis"
ht$Species[ht$Species == "Lithocarpus glabra"] <- "Lithocarpus glaber"
ht$Species[ht$Species == "Magnoliaceae glanca"] <- "Magnoliaceae blanca"
ht$Species[ht$Species == "Metasequoia glytpstroboides"] <- "Metasequoia glyptostroboides"
ht$Species[ht$Species == "Neea steinbachii"] <- "Weingartia steinbachii"
ht$Species[ht$Species == "Ocotea sp.2"] <- "Ocotea sp."
ht$Species[ht$Species == "Picea sitkensis"] <- "Picea sitchensis"
ht$Species[ht$Species == "Pinus flexibilis"] <- "Pinus flexilis"
ht$Species[ht$Species == "Pinus tabulaeformis"] <- "Pinus tabuliformis"
ht$Species[ht$Species == "Populus deltoides_petrowskyana1"] <- "Populus deltoides petrowskyana"
ht$Species[ht$Species == "Populus deltoides_petrowskyana2"] <- "Populus deltoides petrowskyana"
ht$Species[ht$Species == "Populus deltoides_trichocarpa1"] <- "Populus deltoides trichocarpa"
ht$Species[ht$Species == "Populus deltoides_trichocarpa2"] <- "Populus deltoides trichocarpa"
ht$Species[ht$Species == "Populus laurifolia_nigra"] <- "Populus laurifolia nigra"
ht$Species[ht$Species == "Salix spp"] <- "Salix sp."
ht$Species[ht$Species == "Scyadopytis verticellata"] <- "Sciadopitys verticillata"
ht$Species[ht$Species == "Taxus breviflora"] <- "Taxus brevifolia"
ht$Species[ht$Species == "Thuya plicata"] <- "Thuja plicata"
ht$Species[ht$Species == "Ulmus davidiana_japonica"] <- "Ulmus davidiana japonica"
ht$Species[ht$Species == "Vochysia tyrsoideae"] <- "Vochysia thyrsoidea"

ht$Species[ht$Species == "Pinus tabulaeformis"] <-  "Pinus tabulaeformis"
ht$Species[ht$Species == "Pterospermum lanceaefolium"] <-  "Pterospermum lanceifolium"
ht$Species[ht$Species == "Picea sitkaensis"] <-  "Picea sitchensis"
ht$Species[ht$Species == "Thuya plicata"] <-  "Thuja plicata"
ht$Species[ht$Species == "Oxydendron arboreum"] <-  "Oxydendrum arboreum"
ht$Species[ht$Species == "Fagus grandiflora"] <-  "Fagus grandifolia"
ht$Species[ht$Species == "Launea arborescens"] <-  "Launaea arborescens"
ht$Species[ht$Species == "Retrophyllum minor"] <-  "Retrophyllum minus"
ht$Species[ht$Species == "Beureria cumanensis"] <-  "Bourreria cumanensis"

ht$Species[ht$Species == "Brachyleana neriifolia"] <-  "Brachylaena neriifolia"
ht$Species[ht$Species == "Metasequoia glytpstroboides"] <- "Metasequoia glyptostroboides"
ht$Species[ht$Species == "Picea sitkensis"] <- "Picea sitchensis"
ht$Species[ht$Species == "Pinus flexibilis"] <- "Pinus flexilis"
ht$Species[ht$Species == "Thuya plicata"] <- "Thuja plicata"
ht$Species[ht$Species == "Populus FRI"] <- "Populus sp."
ht$Species[ht$Species == "Populus GAV"] <- "Populus sp."
ht$Species[ht$Species == "Populus TRI"] <- "Populus sp."
ht$Species[ht$Species == "Salix BJO"] <- "Salix sp."
ht$Species[ht$Species == "Salix DEL"] <- "Salix sp."
ht$Species[ht$Species == "Salix Q83"] <- "Salix sp."
ht$Species[ht$Species == "Eucalyptus grandis_camaldulensis"] <- "Eucalyptus grandis camaldulensis"
ht$Species[ht$Species == "Eucalyptus grandis_urophylla"] <- "Eucalyptus grandis urophylla"
ht$Species[ht$Species == "Dryandra cirisoides"] <- "Dryandra cytisoides"
ht$Species[ht$Species == "Scyadopytis verticellata"] <- "Sciadopitys verticillata"
ht$Species[ht$Species == "Prumnopytis ferruginea"] <- "Prumnopitys ferruginea"
ht$Species[ht$Species == "Retrophyllum minor"] <-  "Retrophyllum minus"
ht$Species[ht$Species == "Alnus viridis_sinuata"] <-  "Alnus viridis_sinuata"
ht$Species[ht$Species == "Vochysia tyrsoideae"] <- "Vochysia thyrsoidea"
ht$Species[ht$Species == "Beureria cumanensis"] <- "Bourreria cumanensis"
ht$Species[ht$Species == "Capparis aristiguetae"] <- "Capparis aristigueta"
ht$Species[ht$Species == "Taxus breviflora"] <- "Taxus brevifolia"
ht$Species[ht$Species == "Homolanthus novoguinensis"] <- "Homalanthus novoguineensis"
ht$Species[ht$Species == "Fagus grandiflora"] <-  "Fagus grandifolia"
ht$Species[ht$Species == "Neea steinbachii"] <- "Weingartia steinbachii"
ht$Species[ht$Species == "Populus laurifolia_nigra"] <- "Populus laurifolia nigra"
ht$Species[ht$Species == "Populus deltoides_petrowskyana1"] <- "Populus deltoides petrowskyana"
ht$Species[ht$Species == "Populus deltoides_petrowskyana2"] <- "Populus deltoides petrowskyana"
ht$Species[ht$Species == "Populus deltoides_trichocarpa1"] <- "Populus deltoides trichocarpa"
ht$Species[ht$Species == "Populus deltoides_trichocarpa2"] <- "Populus deltoides trichocarpa"
ht$Species[ht$Species == "Baccaurea ramilfora"] <-  "Baccaurea ramiflora"
ht$Species[ht$Species == "Craibiodendron scleranthum_kwangtungense"] <- "Caryodendron scleranthum_kwangtungense"
ht$Species[ht$Species == "Pterospermum lanceaefolium"] <-  "Pterospermum lanceifolium"
ht$Species[ht$Species == "Ocotea sp.2"] <- "Ocotea sp."
ht$Species[ht$Species == "Ulmus davidiana_japonica"] <- "Ulmus davidiana japonica"
ht$Species[ht$Species == "Magnoliaceae glanca"] <- "Magnoliaceae blanca"
ht$Species[ht$Species == "Cyclobalanopsis myrsinaefolia"] <- "Cyclobalanopsis myrsinifolia"
ht$Species[ht$Species == "Quercus phillyraeoides"] <- "Quercus phillyreoides"
ht$Species[ht$Species == "Lithocarpus glabra"] <- "Lithocarpus glaber"
ht$Species[ht$Species == "Pinus tabulaeformis"] <-  "Pinus tabulaeformis"

ht$Species[ht$Species == "Psuedotsuga menziesii"] <-  "Pseudotsuga menziesii"
ht$Species[ht$Species == "Nectandra purpurascens"] <-  "Nectandra purpurea"
ht$Species[ht$Species == "Brachyleana neriifolia"] <-  "Brachylaena neriifolia"
ht$Species[ht$Species == "Nageia fleureii"] <-  "Nageia fleurety"
ht$Species[ht$Species == "Metasequoia glytpstroboides"] <- "Metasequoia glyptostroboides"
ht$Species[ht$Species == "Picea sitkensis"] <- "Picea sitchensis"
ht$Species[ht$Species == "Pinus flexibilis"] <- "Pinus flexilis"
ht$Species[ht$Species == "Thuya plicata"] <- "Thuja plicata"
ht$Species[ht$Species == "Populus X"] <- "Populus sp."
ht$Species[ht$Species == "Salix X"] <- "Salix sp."
ht$Species[ht$Species == "Salix schwereniixviminalis"] <- "Salix viminalis"
ht$Species[ht$Species == "Populus trichocarpaxdeltoides"] <- "Populus trichocarpa deltoides"
ht$Species[ht$Species == "Eucalyptus grandis_camaldulensis"] <- "Eucalyptus grandis camaldulensis"
ht$Species[ht$Species == "Eucalyptus grandis_urophylla"] <- "Eucalyptus grandis urophylla"
ht$Species[ht$Species == "Dryandra cirisoides"] <- "Dryandra cytisoides"
ht$Species[ht$Species == "Scyadopytis verticellata"] <- "Sciadopitys verticillata"
ht$Species[ht$Species == "Juniperus virgiana"] <- "Juniperus virginiana"
ht$Species[ht$Species == "Prumnopytis ferruginea"] <- "Prumnopitys ferruginea"
ht$Species[ht$Species == "Retrophyllum minor"] <-  "Retrophyllum minus"
ht$Species[ht$Species == "Salix spp."] <-  "Salix sp."
ht$Species[ht$Species == "Vochysia tyrsoideae"] <- "Vochysia thyrsoidea"
ht$Species[ht$Species == "Protium sp"] <- "Protium sp."
ht$Species[ht$Species == "Beureria cumanensis"] <- "Bourreria cumanensis"
ht$Species[ht$Species == "Capparis aristiguetae"] <- "Capparis aristigueta"
ht$Species[ht$Species == "Lonchocarpus dipteronereus"] <- "Lonchocarpus dipteroneurus"
ht$Species[ht$Species == "Taxus breviflora"] <- "Taxus brevifolia"
ht$Species[ht$Species == "Homolanthus novoguinensis"] <- "Homalanthus novoguineensis"
ht$Species[ht$Species == "Oxydendron arboreum"] <-  "Oxydendrum arboreum"
ht$Species[ht$Species == "Pseudostuga menziesii"] <-  "Pseudotsuga menziesii"
ht$Species[ht$Species == "Fagus grandiflora"] <-  "Fagus grandifolia"
ht$Species[ht$Species == "Aspidosperma cylindrocarpum"] <-  "Aspidosperma cylindrocarpon"
ht$Species[ht$Species == "Neea steinbachii"] <- "Weingartia steinbachii"
ht$Species[ht$Species == "Launea arborescens"] <-  "Launaea arborescens"
ht$Species[ht$Species == "Populus laurifoliaxnigra"] <- "Populus laurifolia nigra"
ht$Species[ht$Species == "Populus balsamiferaxeuramericana"] <- "Populus balsamifera euroamericana"
ht$Species[ht$Species == "Populus balsamiferaxsimonii"] <- "Populus balsamifera simonii"
ht$Species[ht$Species == "Populus deltoidesxpetrowskyana"] <- "Populus deltoides petrowskyana"
ht$Species[ht$Species == "Baccaurea ramilfora"] <-  "Baccaurea ramiflora"


### DATABASE ------------------------------------------------------------------------------------------------- ####

# How much it changes?
dif <- ht %>% filter(sp != ht$Species) %>% dplyr::select(Species)
length(dif$Species) # 486 changes (77 comparing from hammond dataset)

# Duplicated species
ht_unique <- ht[!duplicated(ht$Species), ]
length(ht$Species) - length(ht_unique$Species) # 175 species duplicated

duplSpp <- duplicated(ht$Species)
dupl.df <- ht[duplicated(ht$Species), ]
dupl.df <- arrange(dupl.df, Species)

# Aggregate data
names(ht)
catVar <- c("order", "family", "genus", "group", "Growth.form", "life.form", "leaf.form", "biome")
numVar <- names(ht)[c(10:length(ht))]

for(var in numVar){
  ht[, var] <- as.numeric(ht[, var])
}

for(var in catVar){
  ht[, var] <- as.character(ht[, var])
}

ht <- ht[order(ht$Species), ]
head(ht)

ht_aggr_mean <- stats::aggregate(ht[, numVar], by = list(ht$Species), FUN = mean.fun)
ht_aggr_mode <- stats::aggregate(ht[, catVar], by = list(ht$Species), FUN = mode.fun)

ht_aggr <- merge(ht_aggr_mean, ht_aggr_mode, by = "Group.1")

names(ht_aggr)[names(ht_aggr) == "Group.1"] <- "Species"
length(ht_aggr$Species) # 47540

## Checking
head(ht_aggr)

ht[which(ht$Species == "Acer saccharum"), "MinWP_pd"]
ht_aggr[which(ht_aggr$Species == "Acer saccharum"),  "MinWP_pd"]

ht_tax <- arrange(ht_aggr, Species)
head(ht_tax)

### TAXONLOOKUP TO FILL TAXONOMIC GROUPS GAPS ---------------------------------------------------------------- ####

# New families (avoiding empty new data)

taxVars <- c("family", "order", "group")

# Genus

for(i in 1:length(ht_tax$Species)){
  ht_tax$genus[i] <- unlist(str_split(ht_tax$Species[i], " "))[1]
}

# Taxonlookup
species <- as.character(ht_tax$Species)
taxonomy <- lookup_table(species, by_species = TRUE, family.tax = "plantlist")
taxonomy <- cbind("Species" = rownames(taxonomy), data.frame(taxonomy, row.names=NULL))
taxonomy <- taxonomy %>% dplyr::select(! genus)

# Character
for(i in taxVars){
  ht_tax[, i] <- as.character(ht_tax[, i])
}

# Merge
ht_tax_compl <- merge(ht_tax, taxonomy, by = "Species", all.x = T) # priorizing what we already have
head(ht_tax_compl)

# Combinate
for(i in taxVars){
  ht_tax_compl <- comb.fun(df = ht_tax_compl, var =  i)
}

# Factors
for(i in c("genus", taxVars)){
  ht_tax_compl[, i] <- as.factor(ht_tax_compl[, i])
}

# Checking
summary(ht_tax_compl$genus) # 0 NAs
summary(ht_tax_compl$family) # 93 NAs
summary(ht_tax_compl$order) # 99 NAs
summary(ht_tax_compl$group) # 93 NAs

# Manual edit

ht_tax_compl$genus[ht_tax_compl$genus == "Dypterix"] <-  "Dipteryx"

# Order

ht_tax_compl <- ht_tax_compl[, names(ht)]
head(ht_tax_compl)
length(ht_tax_compl$Species) # (47197)

# Save

write.csv(ht_tax_compl, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_2_taxonomy_scrub.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_2_taxonomy_scrub.csv"))


#### SPECIES DISTRIBUTIONS ---------------------------------------------------------------------------------------------- ####

# New names (corrected)
distNewSpecies <- paste(distSpp$New.Genus, distSpp$New.Species, sep="_") # New names

length((which(!dist_sp %in% distNewSpecies))) # 473 differences
dist_unique <- distNewSpecies[!duplicated(distNewSpecies)]
length(distNewSpecies) - length(dist_unique) # 148 rows duplicated

# Aggregate data
all_species$Species <- distNewSpecies

agr_all_species <- raster::aggregate(all_species, by = "Species", dissolve = T)

length(all_species$Species) - length(agr_all_species$Species)

agr_all_species$Species[which(agr_all_species$Species == "Blumeodendron_kurzii")] # We know that this specie was duplicated, so we are checking here if it still in the database
new_all_species <- agr_all_species

# Save

save(new_all_species, file = paste0(processed.data, 'all_species_2020.r'))
print(paste0(processed.data, 'all_species_2020.r'))

save(new_all_species, file = paste0(processed.data, 'species_distributions/all_species_2020.r'))
print(paste0(processed.data, 'species_distributions/all_species_2020.r'))


#### PHYLOGENY ---------------------------------------------------------------------------------------------- ####

# New names (corrected)
head(trSpp)
trNewSpecies <- paste(trSpp$New.Genus, trSpp$New.Species, sep=" ") # New names

length((which(!tr_sp %in% trNewSpecies))) # 473 differences
tr_unique <- trNewSpecies[!duplicated(trNewSpecies)]
length(trNewSpecies) - length(tr_unique) # 148 rows duplicated

# Dataframe of old and new species
tr_species <- data.frame("Old_sp" = str_replace_all(tr_sp, " ", "_"),
                         "New_sp" = str_replace_all(trNewSpecies, " ", "_"))

# Species with changes in its specific epitetum (so they can still part of the phylogeny only changing its name)
genus <- str_extract(trSpp$Taxon, "(.*)(?=\\s)") # Better to use this compared to the genus column of the dataframe
trSpp$genus <- genus
trSpEpi <- trSpp[which(trSpp$genus == trSpp$New.Genus), ]
length(trSpEpi$Taxon) # 46434 to keep
trSpEpi <- paste(trSpEpi$New.Genus, trSpEpi$New.Species, sep = "_") # Species to keep in the phyogeny just changing names


### Keep only standard names not duplicated and from which changes at genus level has not been made

# Remove duplicates
new_tr <- tr
unique <- tr_species %>% filter(tr_species$New_sp %in% trSpEpi) %>% distinct(New_sp, .keep_all = T) # We want to keep those duplicates that has not changed their genus
length(unique$New_sp)

# Dropping duplicated species
new_tr <- drop.tip(new_tr, which(!new_tr$tip.label %in% unique$Old_sp))

# Bring tip names to a standard nomenclature dropping duplicates (species already in the phylogeny)

for(i in 1:length(new_tr$tip.label)){
  old_sp <- new_tr$tip.label[i]
  new_sp <- tr_species[which(tr_species$Old_sp == old_sp), "New_sp"]
  new_tr$tip.label[new_tr$tip.label == old_sp] <- new_sp
}

new_tr$tip.label[which(new_tr$tip.label == "Blumeodendron_kurzii")] # We know that this specie was duplicated, so we are checking here if it still in the phylogeny

# Remove species which genus has changed

new_tr <- drop.tip(new_tr, which(!new_tr$tip.label %in% str_replace_all(trSpEpi, " ", "_")))
new_tr # 46306
length(unique(new_tr$tip.label)) # No duplicates, god to go

write.tree(new_tr, "data/phylogeny/treevol_2021.tre")
print("data/phylogeny/treevol_2021.tre")

write.tree(new_tr, "data/treevol_2021.tre")
print("data/treevol_2021.tre")

