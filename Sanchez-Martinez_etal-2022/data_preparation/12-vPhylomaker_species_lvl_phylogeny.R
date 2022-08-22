#### VPHYLOMAKER SPECIES LEVEL PHYLOGENY --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #### 

print("Species level phylogeny for plotting...")

source("code/data_preparation/0-data_preparation_init.R")

if(!file.exists())
#### Vphylomaker ####

load(paste0(processed.data, "treevol_globalDatabase_2021_georeferenced.RData")) # object called te

# phylomaker phylgeny 
splist <- data.frame("species" = te$Species, "genus" = te$genus, "family" = te$family)

for(i in 1:length(splist$species)){
  splist$genus[i] <- unlist(str_split(splist$species[i], " "))[1]
}

length(splist$species)
splist <- splist %>% filter(!is.na(splist$species)) # we have some species distributions not present in ht dataset (so we don't really need them for now, but I included them for future analyses)

ht.tree <- phylo.maker(splist, tree = GBOTB.extended, scenarios = "S1", r = 1)
sp_pm.tr <- ht.tree$scenario.1

sp_pm.tr$Nnode # 14946 nodes
length(sp_pm.tr$tip.label) # 45455 tips, which means a lot of polytomies. But we are interested in the topology at least, at the genus level, so it will work.

sp_pm.tr$node.label <- NULL


write.tree(sp_pm.tr, paste0(processed.data, "phylogeny/treevol_phylogeny_species_lvl_2021.tre"))
print(paste0(processed.data, "phylogeny/treevol_phylogeny_species_lvl_2021.tre"))
