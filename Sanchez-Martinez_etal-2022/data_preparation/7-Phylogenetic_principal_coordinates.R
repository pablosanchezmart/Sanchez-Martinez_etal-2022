### PHYLOGENETIC DISTANCES MATRIX DECOMPOSITION ------------------------------------------------------------------------------------------ ####

print("Calculating phylogenetic coordinates...")

source("code/data_preparation/0-data_preparation_init.R")

# P Sanchez-Marinez

ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_4_env_PCs.csv"), header = T)

### DATA NOT IN PHYLO

# Genus tree
# gen.tr <- read.tree("data/phylogeny/treevol_phylogeny_genus_lvl_2020.tre") # old Genus phylo used
gen.tr <- read.tree("data/phylogeny/treevol_phylogeny_genus_lvl_2021.tre") # Genus

genus <- levels(as.factor(ht$genus))
length(genus) # 4229 genera in data
length(gen.tr$tip.label) # 3449 genera in phylo (3488)

genus_notInPhylo <- genus[which(!genus %in% gen.tr$tip.label)]
length(genus_notInPhylo) # 784 genus in data and not in phylo (741)

species_not_in_genusPhylo <- ht[which(ht$genus %in% genus_notInPhylo), "Species"]
length(species_not_in_genusPhylo) # 1804 species in data which genus is not in genusPhylo (1745)

a <- ht[which(ht$genus %in% genus_notInPhylo), ]
lost_P50 <- a %>% filter(!is.na(P50)) %>% dplyr::select(Species)
length(lost_P50$Species) #146 P50 values (114)
lost_Pmin <- a %>% filter(!is.na(MinWP_md)) %>% dplyr::select(Species)
length(lost_Pmin$Species) #44 P50 values (36)

# write.csv(genus, "data/phylogeny/genus_list_23-04-2021.csv")
# print("data/phylogeny/genus_list_23-04-2021.csv")

# write.csv(genus_notInPhylo, "data/phylogeny/genus_list_not_in_phylo_23-04-2021.csv")
# print("data/phylogeny/genus_list_not_in_phylo_23-04-2021.csv")

# Genus level phylogeny (run only once)
if(!file.exists(paste0(processed.data, "treevol_workflow/phylo_eigenvectors.csv"))){
    mat <- cophenetic.phylo(gen.tr)
    mat <- as.data.frame(mat)
    pCordA <- pcoa(mat)
    
    phyloEigenV <- as.data.frame(pCordA$vectors)
    phyloEigenV$genus <- rownames(phyloEigenV)
    
    write.csv(phyloEigenV, paste0(processed.data, "treevol_workflow/phylo_eigenvectors.csv"), row.names = F)
    print(paste0(processed.data, "treevol_workflow/phylo_eigenvectors.csv"))
} else {
    # Phylogenetic principal coordinates merge
    phyloEigenV <-  read.csv(paste0(processed.data, "treevol_workflow/phylo_eigenvectors.csv"), header = T)   
}

ht_phylo <- merge(ht, phyloEigenV[, c(1:100, length(phyloEigenV))], by = "genus", all.x = T)
names(ht_phylo) <- str_replace_all(names(ht_phylo), "Axis.", "Phylo_axis_")

ht_phylo <- ht_phylo[, c("Species", names(ht_phylo)[-2])]
head(ht_phylo)
length(ht_phylo$Species) # (47197)


write.csv(ht_phylo, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_5_phylo_PCs.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_5_phylo_PCs.csv"))
