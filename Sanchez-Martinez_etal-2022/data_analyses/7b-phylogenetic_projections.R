##### PHYLOGENETIC PROJECTION ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #####

remove(list = ls())
gc()
# install.packages("ggplot2", dependencies = T)

#### PACKAGES ---------------------------------------------------------------------------------------------------------------------- ####

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

#### FUNCTIONS ------------------------------------------------------------------------------------------------------------------- ####

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

display.brewer.all(colorblindFriendly = TRUE)

#### DATA AND PARAMETERS --------------------------------------------------------------------------------------------------------- ####

h <- 7.87
w <- 7.87

h <- 7.87 * 1.5
w <- 7.87

predictors <-  c(paste0("PC", 1:5), paste0("Phylo_axis_", 1:5))

# Mortality data species-wise

md_species <- read.csv("data/mortality/Hammond_etal_database_withspecies_transformed_species_wise.csv", header = T)
md_species$species <- as.factor(md_species$species)
length(unique(md_species$species))
md_species_count <- as.data.frame(table(md_species$species)) 
md_species_count
names(md_species_count) <- c("Species", "mortality_events")

#### OBSERVED DATA ---------------------------------------------------------------------------------------------------------------- ####

obsTraits <- read.csv(paste0(output.dir, "RF_imputations/RFimputation_obsData_100_bi_climate_5_phylo_5_group.csv"), header = T) # 1051 species
obsTraits$sp <- str_replace_all(obsTraits$Species, " ", "_")

obsTraits <- join(obsTraits, md_species_count, by = "Species", type = "left")
obsTraits[which(is.na(obsTraits$mortality_events)), "mortality_events"] <- 0
obsTraits[which(obsTraits$mortality_events > 0), "mortality_binomial"] <- 1
obsTraits[which(is.na(obsTraits$mortality_binomial)), "mortality_binomial"] <- 0

sum(obsTraits$mortality_binomial) # 122 (145) species with mortality events for observed data

#### Observed Data genus-level data ----------------------------------------------------------------------------------------------------------------- ####

# Data intersection

tstGen.df <- obsTraits[obsTraits$genus %in% gen.tr$tip.label, ]
length(tstGen.df$genus) # 1051 species in common

# Genus level aggregation

names(tstGen.df)
tstGen1.df <- aggregate(tstGen.df[, c("MinWP_md", "MinWP_md_sd", "MinWP_md_cv", "P50","P50_sd", "P50_cv", "HSM", "HSM_sd",
                                      "HSM_cv", "MinWP_soil_SX", predictors)], 
                        by = list(tstGen.df$genus), FUN = mean)
tstGen2.df <- aggregate(tstGen.df[, c("mortality_events")], by = list(tstGen.df$genus), FUN = sum)
names(tstGen2.df)[2] <- "mortality_events"
tstGen3.df <- aggregate(tstGen.df[, c("mortality_binomial")], by = list(tstGen.df$genus), FUN = max)
names(tstGen3.df)[2] <- "mortality_binomial"
tstGen4.df <- aggregate(tstGen.df[, c("family", "order", "group")], by = list(tstGen.df$genus), FUN = mode.fun)
tstGen.df <- join(tstGen1.df, tstGen2.df, by = "Group.1")

tstGen.df <- join(tstGen.df, tstGen3.df, by = "Group.1")
tstGen.df <- join(tstGen.df, tstGen4.df, by = "Group.1")
head(tstGen.df)
tstGen.df <- dplyr::rename(tstGen.df, "genus" = "Group.1")

length(tstGen.df$genus) # 498

tstGen.df$genus <- as.factor(tstGen.df$genus)
tstGen.df$order <- as.factor(tstGen.df$order)
tstGen.df$fmly <- as.factor(tstGen.df$family)
tstGen.df$group <- as.factor(tstGen.df$group)
summary(tstGen.df)

tstGen.tr <- drop.tip(gen.tr, gen.tr$tip.label[!gen.tr$tip.label %in% tstGen.df$genus])
rownames(tstGen.df) <- tstGen.df$genus
tstGen.df <- tstGen.df[tstGen.tr$tip.label,]
tstGen.df <- tstGen.df[, c("genus", names(tstGen.df))] # we need it as a first column

#### OD genus level trait vertical plot ------------------------------------------------------------------------- ####

mean.val <- tstGen.df %>% filter(group == "Angiosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>% 
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
             group_by(mortality_binomial) %>%
               summarize_all(., mean)
mean.val

mean.val <- tstGen.df %>% filter(group == "Gymnosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>% 
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., mean)
mean.val

pdf(paste0(rslts.dir, "figures/phylogeny/obsData_genus_mean_values_rf.pdf"), height = h, width = w)
plotPhyloFun(phylo = tstGen.tr, traits = tstGen.df, brwPal = "BrBG", FontsizeFactor = 0.1)
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/obsData_genus_mean_values_rf.pdf"))


pdf(paste0(rslts.dir, "figures/phylogeny/obsData_genus_sd_values_rf.pdf"), height = h, width = w)
plotPhyloSDFun(phylo = tstGen.tr, traits = tstGen.df, brwPal = "BrBG", FontsizeFactor = 0.1)
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/obsData_genus_sd_values_rf.pdf"))

### Species level ###

#### OD species-level data ----------------------------------------------------------------------------------------------------------------- ####

# Data intersection

tstSp.df <- obsTraits[obsTraits$sp %in% sp.tr$tip.label, ]
length(tstSp.df$Species) # 789 species in common

tstSp.df$genus <- as.factor(tstSp.df$genus)
tstSp.df$order <- as.factor(tstSp.df$order)
tstSp.df$fmly <- as.factor(tstSp.df$family)
tstSp.df$group <- as.factor(tstSp.df$group)
summary(tstSp.df)

tstSp.tr <- drop.tip(sp.tr, sp.tr$tip.label[!sp.tr$tip.label %in% tstSp.df$sp])
rownames(tstSp.df) <- tstSp.df$sp
tstSp.df <- tstSp.df[tstSp.tr$tip.label,]
tstSp.df <- tstSp.df[, c("sp", names(tstSp.df))] # we need it as a first column


#### OD species-level trait vertical plot ------------------------------------------------------------------------------------------ ####

pdf(paste0(rslts.dir, "figures/phylogeny/obsData_species_values_rf.pdf"), height = h, width = w)
plotPhyloFun(phylo = tstSp.tr, traits = tstSp.df, brwPal = "BrBG", FontsizeFactor = 0.05)
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/obsData_species_values_rf.pdf"))

wholeMean <- tstSp.df %>% filter(group == "Angiosperms")
mean(wholeMean$HSM)

tstSp.df$N_neg_HSM <- ifelse(tstSp.df$HSM < 0, "negative HSM", "positive HSM")
tstSp.df$N_low05_HSM <- ifelse(tstSp.df$HSM < 0.5, "HSM < 0.5", "HSM > 0.5")
tstSp.df$N_low1_HSM <- ifelse(tstSp.df$HSM < 1, "HSM < 1", "HSM > 1")

mean.val <- tstSp.df %>% filter(group == "Angiosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>%
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., mean)
mean.val

sd.val <- tstSp.df %>% filter(group == "Angiosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>%
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., sd)
sd.val

mean.val <- tstSp.df %>% filter(group == "Gymnosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>%
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., mean)
mean.val

sd.val <- tstSp.df %>% filter(group == "Gymnosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>%
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., sd)
sd.val

neg_HSM <- tstSp.df %>% dplyr::select(N_neg_HSM, mortality_binomial) %>%
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., count)

N_by_order <- table(tstSp.df$order, tstSp.df$mortality_binomial)
N_by_order

N_by_group <- table(tstSp.df$group, tstSp.df$mortality_binomial)
N_by_group

# Standard deviation
pdf(paste0(rslts.dir, "figures/phylogeny/obsData_species_sd_rf.pdf"), height = h, width = w)
plotPhyloSDFun(phylo = tstSp.tr, traits = tstSp.df, brwPal = "BrBG", FontsizeFactor = 0.05, 
             problematicSpp = c("Fagus_sylvatica", "Ficus_variegata"))
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/obsData_species_sd_rf.pdf"))


#### IMPUTED DATA ---------------------------------------------------------------------------------------------------------------- ####

allTraits <- read.csv(paste0(output.dir, "RF_imputations/RFimputation_100_bi_climate_5_phylo_5_group.csv"), header = T) # 45345 species
allTraits$sp <- str_replace_all(allTraits$Species, " ", "_")

allTraits <- join(allTraits, md_species_count, by = "Species", type = "left")
allTraits[which(is.na(allTraits$mortality_events)), "mortality_events"] <- 0
allTraits[which(allTraits$mortality_events > 0), "mortality_binomial"] <- 1
allTraits[which(is.na(allTraits$mortality_binomial)), "mortality_binomial"] <- 0
sum(allTraits$mortality_events)

#### ID genus-level data --------------------------------------------------------------------------------------------------------- ####

# Data intersection

tstGen.df <- allTraits[allTraits$genus %in% gen.tr$tip.label, ]
length(tstGen.df$Species) # 44697 species in common

# Genus level aggregation

names(tstGen.df)
tstGen1.df <- aggregate(tstGen.df[, c("MinWP_md", "P50", "HSM", "MinWP_soil_SX")], 
                        by = list(tstGen.df$genus), FUN = min)
tstGen2.df <- aggregate(tstGen.df[, c("MinWP_md_sd", "MinWP_md_cv","P50_sd", "P50_cv", "HSM_sd",  "HSM_cv", predictors)], 
                        by = list(tstGen.df$genus), FUN = mean)
tstGen3.df <- aggregate(tstGen.df[, c("mortality_events")], by = list(tstGen.df$genus), FUN = sum)
names(tstGen3.df)[2] <- "mortality_events"
tstGen4.df <- aggregate(tstGen.df[, c("mortality_binomial")], by = list(tstGen.df$genus), FUN = max)
names(tstGen4.df)[2] <- "mortality_binomial"
tstGen5.df <- aggregate(tstGen.df[, c("family", "order", "group")], by = list(tstGen.df$genus), FUN = mode.fun)

tstGen.df <- join(tstGen1.df, tstGen2.df, by = "Group.1")
tstGen.df <- join(tstGen.df, tstGen3.df, by = "Group.1")
tstGen.df <- join(tstGen.df, tstGen4.df, by = "Group.1")
tstGen.df <- join(tstGen.df, tstGen5.df, by = "Group.1")
names(tstGen.df)

tstGen.df <- dplyr::rename(tstGen.df, genus = Group.1)

length(tstGen.df$genus) # 3359

tstGen.df$order <- as.factor(tstGen.df$order)
tstGen.df$fmly <- as.factor(tstGen.df$family)
head(tstGen.df)

tstGen.tr <- drop.tip(gen.tr, gen.tr$tip.label[!gen.tr$tip.label %in% tstGen.df$genus])
rownames(tstGen.df) <- tstGen.df$genus
tstGen.df <- tstGen.df[tstGen.tr$tip.label,]

#### ID trait vertical plot ---------------------------------------------------------------------------------- ####

pdf(paste0(rslts.dir, "figures/phylogeny/ImpData_genus_mean_values_rf.pdf"), height = h, width = w)
plotPhyloFun(phylo = tstGen.tr, traits = tstGen.df, brwPal = "BrBG", FontsizeFactor = 0.01)
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/ImpData_genus_mean_values_rf.pdf"))


pdf(paste0(rslts.dir, "figures/phylogeny/ImpData_genus_SD_values_rf.pdf"), height = h, width = w)
plotPhyloSDFun(phylo = tstGen.tr, traits = tstGen.df, brwPal = "BrBG", FontsizeFactor = 0.01)
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/ImpData_genus_SD_values_rf.pdf"))

mean.val <- tstGen.df %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>% 
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., mean)
mean.val

mean.val <- tstGen.df %>% filter(group == "Angiosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>% 
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., mean)
mean.val

mean.val <- tstGen.df %>% filter(group == "Gymnosperms") %>% dplyr::select(MinWP_md, P50, HSM, mortality_binomial) %>% 
  mutate(mortality_binomial = factor(mortality_binomial)) %>%
  group_by(mortality_binomial) %>%
  summarize_all(., mean)
mean.val


#### ID species-level data only with  mortality ---------------------------------------------------------------------------------- ####

tstSp.df <- allTraits[allTraits$sp %in% sp.tr$tip.label, ]
length(tstSp.df$Species) # 45043 species in common

tstSp.df$genus <- as.factor(tstSp.df$genus)
tstSp.df$order <- as.factor(tstSp.df$order)
tstSp.df$fmly <- as.factor(tstSp.df$family)
tstSp.df$group <- as.factor(tstSp.df$group)
summary(tstSp.df)

tstSp.tr <- drop.tip(sp.tr, sp.tr$tip.label[!sp.tr$tip.label %in% tstSp.df$sp])
rownames(tstSp.df) <- tstSp.df$sp
tstSp.df <- tstSp.df[tstSp.tr$tip.label,]
tstSp.df <- tstSp.df[, c("sp", names(tstSp.df))] # we need it as a first column

# only with observed mortality
tstSp.df_md <- tstSp.df[which(tstSp.df$mortality_binomial == 1), ]
length(tstSp.df_md$genus) # 139 (478) species
tstSp.df_md.tr <- drop.tip(sp.tr, sp.tr$tip.label[!sp.tr$tip.label %in% tstSp.df_md$sp])
rownames(tstSp.df_md) <- tstSp.df_md$sp

tstSp.df_md <- tstSp.df_md[tstSp.df_md.tr$tip.label,]

#### ID trait vertical plot (only species with  mortality) ------------------------------------------------------------------------- ####

pdf(paste0(rslts.dir, "figures/phylogeny/mortality_data_species_mean_values_rf.pdf"), height = h/1.5, width = w/1.5)
plotMortPhyloFun(phylo = tstSp.df_md.tr, traits = tstSp.df_md, brwPal = "BrBG", FontsizeFactor = 0.1)
dev.off()
print(paste0(rslts.dir, "figures/phylogeny/mortality_data_species_mean_values_rf.pdf"))

## HSM and mortality proportions

# Negative HSM
neg_HSM_all_species_traits_polygons <- allTraits %>% filter(HSM < 0)
length(neg_HSM_all_species_traits_polygons$HSM) # 8899 (3434) species with negative HSM
length(allTraits$HSM)

length(neg_HSM_all_species_traits_polygons$HSM)/ length(allTraits$HSM) # (8.09%)


# low HSM (<0.5)

low_HSM_all_species_traits_polygons <- allTraits %>% filter(HSM < 0.5)
length(low_HSM_all_species_traits_polygons$HSM) # 35718 (24101) species with HSM < 0.5
length(low_HSM_all_species_traits_polygons$HSM)/length(allTraits$HSM)

# low HSM (<1)

low_HSM_all_species_traits_polygons <- allTraits %>% filter(HSM < 1)
length(low_HSM_all_species_traits_polygons$HSM) # 44410 (42203) species with HSM < 1

length(low_HSM_all_species_traits_polygons$HSM)/length(allTraits$HSM) # 98.1 (92.9) species with HSM < 1


#### LOOK CLOSER TO SOME GROUPS -------------------------------------------------------------------------------------------------- ####

plotPhyloFamFun(phylo = tstSp.tr, traits = tstSp.df, brwPal = "BrBG", panelSpace = 4, FontsizeFactor = 0.1,
                            orderName = "Pinales")
