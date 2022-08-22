#### FINAL DATABASE ARRANGEMENTS ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ######


print("Final database arrengements and ready to go!")

source("code/data_preparation/0-data_preparation_init.R")

#### DATA ----------------------------------------------------------------------------------------------------------------------- ####

# Traits and predictors
ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_6_mean_MinWP_soil.csv"), header = T)

summary(as.factor(ht$genus))
### Join Agiosperms/gymnosperms and leaf form ###
ht$group_leaf <- paste0(ht$group, "_", ht$leaf.form)

ht$group_leaf <- as.character(ht$group_leaf)
ht[which(ht$group_leaf == "Gymnosperms_D"), "group_leaf"] <- "Gymnosperms"
ht[which(ht$group_leaf == "Gymnosperms_E"), "group_leaf"] <- "Gymnosperms"
ht[which(ht$group_leaf == "Gymnosperms_W"), "group_leaf"] <- "Gymnosperms"
ht[which(ht$group_leaf == "Gymnosperms_NA"), "group_leaf"] <- "Gymnosperms"

for(i in 1:length(ht$group_leaf)){
  if(!is.na(str_extract(ht$group_leaf[i], "NA"))){
    ht$group_leaf[i] <- NA
  }
}
ht$group_leaf <- as.factor(ht$group_leaf)


# Variables order
ht$sp <- str_replace_all(ht$Species, " ", "_")
length(names(ht))
ht <- ht[, c("Species", "sp", "genus", "family", "order", "group", "Growth.form", "life.form", "leaf.form", "group_leaf", "biome", 
             "Hmax", "Hact", "Ks", "Kl", "P50", "P88", "MinWP_pd", "MinWP_md","HSM", "AlAs", "WD", "Vdia", "VD", "VLmax", "Hv", "PItlp", "WUE", "SLA", "Aarea", "Amass", "LL", "pLL", "N", "Narea", "Rd", "gmax",
             "MAT", "TMean_Warmest", "TMean_Coldest", "AP", "Prec_Wettest", "Prec_Driest", "Prec_Seasonality", "Prec_WetQ", "Prec_DryQ", "Prec_WarmQ", "Prec_ColdQ", "Diurnal_TRange", "Isotermality", "TSeasonality", "TMax.Warmest", "TMin.Coldest", "TRange", "TMean_Wettest", "TMean_Driest", "srad", "sradMin", "tmax", "tmin", "vapr", "windMax", "wind", "ABDRock", "SWC_200cm", "AI", "VPDMax",
             "MinWP_soil_SX","CEC_30cm", "clay_30cm", "OCarbon_30cm", "ph_30cm", "silt_30cm", "sand_30cm", 
             "log_negP50", "log_negMinWP_md", "log_Ks", "log_Hv", "log_Kl", "sqrt_AP", "sqrt_Prec_Wettest", "sqrt_Prec_Driest", "sqrt_Prec_Seasonality", "sqrt_Prec_WetQ", "sqrt_Prec_DryQ", "sqrt_Prec_WarmQ", "log_Prec_ColdQ", "log_TSeasonality", "log_TRange", "log_windMax", "log_wind", "log_ABDRock", "sqrt_AI", "log_VPDMax", "log_CEC_30cm", "log_OCarbon_30cm", 
             "PC1", "PC2", "PC3", "PC4",  "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", paste0("Phylo_axis_", 1:100))]

length(names(ht)) # 206
length(ht$Species) # (47197)


### Filtering extremely high values ###

length(ht$P50[which(ht$P50 > -0.5)]) # 39 P50 values extremely high
ht$P50[which(ht$P50 > -0.5)] <- NA

length(ht$MinWP_md[which(ht$MinWP_md > -0.5)]) # 15 Pmin values extremely high
ht$MinWP_md[which(ht$MinWP_md > -0.5)] <- NA

length(names(ht))
length(ht$Species) # (47197)

write.csv(ht, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2021.csv"))


# Number of observations for hydraulic traits Sanchez-Martinez et al. 2021
length(ht$Species)
a <- ht %>% filter(!is.na(MinWP_md))
length(a$Species) # 536 (819 with hammond data, 685 after filtering high values)
a <- ht %>% filter(!is.na(P50))
length(a$Species) # 863 (1678 with hammond data, 1376 after filtering high values) (1073 using P88 angiosperms)

a <- ht %>% filter(!is.na(P50) | !is.na(MinWP_md))
length(a$Species) # (1832 with hammond data, 1531 after filtering high values) (1389 using P88 angiosperms)
