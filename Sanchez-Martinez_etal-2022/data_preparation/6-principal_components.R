#### PRINCIPAL COMPONENTS CALCULATION ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

## P Sanchez-Martinez

print("Running PCA computation...")

source("code/data_preparation/0-data_preparation_init.R")

### FUNCTIONS -----------------------------------------------------------------------------------------------------------#### 

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


### VARIABLES AND DATA --------------------------------------------------------------------------------------------------####

ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_3_varTransf.csv"), header = T)

names(ht)

hydVars <- c("log_negP50", "log_negMinWP_md", "HSM")
envVars <- c("MAT", "TMean_Warmest", "TMean_Coldest",
             "sqrt_AP", "sqrt_Prec_Wettest", "sqrt_Prec_Driest", "sqrt_Prec_Seasonality", "sqrt_Prec_WetQ", "sqrt_Prec_DryQ", "sqrt_Prec_WarmQ", "Prec_ColdQ", 
             "Diurnal_TRange", "Isotermality", "log_TSeasonality", "TMax.Warmest", "TMin.Coldest", "log_TRange", "TMean_Wettest", 
             "TMean_Driest", "srad", "sradMin", "tmax", "tmin", "vapr", "log_windMax", "log_wind", "log_ABDRock", "SWC_200cm", "sqrt_AI", 
             "log_VPDMax", "log_CEC_30cm", "clay_30cm", "log_OCarbon_30cm", "ph_30cm", "silt_30cm", "sand_30cm")

envVars <- c("MAT",
             "sqrt_AP", "sqrt_Prec_Wettest", "sqrt_Prec_Driest", "sqrt_Prec_Seasonality", "sqrt_Prec_WarmQ", "Prec_ColdQ", "Isotermality",
             "log_TSeasonality", "TMin.Coldest", "TMean_Wettest", 
             "TMean_Driest", "srad", "vapr", "log_wind", "log_ABDRock", "SWC_200cm", 
             "log_VPDMax", "log_CEC_30cm", "clay_30cm", "log_OCarbon_30cm", "ph_30cm", "sand_30cm")


## Complete data
ht_pca <- completeFun(ht, envVars)
length(ht_pca$Species)
summary(ht_pca)
# Auxiliar database used later to merge PCs with HydraTRY
ht_pc <- ht_pca

# Environmental variables
ht_pca <- ht_pc[, envVars]

# Hydraulic traits
ht_pca_hyd <- ht_pc[, hydVars]

# For group plot
ht_pca_gr <- ht_pc$group
# For leaf form plot
ht_pca_de <- ht_pc$leaf.form

### PCA results ----------------------------------------------------------------------------------------------------####

pca <- prcomp(ht_pca, scale = T)
summary(pca)


#Eighenvalues
get_eig(pca)
fviz_eig(pca, addlabels = T, hjust=0.3)
dev.off()

# Results for Variables
res.var <- get_pca_var(pca)
round(res.var$cor, 3)            # Coordinates
round(res.var$contrib, 3)        # Contributions to the PCs
round(res.var$cos2, 3)           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation


### HT - PCs PEARSON CORRELATIONS ------------------------------------------------------------------------------------------------ ####

## Predict coordinates and compute cos2
hyd.coord <- cor(ht_pca_hyd, (pca$x)*-1, use = "pairwise.complete.obs")
hyd.cos2 <- hyd.coord^2

## Correlations
hyd.coord[, 1:3]


### MERGE WITH HydraTRY DATASET ------------------------------------------------------------------------------------------#### 

## Save 3 PC

pca_merge <- pca$x * -1

length(pca_merge[, 1])

length(ht_pc$Species)
ht_pc$PC1 <- pca_merge[, 1]
ht_pc$PC2 <- pca_merge[, 2]
ht_pc$PC3 <- pca_merge[, 3]
ht_pc$PC4 <- pca_merge[, 4]
ht_pc$PC5 <- pca_merge[, 5]
ht_pc$PC6 <- pca_merge[, 6]
ht_pc$PC7 <- pca_merge[, 7]
ht_pc$PC8 <- pca_merge[, 8]
ht_pc$PC9 <- pca_merge[, 9]
ht_pc$PC10 <- pca_merge[, 10]

summary(ht_pc)
length(ht_pc$Species)

## Join with Ht dataset
names(ht_pc)
ht_pca <- ht_pc[, c("Species", paste0("PC", 1:10))]
ht_pca <- join(ht, ht_pca, by = "Species")
summary(ht_pca)
length(ht_pca$Species) # (47197)

write.csv(ht_pca, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_4_env_PCs.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_4_env_PCs.csv"))
