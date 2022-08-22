#### DATASET PREPARATION --------------------------------------------------------------------------------------------------------------- ####

## P Sanchez-Martinez

print("Dataset preparation ...")

source("code/data_preparation/0-data_preparation_init.R")


#### FUNCTIONS ------------------------------------------------------------------------------------ ####


## VPDmax warmest month calculation
vpdMaxFun <- function(df = env, tmaxVars = c(30:41)){
  df$VPDMax <- numeric(length(df[,1]))
  print("Calculating species VPDmax warmest month ...")
  for(i in 1:length(df$species)){
    if(all(is.na(df[i, tmaxVars]))){
      df$VPDMax[i] <- NA
    } else {
      tmax.month <- which.max(df[i , tmaxVars])                                                     # Warmest month
      if(tmax.month < 10){
        tmax.month <- paste0("0", as.character(tmax.month))
      }
      vapr.month <- paste0("vapr_", tmax.month)                                                  # vapr (kPa) of the warmest month           
      tmax.month <- paste0("tmax_", tmax.month)                                                  # Tmax value warmest month
      actVapr <- df[i, vapr.month]                                                               # vapr value warmest month
      satVp <- SVP(df[i, tmax.month], isK = F, formula = c("Clausius-Clapeyron", "Murray"))      # SatVP (hPa) (Tmax warmest month)
      df$VPDMax[i] <- (satVp/10) - actVapr                                                       # VPDmax = vapr(warmest month) - Tmax (warmest month)
    }
  }
  print("Removing monthly tmax and vapr variables...")
  df <- df[, -tmaxVars]
  vaprVars <- paste0("vapr_0", 1:9)
  vaprVars <- c(vaprVars, paste0("vapr_", 10:12))
  vaprVars <- which(names(df) %in% vaprVars)
  df <- df[, -vaprVars]
  return(df)
}# end of function

## Combination function
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

#### DATA ------------------------------------------------------------------------------------------------------------------------------ ####

# Hydra TRY

ht <- read.csv("data/traits/HydraTRY_Pablo_2020.csv", header = T)

### Merge environmental data

# Aridity index
ai <- read.csv("data/species_data/all_species_env_ai.csv", header = T)
summary(ai)
# Bioclim
bio <- read.csv("data/species_data/all_species_env_bio.csv", header = T)
summary(bio)
# Soilgrids
sg <- read.csv("data/species_data/all_species_env_sg.csv", header = T)
summary(sg)
# Worldclim
tmax <- read.csv("data/species_data/all_species_env_tmax.csv", header = T)
summary(tmax)

vapr <- read.csv("data/species_data/all_species_env_vapr.csv", header = T)
summary(vapr)

climAnnual <- read.csv("data/species_data/all_species_env_annual.csv", header = T)
summary(climAnnual)

env <- merge(ai, bio, by = "species")
env <- merge(env, sg, by = "species")
env <- merge(env, tmax, by = "species")
env <- merge(env, vapr, by = "species")
env <- merge(env, climAnnual, by = "species")

### Change names

envNames <- str_remove_all(names(env), "wc2.1_2.5m_")
envNames <- str_remove_all(envNames, "sg_2.5m_")
names(env) <- envNames

### VPDMax calculation --------------------------------------------------------------------------------------------------------- ####

names(env)
env.vpd <- vpdMaxFun(df = env,  tmaxVars = c(30:41))

### Rename final dataset ------------------------------------------------------------------------------------------------------ ####

names(ht)
env.ht <- dplyr::rename(env.vpd, 
                 # Worldclim
                 MAT = bio_1, Diurnal_TRange = bio_2, Isotermality = bio_3, TSeasonality = bio_4, TMax.Warmest = bio_5,
                 TMin.Coldest = bio_6, TRange = bio_7, TMean_Wettest = bio_8, TMean_Driest = bio_9, TMean_Warmest = bio_10, TMean_Coldest = bio_11,
                 AP = bio_12, Prec_Wettest = bio_13, Prec_Driest = bio_14, Prec_Seasonality = bio_15, Prec_WetQ = bio_16, Prec_DryQ = bio_17, Prec_WarmQ = bio_18,
                 Prec_ColdQ = bio_19, windMax = wind_max_annual, wind = wind_mean_annual, vapr = vapr_mean_annual, srad = srad_mean_annual, sradMin = srad_min_annual, 
                 tmax = tmax_mean_annual, tmin = tmin_mean_annual,
                 # Soilgrids
                 ABDRock = BDTICM_M_1km_ll, CEC_30cm = CECSOL_M_sl4_1km_ll, clay_30cm = CLYPPT_M_sl4_1km_ll, OCarbon_30cm = ORCDRC_M_sl4_1km_ll, ph_30cm = PHIHOX_M_sl4_1km_ll,
                 silt_30cm = SLTPPT_M_sl4_1km_ll, sand_30cm = SNDPPT_M_sl4_1km_ll, SWC_200cm = WWP_M_sl7_1km_ll, 
                 # aridity index
                 AI = ai_et0
                 )

### FINAL DATASET

env.comb <- env.ht[, which(names(env.ht) %in% names(ht))]
names(env.comb)

ht.comb <- ht[, which(names(ht) %in% names(env.comb))]
names(ht.comb)

ht.trait <- ht[, c(95, 88, 1:33)]
names(ht.trait)

# Traits and common environmental variables
ht.comb <- merge(ht.trait, ht.comb, by = "species", all = T)
names(ht.comb)
length(ht.comb$species)

# Environmental variables for the whole dataset
comb.df <- merge(env.comb, ht.comb, by = "species", all = T)
names(comb.df)
length(comb.df$species) # 47264

### XFT dataset (Hammond et al. 2021) ------------------------------------------------------------------------------------------------------ ####

if(isTRUE(hammond_xft_data | isTRUE(hammond_xft_data_P88))){
  
xft <- read.csv("data/traits/XFT_2022_hammond_Sanchez-Martinez.csv", header = T)

## Filter curves

xft <- xft %>% filter(Curve == "S", Developmental.stage %in% c("a", "A"))

summary(as.factor(xft$Curve))
summary(as.factor(xft$Developmental.stage))

xft <-  xft %>% dplyr::select(binomial_gnr, Group, Family, Growth.form, lat, long, Height.max..m., Height.actual..m., P50, P88, Ks..kg.m.1.MPa.1.s.1., KL..kg.m.1.MPa.1.s.1., psi.min.predawn..MPa., psi..min.midday..MPa., Huber.value..m2.m.2.,
                       Leaf.to.sapwood.area.ratio..m2.m.2., SLA..cm2.g.1., Gs..mol.m.2.s.1., TLP..MPa.) %>% 
  dplyr::rename(Species = binomial_gnr, group = Group, family = Family, latitude_mean = lat, longitude_mean = long, Hmax = Height.max..m., Hact = Height.actual..m., Ks = Ks..kg.m.1.MPa.1.s.1., Kl = KL..kg.m.1.MPa.1.s.1., MinWP_pd = psi.min.predawn..MPa., MinWP_md = psi..min.midday..MPa.,
                Hv = Huber.value..m2.m.2., SLA = SLA..cm2.g.1., gmax = Gs..mol.m.2.s.1., PItlp = TLP..MPa.)


xft$species <- str_replace(xft$Species, " ", "_")

# Aggregate data
names(xft)
catVar <- c("family", "group", "Growth.form", "species")
numVar <- c("latitude_mean", "longitude_mean", "Hmax", "Hact", "P50", "P88", "Ks", "Kl", "MinWP_pd", "MinWP_md", "Hv", "Leaf.to.sapwood.area.ratio..m2.m.2.", "SLA", "gmax", "PItlp")

for(var in numVar){
  xft[, var] <- as.numeric(xft[, var])
}

for(var in catVar){
  xft[, var] <- as.character(xft[, var])
}

xft_aggr_mean <- aggregate(xft[, numVar], by = list(xft$Species), FUN = mean.fun)
xft_aggr_mode <- aggregate(xft[, catVar], by = list(xft$Species), FUN = mode.fun)

xft_aggr <- merge(xft_aggr_mean, xft_aggr_mode, by = "Group.1")

names(xft_aggr)[names(xft_aggr) == "Group.1"] <- "Species"

xft[xft$Species == "Abies alba",]
head(xft_aggr)

levels(as.factor(ht$group))
levels(as.factor(xft_aggr$group))
xft_aggr$group <-  dplyr::recode(xft_aggr$group, angiosperm = "Angiosperms", gymnosperm = "Gymnosperms")
xft_aggr <- xft_aggr %>% filter(group %in% c("Angiosperms", "Gymnosperms"))
length(xft_aggr$Species) # 1553 (999 after filtering)

# Combine

length(comb.df$species)

comb.df <- merge(comb.df, xft_aggr, by = "species", all = T)
length(comb.df$species)
}

### Final dataset combination ------------------------------------------------------------------------------------------------------ ####

for(var in names(comb.df)){
  if(!is.na(str_extract(var, "\\.x"))){
    var <- str_remove(var, "\\.x")
    comb.df <- comb.fun(comb.df, var)
  }
}

# Change P50 for P88 in angiosperms
if(isTRUE(hammond_xft_data_P88)){
  comb.df[which(comb.df$group == "Angiosperms"), "P50"] <- comb.df[which(comb.df$group == "Angiosperms"), "P88"]
  summary(comb.df$P50)
}

### Save dataset ----------------------------------------------------------------------------------------- ####

comb.df$Species <- comb.df$species
comb.df$species <- NULL
names(comb.df)

comb.df <- comb.df[, c("Species","genus", "family", "order", "group", "Growth.form", "life.form", "leaf.form", "biome", "Hmax", "Hact", 
             "Ks", "Kl", "P50","P88",  "MinWP_pd",  "MinWP_md",  "HSM", "AlAs", "WD", "Vdia", "VD", "VLmax", "Hv", "PItlp", "WUE", "SLA", 
             "Aarea", "Amass", "LL", "pLL", "N", "Narea", "Rd", "gmax", "MAT", "TMean_Warmest", "TMean_Coldest",
             "AP", "Prec_Wettest", "Prec_Driest", "Prec_Seasonality", "Prec_WetQ", "Prec_DryQ", "Prec_WarmQ", "Prec_ColdQ", 
             "Diurnal_TRange", "Isotermality", "TSeasonality", "TMax.Warmest", "TMin.Coldest", "TRange", "TMean_Wettest", 
             "TMean_Driest", "srad", "sradMin", "tmax", "tmin", "vapr", "windMax", "wind", "ABDRock", "SWC_200cm", "AI", 
             "VPDMax", "CEC_30cm", "clay_30cm", "OCarbon_30cm", "ph_30cm", "silt_30cm", "sand_30cm")]


names(comb.df)
summary(comb.df)
length(comb.df$Species) # 47264 (47598 with hammond data)

a <- comb.df %>% filter(!is.na(MinWP_md))
length(a$Species) # 536 (831 with xft data, 720 after filtering)
a <- comb.df %>% filter(!is.na(P50))
length(a$Species) # 863 (1716 with xft, 1470 after filtering) (1087 with P88 xft)

a <- comb.df %>% filter(!is.na(P88))
length(a$Species) # 863 (1716 with xft, 735 after filtering) (1087 with P88 xft)

write.csv(comb.df, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_1.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_1.csv"))
