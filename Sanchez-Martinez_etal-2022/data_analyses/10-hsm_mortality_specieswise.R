##### PHYLOGENETIC PROJECTION -------------------------------------------------- #####

remove(list = ls())
gc()
# install.packages("ggplot2", dependencies = T)

#### PACKAGES ------------------------------------------------------------------ ####

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

h <- 2.61
w <- 7.87

tables.dir <- paste0(rslts.dir, "tables/speies_lvl/")
rslts.dir <- paste0(rslts.dir, "figures/HSM_mortality/species_lvl/")

#### FUNCTIONS ----------------------------------------------------------------- ####

logitModSp <- function(frml = "p ~ HSM", dta = hsm_md_sp, nIter = 5, plotCurves = T, fmly = "binomial"){
  
  ## Modelling data preparation
  rslts.df <- data.frame()
  r2.vec <- numeric()
  aic.vec <- numeric()
  auc.vec <- numeric()
  fitted.df <- data.frame()
  
  # Variables (carefull if generalizing that, very code-specific)
  
  vars <- str_replace_all(frml, "\\+", " ")
  vars <- str_replace_all(vars, "\\*", " ")
  vars <- str_replace_all(vars, "\\~", " ")
  
  vars <- unlist(str_split(vars, pattern = " "))
  vars <- vars[-which(vars == "")]
  
  dta  <- dta %>% dplyr::select(vars)
  
  if("mortality_events" %in% vars){
    dta$p <- dta$mortality_events
  }
  
  # Modelling
  for(i in 1:nIter){
    print(i)
    # Random points
    dta <- na.omit(dta)
    mort.dta <- dta %>% dplyr::filter(p > 0)
    nbgr <- length(mort.dta[ ,1])
    bg.dta <- dta %>% dplyr::filter(p == 0)
    bg.dta <- slice_sample(bg.dta, n = nbgr) # Same amount of mortality data, but we also used 10000 random points with similar results
    mod.dta <- rbind(mort.dta, bg.dta)
    
    # Model
    mod <- glm(as.formula(frml), data = mod.dta, family = fmly)
    smry <- summary(mod)
    
    # Coefficients
    coef.df <- as.data.frame(smry$coefficients[, c("Estimate", "Std. Error", "Pr(>|z|)")])
    coef.df$Id <- row.names(coef.df)
    rslts.df <- rbind(rslts.df, coef.df)
    
    if(isTRUE(plotCurves)){
      # Predict (for curve response)
      fitted <- data.frame("fit" = predict(mod, mod.dta, type = "response"), "iter" = i)
      fitted <- cbind(fitted, mod.dta[, vars])
      fitted.df <- rbind(fitted.df, fitted) 
    }
    
    # R2
    r2 <- lrm(as.formula(frml), data = mod.dta)
    r2 <- r2$stats["R2"]
    r2.vec[i] <- r2
    # AIc
    aic <- smry$aic
    aic.vec[i] <- aic
    
    # AUC
    if(fmly == "binomial"){
    # test and train
    sample = sample.split(mod.dta$p, SplitRatio = 0.80)
    train = subset(mod.dta, sample == TRUE)
    test  = subset(mod.dta, sample == FALSE)
    
    aucMod <- glm(as.formula(frml), data = train, family = fmly)
    eval <- dismo::evaluate(p=test[which(test$p > 0), ], a=test[which(test$p == 0), ], model = aucMod)
    auc <- eval@auc
    auc.vec[i] <- auc

  } else {
    auc.vec <- NA
  }
  }
  
  # Summarize coefficients
  coefficients.mean <- rslts.df %>% group_by(Id) %>% summarize_all(., mean)
  coefficients.sd <- rslts.df %>% group_by(Id) %>% summarize_all(., sd)
  coefficients.rslts <- cbind("model" = frml, coefficients.mean, coefficients.sd[, -1])
  names(coefficients.rslts) <- c("model", "parameter", "estimate_mean", "st_error_mean", "p_value_mean", "estimate_sd", "st_error_sd", "p_value_sd")
  r2.rslts <- data.frame("model" = frml,"R2_mean" = mean(r2.vec), "R2_sd" = sd(r2.vec))
  auc.rslts <- data.frame("model" = frml,"auc_mean" = mean(auc.vec), "auc_sd" = sd(auc.vec))
  
  
  # Plots
  
  mod.dta[, vars[1]] <- as.factor(mod.dta[, vars[1]])
  
  if(length(vars) == 2){
    p <- ggplot(mod.dta, aes(x = mod.dta[, vars[2]], y = mod.dta[, vars[1]], colour = mod.dta[, vars[1]])) + geom_boxplot() + 
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.text = element_blank(), legend.title = element_blank())  +
      scale_color_manual(values=c("1" ="#e34a33", "0" = "#2ca25f"))
    boxPlot_leg <- get_legend(p)
    p <-  p + theme(legend.position = "none")
    dev.off()
    # Save
    pdf(paste0(rslts.dir, "boxplots/", vars[2], "_boxplot.pdf"), height = w/3, width = w/3)
    print(p)
    dev.off()
    print(paste0(rslts.dir, "boxplots/", vars[2], "_boxplot.pdf"))
    
    if(isTRUE(plotCurves)){
      curve <- ggplot(fitted.df, aes(y = fit, x =  fitted.df[, vars[2]]))
      curve <- curve + geom_smooth(method = "gam", color = "#e41a1c")  +
        scale_color_brewer(palette = "Spectral") + theme(legend.position = "none") + theme(axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.title = element_blank())
      
      # Save curves
      pdf(paste0(rslts.dir, "curves/", vars[2], "_curves.pdf"), height = w/3, width = w/3)
      # png(paste0("results/figures/HSM_mortality/curves/", vars[2], "_curves.png"), height = w/2, width = w/2, units = "in", res = 600)
      print(curve)
      dev.off()
      print(paste0(rslts.dir, "curves/", vars[2], "_curves.pdf"))
    }
  }
  if(length(vars) == 3){
    
    p <- ggplot(mod.dta, aes(x = mod.dta[, vars[3]], y =  mod.dta[, vars[2]], colour = mod.dta[, vars[1]])) +
      geom_boxplot() + 
      theme(axis.text.x = element_blank(), axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.text = element_blank(), legend.title = element_blank()) +
      scale_color_manual(values=c("1" ="#e34a33", "0" = "#2ca25f"))
    boxPlot_leg <- get_legend(p)
    p <-  p + theme(legend.position = "none")
    
    # Save boxplot
    pdf(paste0(rslts.dir, "boxplots/", vars[2], "_", vars[3], "_boxplot.pdf"), height = w/3, width = w/3)
    print(p)
    dev.off()
    print(paste0(rslts.dir, "boxplots/", vars[2], "_", vars[3], "_boxplot.pdf"))
    
    if(isTRUE(plotCurves)){
      curve <- ggplot(fitted.df, aes(y = fit, x =  fitted.df[, vars[2]], col =  fitted.df[, vars[3]]))
      # curve <- curve +  geom_point(alpha = 0.1, size = 0.5)
      curve <- curve + geom_smooth(aes(y = fit, x =  fitted.df[, vars[2]], col =  fitted.df[, vars[3]]), method = "gam")  + theme(legend.position = "none")
      curve <- curve +  theme(axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.title = element_blank()) +
        scale_color_brewer(palette = "Spectral")
      
      # Save curves
      pdf(paste0(rslts.dir, "curves/", vars[2], "_", vars[3], "_curves.pdf"), height = w/3, width = w/3)
      # png(paste0("results/figures/HSM_mortality/curves/", vars[2], "_", vars[3], "_curves.png"), height = w/2, width = w/2, units = "in", res = 600)
      print(curve)
      dev.off()
      print(paste0(rslts.dir, "curves/", vars[2], "_", vars[3], "_curves.pdf"))
    }
  }
  # Results
  rslts <- list()
  rslts$coefficients <- coefficients.rslts
  rslts$R2 <- r2.rslts
  rslts$aic <- mean(aic.vec)
  rslts$auc <- auc.rslts
  
  
  # rslts$model_plot <- mod
  # rslts$boxPlot <- p
  # if(isTRUE(plotCrvs)){
  #   rslts$curve <- curve 
  # }
  if(length(vars) > 2){
    rslts$boxPlot_leg <- boxPlot_leg
  }
  # rslts$mod.data <- mod.dta
  return(rslts)
}

#### DATA ---------------------------------------------------------------------- ####

# Mortality data species-wise

md_species <- read.csv("data/mortality/Hammond_etal_database_withspecies_transformed_species_wise.csv", header = T)
md_species$species <- as.factor(md_species$species)
length(unique(md_species$species))
md_species_count <- as.data.frame(table(md_species$species)) 
md_species_count
names(md_species_count) <- c("Species", "mortality_events")
summary(md_species_count)

# Imputed traits data

allTraits_imp <- read.csv(paste0(output.dir, "RF_imputations/RFimputation_100_bi_climate_5_phylo_5_group.csv"), header = T) # 45345 species
allTraits_imp$sp <- str_replace_all(allTraits_imp$Species, " ", "_")

allTraits <- join(allTraits_imp, md_species_count, by = "Species", type = "left")
allTraits[which(is.na(allTraits$mortality_events)), "mortality_events"] <- 0
allTraits[which(allTraits$mortality_events > 0), "mortality_binomial"] <- 1
allTraits[which(is.na(allTraits$mortality_binomial)), "mortality_binomial"] <- 0

allTraits <- dplyr::rename(allTraits, p = mortality_binomial)

# Observed traits data

ht_md <- join(ht, md_species_count, by = "Species", type = "left")
ht_md[which(is.na(ht$mortality_events)), "mortality_events"] <- 0
ht_md[which(ht_md$mortality_events > 0), "mortality_binomial"] <- 1
ht_md[which(is.na(ht_md$mortality_binomial)), "mortality_binomial"] <- 0

ht_md <- dplyr::rename(ht_md, p = mortality_binomial)

#### MODELS  ------------------------------------------------------------------- ####

# Using imputed data

a <- allTraits %>% filter(!is.na(HSM),!is.na(p))
count(a$p)

hsm_spMod <- logitModSp(frml = "p ~ HSM", dta = allTraits, nIter = 100, plotCurves = T, fmly = "binomial")
hsm_spMod$coefficients
hsm_spMod$R2
hsm_spMod$auc

hsm_spMod_group <- logitModSp(frml = "p ~ HSM * group", dta = allTraits, nIter = 100, plotCurves = T)
hsm_spMod_group$coefficients
hsm_spMod_group$R2
hsm_spMod_group$auc


hsm_splm <- logitModSp(frml = "mortality_events ~ HSM", dta = allTraits, nIter = 100, plotCurves = T,  fmly = "poisson")
hsm_splm$coefficients
hsm_splm$R2

hsm_splm_group <- logitModSp(frml = "mortality_events ~ HSM * group", dta = allTraits, nIter = 100, plotCurves = T,  fmly = "poisson")
hsm_splm_group$coefficients
hsm_splm_group$R2
hsm_splm_group$auc

# Using observed data

a <- ht_md %>% filter(!is.na(HSM),!is.na(p))
count(a$p)

hsm_spMod_obs <- logitModSp(frml = "p ~ HSM", dta = ht_md, nIter = 100, plotCurves = T, fmly = "binomial")
hsm_spMod_obs$coefficients
hsm_spMod_obs$R2
hsm_spMod_obs$auc

hsm_spMod_group_obs <- logitModSp(frml = "p ~ HSM * group", dta = ht_md, nIter = 100, plotCurves = T)
hsm_spMod_group_obs$coefficients
hsm_spMod_group_obs$R2
hsm_spMod_group_obs$auc


hsm_splm_obs <- logitModSp(frml = "mortality_events ~ HSM", dta = ht_md, nIter = 100, plotCurves = T,  fmly = "poisson")
hsm_splm_obs$coefficients
hsm_splm_obs$R2

hsm_splm_group_obs <- logitModSp(frml = "mortality_events ~ HSM * group", dta = ht_md, nIter = 100, plotCurves = T,  fmly = "poisson")
hsm_splm_group_obs$coefficients
hsm_splm_group_obs$R2
hsm_splm_group_obs$auc


#### SPECIES WITH MORTALITY VALUES RESPECT MEAN HSM PER PIXEL ------------------ ####

# Aggregate data each 10km

md_species_aggr <- md_species %>% dplyr::select(species, x, y)

md_species_aggr$x <- round(md_species_aggr$x/10) * 10
md_species_aggr$y <- round(md_species_aggr$y/10) * 10
md_species_aggr$x <- round(md_species_aggr$x)
md_species_aggr$y <- round(md_species_aggr$y)

md_species_aggr <- unique(md_species_aggr)
length(md_species_aggr$x)

HSM_mean <- raster(paste0(output.dir, "HSM/HSM_mean.tif"))
HSM_mean

md_species_aggr <- dplyr::rename(md_species_aggr, Species = species)

allTraits_md <- join(md_species_aggr, allTraits_imp, by = "Species", type = "left")
head(allTraits_md)
length(allTraits_md$Species)

# Extract mean values per mortality point

# plot(HSM_mean)
# points(md_species[, c("x", "y")])

allTraits_md$Mean_HSM_speciesPool <- raster::extract(HSM_mean, allTraits_md[, c("x", "y")])

allTraits_md <- na.omit(allTraits_md)
head(allTraits_md)
length(md_species_aggr$Species)

### T test comparison ####

# Compare HSM of species with mortality to mean species pool HSM values
allTraits_md_hsm_comp <- allTraits_md %>% dplyr::select(Species, group, HSM, Mean_HSM_speciesPool, x, y) %>% dplyr::rename(species_HSM = HSM)
allTraits_md_hsm_comp$group  <- as.factor(allTraits_md_hsm_comp$group)
summary(allTraits_md_hsm_comp)

allTraits_md_hsm_comp$isInferior <- ifelse(allTraits_md_hsm_comp$species_HSM < allTraits_md_hsm_comp$Mean_HSM_speciesPool, 1, 0)
sum(allTraits_md_hsm_comp$isInferior)
length(allTraits_md_hsm_comp$Species)

lowerPoints <- allTraits_md_hsm_comp %>% filter(isInferior == 1)

plot(HSM_mean)
points(allTraits_md_hsm_comp[, c("x", "y")], col = "green")
points(lowerPoints[, c("x", "y")], col = "red")

allTraits_md_hsm_comp$isNegative <- ifelse(allTraits_md_hsm_comp$species_HSM < 0, 1, 0)
sum(allTraits_md_hsm_comp$isNegative)

ggplot(allTraits_md_hsm_comp, aes(x = Mean_HSM_speciesPool, y = species_HSM, color = group)) + geom_point()  + geom_abline()
ggplot(allTraits_md_hsm_comp) + geom_density(aes(x = Mean_HSM_speciesPool)) + geom_density(aes(x = species_HSM, color = "red"))

t.test <- t.test(x = allTraits_md_hsm_comp$species_HSM, y = allTraits_md_hsm_comp$Mean_HSM_speciesPool)      
t.test

## Angiosperms

allTraits_md_hsm_comp_ang <- allTraits_md_hsm_comp %>% filter(group == "Angiosperms")

t.test_ang <- t.test(x = allTraits_md_hsm_comp_ang$species_HSM, y = allTraits_md_hsm_comp_ang$Mean_HSM_speciesPool)      
t.test_ang

allTraits_md_hsm_comp_gym <- allTraits_md_hsm_comp %>% filter(group == "Gymnosperms")

t.test_gym <- t.test(x = allTraits_md_hsm_comp_gym$species_HSM, y = allTraits_md_hsm_comp_gym$Mean_HSM_speciesPool)      
t.test_gym

## Gymnosperms

### Linear model ####

species_HSM.df <- data.frame("variable" = "species_HSM", "value" = allTraits_md_hsm_comp$species_HSM, "group" = allTraits_md_hsm_comp$group)
pools_HSM.df <- data.frame("variable" = "pools_HSM", "value" = allTraits_md_hsm_comp$Mean_HSM_speciesPool, "group" = allTraits_md_hsm_comp$group)

allTraits_md_hsm_comp_lm <- rbind(species_HSM.df, pools_HSM.df)

# By species

lmod <- lm(value ~ variable, data = allTraits_md_hsm_comp_lm)
summary(lmod)

lmod_group <- lm(value ~ variable * group, data = allTraits_md_hsm_comp_lm)
summary(lmod_group)
