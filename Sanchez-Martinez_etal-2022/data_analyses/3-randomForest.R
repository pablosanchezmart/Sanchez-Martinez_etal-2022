### RANDOM FORESTS IMPUTATION------------------------------------------------------------------------------------- ####


# P Sanchez-Martinez

remove(list = ls())
source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

forceRun <- T

#### CROSS VALIDATION --------------------------------------------------------------------------------------------------------------------- ####

# Parallelization options
showConnections()
cl <- makeCluster(nClusters)
registerDoParallel(cl)
for(propNA in c(0.1, 0.2, 0.3, 0.5, 0.7, 0.8)){
  print(propNA)
  if(file.exists(paste0(output.dir, "RF_imputations/models_performance_",  as.character(gsub("\\.", "", propNA)), "NA.csv")) && forceRun == F){
    r2.rslts <- read.csv(paste0(output.dir, "RF_imputations/models_performance_",  as.character(gsub("\\.", "", propNA)), "NA.csv"), header = T) # Previous results 
  } else{
    r2.rslts <- data.frame()
  }
mod.lst <- list()
for(i in 1:length(RFmodels$Model)){
  mdl.name <- as.character(RFmodels[i, "Model"])
  response <- unlist(strsplit(RFmodels[i, "Response"], split=", "))
  predictors <- unlist(strsplit(RFmodels[i, "Predictors"], split=", "))
  mod.rslts <- RFCrosValFun(name = mdl.name, data = rf.df, resp = response, pred = predictors, NasP = propNA, nIter = CVnIterations)
  mod.lst[[mdl.name]] <- mod.rslts
  r2.rslts <- rbind(r2.rslts, mod.rslts$sumR2)
}

write.csv(r2.rslts, paste0(output.dir, "RF_imputations/models_performance_",  as.character(gsub("\\.", "", propNA)), "NA.csv"), row.names = F)
print(paste0(output.dir, "RF_imputations/models_performance_",  as.character(gsub("\\.", "", propNA))))
}
stopCluster(cl)

### Unify model performance results ###

r2.files <- list.files(paste0(output.dir, "RF_imputations/"), pattern = "models_performance_")
r2.all <- data.frame("Model" = read.csv(paste0(output.dir, "RF_imputations/", r2.files[1]))[, "Model"])
r2.all$Variable <- read.csv(paste0(output.dir, "RF_imputations/", r2.files[1]))[, "Variable"]

for(i in 1:length(r2.files)){
  r2.df <- read.csv(paste0(output.dir, "RF_imputations/", r2.files[i]))
  r2.df <- r2.df[, c("Model", "Variable", "R2", "R2_sd")]
  NAprop <- parse_number(r2.files[i])
  names(r2.df)[c(3, 4)] <- paste0(names(r2.df)[c(3, 4)], "_0",NAprop)
  r2.all <- cbind(r2.all, round(r2.df[, c(3, 4)], 3))
}
head(r2.all)
write.csv(r2.all, paste0(output.dir, "RF_imputations/models_performance.csv"), row.names = F)
print(paste0(output.dir, "RF_imputations/models_performance.csv"))

#### PREDICTING USING BEST MISSING VALUES MULTIVARIATE RANDOM FOREST ------------------------------------------------------------------------------------------- ####

# Parallelization options
showConnections()
cl <- makeCluster(nClusters)
registerDoParallel(cl)

# set.seed(123)
imp.lst <- list()
OBBE.lst <- list()
model.name <- RFmodels$Model[24]
for(i in 1:impnIterations){
  print(i)
  response <- unlist(strsplit(RFmodels[RFmodels$Model == model.name, "Response"], split=", "))
  predictors <- unlist(strsplit(RFmodels[RFmodels$Model == model.name, "Predictors"], split=", "))
  dta <- rf.df[complete.cases(rf.df[, c(predictors)]), ]
  Species <- dta$Species
  dta <- dta %>% dplyr::select(c(response, predictors))
  rf <-  missForest(dta, verbose = F, variablewise = T, parallelize = "forests", maxiter = 50)
  names(rf$OOBerror) <- names(dta)
  OBBE.lst[[i]] <- list(rf$OOBerror)
  rfimp <- cbind(Species, rf$ximp)
  imp.lst[[i]] <- rfimp
  
}
stopCluster(cl)
save(imp.lst, file = paste0(output.dir, "RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))
print(paste0(output.dir, "RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))

#### RESULTS ----------------------------------------------------------------------------------------------------------------------------- ####

# Errors summarization
meanOBBE <- OBBE.lst %>%
  bind_rows() %>%
  summarise_all(., mean)

sdOBBE <- OBBE.lst %>%
  bind_rows() %>%
  summarise_all(., sd)

names(sdOBBE) <- paste0(names(sdOBBE), "_sd")

OBBE.rslts <- cbind(meanOBBE, sdOBBE)

write.csv(OBBE.rslts, paste0(output.dir, "RF_imputations/OBBE_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0(output.dir, "RF_imputations/OBBE_", impnIterations, "_" , model.name ,".csv"))

# HSM calculation and imputation summarisation

load(paste0(output.dir, "RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))

# Results in dataframe
imp.df <- imp.lst %>%
  bind_rows() %>%
  group_by(., Species) %>%
  dplyr::select(response)

# HSM calculation
imp.df$HSM <- imp.df$MinWP_md - imp.df$P50

# Mean imputation values
meanImp <- imp.df %>%
  summarise_all(., mean)

# Sd of imputation values
sdImp <- imp.df %>%
  summarise_all(., sd)

sdImp <- dplyr::rename(sdImp, P50_sd = P50, MinWP_md_sd = MinWP_md, HSM_sd = HSM) %>%
  dplyr::select(MinWP_md_sd, P50_sd, HSM_sd)

# Cv of imputation values
cvImp <- imp.df %>%
  summarise_all(., cv) %>% dplyr::select(MinWP_md, P50, HSM)
cvImp <- abs(cvImp)

cvImp <- dplyr::rename(cvImp, P50_cv = P50, MinWP_md_cv = MinWP_md, HSM_cv = HSM) %>%
  dplyr::select(MinWP_md_cv, P50_cv, HSM_cv)

imp.rslts <- cbind(meanImp, sdImp, cvImp, rf.df[complete.cases(rf.df[, c(predictors)]), c("genus", "family", "order", predictors)])
if(!"MinWP_soil_SX" %in% predictors){
  imp.rslts <- cbind(imp.rslts, "MinWP_soil_SX" = rf.df[complete.cases(rf.df[, c(predictors)]), "MinWP_soil_SX"])
}
head(imp.rslts)

write.csv(imp.rslts, paste0(output.dir, "RF_imputations/RFimputation_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0(output.dir, "RF_imputations/RFimputation_", impnIterations, "_" , model.name ,".csv"))


## HSM R2

hsmR2 <- rf.df %>% filter(!is.na(P50) & !is.na(MinWP_md))

impHsm <- read.csv(paste0(output.dir, "RF_imputations/RFimputation_", impnIterations, "_" , RFmodels$Model[24] ,".csv"), header = T) %>% 
  dplyr::select(Species, HSM) %>%
  dplyr::rename(impHSM = HSM)

hsmR2 <- merge(hsmR2, impHsm, by = "Species", all.x = T)

hsmR2.dif <- hsmR2 %>% filter(HSM != impHSM)

mod <- lm(HSM ~ impHSM, data = hsmR2.dif)
summary(mod)

#### PREDICTING USING BEST MISSING VALUES MULTIVARIATE RANDOM FOREST ONLY FOR OBSERVED DATA (OBSERVATIONS WITH SOME TRAIT RECORD) ------------------------------------------------------------------------------------------- ####

# Parallelization options
showConnections()
cl <- makeCluster(nClusters)
registerDoParallel(cl)

# set.seed(123)
imp.lst <- list()
OBBE.lst <- list()
model.name <- RFmodels$Model[24]
for(i in 1:impnIterations){
  print(i)
  response <- unlist(strsplit(RFmodels[RFmodels$Model == model.name, "Response"], split=", "))
  predictors <- unlist(strsplit(RFmodels[RFmodels$Model == model.name, "Predictors"], split=", "))
  dta <- obs.rf.df[complete.cases(obs.rf.df[, c(predictors)]), ]
  Species <- dta$Species
  dta <- dta %>% dplyr::select(c(response, predictors))
  rf <-  missForest(dta, verbose = F, variablewise = T, parallelize = "forests", maxiter = 50)
  names(rf$OOBerror) <- names(dta)
  OBBE.lst[[i]] <- list(rf$OOBerror)
  rfimp <- cbind(Species, rf$ximp)
  imp.lst[[i]] <- rfimp
  
}
stopCluster(cl)
save(imp.lst, file = paste0(output.dir, "RF_imputations/RF_obsData_results_", impnIterations, "_" , model.name ,".RData"))
print(paste0(output.dir, "RF_imputations/RF_obsData_results_", impnIterations, "_" , model.name ,".RData"))

#### RESULTS ----------------------------------------------------------------------------------------------------------------------------- ####

# Errors summarization
meanOBBE <- OBBE.lst %>%
  bind_rows() %>%
  summarise_all(., mean)

sdOBBE <- OBBE.lst %>%
  bind_rows() %>%
  summarise_all(., sd)

names(sdOBBE) <- paste0(names(sdOBBE), "_sd")

OBBE.rslts <- cbind(meanOBBE, sdOBBE)

write.csv(OBBE.rslts, paste0(output.dir, "RF_imputations/OBBE_obsData_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0(output.dir, "RF_imputations/OBBE_obsData_", impnIterations, "_" , model.name ,".csv"))

# HSM calculation and imputation summarisation

load(paste0(output.dir, "RF_imputations/RF_obsData_results_", impnIterations, "_" , model.name ,".RData"))

# Results in dataframe
imp.df <- imp.lst %>%
  bind_rows() %>%
  group_by(., Species) %>%
  dplyr::select(response)

# HSM calculation
imp.df$HSM <- imp.df$MinWP_md - imp.df$P50

# Mean imputation values
meanImp <- imp.df %>%
  summarise_all(., mean)

# Sd of imputation values
sdImp <- imp.df %>%
  summarise_all(., sd)

sdImp <- dplyr::rename(sdImp, P50_sd = P50, MinWP_md_sd = MinWP_md, HSM_sd = HSM) %>%
  dplyr::select(MinWP_md_sd, P50_sd, HSM_sd)

# Cv of imputation values
cvImp <- imp.df %>%
  summarise_all(., cv) %>% dplyr::select(MinWP_md, P50, HSM)
cvImp <- abs(cvImp)

cvImp <- dplyr::rename(cvImp, P50_cv = P50, MinWP_md_cv = MinWP_md, HSM_cv = HSM) %>%
  dplyr::select(MinWP_md_cv, P50_cv, HSM_cv)

imp.rslts <- cbind(meanImp, sdImp, cvImp, obs.rf.df[ complete.cases(obs.rf.df[, c(predictors)]), c("genus", "family", "order", predictors)])
if(!"MinWP_soil_SX" %in% predictors){
  imp.rslts <- cbind(imp.rslts, "MinWP_soil_SX" = obs.rf.df[complete.cases(obs.rf.df[, c(predictors)]), "MinWP_soil_SX"])
}
head(imp.rslts)

write.csv(imp.rslts, paste0(output.dir, "RF_imputations/RFimputation_obsData_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0(output.dir, "RF_imputations/RFimputation_obsData_", impnIterations, "_" , model.name ,".csv"))
