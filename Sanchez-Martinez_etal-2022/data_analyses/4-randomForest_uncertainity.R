### RANDOM FORESTS IMPUTATION------------------------------------------------------------------------------------- ####


# P Sanchez-Martinez

remove(list = ls())
source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

#### PREDICTING USING BEST MISSING VALUES MULTIVARIATE RANDOM FOREST (VARYING SPECIES USED) ------------------------------------------------------------------------------------------- ####

# Parallelization options
showConnections()
cl <- makeCluster(nClusters)
registerDoParallel(cl)

model.name <- RFmodels$Model[24]
response <- unlist(strsplit(RFmodels[RFmodels$Model == model.name, "Response"], split=", "))
predictors <- unlist(strsplit(RFmodels[RFmodels$Model == model.name, "Predictors"], split=", "))

imp.lst <- list() # final list of results
for(n in 1:impnIterations){
  # NA creation
  train.df <- rf.df %>% dplyr::select(c("Species", response, predictors))
  train.df <- train.df[complete.cases(train.df[, c(predictors)]), ]
  species <- train.df$Species
  
  respNa.df <- prodNA(train.df[, c("P50", "MinWP_md")], noNA = 0.2)
  train.df$P50 <- respNa.df$P50
  train.df$MinWP_md <- respNa.df$MinWP_md
  
  # Predictive model
  train.df$Species <- NULL
  train.df$Species <- NULL
  
  inteRslts <- data.frame() # each iteration results
  rfImpCv <- missForest(train.df, parallelize = "forests", maxiter = 50, verbose = F, variablewise = T)
  # Imputed values using only a proportion of the species with data
  ximp <- rfImpCv$ximp
  ximp$Species <- species
  imp.lst[[n]] <- ximp
}
stopCluster(cl)


save(imp.lst, file = paste0("outputs/RF_imputations/RF_uncertainity_80_", impnIterations, "_" , model.name ,".RData"))
print(paste0("outputs/RF_imputations/RF_uncertainity_80", impnIterations, "_" , model.name ,".RData"))

#### RESULTS ----------------------------------------------------------------------------------------------------------------------------- ####

# HSM calculation and imputation summarisation

load(paste0("outputs/RF_imputations/RF_uncertainity_80_", impnIterations, "_" , model.name ,".RData"))

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
  summarise_all(., cv)
cvImp[,-1] <- abs(cvImp[,-1])

cvImp <- dplyr::rename(cvImp, P50_cv = P50, MinWP_md_cv = MinWP_md, HSM_cv = HSM) %>%
  dplyr::select(MinWP_md_cv, P50_cv, HSM_cv)

imp.rslts <- cbind(meanImp, sdImp, cvImp, rf.df[complete.cases(rf.df[, c(predictors)]), c("genus", "family", "order", predictors)])

if(!"MinWP_soil_SX" %in% predictors){
  imp.rslts <- cbind(imp.rslts, "MinWP_soil_SX" = rf.df[complete.cases(rf.df[, c(predictors)]), "MinWP_soil_SX"])
}
head(imp.rslts)

write.csv(imp.rslts, paste0("outputs/RF_imputations/RFimputation_uncertainity_80_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0("outputs/RF_imputations/RFimputation_uncertainity_80_", impnIterations, "_" , model.name ,".csv"))


#### PREDICTING USING BEST MISSING VALUES MULTIVARIATE RANDOM FOREST (ONLY PHYLOGENY) ------------------------------------------------------------------------------------------- ####

# Parallelization options
showConnections()
cl <- makeCluster(detectCores()/4)
registerDoParallel(cl)

# set.seed(123)
imp.lst <- list()
OBBE.lst <- list()
model.name <- RFmodels$Model[12]
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
save(imp.lst, file = paste0("outputs/RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))
print(paste0("outputs/RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))

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

write.csv(OBBE.rslts, paste0("outputs/RF_imputations/OBBE_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0("outputs/RF_imputations/OBBE_", impnIterations, "_" , model.name ,".csv"))

# HSM calculation and imputation summarisation

load(paste0("outputs/RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))

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
  summarise_all(., cv)
cvImp[,-1] <- abs(cvImp[,-1])

cvImp <- dplyr::rename(cvImp, P50_cv = P50, MinWP_md_cv = MinWP_md, HSM_cv = HSM) %>%
  dplyr::select(MinWP_md_cv, P50_cv, HSM_cv)

imp.rslts <- cbind(meanImp, sdImp, cvImp, rf.df[complete.cases(rf.df[, c(predictors)]), c("genus", "fmly", "order", predictors)])
if(!"MinWP_soil_SX" %in% predictors){
  imp.rslts <- cbind(imp.rslts, "MinWP_soil_SX" = rf.df[complete.cases(rf.df[, c(predictors)]), "MinWP_soil_SX"])
}
head(imp.rslts)

write.csv(imp.rslts, paste0("outputs/RF_imputations/RFimputation_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0("outputs/RF_imputations/RFimputation_", impnIterations, "_" , model.name ,".csv"))



#### PREDICTING USING BEST MISSING VALUES MULTIVARIATE RANDOM FOREST (ONLY CLIMATE) ------------------------------------------------------------------------------------------- ####

# Parallelization options
showConnections()
cl <- makeCluster(detectCores()/4)
registerDoParallel(cl)

# set.seed(123)
imp.lst <- list()
OBBE.lst <- list()
model.name <- RFmodels$Model[6]
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
save(imp.lst, file = paste0("outputs/RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))
print(paste0("outputs/RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))

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

write.csv(OBBE.rslts, paste0("outputs/RF_imputations/OBBE_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0("outputs/RF_imputations/OBBE_", impnIterations, "_" , model.name ,".csv"))

# HSM calculation and imputation summarisation

load(paste0("outputs/RF_imputations/RF_results_", impnIterations, "_" , model.name ,".RData"))

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
  summarise_all(., cv)
cvImp[,-1] <- abs(cvImp[,-1])

cvImp <- dplyr::rename(cvImp, P50_cv = P50, MinWP_md_cv = MinWP_md, HSM_cv = HSM) %>%
  dplyr::select(MinWP_md_cv, P50_cv, HSM_cv)

imp.rslts <- cbind(meanImp, sdImp, cvImp, rf.df[complete.cases(rf.df[, c(predictors)]), c("genus", "fmly", "order", predictors)])
if(!"MinWP_soil_SX" %in% predictors){
  imp.rslts <- cbind(imp.rslts, "MinWP_soil_SX" = rf.df[complete.cases(rf.df[, c(predictors)]), "MinWP_soil_SX"])
}
head(imp.rslts)

write.csv(imp.rslts, paste0("outputs/RF_imputations/RFimputation_", impnIterations, "_" , model.name ,".csv"), row.names = F)
print(paste0("outputs/RF_imputations/RFimputation_", impnIterations, "_" , model.name ,".csv"))

