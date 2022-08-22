#### VARIABLES TRANSFORMATION AND CORRELATION ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

## P Sanchez-Martinez

print("Checking variables normality and transforming some of them to approach normality")

source("code/data_preparation/0-data_preparation_init.R")

plotDistributions <- F

#### FUNCTIONS---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

transfPlotFun <- function(df = ht, var = "MAT"){
  dens <- ggdensity(df[, var], main = paste0("Density plot of ", var), xlab = var)
  qq <- ggqqplot(df[, var], main = paste0("QQ plot of ", var))
  
  ## Log transformation
  if(all(na.omit(df[, var] < 0))){
    df[, var] <- -df[, var]
  }
  log.vec <- log(df[, var])
  dens_log <- ggdensity(log.vec, main = paste0("Density plot of log(", var, ")"), xlab = paste0("log(", var, ")"))
  qq_log <- ggqqplot(log.vec, main = paste0("QQ plot of log(", var, ")"))
  
  ## Sqrt transformation
  sqrt.vec <- sqrt(df[, var])
  dens_sqrt <- ggdensity(sqrt.vec, main = paste0("Density plot of sqrt(", var, ")"), xlab = paste0("sqrt(", var, ")"))
  qq_sqrt <- ggqqplot(sqrt.vec, main = paste0("QQ plot of sqrt(", var, ")"))
  
  p <- ggarrange(dens, qq, dens_log, qq_log, dens_sqrt, qq_sqrt)
  print(p)
}

checkPredictorCorrelations <- function(dta = ht, vars = c("log_nP50", "log_nYmin")){
  cat("Checking predictors correlations ...\n")
  corData <- dta[complete.cases(dta[, vars]), ]
  cor_Vsp <- round(cor(corData[, vars]), 3)
  f.out <- paste0( output.dir, "vars_colinearity.csv")
  write.csv(cor_Vsp, f.out, row.names=T)
  cat("  ==>", f.out, "\n")
  return(cor_Vsp)
}

### DATA ---------------------------------------------------------------------------------------------- ####

ht <- read.csv(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_2_taxonomy_scrub.csv"), header = T)


#### HYDRAULIC TRAITS --------------------------------------------------------------------------------- ####

if(isTRUE(plotDistributions)){
# P50
transfPlotFun(df = ht, var = "P50")

# MinWP_md
transfPlotFun(df = ht, var = "MinWP_md")

# HSM
transfPlotFun(df = ht, var = "HSM")

#### ENVIRONMENTAL VARS ----------------------------------------------------------------------------------- ####

names(ht)
vars <- names(ht)[35:length(ht)]
for(var in vars){
  print(var)
  transfPlotFun(df = ht, var = var)
}

}

#### VARIABLE TRANSFORMATIONS ---------------------------------------------------------------------------- ####

ht$log_negP50 <- log(-ht$P50)
ht$log_negP88 <- log(-ht$P88)
ht$log_negMinWP_md <- log(-ht$MinWP_md)
ht$log_Ks <- log(ht$Ks)
ht$log_Hv <- log(ht$Hv)
ht$log_Kl <- log(ht$Kl)

ht$sqrt_AP <- sqrt(ht$AP)
ht$sqrt_Prec_Wettest <- sqrt(ht$Prec_Wettest)
ht$sqrt_Prec_Driest <- sqrt(ht$Prec_Driest)
ht$sqrt_Prec_Seasonality <- sqrt(ht$Prec_Seasonality)
ht$sqrt_Prec_WetQ <- sqrt(ht$Prec_WetQ)
ht$sqrt_Prec_DryQ <- sqrt(ht$Prec_DryQ)
ht$sqrt_Prec_WarmQ <- sqrt(ht$Prec_WarmQ)
ht$log_Prec_ColdQ <- log(ht$Prec_ColdQ)
ht$log_TSeasonality <- log(ht$TSeasonality)
ht$log_TRange <- log(ht$TRange)
ht$log_windMax <- log(ht$windMax)
ht$log_wind <- log(ht$wind)
ht$log_ABDRock <- log(ht$ABDRock)
ht$sqrt_AI <- sqrt(ht$AI)
ht$log_VPDMax <- log(ht$VPDMax)
ht$log_CEC_30cm <- log(ht$CEC_30cm)
ht$log_OCarbon_30cm <- log(ht$OCarbon_30cm)

# Factors
names(ht)
facVars <- names(ht)[1:9]
# Factors
for(i in facVars){
  ht[, i] <- as.factor(ht[, i])
}

# Save final dataset

names(ht)
length(ht$Species) # (47581)

write.csv(ht, paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_3_varTransf.csv"), row.names = F)
print(paste0(processed.data, "treevol_workflow/treevol_globalDatabase_2020_3_varTransf.csv"))


#### VARIABLES COLINEALITY ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

## Correlation between environmental variables

envVarsCor <- checkPredictorCorrelations(dta = ht, vars = c(35:70))
envVarsCor
