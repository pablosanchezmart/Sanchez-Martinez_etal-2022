##### HSM AND MORTALITY ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #####

remove(list = ls())
gc()

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")


options(contrasts=c("contr.treatment","contr.poly")) #Unordered and Ordered categorical variables
# options(contrasts=c("contr.sum","contr.poly")) # Diferencia de cada grupo respecto del promedio.


# Figures size

h <- 2.61
w <- 7.87

maxent.reslts.dir <- paste0(output.dir, "maxent")
plotCrvs <- F # To plot curves is slower and require memory, so we'll only set this true the first time we run the code or when plotting is needed
ForceRun <- T  # needed to plot predictions properly
AggregatedPoints <- T
memory.limit(50000) # We'll need more memory here


results.dir <- rslts.dir
tables.dir <- NULL

if(isTRUE(AggregatedPoints)){
  tables.dir <- paste0(results.dir, "tables/aggregated_points/")
  rslts.dir <- paste0(results.dir, "figures/HSM_mortality/aggregated_points/")
} else{
  tables.dir <- paste0(results.dir, "tables/")
  rslts.dir <- paste0(results.dir, "figures/HSM_mortality/")
}


# Number of background points for maxent model
numberBg <- 10000

#### FUNCTIONS -------------------------------------------------------------------------------------------------------------------------------------- ####

logitMod <- function(frml = "p ~ HSM_min * biome.synth", dta = hsm_md.df, nIter = 2, plotCurves = plotCrvs){
  
  ## Modelling data preparation
  rslts.df <- data.frame()
  r2.vec <- numeric()
  aic.vec <- numeric()
  auc.vec <- numeric()
  trends.df <- data.frame()
  anova.rslts <- data.frame()
  fitted.df <- data.frame()
  mantelTestAutoc.vec <- numeric()
  
  # Variables (carefull if generalizing that, very code-specific)

  vars <- str_replace_all(frml, "\\+", " ")
  vars <- str_replace_all(vars, "\\*", " ")
  vars <- str_replace_all(vars, "\\~", " ")
  vars <- str_replace_all(vars, "0", " ")
  
  vars <- unlist(str_split(vars, pattern = " "))
  vars <- vars[-which(vars == "")]
  
  # Modelling
  for(i in 1:nIter){
    print(i)
    # Random points
    dta <- na.omit(dta)
    mort.dta <- dta %>% dplyr::filter(p == 1)
    nbgr <- length(mort.dta[ ,1])
    bg.dta <- dta %>% dplyr::filter(p == 0)
    bg.dta <- slice_sample(bg.dta, n = nbgr) # Same amount of mortality data, but we also used 10000 random points with similar results
    mod.dta <- rbind(mort.dta, bg.dta)
    
    # Model
    mod <- glm(as.formula(frml), data = mod.dta, family = "binomial")
    smry <- summary(mod)
    
    if(isTRUE(plotCurves)){
      # Predict (for curve response)
      fitted <- data.frame("fit" = predict(mod, mod.dta, type = "response"), "iter" = i)
      fitted <- cbind(fitted, mod.dta[, vars])
      fitted.df <- rbind(fitted.df, fitted) 
    }
    
    # Coefficients
    coef.df <- as.data.frame(smry$coefficients[, c("Estimate", "Std. Error", "Pr(>|z|)")])
    coef.df$Id <- row.names(coef.df)
    rslts.df <- rbind(rslts.df, coef.df)
    
    # Predictor significance
    anova.df <- as.data.frame(anova(mod, test = "Chisq"))
    anova.df <- cbind("variable" = row.names(anova.df), anova.df)
    anova.rslts <- rbind(anova.rslts, anova.df)
    
    # R2
    frml_r2 <- str_replace_all(frml, "0 +", " ")
    frml_r2 <- str_replace_all(frml, " 1 + ", "")
    
    r2 <- lrm(as.formula(frml_r2), data = mod.dta)
    r2 <- r2$stats["R2"]
    r2.vec[i] <- r2
    
    # AIc
    aic <- smry$aic
    aic.vec[i] <- aic
    
    # AUC

    # test and train
    sample = sample.split(mod.dta$p, SplitRatio = 0.80)
    train = subset(mod.dta, sample == TRUE)
    test  = subset(mod.dta, sample == FALSE)
    
    aucMod <- glm(as.formula(frml), data = train, family = "binomial")
    eval <- dismo::evaluate(p=test[which(test$p == 1), ], a=test[which(test$p == 0), ], model = aucMod)
    auc <- eval@auc
    auc.vec[i] <- auc
    
    # Significance of interaction by its own
    
    if(length(vars) > 2 && !is.na(str_extract(frml, "biome.synth")) | !is.na(str_extract(frml, "ft"))){
      if(!is.na(str_extract(frml, "biome.synth"))){
        trends.mod <- test(emtrends(mod, pairwise ~ biome.synth, var= vars[2]))
        trends.mod <- as.data.frame(trends.mod$emtrends)
        trends.df <- rbind(trends.df, trends.mod)
      }
      if(!is.na(str_extract(frml, "ft"))){
        trends.mod <- test(emtrends(mod, pairwise ~ ft, var= vars[2]))
        trends.mod <- as.data.frame(trends.mod$emtrends)
        trends.df <- rbind(trends.df, trends.mod)
      }
    }  
    # Spatial autocorrelation of residuals
    
    mod.dta$residuals <- mod$residuals
    
    geoDist <- dist(cbind(mod.dta$x, mod.dta$y))
    resDist <- dist(mod.dta$residuals)
    
    mantelTestAutocor <- mantel.rtest(geoDist, resDist, nrepet = 10)
    
    mantelTestAutoc.vec[i] <- mantelTestAutocor$obs
  }
  
  # Summarize coefficients
  coefficients.mean <- rslts.df %>% group_by(Id) %>% summarize_all(., mean)
  coefficients.sd <- rslts.df %>% group_by(Id) %>% summarize_all(., sd)
  coefficients.rslts <- cbind("model" = frml, coefficients.mean, coefficients.sd[, -1])
  names(coefficients.rslts) <- c("model", "parameter", "estimate_mean", "st_error_mean", "p_value_mean", "estimate_sd", "st_error_sd", "p_value_sd")
  r2.rslts <- data.frame("model" = frml,"R2_mean" = mean(r2.vec), "R2_sd" = sd(r2.vec))
  auc.rslts <- data.frame("model" = frml,"auc_mean" = mean(auc.vec), "auc_sd" = sd(auc.vec))
  
  # Significative effects for plot
  
  coefficients.rslts[which(coefficients.rslts$p_value_mean > 0.05), "p_value_mean_significance"] <- "NS"
  coefficients.rslts[which(coefficients.rslts$p_value_mean < 0.05), "p_value_mean_significance"] <- "."
  coefficients.rslts[which(coefficients.rslts$p_value_mean < 0.01), "p_value_mean_significance"] <- "*"
  coefficients.rslts[which(coefficients.rslts$p_value_mean < 0.001), "p_value_mean_significance"] <- "**"
  coefficients.rslts[which(coefficients.rslts$p_value_mean < 0.0001), "p_value_mean_significance"] <- "***"
  
  # Summarize predictors significance
  anova.mean <- anova.rslts %>% group_by(variable) %>% summarize_all(., mean)
  anova.sd <- anova.rslts %>% group_by(variable) %>% summarize_all(., sd)
  anova.rslts <-   cbind("model" = frml, anova.mean, anova.sd[, -1])
  
  names(anova.rslts) <- c("model", "variable", "df_mean", "deviance_mean", "resid_df_mean", "resid_dev_mean", "pr_Chi_mean", "df_sd", "deviance.sd", "resid_df_sd", "resid_dev_sd", "pr_Chi_sd")
  
  anova.rslts[which(anova.rslts$pr_Chi_mean > 0.05), "pr_Chi_mean_significance"] <- "NS"
  anova.rslts[which(anova.rslts$pr_Chi_mean < 0.05), "pr_Chi_mean_significance"] <- "."
  anova.rslts[which(anova.rslts$pr_Chi_mean < 0.01), "pr_Chi_mean_significance"] <- "*"
  anova.rslts[which(anova.rslts$pr_Chi_mean < 0.001), "pr_Chi_mean_significance"] <- "**"
  anova.rslts[which(anova.rslts$pr_Chi_mean < 0.0001), "pr_Chi_mean_significance"] <- "***"
  
  # Summarize trends (interaction significances by its own)
  if(length(vars) > 2 && !is.na(str_extract(frml, "biome.synth")) | !is.na(str_extract(frml, "ft"))){
    if(!is.na(str_extract(frml, "biome.synth"))){
      trends.mean <- trends.df %>% group_by(biome.synth) %>% summarize_all(., mean)
      trends.sd <- trends.df %>% group_by(biome.synth) %>% summarize_all(., sd)
      trends.rslts <- cbind("model" = frml, trends.mean, trends.sd[, -1])
      
      trends.rslts <- trends.rslts[, c(-5, -6, -10, -11)]
      
      # Significative effects for plot
      names(trends.rslts) <- c("model", "factor", "trend_mean", "st_error_mean", "p_value_mean",
                               "trend_sd", "st_error_sd", "p_value_sd")
      
      trends.rslts[which(trends.rslts$p_value_mean > 0.05), "p_value_mean_significance"] <- "NS"
      trends.rslts[which(trends.rslts$p_value_mean < 0.05), "p_value_mean_significance"] <- "."
      trends.rslts[which(trends.rslts$p_value_mean < 0.01), "p_value_mean_significance"] <- "*"
      trends.rslts[which(trends.rslts$p_value_mean < 0.001), "p_value_mean_significance"] <- "**"
      trends.rslts[which(trends.rslts$p_value_mean < 0.0001), "p_value_mean_significance"] <- "***"
    }
    
    if(!is.na(str_extract(frml, "ft"))){
      trends.mean <- trends.df %>% group_by(ft) %>% summarize_all(., mean)
      trends.sd <- trends.df %>% group_by(ft) %>% summarize_all(., sd)
      trends.rslts <- cbind("model" = frml, trends.mean, trends.sd[, -1])
      
      trends.rslts <- trends.rslts[, c(-5, -6, -10, -11)]
      
      # Significative effects for plot
      names(trends.rslts) <- c("model", "factor", "trend_mean", "st_error_mean", "p_value_mean",
                               "trend_sd", "st_error_sd", "p_value_sd")
      
      trends.rslts[which(trends.rslts$p_value_mean > 0.05), "p_value_mean_significance"] <- "NS"
      trends.rslts[which(trends.rslts$p_value_mean < 0.05), "p_value_mean_significance"] <- "."
      trends.rslts[which(trends.rslts$p_value_mean < 0.01), "p_value_mean_significance"] <- "*"
      trends.rslts[which(trends.rslts$p_value_mean < 0.001), "p_value_mean_significance"] <- "**"
      trends.rslts[which(trends.rslts$p_value_mean < 0.0001), "p_value_mean_significance"] <- "***"
    }
  }
  
  # Plots
  
  mod.dta[, vars[1]] <- as.factor(mod.dta[, vars[1]])
  
  if(!is.na(str_extract(frml, "log_ai")) && length(vars) > 2){
    vars <- vars[-which(vars == "log_ai")]
  }
  
  if(length(vars) == 2){
    p <- ggplot(mod.dta, aes(x = mod.dta[, vars[2]], y = mod.dta[, vars[1]], colour = mod.dta[, vars[1]]), fill = mod.dta[, vars[1]] %in% sign_coefs$parameter) + geom_boxplot() + 
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.text = element_blank(), legend.title = element_blank())  +
      scale_color_manual(values=c("1" ="#e34a33", "0" = "#2ca25f")) +  xlim(floor(min( mod.dta[, vars[2]])), ceiling(max( mod.dta[, vars[2]])))
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
      # curve <- curve + geom_point(alpha = 0.1, size = 0.5)
      curve <- curve + geom_smooth(method = "gam", color = "#e41a1c")  +
        scale_color_brewer(palette = "Spectral") + theme(legend.position = "none") + theme(axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.title = element_blank()) +
        xlim(floor(min( mod.dta[, vars[2]])), ceiling(max( mod.dta[, vars[2]])))
      
      # Save curves
      pdf(paste0(rslts.dir, "curves/", vars[2], "_curves.pdf"), height = w/3, width = w/3)
      # png(paste0("results/figures/HSM_mortality/curves/", vars[2], "_curves.png"), height = w/2, width = w/2, units = "in", res = 600)
      print(curve)
      dev.off()
      print(paste0(rslts.dir, "curves/", vars[2], "_curves.pdf"))
    }
  }
  if(length(vars) > 2 && !is.na(str_extract(frml, "biome.synth")) | !is.na(str_extract(frml, "ft"))){
    
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
      curve <- curve +  theme(axis.title.y = element_blank(),  axis.title.x = element_blank(), legend.title = element_blank(), axis.text.y = element_blank()) +
        scale_color_brewer(palette = "Spectral")
      curve <-  curve + theme(legend.position = "none") +   xlim(floor(min( mod.dta[, vars[2]])), ceiling(max( mod.dta[, vars[2]])))
      
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
  rslts$anova <- anova.rslts
  rslts$R2 <- r2.rslts
  rslts$aic <- mean(aic.vec)
  rslts$auc <- auc.rslts
  rslts$mantelTest <- mean(mantelTestAutoc.vec)
  
  # rslts$model_plot <- mod
  # rslts$boxPlot <- p
  # if(isTRUE(plotCrvs)){
  #   rslts$curve <- curve 
  # }
  if(length(vars) > 2 && !is.na(str_extract(frml, "biome.synth")) | !is.na(str_extract(frml, "ft"))){
    rslts$interaction_trends <- trends.rslts
    rslts$boxPlot_leg <- boxPlot_leg
  }
  # rslts$mod.data <- mod.dta
  return(rslts)
}

predMod <- function(mod.var = c("p", "HSM_mean", "negative_HSM_count",  "negative_HSM_prop", "HSM_variance", "biome.synth", "ft"),
                    dta = hsm_md.df, nIter = 2, nbg = 10000, maxent.args = c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                                                                             "hinge=TRUE", "threshold=FALSE", "addSamplesToBackground=FALSE",
                                                                             "responsecurves=TRUE", "writeplotdata","pictures=TRUE", "jackknife=TRUE"), 
                    maxent.results.dir){
  ## Modelling data preparation
  rslts <- list()
  testAucs <- numeric()
  testAics <- numeric()
  pred_crosval.lst <- list()
  
  # Data preparation
  dta <- na.omit(dta)
  dta$biome.synth <- as.factor(as.numeric(as.factor(dta$biome.synth)))
  dta$ft <- as.factor(as.numeric(as.factor(dta$ft)))
  
  # prediction database (keeping only background points, modelling variables and coordinates to project)
  pred.df <- dta %>% filter(p == 0)  # this is the background (which is the couverage of the environmental variables)
  pred.df$p <- NULL
  pred.df$predicted_p <- 0
  
  pred_crosval.df <- pred.df
  
  # Keep only modelling variables
  occ <- dta %>% filter( p == 1) %>% dplyr::select(x, y)
  dta <- dta[, mod.var]
  
  # Modelling
  for(i in 1:nIter){
    print(i)
    
    # Random points
    bg.dta <- dta %>% dplyr::filter(p == 0)
    bg.dta <- slice_sample(bg.dta, n = nbg)
    mort.dta <- dta %>% dplyr::filter(p == 1)
    mod.dta <- rbind(mort.dta, bg.dta)
    
    # test and train
    sample = sample.split(mod.dta$p, SplitRatio = 0.80)
    train = subset(mod.dta, sample == TRUE)
    test  = subset(mod.dta, sample == FALSE)
    
    # Model to evaluate
    mod <- maxent(x=train[, -1], p = train[, 1], args = maxent.args)
    
    # Evaluate
    df.sp.test <- test %>% filter(p == 1)
    df.bg.test <- test %>% filter(p == 0)
    # df.bg.test <- slice_sample(df.bg.test, n = length(df.sp.test$p))
    
    df.sp.test.predictions <- predict(mod, df.sp.test)
    df.bg.test.predictions <- predict(mod, df.bg.test)
    
    # eval <- dismo::evaluate(p=df.sp.test.predictions, a=df.bg.test.predictions, model = mod)
    # eval <- dismo::evaluate(p=df.sp.test.predictions, a=df.bg.test.predictions)
    eval <- dismo::evaluate(p=df.sp.test, a=df.bg.test, model = mod)
    # plot(eval, "ROC")
    
    testAuc <- eval@auc
    testAucs[i] <- testAuc
    
    pred_crosval.df[, paste0("predicted_p_", i)] <- predict(mod, pred.df[, mod.var[-1]]) # predict using the whole world dataframe (excluding p)
    # pred_crosval.df <-  pred_crosval.df[, c("x", "y", "predicted_p")]
    # head(pred_crosval.df)
    # coordinates(pred_crosval.df) <-pred_crosval.df[, c("x", "y")]
    # proj4string(pred_crosval.df) <- behrmann
    # length(pred_crosval.df$ID)
    # pred_crosval.df <- raster::aggregate(pred_crosval.df, FUN = mean, factor = 2)
    gc()
    removeTmpFiles()
  }
  
  # summarize test AUCs
  testAucs.df <- data.frame("testAUC_mean" = mean(testAucs), "testAUC_sd" = sd(testAucs))
  
  ## Crosvalidation standard deviation projection
  
  # Results in dataframe
  pred_crosval.df$predicted_p_mean <- apply(pred_crosval.df[, paste0("predicted_p_", 1:nIter)], 1, mean)
  pred_crosval.df$predicted_p_sd <- apply(pred_crosval.df[, paste0("predicted_p_", 1:nIter)], 1, sd)

  # plots 
  cols <- rev(brewer.pal(10, "Spectral"))
  
  p_unc_sd <- ggplot(pred_crosval.df) + geom_raster(aes(x=x, y=y, fill=predicted_p_sd))
  p_unc_sd <- p_unc_sd + scale_fill_gradientn(colours= cols)
  p_unc_sd <- p_unc_sd + coord_equal()
  p_unc_sd <- p_unc_sd + theme_light() +  theme(panel.grid.minor.x = element_blank(),
                                  panel.grid.minor.y = element_blank(),
                                  axis.title = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  legend.title = element_blank(),
                                  legend.position = "left", #legend.position = c(0.21, 0.19)
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  p_unc_mean <- ggplot(pred_crosval.df) + geom_raster(aes(x=x, y=y, fill=predicted_p_mean))
  p_unc_mean <- p_unc_mean + scale_fill_gradientn(colours= cols)
  p_unc_mean <- p_unc_mean + coord_equal()
  p_unc_mean <- p_unc_mean + theme_light() +  theme(panel.grid.minor.x = element_blank(),
                                                panel.grid.minor.y = element_blank(),
                                                axis.title = element_blank(),
                                                axis.text.x = element_blank(),
                                                axis.text.y = element_blank(),
                                                axis.ticks.x = element_blank(),
                                                axis.ticks.y = element_blank(),
                                                legend.title = element_blank(),
                                                legend.position = "left", #legend.position = c(0.21, 0.19)
                                                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # # Predict
  mod <- maxent(x=mod.dta[, -1], p = mod.dta[, 1], args = maxent.args, path = maxent.results.dir)
  
  pred.df$predicted_p <- predict(mod, pred.df[, mod.var[-1]]) # predict using the whole world dataframe (excluding p)
  
  gc()
  removeTmpFiles()
  
  # plot projection
  
  p <- ggplot(pred.df) + geom_raster(aes(x=x, y=y, fill=predicted_p))
  p <- p + scale_fill_gradientn(colours= cols)
  p <- p + coord_equal()
  p <- p + theme_light() +  theme(panel.grid.minor.x = element_blank(),
                                  panel.grid.minor.y = element_blank(),
                                  axis.title = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  legend.title = element_blank(),
                                  legend.position = "left", #legend.position = c(0.21, 0.19)
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  print(p)
  
  # Results
  
  rslts$model <- mod
  rslts$testAucs <- testAucs
  rslts$testAUC <- testAucs.df
  rslts$pred.df <- pred.df
  rslts$projection <- p
  rslts$pred_crosval_sd_projection <- p_unc_sd
  rslts$pred_crosval_mean_projection <- p_unc_mean
  
  return(rslts)
}


#### DATA ------------------------------------------------------------------------------------------------------------------------------------------- ####

if(!file.exists(paste0(processed.data, "mortality/mortality_logistic_worldwideProjection_dataset.csv")) | isTRUE(ForceRun)){
  # Mortality data
  load(file = "outputs/mortality/spatial_mortality_traits.RData") # object name: md_all_species_traits
  
  md <- read.csv("data/mortality/Hammond_etal_database_withspecies_transformed.csv", header = T)
  load("data/mortality/Hammond_etal_database_withspecies_transformed_geo.RData") # object named md.sf_trans
  # Biomes and functional types
  
  # # Land cover data (functional types)
  lc.r <- raster("data/species_data/ecoregions/land_cover_copernicus/land_cover.tif")
  
  # Load species distributions with imputed traits
  load(paste0(output.dir, "species_traits/all_species_traits_RFimp.RData"))
  all_species_traits_polygons <- sf::st_cast(all_species_traits, 'MULTIPOLYGON')
  
  # # Biomes
  biomes.r <- raster("outputs/biomes_ecoregions/biomes.tif")
  
  #### HSM and environment -------------------------------------------------------------------------------------------------------------------------------------------- ####
  
  ### Whole database preparation ####
  
  # HSM
  hsm.st <- stack(list.files(paste0(output.dir, "HSM/"), full.names = T))
  P50_st <- stack(list.files(paste0(output.dir, "P50/"), full.names = T))
  Pmin_st <- stack(list.files(paste0(output.dir, "Pmin/"), full.names = T))
  
  # Aridity index (moisture index)
  fnames <- list.files(path= "data/species_data/env_data/2.5m/", full.names=TRUE)
  fnames <- fnames[c(1, 29, 30, 44, 4, 8)]
  env.st <- stack(fnames)
  env.st$ai_et0 <- env.st$ai_et0/10000
  names(env.st) <- c("ai", "AP", "Srad", "Tmax", "Clay", "Sand")
  
  beginCluster(n = 2)
  env.st <- projectRaster(env.st, hsm.st)
  endCluster()
  
  # Biomes
  beginCluster(n = 2)
  biomes.r <- projectRaster(biomes.r, hsm.st, method = "ngb")
  endCluster()
  ecoRegions.sf <- st_read("data/species_data/ecoregions/Ecoregions2017.shp")
  
  beginCluster(n = 2)
  lc.r <- projectRaster(lc.r, hsm.st, method = "ngb")
  endCluster()
  
  # Stack HSM, environmental variables, biomes and land cover
  
  hsm.st$biomes.r <- biomes.r
  hsm.st$land_cover <- lc.r
  
  hsm.st <- stack(hsm.st, env.st, P50_st, Pmin_st)
  
  bg.df <- data.frame(as(hsm.st, "SpatialPixelsDataFrame"))
  bg.df$p <- 0
  bg.df$p[which(is.na(bg.df$HSM_mean))] <- NA  # to avoid background in the ocean (or places where we do not have HSM data)
  
  # Add mortality ocurrence points
  md.coord <- md[, c("x", "y")]
  # keep one mortality point each 10 km2?
  if(isTRUE(AggregatedPoints)){
    md.coord <- round(md.coord/10) * 10 
  } else{
    md.coord <- round(md.coord) 
  }
  md.coord <- unique(md.coord[, c("x", "y")])

  md_extr <- data.frame(raster::extract(hsm.st, md.coord))
  md.coord$p <- 1
  md.df <- cbind(md_extr, md.coord)
  
  # For those with na, add buffer
  md.df[which(is.na(md.df$biomes.r)), "biomes.r"] <- raster::extract(hsm.st$biomes.r, md.df[which(is.na(md.df$biomes.r)), c("x", "y")],
                                                                                   buffer = 20, fun = modal)
  summary(md.df)
  
  # Presence-backgroud database
  hsm_md.df <- rbind(md.df, bg.df) # to keep presences and when not present, absences (and not presences and absences where present)
  hsm_md.df$biomes.r <- as.factor(hsm_md.df$biomes.r)
  hsm_md.df$land_cover <- as.factor(hsm_md.df$land_cover)
  
  # BIOME AND FUNCTIONAL TYPE
  
  # Biome reclassification
  
  biome.recl <- data.frame("BIOME_NAME" = unique(as.factor(ecoRegions.sf$BIOME_NAME))[-15],
                           "BIOME_NUM" = unique(as.factor(ecoRegions.sf$BIOME_NUM)))
  
  biome.recl$BIOME_SYNTH <- c("Boreal", "Tropical and subtropical moist", "Mediterranean", "Desert and xeric", "Temperate",
                              "Boreal", "Temperate",
                              "Temperate", "Others", "Tropical and subtropical dry",
                              "Others", "Tropical and subtropical dry", "Tropical and subtropical dry", "Tropical and subtropical dry")
  
  for(biome_num in levels(hsm_md.df$biomes.r)){
    biome_name <- as.character(biome.recl[biome.recl$BIOME_NUM == biome_num, "BIOME_NAME"])
    hsm_md.df[which(hsm_md.df$biomes.r == biome_num), "biome.name"] <- biome_name
  
    biome_synth <- as.character(biome.recl[biome.recl$BIOME_NUM == biome_num, "BIOME_SYNTH"])
    hsm_md.df[which(hsm_md.df$biomes.r == biome_num), "biome.synth"] <- biome_synth
  
  }
  
  hsm_md.df$biome.name <- as.factor(hsm_md.df$biome.name)
  hsm_md.df$biome.synth <- as.factor(hsm_md.df$biome.synth)
  summary(hsm_md.df$biome.synth)
  
  # Functional type map (from land cover copernicus)
  
  ft.recl <- data.frame("FT_NAME" = c("Crops and Grassland", "Crops and Grassland", "Crops and Grassland", "Crops and Grassland", "Mosaic", "Mosaic", "Broadleaved evergreen",  "Broadleaved deciduous", "Broadleaved deciduous", "Broadleaved deciduous", "Needleleaved",  "Needleleaved", "Needleleaved",  "Needleleaved", "Needleleaved", "Mixed forest", "Mosaic", "Mosaic", "Shrubland",  "Shrubland",  "Shrubland", "Crops and Grassland", "Others", "Others", "Others",   "Others", "Others", "Others", "Others", "Others", "Others", "Others", "Others","Others",  "Others"),
                        "LC_NUM" = c("10",                         "11",                       "12",              "20",          "30",     "40",            "50",                      "60",                  "61",                   "62",                  "70",             "71",              "72",               "80",       "81",              "90",    "100",   "110",      "120",       "121",      "122",                 "130",     "140",    "150",   "152",      "153",     "160",   "170",     "180",     "190", "200",  "201", "202",      "210",   "220")
  )
  
  hsm_md.df$ft <- character(length = length(hsm_md.df$HSM_mean))
  for(lc in levels(hsm_md.df$land_cover)){
    ft_name <- as.character(ft.recl[ft.recl$LC_NUM == lc, "FT_NAME"])
    hsm_md.df[which(hsm_md.df$land_cover == lc), "ft"] <- ft_name
  }
  hsm_md.df$ft[hsm_md.df$ft == ""] <- NA
  hsm_md.df$ft <- as.factor(hsm_md.df$ft)
  summary(hsm_md.df$ft)
  
  ## Checking normality of predictors
  # densityplot(hsm_md.df$HSM_variance)
  
  hsm_md.df$log_negative_HSM_count <- log(hsm_md.df$negative_HSM_count)
  hsm_md.df$log_negative_HSM_count[hsm_md.df$log_negative_HSM_count == -Inf] <- 0
  
  hsm_md.df$sqrt_negative_HSM_count <- sqrt(hsm_md.df$negative_HSM_count)
  hsm_md.df$log_HSM_variance <- log(hsm_md.df$HSM_variance)
  hsm_md.df$log_HSM_variance[hsm_md.df$log_HSM_variance == -Inf] <- 0
  
  hsm_md.df$log_ai <- log(hsm_md.df$ai)
  hsm_md.df$log_ai[hsm_md.df$log_ai == -Inf] <- 0
  
  write.csv(hsm_md.df, paste0(processed.data, "mortality/mortality_logistic_worldwideProjection_dataset.csv"), row.names =  F)
  print(paste0(processed.data, "mortality/mortality_logistic_worldwideProjection_dataset.csv"))

} else {
  hsm_md.df <- read.csv(paste0(processed.data, "mortality/mortality_logistic_worldwideProjection_dataset.csv"), header = T)
}

# Number of observations per biome and functional type
N_biome_ft.df <- hsm_md.df %>% filter(p ==1) %>% dplyr::select(biome.synth, ft)
length(N_biome_ft.df$biome.synth)

N_biome.df <- count(N_biome_ft.df$biome.synth)
N_ft.df <- count(N_biome_ft.df$ft)

write.csv(N_biome.df, paste0(tables.dir, "number_mortalityObs_biome.csv"))
print(paste0(tables.dir, "number_mortalityObs_biome.csv"))


#### UNIVARIATE AND BIVARIATE MODELS ---------------------------------------------------------------------------------------------------------------- ####

if(!dir.exists((paste0(output.dir, "HSM_mortality_logit_models")))){
  dir.create(paste0(output.dir, "HSM_mortality_logit_models"))
}

### Previous reesults ####

if(file.exists(paste0(tables.dir, "glm_auc_100ite.csv"))){
  all_auc <- read.csv(paste0(tables.dir, "glm_auc_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_auc_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_aic_100ite.csv"))){
  all_aics <- read.csv(paste0(tables.dir, "glm_aic_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_aic_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_coefficients_100ite.csv"))){
  all_coefs <- read.csv(paste0(tables.dir, "glm_coefficients_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_coefficients_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_anova_100ite.csv"))){
  all_anova <- read.csv(paste0(tables.dir, "glm_anova_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_anova_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_anova_ai_100ite.csv"))){
  all_anova_ai <- read.csv(paste0(tables.dir, "glm_anova_ai_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_anova_ai_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_r2_100ite.csv"))){
  all_r2 <- read.csv(paste0(tables.dir, "glm_r2_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_r2_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_r2_ai_100ite.csv"))){
  all_r2_ai <- read.csv(paste0(tables.dir, "glm_r2_ai_100ite.csv"), header = T)
  print(paste0("reading previous results ", "glm_r2_ai_100ite.csv"))
}

if(file.exists(paste0(tables.dir, "glm_interactions_trends.csv"))){
  all_trends <- read.csv(paste0(tables.dir, "glm_interactions_trends.csv"), header = T)
  print(paste0("reading previous results ", "glm_interactions_trends.csv"))
}

if(file.exists(paste0(tables.dir, "glm_interactions_trends_ai.csv"))){
  all_trends_ai <- read.csv(paste0(tables.dir, "glm_interactions_trends_ai.csv"), header = T)
  print(paste0("reading previous results ", "glm_interactions_trends_ai.csv"))
}


### Mortality vs. aridity (moisture) index (climate) ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/aridity_idex.RData")) && isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/aridity_idex.RData"))
} else {
  md_ai.mod <- logitMod(frml = "p ~ log_ai", dta = hsm_md.df, nIter = 100)
  md_ai.mod$mantelTest
  save(md_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/aridity_idex.RData"))
  cat(output.dir, "HSM_mortality_logit_models/aridity_idex.RData\n")
}


### Mortality vs. aridity (moisture) index (climate) * biome ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/aridity_idex_biome.RData")) && isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/aridity_idex_biome.RData"))
} else {
  md_ai_biome.mod <- logitMod(frml = "p ~ log_ai * biome.synth", dta = hsm_md.df, nIter = 100)
  
  save(md_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/aridity_idex_biome.RData"))
  cat(output.dir, "HSM_mortality_logit_models/aridity_idex_biome.RData\n")
}

### Mortality vs. tmax (climate) ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/tmax.RData")) && isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/tmax.RData"))
} else {
  md_tmax.mod <- logitMod(frml = "p ~ Tmax", dta = hsm_md.df, nIter = 100)
  
  save(md_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/tmax.RData"))
  cat(output.dir, "HSM_mortality_logit_models/tmax.RData\n")
}

### Mortality vs. AP (climate) ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/AP.RData")) && isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/AP.RData"))
} else {
  md_ap.mod <- logitMod(frml = "p ~ AP", dta = hsm_md.df, nIter = 100)
  
  save(md_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/AP.RData"))
  cat(output.dir, "HSM_mortality_logit_models/AP.RData\n")
}


### Mortality vs. HSM mean (and randompoints plot and boxplot legend) ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean.RData"))
} else {
  md_hsm_mean.mod <- logitMod(frml = "p ~ HSM_mean", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_mean.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_mean.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_mean.RData\n")
}

# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ai.RData"))
} else {
  md_hsm_mean_ai.mod <- logitMod(frml = "p ~ HSM_mean + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_mean_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_mean_ai.RData\n")
}

## plot random points (one set of background)

p <- plotMapFun(hsm.st$HSM_mean, var = "HSM_mean")
p <- p + geom_point(data = md_hsm_mean.mod$mod.data[md_hsm_mean.mod$mod.data$p == 0, ], aes(x = x, y = y), colour = "#2ca25f", shape = 3, size = 0.2)
p <- p + geom_point(data = md_hsm_mean.mod$mod.data[md_hsm_mean.mod$mod.data$p == 1, ], aes(x = x, y = y), colour = "#e34a33", shape = 3, size = 0.2) + theme(legend.position = "right")
p

png(paste0(rslts.dir, "mortality_background_points.png"), height = h, width = w, units = "in", res = 1000)
p
dev.off()
print(paste0(rslts.dir, "mortality_background_points.png"))

  
### Mortality vs. HSM mean * biome ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_biome.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_biome.RData"))
} else {
  md_hsm_mean_biome.mod <- logitMod(frml = "p ~ HSM_mean * biome.synth", dta = hsm_md.df, nIter = 10)
  gc()
  removeTmpFiles()
  save(md_hsm_mean_biome.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_biome.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_mean_biome.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_biome_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_biome_ai.RData"))
} else {
  md_hsm_mean_biome_ai.mod <- logitMod(frml = "p ~ HSM_mean * biome.synth + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_mean_biome_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_biome_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_mean_biome_ai.RData\n")
}

# Biome boxplot legend
pdf("results/figures/HSM_mortality/HSM__mortality_biome_boxplot_legend.pdf", height = w/15, width = w/18)
plot(md_hsm_mean_biome.mod$boxPlot_leg)
dev.off()
print("results/figures/HSM_mortality/HSM__mortality_biome_boxplot_legend.pdf")


### Mortality vs. HSM mean * functional type ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ft.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ft.RData"))
} else {
  md_hsm_mean_ft.mod <- logitMod(frml = "p ~ HSM_mean * ft", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_mean_ft.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ft.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_mean_ft.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ft_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ft_ai.RData"))
} else {
  md_hsm_mean_ft_ai.mod <- logitMod(frml = "p ~ HSM_mean * ft + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_mean_ft_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_mean_ft_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_mean_ft_ai.RData\n")
}


### Mortality vs. HSM min ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_min.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_min.RData"))
} else {
  md_hsm_min.mod <- logitMod(frml = "p ~ HSM_min", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_min.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_min.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_min.RData\n")
}

# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ai.RData"))
} else {
  md_hsm_min_ai.mod <- logitMod(frml = "p ~ HSM_min + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_min_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_min_ai.RData\n")
}


### Mortality vs. HSM min * biome ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_biome.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_biome.RData"))
} else {
  md_hsm_min_biome.mod <- logitMod(frml = "p ~ HSM_min * biome.synth", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_min_biome.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_min_biome.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_min_biome.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_biome_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_biome_ai.RData"))
} else {
  md_hsm_min_biome_ai.mod <- logitMod(frml = "p ~ HSM_min * biome.synth + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_min_biome_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_min_biome_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_min_biome_ai.RData\n")
}


### Mortality vs. HSM min * functional type ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ft.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ft.RData"))
} else {
  md_hsm_min_ft.mod <- logitMod(frml = "p ~ HSM_min * ft", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_min_ft.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ft.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_min_ft.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ft_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ft_ai.RData"))
} else {
  md_hsm_min_ft_ai.mod <- logitMod(frml = "p ~ HSM_min * ft + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_min_ft_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_min_ft_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_min_ft_ai.RData\n")
}


### Mortality vs. HSM diversity ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance.RData"))
} else {
  md_hsm_variance.mod <- logitMod(frml = "p ~ log_HSM_variance", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ai.RData"))&& isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ai.RData"))
} else {
  md_hsm_variance_ai.mod <- logitMod(frml = "p ~ log_HSM_variance + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_ai.RData\n")
}


### Mortality vs. HSM diversity * biome ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_biome.RData")) && isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_biome.RData"))
} else {
  md_hsm_variance_biome.mod <- logitMod(frml = "p ~ log_HSM_variance * biome.synth", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_biome.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_biome.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_biome.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_biome_ai.RData")) && isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_biome_ai.RData"))
} else {
  md_hsm_variance_biome_ai.mod <- logitMod(frml = "p ~ log_HSM_variance * biome.synth + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_biome_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_biome_ai.RData\n")
}


### Mortality vs. HSM diversity * functional type ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ft.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ft.RData"))
} else {
  md_hsm_variance_ft.mod <- logitMod(frml = "p ~ log_HSM_variance * ft", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_ft.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ft.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_ft.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ft_ai.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ft_ai.RData"))
} else {
  md_hsm_variance_ft_ai.mod <- logitMod(frml = "p ~ log_HSM_variance * ft + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_ft_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_ft_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_ft.RData\n")
}


### Mortality vs. (sqrt) number of species with negative HSM ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount.RData"))
} else {
  md_hsm_sqrt_negCount.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ai.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ai.RData"))
} else {
  md_hsm_sqrt_negCount_ai.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ai.RData\n")
}


# Including N-species

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_NSpecies.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_NSpecies.RData"))
} else {
  md_hsm_sqrt_negCount_Nspp.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count + N_species", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_NSpecies.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_NSpecies.RData\n")
}


### Mortality vs. (sqrt) number of species with negative HSM * biome ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome.RData"))
} else {
  md_hsm_sqrt_negCount_biome.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count * biome.synth", dta = hsm_md.df, nIter = 100) 
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount_biome.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome_ai.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome_ai.RData"))
} else {
  md_hsm_sqrt_negCount_biome_ai.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count * biome.synth + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount_biome_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_biome_ai.RData\n")
}


### Mortality vs. (sqrt) number of species with negative HSM * functional type ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft.RData"))
} else {
  md_hsm_sqrt_negCount_ft.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count * ft", dta = hsm_md.df, nIter = 100) 
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount_ft.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft.RData\n")
}


# Including AI

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft_ai.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft_ai.RData"))
} else {
  md_hsm_sqrt_negCount_ft_ai.mod <- logitMod(frml = "p ~ sqrt_negative_HSM_count * ft + log_ai", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_sqrt_negCount_ft_ai.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft_ai.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_sqrt_negCount_ft_ai.RData\n")
}


### Mortality vs. HSM diversity * HSM min ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_min.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_min.RData"))
} else {
  md_hsm_variance_min.mod <- logitMod(frml = "p ~ HSM_min + log_HSM_variance", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_min.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_min.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_min.RData\n")
}


if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_min_int.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_min_int.RData"))
} else {
  md_hsm_variance_min_int.mod <- logitMod(frml = "p ~ HSM_min * log_HSM_variance", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_variance_min_int.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_variance_min_int.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_variance_min_int.RData\n")
}
md_hsm_variance_min_int.mod$coefficients

### Mortality vs. HSM max and HSM min ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min.RData"))
} else {
  md_hsm_max_min.mod <- logitMod(frml = "p ~ HSM_min + HSM_max", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_max_min.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_max_min.RData\n")
}


if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min_int.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min_int.RData"))
} else {
  md_hsm_max_min_int.mod <- logitMod(frml = "p ~ HSM_min * HSM_max", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_max_min_int.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min_int.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_max_min_int.RData\n")
}


### Mortality vs. HSM max, HSM min and HSM richness ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min_variance.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min_variance.RData"))
} else {
  md_hsm_max_min_variance.mod <- logitMod(frml = "p ~ HSM_min + HSM_max + log_HSM_variance", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(md_hsm_max_min_variance.mod, file = paste0(output.dir, "HSM_mortality_logit_models/hsm_max_min_variance.RData"))
  cat(output.dir, "HSM_mortality_logit_models/hsm_max_min_variance.RData\n")
}
md_hsm_max_min_variance.mod$coefficients

### Mortality vs. all HSM variables ####

if(file.exists(paste0(output.dir, "HSM_mortality_logit_models/maxent_vars_model.RData"))| isFALSE(ForceRun)){
  load(paste0(output.dir, "HSM_mortality_logit_models/all_hsm.RData"))
} else {
  maxent_vars_model.mod <- logitMod(frml = "p ~ Tmax + ai + log_HSM_variance + sqrt_negative_HSM_count + biome.synth + ft", dta = hsm_md.df, nIter = 100)
  gc()
  removeTmpFiles()
  save(maxent_vars_model.mod, file = paste0(output.dir, "HSM_mortality_logit_models/maxent_vars_model.RData"))
  cat(output.dir, "HSM_mortality_logit_models/maxent_vars_model.RData\n")
}

### Results all models ####

# all models autocorrelation

 min(  md_ai.mod$mantelTest, 
                  md_ap.mod$mantelTest,
                  md_tmax.mod$mantelTest, 
                  md_hsm_mean.mod$mantelTest, md_hsm_mean_biome.mod$mantelTest, md_hsm_mean_ft.mod$mantelTest,
                  md_hsm_min.mod$mantelTest, md_hsm_min_biome.mod$mantelTest, md_hsm_min_ft.mod$mantelTest,
                  md_hsm_sqrt_negCount.mod$mantelTest, md_hsm_sqrt_negCount_biome.mod$mantelTest, md_hsm_sqrt_negCount_ft.mod$mantelTest,
                  md_hsm_variance.mod$mantelTest, md_hsm_variance_biome.mod$mantelTest, md_hsm_variance_ft.mod$mantelTest)


# Model auc

all_auc <- rbind(  md_ai.mod$auc, 
                   md_ap.mod$auc,
                   md_tmax.mod$auc, 
                   md_hsm_mean.mod$auc, md_hsm_mean_biome.mod$auc, md_hsm_mean_ft.mod$auc,
                   md_hsm_min.mod$auc, md_hsm_min_biome.mod$auc, md_hsm_min_ft.mod$auc,
                   # md_hsm_max.mod$auc, md_hsm_max_biome.mod$auc, md_hsm_max_ft.mod$auc,
                   md_hsm_sqrt_negCount.mod$auc, md_hsm_sqrt_negCount_biome.mod$auc, md_hsm_sqrt_negCount_ft.mod$auc,
                   md_hsm_variance.mod$auc, md_hsm_variance_biome.mod$auc, md_hsm_variance_ft.mod$auc)

write.csv(all_auc, paste0(tables.dir, "glm_auc_100ite.csv"), row.names = F)
print(paste0(tables.dir, "glm_auc_100ite.csv"))

# Model aic
all_aics <- rbind("md_ai.mod" = md_ai.mod$aic, "md_ap.mod" = md_ap.mod$aic, "md_tmax.mod" = md_tmax.mod$aic,
  "md_hsm_mean.mod" = md_hsm_mean.mod$aic, "md_hsm_mean_biome.mod" = md_hsm_mean_biome.mod$aic, "md_hsm_mean_ft.mod" = md_hsm_mean_ft.mod$aic,
                  "md_hsm_min.mod" = md_hsm_min.mod$aic, "md_hsm_min_biome.mod" = md_hsm_min_biome.mod$aic, "md_hsm_min_ft.mod" = md_hsm_min_ft.mod$aic,
                  # "md_hsm_max.mod" = md_hsm_max.mod$aic, "md_hsm_max_biome.mod" = md_hsm_max_biome.mod$aic, "md_hsm_max_ft.mod" = md_hsm_max_ft.mod$aic,
                  "md_hsm_sqrt_negCount.mod" = md_hsm_sqrt_negCount.mod$aic, "md_hsm_sqrt_negCount_biome.mod" = md_hsm_sqrt_negCount_biome.mod$aic, "md_hsm_sqrt_negCount_ft.mod" = md_hsm_sqrt_negCount_ft.mod$aic,
                  "md_hsm_variance.mod" = md_hsm_variance.mod$aic, "md_hsm_variance_biome.mod" = md_hsm_variance_biome.mod$aic, "md_hsm_variance_ft.mod" = md_hsm_variance_ft.mod$aic)
colnames(all_aics) <- "aic"

write.csv(all_aics, paste0(tables.dir, "glm_aic_100ite.csv"), row.names = T)
print(paste0(tables.dir, "glm_aic_100ite.csv"))

# Model coefficients
all_coefs <- rbind(md_hsm_mean.mod$coefficients, md_hsm_mean_biome.mod$coefficients, md_hsm_mean_ft.mod$coefficients,
                   md_hsm_min.mod$coefficients, md_hsm_min_biome.mod$coefficients, md_hsm_min_ft.mod$coefficients,
                   # md_hsm_max.mod$coefficients, md_hsm_max_biome.mod$coefficients, md_hsm_max_ft.mod$coefficients,
                   md_hsm_sqrt_negCount.mod$coefficients, md_hsm_sqrt_negCount_biome.mod$coefficients, md_hsm_sqrt_negCount_ft.mod$coefficients,
                   md_hsm_variance.mod$coefficients, md_hsm_variance_biome.mod$coefficients, md_hsm_variance_ft.mod$coefficients)

write.csv(all_coefs, paste0(tables.dir, "glm_coefficients_100ite.csv"), row.names = F)
print(paste0(tables.dir, "glm_coefficients_100ite.csv"))


# all_coefs_ai <- rbind(md_hsm_mean_ai.mod$coefficients, md_hsm_mean_biome_ai.mod$coefficients, md_hsm_mean_ft_ai.mod$coefficients,
#                       md_hsm_variance_ai.mod$coefficients, md_hsm_variance_biome_ai.mod$coefficients, md_hsm_variance_ft_ai.mod$coefficients,
#                       md_hsm_min_ai.mod$coefficients, md_hsm_min_biome_ai.mod$coefficients, md_hsm_min_ft_ai.mod$coefficients,
#                       # md_hsm_max_ai.mod$coefficients, md_hsm_max_biome_ai.mod$coefficients, md_hsm_max_ft_ai.mod$coefficients,
#                       md_hsm_sqrt_negCount_ai.mod$coefficients, md_hsm_sqrt_negCount_biome_ai.mod$coefficients, md_hsm_sqrt_negCount_ft_ai.mod$coefficients)
# 
# write.csv(all_coefs_ai, paste0(tables.dir, "glm_coefficients_ai_100ite.csv"), row.names = F)
# print(paste0(tables.dir, "glm_coefficients_ai_100ite.csv"))

# Variable significance (anova)
all_anova <- rbind(md_hsm_mean.mod$anova, md_hsm_mean_biome.mod$anova, md_hsm_mean_ft.mod$anova,
                   md_hsm_min.mod$anova, md_hsm_min_biome.mod$anova, md_hsm_min_ft.mod$anova,
                   # md_hsm_max.mod$anova, md_hsm_max_biome.mod$anova, md_hsm_max_ft.mod$anova,
                   md_hsm_sqrt_negCount.mod$anova, md_hsm_sqrt_negCount_biome.mod$anova, md_hsm_sqrt_negCount_ft.mod$anova,
                   md_hsm_variance.mod$anova, md_hsm_variance_biome.mod$anova, md_hsm_variance_ft.mod$anova)

head(all_anova)

all_anova[, c(-1, -2, -13)] <- round(all_anova[, c(-1, -2, -13)], 3)

all_anova <- as_tibble(all_anova)

write.csv(all_anova, paste0(tables.dir, "glm_anova_100ite.csv"), row.names = F)
print(paste0(tables.dir, "glm_anova_100ite.csv"))

all_anova_ai <- rbind(md_hsm_mean_ai.mod$anova, md_hsm_mean_biome_ai.mod$anova, md_hsm_mean_ft_ai.mod$anova,
                   md_hsm_variance_ai.mod$anova, md_hsm_variance_biome_ai.mod$anova, md_hsm_variance_ft_ai.mod$anova,
                   md_hsm_min_ai.mod$anova, md_hsm_min_biome_ai.mod$anova, md_hsm_min_ft_ai.mod$anova,
                   # md_hsm_max_ai.mod$anova, md_hsm_max_biome_ai.mod$anova, md_hsm_max_ft_ai.mod$anova,
                   md_hsm_sqrt_negCount_ai.mod$anova, md_hsm_sqrt_negCount_biome_ai.mod$anova, md_hsm_sqrt_negCount_ft_ai.mod$anova)

write.csv(all_anova_ai, paste0(tables.dir, "glm_anova_ai_100ite.csv"), row.names = F)
print(paste0(tables.dir, "glm_anova_ai_100ite.csv"))

# Models R2
all_r2 <- rbind(   md_ai.mod$R2, 
                   md_ap.mod$R2,
                   md_tmax.mod$R2,
                   md_hsm_mean.mod$R2, md_hsm_mean_biome.mod$R2, md_hsm_mean_ft.mod$R2,
                   md_hsm_min.mod$R2, md_hsm_min_biome.mod$R2, md_hsm_min_ft.mod$R2,
                   # md_hsm_max.mod$R2, md_hsm_max_biome.mod$R2, md_hsm_max_ft.mod$R2,
                   md_hsm_sqrt_negCount.mod$R2, md_hsm_sqrt_negCount_biome.mod$R2, md_hsm_sqrt_negCount_ft.mod$R2,
                   md_hsm_variance.mod$R2, md_hsm_variance_biome.mod$R2, md_hsm_variance_ft.mod$R2)

write.csv(all_r2, paste0(tables.dir, "glm_r2_100ite.csv"), row.names = F)
print(paste0(tables.dir, "glm_r2_100ite.csv"))

all_r2_ai <- rbind(md_hsm_mean_ai.mod$R2, md_hsm_mean_biome_ai.mod$R2, md_hsm_mean_ft_ai.mod$R2,
                      md_hsm_variance_ai.mod$R2, md_hsm_variance_biome_ai.mod$R2, md_hsm_variance_ft_ai.mod$R2,
                      md_hsm_min_ai.mod$R2, md_hsm_min_biome_ai.mod$R2, md_hsm_min_ft_ai.mod$R2,
                      # md_hsm_max_ai.mod$R2, md_hsm_max_biome_ai.mod$R2, md_hsm_max_ft_ai.mod$R2,
                      md_hsm_sqrt_negCount_ai.mod$R2, md_hsm_sqrt_negCount_biome_ai.mod$R2, md_hsm_sqrt_negCount_ft_ai.mod$R2)

write.csv(all_r2_ai, paste0(tables.dir, "glm_r2_ai_100ite.csv"), row.names = F)
print(paste0(tables.dir, "glm_r2_ai_100ite.csv"))

# Interaction trends significance (for boxplot)

all_trends <- rbind(md_hsm_mean_biome.mod$interaction_trends, md_hsm_mean_ft.mod$interaction_trends,
                    md_hsm_min_biome.mod$interaction_trends, md_hsm_min_ft.mod$interaction_trends,
                    # md_hsm_max_biome.mod$interaction_trends, md_hsm_max_ft.mod$interaction_trends,
                    md_hsm_sqrt_negCount_biome.mod$interaction_trends, md_hsm_sqrt_negCount_ft.mod$interaction_trends,
                    md_hsm_variance_biome.mod$interaction_trends, md_hsm_variance_ft.mod$interaction_trends)

write.csv(all_trends, paste0(tables.dir, "glm_interactions_trends.csv"), row.names = F)
print(paste0(tables.dir, "glm_interactions_trends.csv"))

all_trends_ai <- rbind(md_hsm_mean_biome_ai.mod$interaction_trends, md_hsm_mean_ft_ai.mod$interaction_trends,
                    md_hsm_variance_biome_ai.mod$interaction_trends, md_hsm_variance_ft_ai.mod$interaction_trends,
                    md_hsm_min_biome_ai.mod$interaction_trends, md_hsm_min_ft_ai.mod$interaction_trends,
                    # md_hsm_max_biome_ai.mod$interaction_trends, md_hsm_max_ft_ai.mod$interaction_trends,
                    md_hsm_sqrt_negCount_biome_ai.mod$interaction_trends, md_hsm_sqrt_negCount_ft_ai.mod$interaction_trends)

write.csv(all_trends_ai, paste0(tables.dir, "glm_interactions_trends_ai.csv"), row.names = F)
print(paste0(tables.dir, "glm_interactions_trends_ai.csv"))


#### PREDITIVE MODELS (MAXENT) ---------------------------------------------------------------------------------------------------------------- ####

tables.dir <- paste0(results.dir, "tables/")
rslts.dir <- paste0(results.dir, "figures/HSM_mortality/")

# Correlation between variables (variables included in models will have a correlation loewr than 0)
source("http://www.sthda.com/upload/rquery_cormat.r")

# pdf(paste0(rslts.dir, "/correlation_matrix.pdf"))
correlation.rslts <- cor(hsm_md.df[, c("AP", "Srad", "Tmax", "Clay", "Sand", "ai", 
                                      "HSM_mean", "HSM_min", "negative_HSM_count", "negative_HSM_prop", "HSM_variance", "N_species")], use = "pairwise.complete.obs") # with p value
# dev.off()
# print(paste0(rslts.dir, "/correlation_matrix.pdf"))
round(correlation.rslts, 2)


### Using FT + biome  ####

if(!file.exists(paste0(maxent.reslts.dir, "/ft_biome_model_results", numberBg, ".RData")) | isTRUE(ForceRun)){
  me_ft_biome.rslts <- predMod(mod.var = c("p", "biome.synth", "ft"),
                               dta = hsm_md.df, nIter = 10, nbg = numberBg, maxent.args = c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                                                                                            "hinge=TRUE", "threshold=FALSE", "addSamplesToBackground=FALSE",
                                                                                            "responsecurves=TRUE", "writeplotdata","pictures=TRUE", "jackknife=TRUE"),
                               maxent.results.dir = paste0(maxent.reslts.dir, "/ft_biome"))
  
  save(me_ft_biome.rslts, file = paste0(maxent.reslts.dir, "/ft_biome_model_results", numberBg, ".RData"))
  print(paste0(maxent.reslts.dir, "/ft_biome_model_results", numberBg, ".RData"))
} else {
  load(paste0(maxent.reslts.dir, "/ft_biome_model_results", numberBg, ".RData")) 
}

# Evaluation
me_ft_biome.rslts$testAUC
paste0(maxent.reslts.dir, "/ft_biome")


# Projection
png(paste0(rslts.dir, "projection/ft_biome_model_mortality_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_ft_biome.rslts$projection + theme(legend.position = "none") 
p
dev.off()
print(paste0(rslts.dir, "projection/ft_biome_model_mortality_projection.png"))


# response curves
response(me_ft_biome.rslts$model)

# variable contribution (for the prediction)
df <- data.frame("variable" = names(me_ft_biome.rslts$model@results[c(7,8), ]), 
                 "contribution" = me_ft_biome.rslts$model@results[c(7,8), ],
                 "permutation.importance" = me_ft_biome.rslts$model@results[c(9, 10), ])
df <- df[order(df$contribution, decreasing = T), ]
df$variable <-  as.factor(df$variable)

p <- ggplot(df, aes(contribution, reorder(variable, order(contribution, decreasing = F)), fill = permutation.importance)) + geom_col() + theme_minimal() + 
  theme(axis.text.y = element_blank(), axis.title = element_blank()) + 
  scale_fill_gradient(low = "#fee0d2", high = "#de2d26")
leg <- get_legend(p)
p <- p + theme(legend.position = "none")

pdf(paste0(rslts.dir, "projection/ft_biome_model_mortality_var_contr.pdf"), height = h, width = h)
p
dev.off()
print(paste0(rslts.dir, "projection/ft_biome_model_mortality_var_contr.pdf"))
df

pdf(paste0(rslts.dir, "projection/ft_biome_model_mortality_var_contr_legend.pdf"))
plot(leg)
dev.off()
print(paste0(rslts.dir, "projection/ft_biome_model_mortality_var_contr_legend.pdf"))

remove(me_ft_biome.rslts)

### Using FT + biome + environmental data ####

if(!file.exists(paste0(maxent.reslts.dir, "/clim_model_results", numberBg, ".RData")) | isTRUE(ForceRun)){
    me_env.rslts <- predMod(mod.var = c("p", "Tmax", "Clay", "Sand", "ai", "biome.synth", "ft"),
                            dta = hsm_md.df, nIter = 10, nbg = numberBg, maxent.args = c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                                                                                       "hinge=TRUE", "threshold=FALSE", "addSamplesToBackground=FALSE",
                                                                                       "responsecurves=TRUE", "writeplotdata","pictures=TRUE", "jackknife=TRUE"),
                            maxent.results.dir = paste0(maxent.reslts.dir, "/clim"))
    save(me_env.rslts, file = paste0(maxent.reslts.dir, "/clim_model_results", numberBg, ".RData"))
    print(paste0(maxent.reslts.dir, "/clim_model_results", numberBg, ".RData"))
} else {
  load(file = paste0(maxent.reslts.dir, "/clim_model_results", numberBg, ".RData"))
}

# Evaluation
me_env.rslts$testAUC

# Projection
png(paste0(rslts.dir, "projection/clim_model_mortality_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_env.rslts$projection + theme(legend.position = "none")
p
dev.off()
print(paste0(rslts.dir, "projection/clim_model_mortality_projection.png"))

# variable contribution (for the prediction)
df <- data.frame("variable" = names(me_env.rslts$model@results[7:12, ]), 
                 "contribution" = me_env.rslts$model@results[7:12, ],
                 "permutation.importance" = me_env.rslts$model@results[13:18, ])
df <- df[order(df$contribution, decreasing = T), ]
df$variable <-  as.factor(df$variable)

p <- ggplot(df, aes(contribution, reorder(variable, order(contribution, decreasing = F)), fill = permutation.importance)) + geom_col() + theme_minimal() + 
  theme(axis.text.y = element_blank(), axis.title = element_blank()) + 
  scale_fill_gradient(low = "#fee0d2", high = "#de2d26")
leg <- get_legend(p)
p <- p + theme(legend.position = "none")

pdf(paste0(rslts.dir, "projection/clim_model_mortality_var_contr.pdf"), height = h, width = h)
p
dev.off()
print(paste0(rslts.dir, "projection/clim_model_mortality_var_contr.pdf"))
df

pdf(paste0(rslts.dir, "projection/clim_model_mortality_var_contr_legend.pdf"))
plot(leg)
dev.off()
print(paste0(rslts.dir, "projection/clim_model_mortality_var_contr_legend.pdf"))

# response curves
response(me_env.rslts$model)

# Projection sd crosvalidation (uncertainity)
png(paste0(rslts.dir, "projection/clim_model_mortality_crossvalidation_sd_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_env.rslts$pred_crosval_sd_projection + theme(legend.position = "left") 
p
dev.off()
print(paste0(rslts.dir, "projection/clim_model_mortality_crossvalidation_sd_projection.png"))


# Projection mean crosvalidation (uncertainity)
png(paste0(rslts.dir, "projection/clim_model_mortality_crossvalidation_mean_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_env.rslts$pred_crosval_mean_projection + theme(legend.position = "none") 
p
dev.off()
print(paste0(rslts.dir, "projection/clim_model_mortality_crossvalidation_mean_projection.png"))

remove(me_env.rslts)

### Using FT + biome + environmental data + HSM-based data ####

if(!file.exists(paste0(maxent.reslts.dir, "/hsm_clim_model_results", numberBg, ".RData")) | isTRUE(ForceRun)){
  me_hsm_env.rslts <- predMod(mod.var = c("p", "Tmax", "ai", "negative_HSM_count", "HSM_variance", "biome.synth", "ft"),
                          dta = hsm_md.df, nIter = 10, nbg = numberBg, maxent.args = c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                                                                                     "hinge=TRUE", "threshold=FALSE", "addSamplesToBackground=FALSE",
                                                                                     "responsecurves=TRUE", "writeplotdata","pictures=TRUE", "jackknife=TRUE"),
                          maxent.results.dir = paste0(maxent.reslts.dir, "/hsm_clim"))
  
  save(me_hsm_env.rslts, file = paste0(maxent.reslts.dir, "/hsm_clim_model_results", numberBg, ".RData"))
  print(paste0(maxent.reslts.dir, "/hsm_clim_model_results", numberBg, ".RData"))
} else {
  load(paste0(maxent.reslts.dir, "/hsm_clim_model_results", numberBg, ".RData")) 
}

# Evaluation
me_hsm_env.rslts$testAUC

# Projection
png(paste0(rslts.dir, "projection/hsm_clim_model_mortality_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_hsm_env.rslts$projection + theme(legend.position = "none") 
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_clim_model_mortality_projection.png"))


# response curves
response(me_hsm_env.rslts$model)

# variable contribution (for the prediction)
df <- data.frame("variable" = names(me_hsm_env.rslts$model@results[7:12, ]), 
                 "contribution" = me_hsm_env.rslts$model@results[7:12, ],
                 "permutation.importance" = me_hsm_env.rslts$model@results[13:18, ])
df <- df[order(df$contribution, decreasing = T), ]
df$variable <-  as.factor(df$variable)

p <- ggplot(df, aes(contribution, reorder(variable, order(contribution, decreasing = F)), fill = permutation.importance)) + geom_col() + theme_minimal() + 
  theme(axis.text.y = element_blank(), axis.title = element_blank()) + 
  scale_fill_gradient(low = "#fee0d2", high = "#de2d26")
leg <- get_legend(p)
p <- p + theme(legend.position = "none")

pdf(paste0(rslts.dir, "projection/hsm_clim_model_mortality_var_contr.pdf"), height = h, width = h)
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_clim_model_mortality_var_contr.pdf"))
df

pdf(paste0(rslts.dir, "projection/hsm_clim_model_mortality_var_contr_legend.pdf"))
plot(leg)
dev.off()
print(paste0(rslts.dir, "projection/hsm_clim_model_mortality_var_contr_legend.pdf"))


# Projection sd crosvalidation (uncertainity)
png(paste0(rslts.dir, "projection/hsm_clim_model_mortality_crossvalidation_sd_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_hsm_env.rslts$pred_crosval_sd_projection + theme(legend.position = "left") 
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_clim_model_mortality_crossvalidation_sd_projection.png"))


# Projection mean crosvalidation (uncertainity)
png(paste0(rslts.dir, "projection/hsm_clim_model_mortality_crossvalidation_mean_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_hsm_env.rslts$pred_crosval_mean_projection + theme(legend.position = "none") 
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_clim_model_mortality_crossvalidation_mean_projection.png"))

remove(me_hsm_env.rslts)

### Using HSM-based data ####

if(!file.exists(paste0(maxent.reslts.dir, "/hsm_model_results", numberBg, ".RData")) | isTRUE(ForceRun)){
  
  me_hsm.rslts <- predMod(mod.var = c("p", "HSM_mean", "negative_HSM_count",  "negative_HSM_prop", "HSM_variance", "biome.synth", "ft"),
                          dta = hsm_md.df, nIter = 100, nbg = numberBg, maxent.args = c("linear=FALSE", "product=FALSE", "quadratic=FALSE",
                                                                                        "hinge=TRUE", "threshold=FALSE", "addSamplesToBackground=FALSE",
                                                                                        "responsecurves=TRUE", "writeplotdata","pictures=TRUE", "jackknife=TRUE"),
                          maxent.results.dir = paste0(maxent.reslts.dir, "/hsm"))
  
  save(me_hsm.rslts, file = paste0(maxent.reslts.dir, "/hsm_model_results", numberBg, ".RData"))
  print(paste0(maxent.reslts.dir, "/hsm_model_results", numberBg, ".RData"))
} else {
  load(file = paste0(maxent.reslts.dir, "/hsm_model_results", numberBg, ".RData")) 
}

# Evaluation
me_hsm.rslts$testAUC

# legend (same for all maps)
mort_leg <- get_legend(me_hsm.rslts$projection + theme(legend.position = "bottom", legend.text = element_blank()))
pdf(paste0(rslts.dir, "projection/legend.pdf"), height = 0.3, width = 1.4)
plot(mort_leg)
dev.off()
print(paste0(rslts.dir, "projection/legend.pdf"))

# Projection
png(paste0(rslts.dir, "projection/hsm_model_mortality_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_hsm.rslts$projection + theme(legend.position = "none") 
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_model_mortality_projection.png"))

# variable contribution (for the prediction)
df <- data.frame("variable" = names(me_hsm.rslts$model@results[7:12, ]), 
                 "contribution" = me_hsm.rslts$model@results[7:12, ],
                 "permutation.importance" = me_hsm.rslts$model@results[13:18, ])
df <- df[order(df$contribution, decreasing = T), ]
df$variable <-  as.factor(df$variable)

p <- ggplot(df, aes(contribution, reorder(variable, order(contribution, decreasing = F)), fill = permutation.importance)) + geom_col() + theme_minimal() + 
  theme(axis.text.y = element_blank(), axis.title = element_blank()) + 
  scale_fill_gradient(low = "#fee0d2", high = "#de2d26")
leg <- get_legend(p)
p <- p + theme(legend.position = "none")

pdf(paste0(rslts.dir, "projection/hsm_model_mortality_var_contr.pdf"), height = h, width = h)
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_model_mortality_var_contr.pdf"))
df

pdf(paste0(rslts.dir, "projection/hsm_model_mortality_var_contr_legend.pdf"))
plot(leg)
dev.off()
print(paste0(rslts.dir, "projection/hsm_model_mortality_var_contr_legend.pdf"))

# response curves
response(me_hsm.rslts$model)

# Projection sd crosvalidation (uncertainity)
png(paste0(rslts.dir, "projection/hsm_model_mortality_crossvalidation_sd_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_hsm.rslts$pred_crosval_sd_projection + theme(legend.position = "left") 
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_model_mortality_crossvalidation_sd_projection.png"))


# Projection mean crosvalidation (uncertainity)
png(paste0(rslts.dir, "projection/hsm_model_mortality_crossvalidation_mean_projection.png"), height = h, width = w*0.75, units = "in", res = 1000)
p <- me_hsm.rslts$pred_crosval_mean_projection + theme(legend.position = "none") 
p
dev.off()
print(paste0(rslts.dir, "projection/hsm_model_mortality_crossvalidation_mean_projection.png"))

remove(me_hsm.rslts)

