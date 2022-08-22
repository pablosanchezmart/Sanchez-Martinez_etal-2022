#### FUNCTIONS ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

# P Sanchez-Marinez

#### CROSS-VALIDATION ------------------------------------------------------------------------------------------------------------------------------------------------------- ####

# Cross Validation function
RFCrosValFun <- function(name, data = rf.df, resp = c("P50", "MinWP_md"), pred = c(paste0("PC", 1:10)), NasP = 0.2, nIter = 10){
  rslts <- list() # final list of results
  R2rslts <- data.frame() # final dataframe of r2 results
  for(n in 1:nIter){
    # NA creation
    xTrue <- data %>% dplyr::select(c("Species", resp, pred)) %>% na.omit()
    na.df <- as.data.frame(xTrue[, resp])
    names(na.df) <- resp
    na.df <- prodNA(na.df, noNA = NasP)
    train.df <- xTrue
    train.df[, resp] <- na.df[, resp]
    
    # Observed values
    test.df <- data.frame()
    for(rVar in resp){
      test1 <- train.df %>% filter(is.na(train.df[, rVar]))
      test1 <- xTrue[which(xTrue$Species %in% test1$Species), c("Species", rVar)]
      names(test1)[names(test1) == rVar] <- "Observed_value"
      test1$Variable <- rVar
      test.df <- rbind(test.df, test1)
    }
    
    # Predictive model
    species <- train.df$Species
    train.df$Species <- NULL
    xTrue$Species <- NULL
    
    inteRslts <- data.frame() # each iteration results
    rfImpCv <- missForest(xmis = train.df, parallelize = "forests", maxiter = 50, ntree = 100, 
                          verbose = F, variablewise = T, decreasing = T,
                          xtrue = xTrue)
    rslts[[n]] <- rfImpCv
    # performance calculation
    ximp <- rfImpCv$ximp
    ximp$Species <- species
    for(rVar in resp){
      test1 <- test.df %>% filter(Variable == rVar) 
      pred1 <- ximp %>% filter(Species %in% test1$Species)
      names(pred1)[names(pred1) == rVar] <- "Predicted_value"
      pred1 <- pred1$Predicted_value
      r2.df <- cbind(test1, "predicted_value" = pred1)
      
      rslts1 <- cbind(Iteration = n, Variable = rVar,R2 = R2(pred = r2.df$predicted_value, obs = r2.df$Observed_value), 
                      RMSE = RMSE(pred = r2.df$predicted_value, obs = r2.df$Observed_value), 
                      MAE = MAE(pred = r2.df$predicted_value, obs = r2.df$Observed_value))
      inteRslts <- rbind(inteRslts, rslts1) # variable iteration results
      rslts[[n]][paste0(rVar, "_", "trainSpecies")] <- test.df %>% filter(Variable == rVar) %>% dplyr::select(Species)
    }
    R2rslts <- rbind(R2rslts, inteRslts) # Iteration results
  }
  # Results summarization
  R2rslts$R2 <- as.numeric(R2rslts$R2)
  R2rslts$RMSE <- as.numeric(R2rslts$RMSE)
  R2rslts$MAE <- as.numeric(R2rslts$MAE)
  
  meanR2 <- R2rslts %>% dplyr::select(Variable, R2, RMSE, MAE) %>% group_by(Variable) %>% summarize_all(., mean)
  sdR2 <- R2rslts %>% dplyr::select(Variable, R2, RMSE, MAE) %>% group_by(Variable) %>% summarize_all(., sd)  
  sdR2 <- dplyr::rename(sdR2, R2_sd = R2, RMSE_sd = RMSE, MAE_sd = MAE)
  sumR2 <- cbind("Model" = name, meanR2, sdR2)
  rslts[["R2"]] <- R2rslts
  rslts[["sumR2"]] <- sumR2
  return(rslts)
}


#### IMPUTATION ------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

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

## Imputation keeping real values when having them
dataImpFun <- function(obs.df, fit.df, var){
  # Complete dataset (filling with predictions conditional to the complete model)
  fit.df[, var] <- fit.df$fit
  impData <- merge(x = obs.df[, c("Species","genus", var)], y = fit.df[, c("Species", var)], by = "Species", all.x = F)
  impData[, "obs_fit"] <- if_else(is.na(impData[, paste0(var, ".x")]), "fit", "obs")
  impData <- comb.fun(df = impData, var = var)
  impData <- cbind(impData, fit.df[, c("lwr", "upr")])
  impData[impData$obs_fit == "obs", c("lwr", "upr")] <- NA
  return(impData)
}


# MAPPING ------------------------------------------------------------------------------------------------------------------------------------------------------------------ ####

# Rasterize and stack polygons
st_to_sfFun <- function(polygons.df = all_species_traits){
  sf <- sf::st_as_sf(polygons.df[1, ])
  for(n in 2:length(polygons.df$SOURCE_SHP)){
    print(n)
    sf[n, ] <- sf::st_as_sf(polygons.df[n, ])
  }
  object.size(sf)
  return(sf)
}

# Plot map

plotDataFun <- function(mean.r = NULL, max.r = NULL, min.r = NULL, xText = F){
  
  
  # biome.df <- read.csv("outputs/biomes_ecoregions/biomes_dataframe.csv", header = T)
  # biome.df$biome.synth <- as.factor(biome.df$biome.synth)
  # biome.df$biome.synth_n <- as.numeric(as.factor(biome.df$biome.synth))
  
  # biomes.r <- rasterFromXYZ(biome.df[, c("x", "y", "biome.synth_n")], crs = crs(mean.r))
  # biomes.r <- projectRaster(biomes.r, mean.r, method = "ngb")
  # mean.r <- stack(mean.r, biomes.r)
  mean.df <- data.frame(as(mean.r, "SpatialPixelsDataFrame"))
  names(mean.df)[1] <- "layer"
  
  # mean.df$biome.synth <- dplyr::recode(mean.df$biome.synth_n, "1" = "Boreal", "2" = "Desert and xeric", "3" = "Mediterranean", "4" = "Others", 
                # "5" = "Temperate", "6" = "Tropical and subtropical dry", "7" = "Tropical and subtropical moist")
  # mean.df$biome.synth <- as.factor(mean.df$biome.synth)
  
  
  if(!is.null(min.r)){
    min.df <- data.frame(as(min.r, "SpatialPixelsDataFrame"))
    names(min.df)[1] <- "layer"
  # Minimum by latitude
  min.lat <- min.df %>% group_by(y)  %>% dplyr::select(layer, y)
  min.lat <-aggregate(min.lat, by = list(min.lat$y), FUN = min)  %>% dplyr::select(layer, y)
  }
  if(!is.null(max.r)){
    max.df <- data.frame(as(max.r, "SpatialPixelsDataFrame"))
    names(max.df)[1] <- "layer"
  # Maximum by latitude  
  max.lat <- max.df %>% group_by(y)  %>% dplyr::select(layer, y)
  max.lat <-aggregate(max.lat, by = list(max.lat$y), FUN = max)  %>% dplyr::select(layer, y)
  }
  
  #plot
  #mean values
  p <- ggplot(mean.df) + geom_point(aes(x = y, y = layer), color = "#bdbdbd", alpha = 0.2, size = 0.1) +  #, color = biome.synth
      scale_color_brewer(palette = "Spectral")  + # color = "#bdbdbd" # when not using biomes
    geom_smooth(aes(x = y, y = layer), method = "gam", color = "#636363", size = 0.5) 
  
  if(!is.null(max.r)){
  # max
  p <- p + geom_point(data = max.lat, aes(x = y, y = layer), color = "#9ecae1", shape = 17, size = 0.5, alpha = 0.9) +   
    geom_smooth(data = max.lat, aes(x = y, y = layer), method = "gam", color = "#3182bd", size = 0.5)
  }
  if(!is.null(min.r)){
  # min
  p <- p + geom_point(data = min.lat, aes(x = y, y = layer), color = "#fc9272", shape = 17, size = 0.5, alpha = 0.9) +   
    geom_smooth(data = min.lat, aes(x = y, y = layer), method = "gam", color = "#de2d26", size = 0.5)
  }
  p <- p + theme_light() + theme(axis.title = element_blank(),
                                 legend.title = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 panel.grid.minor.x = element_blank(),
                                 panel.grid.minor.y = element_blank(),
                                 axis.line = element_line(colour = "black", 
                                                          size = 1, linetype = "solid"),
                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  if(isFALSE(xText)){
    p <- p + theme(axis.text.x = element_blank()) 
  }
  p <- p +  coord_flip()
  # print(p)
  return(p)
}

plotMapFun <- function(mean.r = p50_mean, var = "layer"){

  mean.r <- clamp(mean.r, lower = as.numeric(quantile(mean.r, probs = 0.05)), upper = as.numeric(quantile(mean.r, probs = 0.95)), useValues = T)
  
  df <- data.frame(as(mean.r, "SpatialPixelsDataFrame"))
  cols <- brewer.pal(10, "Spectral")
  p <- ggplot(df) + geom_raster(aes(x=x, y=y, fill=df[, var]))
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
  # p
  return(p)
}

plotDataMapFun <- function(mean.r = NULL, max.r = NULL, min.r = NULL, var = "layer", zoomPnts = F){
  
  p <- plotMapFun(mean.r, var)
  legend <- get_legend(p)
  if(isTRUE(zoomPnts)){
    p <- p + geom_point(data = zoomPoints, aes(x = x, y = y), shape = 3) 
  }
  p <- p + theme(legend.position = "none")
  p2 <- plotDataFun(mean.r, max.r, min.r)
  biome_legend <- get_legend(p2)
  p2 <- p2 + theme(legend.position = "none")
  
  p.grob <- set_panel_size(p, width = unit(w-h/2, "in"), height = unit(h, "in"))
  p2.grob <- set_panel_size(p2, width = unit(h/2, "in"), height = unit(h, "in"))
  
  grid.arrange(grobs = list(p.grob, p2.grob), nrow = 1, layout_matrix = rbind(c(1, 1, 1, 2)))
  return(list(legend, biome_legend))
}

plotLogMapFun <- function(data, var = "layer"){
  df <- data.frame(as(data, "SpatialPixelsDataFrame"))
  cols <- rev(brewer.pal(10, "Spectral"))
  p <- ggplot(df) + geom_raster(aes(x=x, y=y, fill=df[, var]))
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
                                  legend.position = c(0.21, 0.19),
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  print(p)
}

discPlotMapFun <- function(data = p50_count){
  reclass_df <- matrix(c(1, 10, 1,
                         10, 50, 2,
                         50, 100, 3,
                         100, 500, 4,
                         500, 1000, 5,
                         1000, 2000, 6,
                         2000, 3000, 7,
                         3000, 6000, 8), ncol = 3, byrow = T
  )
  recl.r <- reclassify(data, reclass_df)
  df <- data.frame(as(recl.r, "SpatialPixelsDataFrame"))
  df$layer <- as.factor(df$layer)
  df$layer <- dplyr::recode(df$layer, "1" = "1-10", "2" = "10-50", "3" = "50-100", "4" = "100-500", 
                     "5" = "500-1,000", "6" = "1,000-2,000", "7" = "2,000-3,000", "8" = "3,000-6,000")
  p <- ggplot(df) + geom_raster(aes(x=x, y=y, fill=layer))
  p <- p + scale_fill_brewer(palette = "Spectral", direction = -1)
  p <- p + coord_equal()
  p <- p + theme(panel.grid = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 legend.title = element_blank(),
  )
  print(p)
}

plotMapPointsFun <- function(mean.r = p50_mean, var = "layer", points = md[, c("x", "y")]){
  df <- data.frame(as(mean.r, "SpatialPixelsDataFrame"))
  cols <- brewer.pal(10, "Spectral")
  p <- ggplot(df) + geom_raster(aes(x=x, y=y, fill=df[, var]))
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
  p <- p + ggplot2::geom_point(data = points, aes(x = x, y = y), shape = 3)
  print(p)
  return(p)
}

# Plot histogram (zoom points)

pointsZoomHistPlot <- function(pointVal = pointsVal, id = "a", var = "P50"){
  
  pointVal <- pointVal %>% filter(point_id == id)
  pointVal <- as.data.frame(pointVal)
  p <- ggplot(pointVal) + geom_histogram(aes(x = pointVal[, var]))
  p <- p + theme_minimal() + theme(axis.title = element_blank(),
                                   legend.title = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.grid.minor.x = element_blank(),
                                   panel.grid.minor.y = element_blank(),
                                   axis.line = element_line(colour = "black", linetype = "solid"),
                                   plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  print(p)
  print(length(pointVal[,1]))
}

# Plot phylo

plotPhyloFun <- function(phylo, traits, brwPal = "BrBG", panelSpace = 4, FontsizeFactor = 0.01,
                         problematicSpp = "none"){
  
  # Some species are not in their place sin the phylogeny (problems within fagales and rosales)
  if(!isTRUE(problematicSpp == "none")){
    traits.clean <- traits[!traits$sp %in% problematicSpp, ]
    phylo.clean <- drop.tip(phylo, problematicSpp)
  } else{
    traits.clean <- traits
    phylo.clean <- phylo
  }

  colsOrd <- colorRampPalette(brewer.pal(10, "BrBG"))(length(unique(traits$order)))
  p <- ggtree(phylo.clean) 
  # name orders
  for(o in unique(traits.clean$order)){
    spp <- rownames(traits.clean[which(traits.clean$order == o), ])
    if(length(spp) == 1){                                               # Some orders have only one species, so MRCA cannot be found
      nd <- which(phylo.clean$tip.label == spp)
    } else{
      nd <- findMRCA(phylo.clean, spp, type = "node") 
    }
    p <- p + geom_cladelabel(node= nd, label= o,
                             fontsize = length(rownames(traits.clean[traits.clean$order == o, ])) * FontsizeFactor + 1) +
      coord_cartesian(clip = 'off')  + theme_tree2(plot.margin= ggtree::margin(6, 6, 6, 6))
  }
  p <- p   %<+% traits.clean + geom_tippoint(aes(color = order))
  
  # Barplots
  # Ymin
  d2 <- data.frame(id = row.names(traits.clean), val = traits.clean$MinWP_md, mortality = as.factor(traits.clean$mortality_binomial))
  p2 <- facet_plot(p, panel = "MinWP_md", data = d2, geom = geom_barh, 
                      mapping =  aes(x = val, fill = mortality), 
                      stat = "identity")
  # P50
  d3 <- data.frame(id = row.names(traits.clean), val = traits.clean$P50, mortality = as.factor(traits.clean$mortality_binomial))
  p3 <- p2 + geom_facet(panel = "P50", data = d3, geom = geom_barh,
                      mapping =  aes(x = val, fill = mortality),  # fill = order
                      stat = "identity")
  # HSM
  d4 <- data.frame(id = row.names(traits.clean), val = traits.clean$HSM, mortality = as.factor(traits.clean$mortality_binomial))
  p4 <- p3 + geom_facet(panel = "HSM", data = d4, geom = ggstance::geom_barh,
                      mapping =  aes(x = val, fill = mortality), stat = "identity")
  p4 <- p4  + scale_fill_manual(values = c("#4daf4a", "#e41a1c")) + scale_color_manual(values = colsOrd) +
    theme(panel.spacing = unit(panelSpace, "lines")) +
    # theme_transparent()  +
    theme(legend.position = "none")
  print(p4)
}

# To look closer if needed
plotPhyloFamFun <- function(phylo, traits, brwPal = "BrBG", panelSpace = 4, FontsizeFactor = 0.01,
                            orderName = "Pinales"){
  
  traits.clean <- traits %>% filter(order == orderName)
  phylo.clean <- keep.tip(phylo, which(phylo$tip.label %in% traits.clean[, 1]))
  
  colsOrd <- colorRampPalette(brewer.pal(10, "BrBG"))(length(unique(traits$order)))
  p <- ggtree(phylo.clean) 
  # name orders
  for(o in unique(traits.clean$family)){
    spp <- rownames(traits.clean[which(traits.clean$family == o), ])
    if(length(spp) == 1){                                               # Some orders have only one species, so MRCA cannot be found
      nd <- which(phylo.clean$tip.label == spp)
    } else{
      nd <- findMRCA(phylo.clean, spp, type = "node") 
    }
    p <- p + geom_cladelabel(node= nd, label= o,
                             fontsize = length(rownames(traits.clean[traits.clean$family == o, ])) * FontsizeFactor + 1) +
      coord_cartesian(clip = 'off')  + theme_tree2(plot.margin= ggtree::margin(6, 6, 6, 6))
  }
  p <- p   %<+% traits.clean + geom_tippoint(aes(color = family))
  
  # Barplots
  # Ymin
  d2 <- data.frame(id = row.names(traits.clean), val = traits.clean$MinWP_md, mortality = as.factor(traits.clean$mortality_binomial))
  p2 <- facet_plot(p, panel = "MinWP_md", data = d2, geom = geom_barh, 
                   mapping =  aes(x = val, fill = mortality), 
                   stat = "identity")
  # P50
  d3 <- data.frame(id = row.names(traits.clean), val = traits.clean$P50, mortality = as.factor(traits.clean$mortality_binomial))
  p3 <- p2 + geom_facet(panel = "P50", data = d3, geom = geom_barh,
                        mapping =  aes(x = val, fill = mortality),  # fill = order
                        stat = "identity")
  # HSM
  d4 <- data.frame(id = row.names(traits.clean), val = traits.clean$HSM, mortality = as.factor(traits.clean$mortality_binomial))
  p4 <- p3 + geom_facet(panel = "HSM", data = d4, geom = ggstance::geom_barh,
                        mapping =  aes(x = val, fill = mortality), stat = "identity")
  p4 <- p4  + scale_fill_manual(values = c("#4daf4a", "#e41a1c")) + scale_color_manual(values = colsOrd) +
    theme(panel.spacing = unit(panelSpace, "lines")) +
    # theme_transparent()  +
    theme(legend.position = "none")
  print(p4)
}

plotMortPhyloFun <- function(phylo, traits, brwPal = "BrBG", panelSpace = 4, FontsizeFactor = 0.01,
                         problematicSpp = "none"){
  
  # Some species are not in their place sin the phylogeny (problems within fagales and rosales)
  if(!isTRUE(problematicSpp == "none")){
    traits.clean <- traits[!traits$sp %in% problematicSpp, ]
    phylo.clean <- drop.tip(phylo, problematicSpp)
  } else{
    traits.clean <- traits
    phylo.clean <- phylo
  }
  
  colsOrd <- colorRampPalette(brewer.pal(10, "BrBG"))(length(unique(traits$order)))
  p <- ggtree(phylo.clean) 
  # name orders
  for(o in unique(traits.clean$order)){
    spp <- rownames(traits.clean[traits.clean$order == o, ])
    if(length(spp) == 1){                                               # Some orders have only one species, so MRCA cannot be found
      nd <- which(phylo.clean$tip.label == spp)
    } else{
      nd <- findMRCA(phylo.clean, spp, type = "node") 
    }
    p <- p + geom_cladelabel(node= nd, label= o,
                             fontsize = length(rownames(traits.clean[traits.clean$order == o, ])) * FontsizeFactor + 1) +
      coord_cartesian(clip = 'off')  + theme_tree2(plot.margin= ggtree::margin(6, 6, 6, 6))
  }
  p <- p   %<+% traits.clean + geom_tippoint(aes(color = order))
  
  # Barplots
  # Ymin
  d2 <- data.frame(id = row.names(traits.clean), val = traits.clean$MinWP_md, mortality = as.factor(traits.clean$mortality_binomial))
  p2 <- facet_plot(p, panel = "MinWP_md", data = d2, geom = geom_barh, 
                   mapping =  aes(x = val, fill = mortality), 
                   stat = "identity")
  # points (option)
  # p2 <- facet_plot(p, panel = "MinWP_md", data = d2, geom = geom_point, mapping =  aes(x = MinWP_md, shape = factor(mortality)))
  
  # P50
  d3 <- data.frame(id = row.names(traits.clean), val = traits.clean$P50, mortality = as.factor(traits.clean$mortality_binomial))
  p3 <- p2 + geom_facet(panel = "P50", data = d3, geom = geom_barh,
                        mapping =  aes(x = val, fill = mortality),  # fill = order
                        stat = "identity")
  # HSM
  d4 <- data.frame(id = row.names(traits.clean), val = traits.clean$HSM, mortality = as.factor(traits.clean$mortality_binomial))
  p4 <- p3 + geom_facet(panel = "HSM", data = d4, geom = ggstance::geom_barh,
                        mapping =  aes(x = val, fill = mortality), stat = "identity")
  p4 <- p4  + scale_fill_manual(values = c("#e41a1c")) + scale_color_manual(values = colsOrd) +
    theme(panel.spacing = unit(panelSpace, "lines")) +
    # theme_transparent()  +
    theme(legend.position = "none")
  print(p4)
}

plotPhyloSDFun <- function(phylo, traits, brwPal = "BrBG", panelSpace = 4, FontsizeFactor = 0.01,
                         problematicSpp = "none"){
  
  # Some species are not in their place sin the phylogeny (problems within fagales and rosales)
  if(!isTRUE(problematicSpp == "none")){
    traits.clean <- traits[!traits$sp %in% problematicSpp, ]
    phylo.clean <- drop.tip(phylo, problematicSpp)
  } else{
    traits.clean <- traits
    phylo.clean <- phylo
  }
  
  colsOrd <- colorRampPalette(brewer.pal(10, "BrBG"))(length(unique(traits$order)))
  p <- ggtree(phylo.clean) 
  # name orders
  for(o in unique(traits.clean$order)){
    spp <- rownames(traits.clean[traits.clean$order == o, ])
    if(length(spp) == 1){                                               # Some orders have only one species, so MRCA cannot be found
      nd <- which(phylo.clean$tip.label == spp)
    } else{
      nd <- findMRCA(phylo.clean, spp, type = "node") 
    }
    p <- p + geom_cladelabel(node= nd, label= o,
                             fontsize = length(rownames(traits.clean[traits.clean$order == o, ])) * FontsizeFactor + 1) +
      coord_cartesian(clip = 'off')  + theme_tree2(plot.margin= ggtree::margin(6, 6, 6, 6))
  }
  p <- p   %<+% traits.clean + geom_tippoint(aes(color = order))
  
  # Barplots
  # Ymin
  d2 <- data.frame(id = row.names(traits.clean), val = traits.clean$MinWP_md_sd, mortality = as.factor(traits.clean$mortality_binomial))
  p2 <- facet_plot(p, panel = "MinWP_md", data = d2, geom = geom_barh, 
                   mapping =  aes(x = val, fill = mortality), 
                   stat = "identity")
  # points (option)
  # p2 <- facet_plot(p, panel = "MinWP_md", data = d2, geom = geom_point, mapping =  aes(x = MinWP_md, shape = factor(mortality)))
  
  # P50
  d3 <- data.frame(id = row.names(traits.clean), val = traits.clean$P50_sd, mortality = as.factor(traits.clean$mortality_binomial))
  p3 <- p2 + geom_facet(panel = "P50", data = d3, geom = geom_barh,
                        mapping =  aes(x = val, fill = mortality),  # fill = order
                        stat = "identity")
  # HSM
  d4 <- data.frame(id = row.names(traits.clean), val = traits.clean$HSM_sd, mortality = as.factor(traits.clean$mortality_binomial))
  p4 <- p3 + geom_facet(panel = "HSM", data = d4, geom = ggstance::geom_barh,
                        mapping =  aes(x = val, fill = mortality), stat = "identity")
  p4 <- p4  + scale_fill_manual(values = c("#4daf4a", "#e41a1c")) + scale_color_manual(values = colsOrd) +
    theme(panel.spacing = unit(panelSpace, "lines")) +
    # theme_transparent()  +
    theme(legend.position = "none")
  print(p4)
}


