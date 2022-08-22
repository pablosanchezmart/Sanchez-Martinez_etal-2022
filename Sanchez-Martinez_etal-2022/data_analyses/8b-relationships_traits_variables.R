#### CORRLATION BETWEEN MINIMUM WATER POTENTIALS ------------------------------------------------------------------- ####

remove(list = ls())
gc()
# install.packages("ggplot2", dependencies = T)

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

#### Pmin plant vs. Pmin soil ------------------------------------------------------------------------------------------------------- ####

cor.df <- ht[complete.cases(ht[, c("log_negMinWP_md", "MinWP_soil_SX")]) ,]

length(cor.df$Species)
cor.df$log_nMinWP_soil_SX <- log(-cor.df$MinWP_soil_SX)

summary(as.factor(cor.df$group)) # 46 gymn, 333 angio
mod_1 <- lm(cor.df$log_negMinWP_md ~ cor.df$log_nMinWP_soil_SX)
summary(mod_1)

mod_2 <- lm(log_negMinWP_md ~ log_nMinWP_soil_SX * group, data = cor.df)
summary(mod_2)

# Relationship (log-log)
png(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_regression.png"),height = 7.87/2, width = 7.87/2, units = "in", res = 1000)
p <- ggplot(cor.df, aes(x = log_nMinWP_soil_SX, y = log_negMinWP_md, color = group)) + geom_point(size = 1.5, alpha = 0.7)
p <- p + geom_smooth(method = "lm",  alpha = 0.3)
p <- p + geom_smooth(method = "lm", color = "black", alpha = 0.3)
p <- p + theme(legend.title = element_blank())
leg <- get_legend(p)
p <- p + theme(legend.position = "none",
               axis.title = element_blank()) 
p
dev.off()
print(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_regression.png"))

png(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_leg.png"),height = 0.5, width = 1.2, units = "in", res = 1000)
plot(leg)
dev.off()
print(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_leg.png"))

# Points
png(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_MinWP_points.png"), height = 7.87/4, width = 7.87/4, units = "in", res = 1000)
p <- ggplot(cor.df) + geom_point(aes(MinWP_soil_SX, MinWP_md, shape = group,  color = group), size = 1, alpha = 0.7) + theme(legend.position = "none", axis.title = element_blank())
p
dev.off()
print(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_MinWP_points.png"))


#### Pmin soil vs. VPDMax ----------------------------------------------------------------------------------------------------------- ####

cor.df <- ht[complete.cases(ht[, c("log_negMinWP_md", "VPDMax")]) ,]

length(cor.df$Species)

summary(as.factor(cor.df$group)) # 46 gymn, 333 angio
mod_1 <- lm(cor.df$log_negMinWP_md ~ cor.df$VPDMax)
summary(mod_1)

mod_2 <- lm(log_negMinWP_md ~ VPDMax * group, data = cor.df)
summary(mod_2)

# Relationship (log-log)
png(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_VPDMax_regression.png"),height = 7.87/2, width = 7.87/2, units = "in", res = 1000)
p <- ggplot(cor.df, aes(x = VPDMax, y = log_negMinWP_md, color = group)) + geom_point(size = 1.5, alpha = 0.7)
p <- p + geom_smooth(method = "lm",  alpha = 0.3)
p <- p + geom_smooth(method = "lm", color = "black", alpha = 0.3)
leg <- get_legend(p)
p <- p + theme(legend.position = "none",
               axis.title = element_blank()) 
p
dev.off()
print(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_VPDMax_regression.png"))

# Points
png(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_MinWP_VPDMax_points.png"), height = 7.87/4, width = 7.87/4, units = "in", res = 1000)
p <- ggplot(cor.df) + geom_point(aes(VPDMax, MinWP_md, shape = group,  color = group), size = 1, alpha = 0.7) + theme(legend.position = "none", axis.title = element_blank())
p
dev.off()
print(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_MinWP_VPDMax_points.png"))


#### Pmin soil vs. VPDmax * Pmin soil ----------------------------------------------------------------------------------------------- ####


cor.df <- ht[complete.cases(ht[, c("log_negMinWP_md", "MinWP_soil_SX", "VPDMax")]) ,]

length(cor.df$Species)
cor.df$log_nMinWP_soil_SX <- log(-cor.df$MinWP_soil_SX)

summary(as.factor(cor.df$group)) # 46 gymn, 333 angio

mod_1 <- lm(cor.df$log_negMinWP_md ~ cor.df$log_nMinWP_soil_SX * cor.df$VPDMax)
summary(mod_1)

mod_2 <- lm(log_negMinWP_md ~ log_nMinWP_soil_SX * VPDMax * group, data = cor.df)
summary(mod_2)

#### Pmin soil vs. P50 -------------------------------------------------------------------------------------------------------------- ####

cor.df <- ht[complete.cases(ht[, c("log_negMinWP_md", "log_negP50")]) ,]

length(cor.df$Species)

summary(as.factor(cor.df$group)) # 46 gymn, 333 angio
mod_1 <- lm(cor.df$log_negP50 ~ cor.df$log_negMinWP_md)
summary(mod_1)

mod_2 <- lm(cor.df$log_negP50 ~ cor.df$log_negMinWP_md * group, data = cor.df)
summary(mod_2)

# Relationship (log-log)
png(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_P50_regression.png"),height = 7.87/2, width = 7.87/2, units = "in", res = 1000)
p <- ggplot(cor.df, aes(x = log_negMinWP_md, y = log_negP50, color = group)) + geom_point(size = 1.5, alpha = 0.7)
p <- p + geom_smooth(method = "lm",  alpha = 0.3)
p <- p + geom_smooth(method = "lm", color = "black", alpha = 0.3)
leg <- get_legend(p)
p <- p + theme(legend.position = "none",
               axis.title = element_blank()) 
p
dev.off()
print(paste0(rslts.dir, "figures/relationship_MinWP/relationship_MinWP_P50_regression.png"))

# Points
png(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_MinWP_P50_points.png"), height = 7.87/4, width = 7.87/4, units = "in", res = 1000)
p <- ggplot(cor.df) + geom_point(aes(MinWP_md, P50, shape = group,  color = group), size = 1, alpha = 0.7) + theme(legend.position = "none", axis.title = element_blank())
p
dev.off()
print(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_MinWP_P50_points.png"))

#### P50 vs. P88 -------------------------------------------------------------------------------------------------------------------- ####

cor.df <- ht[complete.cases(ht[, c("log_negP50", "P88")]) ,]

length(cor.df$Species)
cor.df$log_negP88 <- log(-cor.df$P88)

summary(as.factor(cor.df$group)) # 46 gymn, 333 angio
mod_1 <- lm(cor.df$log_negP50 ~ cor.df$log_negP88)
summary(mod_1)

mod_2 <- lm(cor.df$log_negP50 ~ cor.df$log_negP88 * group, data = cor.df)
summary(mod_2)

# Relationship (log-log)
png(paste0(rslts.dir, "figures/relationship_MinWP/relationship_P50_P88_regression.png"),height = 7.87/2, width = 7.87/2, units = "in", res = 1000)
p <- ggplot(cor.df, aes(x = log_negP88, y = log_negP50, color = group)) + geom_point(size = 1.5, alpha = 0.7)
p <- p + geom_smooth(method = "lm",  alpha = 0.3)
p <- p + geom_smooth(method = "lm", color = "black", alpha = 0.3)
leg <- get_legend(p)
p <- p + theme(legend.position = "none",
               axis.title = element_blank()) 
p
dev.off()
print(paste0(rslts.dir, "figures/relationship_MinWP/relationship_P50_P88_regression.png"))

# Points
png(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_P50_P88_points.png"), height = 7.87/4, width = 7.87/4, units = "in", res = 1000)
p <- ggplot(cor.df) + geom_point(aes(P50, log_negP88, shape = group,  color = group), size = 1, alpha = 0.7) + theme(legend.position = "none", axis.title = element_blank())
p
dev.off()
print(paste0(rslts.dir, "/figures/relationship_MinWP/relationship_P50_P88_points.png"))


#### GEOGRAPHICAL COMPARISON BETWEEN Pmins ---------------------------------------------------------------------------------- ####

Pmins.r <- raster("data/species_data/soil_data/soil_hydraulics/minimum_WP/MinWP_soil_SX.tif")

Pmin.r <- raster(paste0(output.dir, "Pmin/Pmin_mean.tif"))
Pmins.r <- projectRaster(Pmins.r, Pmin.r)

Pmin.st <- stack(Pmins.r, Pmin.r)

plot(Pmin.st)

Pmin.df <- as.data.frame(raster::sampleRandom(Pmin.st, 100000))
# densityplot(Pmin.df$Pmin_mean)

cor(Pmin.df)

mod <- lm(log(-Pmin_mean) ~ log(-MinWP_soil_SX), data = Pmin.df)
summary(mod)

p <- ggplot(Pmin.df, aes(x = MinWP_soil_SX, y = Pmin_mean)) + geom_point(alpha = 0.3, size = 0.1) + geom_smooth(method = "lm", size = 0.5) +
  theme(axis.title = element_blank())
p
png(paste0(rslts.dir, "figures/relationship_MinWP/spatial_relationship_Pmin_Pminsoil.png"), height = h, width = h, units = "in", res = 1000)
p
dev.off()
print("results/figures/comparison_previous_results/spatial_relationship_Pmin_Pminsoil.png")