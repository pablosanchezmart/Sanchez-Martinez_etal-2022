##### PREVIOUS RESULTS  COMPARISON ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #####

remove(list = ls())
gc()

source("code/manuscript/RF/1b-init.R")
source("code/manuscript/RF/2b-functions.R")

library(ncdf4)
library(raster)

rslts.dir <- paste0(rslts.dir, "figures/comparison_previous_results/")

p50_mean <- raster(paste0(output.dir, "/P50/P50_mean.tif"))
HSM_mean <- raster(paste0(output.dir, "/HSM/HSM_mean.tif"))

p50_range <- raster(paste0(output.dir, "/P50/P50_FRic.tif"))
HSM_range <- raster(paste0(output.dir, "/HSM/HSM_FRic.tif"))

h <- 2.61
w <- 7.87

#### LIU DATA --------------------------------------------------------------------------------- ####

p50_mean_liu <- raster("data/Liu_hydraulics_projection/MDF_P50.nc", varname = "P50_m")
plot(p50_mean_liu)
p50_mean_liu <- projectRaster(p50_mean_liu, p50_mean)

p50.st <- stack(p50_mean, p50_mean_liu)
names(p50.st) <- c("P50_mean", "P50_mean_liu")

plot(p50.st)

p50_dif <- p50.st$P50_mean - p50.st$P50_mean_liu
plot(p50_dif)

p50.df <- as.data.frame(raster::sampleRandom(p50.st, 100000))


cor(p50.df)

mod <- lm(log(-P50_mean) ~ log(-P50_mean_liu), data = p50.df)
summary(mod)

p <- ggplot(p50.df, aes(x = log(-P50_mean_liu), y = log(-P50_mean))) + geom_point() + geom_smooth(method = "lm")
p

p50_corr


#### TRUGMAN DATA --------------------------------------------------------------------------------- ####

#### mean P50 ####

nc_data <- nc_open('data/Trugman_hydraulics_projection/CWM_P50_025Deg.nc')


lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)

length(lon)
length(lat)

p50_mean_trug <- ncvar_get(nc_data, "CWM_P50")
dim(p50_mean_trug)

fillvalue <- ncatt_get(nc_data, "CWM_P50", "_FillValue")
fillvalue
p50_mean_trug[p50_mean_trug == fillvalue$value] <- NA

r <- raster(t(p50_mean_trug), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

r_fl <- raster::flip(r, direction = 1)
r_fl <- t(r_fl)

plot(r)
plot(r_fl)
extent(r_fl) <- extent(r)
plot(r_fl)

p50_mean_trug_proj <- projectRaster(r_fl, p50_mean)

p50.st <- stack(p50_mean, p50_mean_trug_proj)
names(p50.st) <- c("P50_mean", "P50_mean_trug")

plot(p50.st)
p50_dif <- p50.st$P50_mean - p50.st$P50_mean_trug
plot(p50_dif)

p50.df <- as.data.frame(raster::sampleRandom(p50.st, 100000))

cor(p50.df)
cor(log(-p50.df))

mod <- lm(log(-P50_mean) ~ log(-P50_mean_trug), data = p50.df)
summary(mod)


mod <- lm(log(-P50_mean) ~ log(-P50_mean_trug), data = p50.df)
summary(mod)


#### mean P50 Plots ####

### Scatterplot

p <- ggplot(p50.df, aes(x = log(-P50_mean_trug), y = log(-P50_mean))) + geom_point(alpha = 0.3, size = 0.1) + geom_smooth(method = "lm", size = 0.5) +
  theme(axis.title = element_blank())
p
png(paste0(rslts.dir, "trugman_scatterplot_log_abs_P50.png"), height = h, width = h, units = "in", res = 1000)
p
dev.off()
print(paste0(rslts.dir, "trugman_scatterplot_log_abs_P50.png"))


p <- ggplot(p50.df, aes(x = P50_mean_trug, y = P50_mean)) + geom_point(alpha = 0.3, size = 0.1)
p
png(paste0(rslts.dir, "trugman_scatterplot_P50.png"), height = h, width = h, units = "in", res = 1000)
p
dev.off()
print(paste0(rslts.dir, "trugman_scatterplot_P50.png"))


### Rasters

p50.spdf <- data.frame(as(p50.st, "SpatialPixelsDataFrame"))

p50.spdf <- p50.spdf[complete.cases(p50.spdf$P50_mean, p50.spdf$P50_mean_trug), ]

p50.spdf$dif <- p50.spdf$P50_mean - p50.spdf$P50_mean_trug

cols <- brewer.pal(10, "Spectral")

p <- ggplot(p50.spdf, aes(x = x, y = y, fill = P50_mean)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "p50_mean_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "p50_mean_raster.png"))


p <- ggplot(p50.spdf, aes(x = x, y = y, fill = P50_mean_trug)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "p50_mean_trug_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "p50_mean_trug_raster.png"))


p <- ggplot(p50.spdf, aes(x = x, y = y, fill = dif)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "p50_mean_difference_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "p50_mean_difference_raster.png"))


#### mean HSM ####

nc_data <- nc_open('data/Trugman_hydraulics_projection/CWM_HSM_025Deg.nc')


lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)

length(lon)
length(lat)

HSM_mean_trug <- ncvar_get(nc_data, "CWM_HSM")
dim(HSM_mean_trug)

fillvalue <- ncatt_get(nc_data, "CWM_HSM", "_FillValue")
fillvalue
HSM_mean_trug[HSM_mean_trug == fillvalue$value] <- NA

r <- raster(t(HSM_mean_trug), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

r_fl <- raster::flip(r, direction = 1)
r_fl <- t(r_fl)

plot(r)
plot(r_fl)
extent(r_fl) <- extent(r)
plot(r_fl)

HSM_mean_trug_proj <- projectRaster(r_fl, HSM_mean)

HSM.st <- stack(HSM_mean, HSM_mean_trug_proj)
names(HSM.st) <- c("HSM_mean", "HSM_mean_trug")

plot(HSM.st)
HSM_dif <- HSM.st$HSM_mean - HSM.st$HSM_mean_trug
plot(HSM.st)

HSM.df <- as.data.frame(raster::sampleRandom(HSM.st, 100000))

cor(HSM.df)

mod <- lm(HSM_mean ~ HSM_mean_trug, data = HSM.df)
summary(mod)

#### mean HSM Plots ####

### Scatterplot

p <- ggplot(HSM.df, aes(x = HSM_mean_trug, y = HSM_mean)) + geom_point(alpha = 0.3, size = 0.1) + geom_smooth(method = "lm", size = 0.5) +
  theme(axis.title = element_blank())
p
png(paste0(rslts.dir, "trugman_scatterplot_HSM.png"), height = h, width = h, units = "in", res = 1000)
p
dev.off()
print(paste0(rslts.dir, "trugman_scatterplot_HSM.png"))


### Rasters

HSM.spdf <- data.frame(as(HSM.st, "SpatialPixelsDataFrame"))

HSM.spdf <- HSM.spdf[complete.cases(HSM.spdf$HSM_mean, HSM.spdf$HSM_mean_trug), ]

HSM.spdf$dif <- HSM.spdf$HSM_mean - HSM.spdf$HSM_mean_trug

cols <- brewer.pal(10, "Spectral")

p <- ggplot(HSM.spdf, aes(x = x, y = y, fill = HSM_mean)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "HSM_mean_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "HSM_mean_raster.png"))


p <- ggplot(HSM.spdf, aes(x = x, y = y, fill = HSM_mean_trug)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "HSM_mean_trug_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "HSM_mean_trug_raster.png"))


p <- ggplot(HSM.spdf, aes(x = x, y = y, fill = dif)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "HSM_mean_difference_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "HSM_mean_difference_raster.png"))


#### range P50 ####

nc_data <- nc_open('data/Trugman_hydraulics_projection/CWM_P50_025Deg.nc')

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)

length(lon)
length(lat)

p50_range_trug <- ncvar_get(nc_data, "max_range_within_grid")
dim(p50_range_trug)

fillvalue <- ncatt_get(nc_data, "max_range_within_grid", "_FillValue")
fillvalue
p50_range_trug[p50_range_trug == fillvalue$value] <- NA

r <- raster(t(p50_range_trug), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

r_fl <- raster::flip(r, direction = 1)
r_fl <- t(r_fl)

plot(r)
plot(r_fl)
extent(r_fl) <- extent(r)
plot(r_fl)

p50_range_trug_proj <- projectRaster(r_fl, p50_range)

p50.st <- stack(p50_range, p50_range_trug_proj)
names(p50.st) <- c("P50_range", "P50_range_trug")

p50_dif <- p50.st$P50_range - p50.st$P50_range_trug

p50.df <- as.data.frame(raster::sampleRandom(p50.st, 100000))

cor(p50.df)

mod <- lm(P50_range ~ P50_range_trug, data = p50.df)
summary(mod)

#### range P50 Plots ####

### Scatterplot

p <- ggplot(p50.df, aes(x = P50_range_trug, y = P50_range)) + geom_point(alpha = 0.3, size = 0.1) + geom_smooth(method = "lm", size = 0.5) +
  theme(axis.title = element_blank())
p
png(paste0(rslts.dir, "trugman_scatterplot_P50_range.png"), height = h, width = h, units = "in", res = 1000)
p
dev.off()
print(paste0(rslts.dir, "trugman_scatterplot_P50_range.png"))

### Rasters

p50.spdf <- data.frame(as(p50.st, "SpatialPixelsDataFrame"))

p50.spdf <- p50.spdf[complete.cases(p50.spdf$P50_range, p50.spdf$P50_range_trug), ]

p50.spdf$dif <- p50.spdf$P50_range - p50.spdf$P50_range_trug

cols <- brewer.pal(10, "Spectral")

p <- ggplot(p50.spdf, aes(x = x, y = y, fill = P50_range)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "p50_range_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "p50_range_raster.png"))


p <- ggplot(p50.spdf, aes(x = x, y = y, fill = P50_range_trug)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "p50_range_trug_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "p50_range_trug_raster.png"))


p <- ggplot(p50.spdf, aes(x = x, y = y, fill = dif)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "p50_range_difference_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "p50_range_difference_raster.png"))


#### range HSM ####

nc_data <- nc_open('data/Trugman_hydraulics_projection/CWM_HSM_025Deg.nc')

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)

length(lon)
length(lat)

HSM_range_trug <- ncvar_get(nc_data, "CWM_range_within_grid")
dim(HSM_range_trug)

fillvalue <- ncatt_get(nc_data, "max_range_within_grid", "_FillValue")
fillvalue
HSM_range_trug[HSM_range_trug == fillvalue$value] <- NA

r <- raster(t(HSM_range_trug), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

r_fl <- raster::flip(r, direction = 1)
r_fl <- t(r_fl)

plot(r)
plot(r_fl)
extent(r_fl) <- extent(r)
plot(r_fl)

HSM_range_trug_proj <- projectRaster(r_fl, HSM_range)

HSM.st <- stack(HSM_range, HSM_range_trug_proj)
names(HSM.st) <- c("HSM_range", "HSM_range_trug")

plot(HSM.st)
HSM_dif <- HSM.st$HSM_range - HSM.st$HSM_range_trug
plot(HSM_dif)

HSM.df <- as.data.frame(raster::sampleRandom(HSM.st, 100000))

cor(HSM.df)

mod <- lm(HSM_range ~ HSM_range_trug, data = HSM.df)
summary(mod)

#### range HSM Plots ####

### Scatterplot

p <- ggplot(HSM.df, aes(x = HSM_range_trug, y = HSM_range)) + geom_point(alpha = 0.3, size = 0.1) + geom_smooth(method = "lm", size = 0.5) +
  theme(axis.title = element_blank())
p
png(paste0(rslts.dir, "trugman_scatterplot_HSM_range.png"), height = h, width = h, units = "in", res = 1000)
p
dev.off()
print(paste0(rslts.dir, "trugman_scatterplot_HSM_range.png"))


### Rasters

HSM.spdf <- data.frame(as(HSM.st, "SpatialPixelsDataFrame"))

HSM.spdf <- HSM.spdf[complete.cases(HSM.spdf$HSM_range, HSM.spdf$HSM_range_trug), ]

HSM.spdf$dif <- HSM.spdf$HSM_range - HSM.spdf$HSM_range_trug

cols <- brewer.pal(10, "Spectral")

p <- ggplot(HSM.spdf, aes(x = x, y = y, fill = HSM_range)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "HSM_range_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "HSM_range_raster.png"))


p <- ggplot(HSM.spdf, aes(x = x, y = y, fill = HSM_range_trug)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "HSM_range_trug_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "HSM_range_trug_raster.png"))


p <- ggplot(HSM.spdf, aes(x = x, y = y, fill = dif)) + geom_raster() + scale_fill_gradientn(colours= cols)
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
p
png(paste0(rslts.dir, "HSM_range_difference_raster.png"))
p
dev.off()
print(paste0(rslts.dir, "HSM_range_difference_raster.png"))
