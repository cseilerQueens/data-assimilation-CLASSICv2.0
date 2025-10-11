.libPaths("/home/cseiler/AMBER/renv/library/R-4.3/x86_64-pc-linux-gnu")

setwd("/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0")

library(raster)
rm(list=ls())

# Number of grid cells to sample:
nGridcells <- 160 # 80 # 60
minFracArea <- 0.99
my.seed <- 7 # 1 for 40 GCs during tests, 7 usually


# Get WWF Biomes
# WWF_BIOME_01: "Tropical & Subtropical Moist Broadleaf Forests"
# WWF_BIOME_02: "Tropical & Subtropical Dry Broadleaf Forests"
# WWF_BIOME_03: "Tropical & Subtropical Coniferous Forests"
# WWF_BIOME_04: "Temperate Broadleaf & Mixed Forests"
# WWF_BIOME_05: "Temperate Conifer Forests"
# WWF_BIOME_06: "Boreal Forests/Taiga"
# WWF_BIOME_07: "Tropical & Subtropical Grasslands, Savannas & Shrublands"
# WWF_BIOME_08: "Temperate Grasslands, Savannas & Shrublands"
# WWF_BIOME_09: "Flooded Grasslands & Savannas"
# WWF_BIOME_10: "Montane Grasslands & Shrublands"
# WWF_BIOME_11: "Tundra"
# WWF_BIOME_12: "Mediterranean Forests, Woodlands & Scrub"
# WWF_BIOME_13: "Deserts & Xeric Shrublands"
# WWF_BIOME_14: "Mangroves"

b01 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_01")
b02 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_02")
b03 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_03")
b04 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_04")
b05 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_05")
b06 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_06")
b07 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_07")
b08 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_08")
b09 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_09")
b10 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_10")
b11 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_11")
b12 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_12")
# b13 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_13")
# b14 <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_T63.nc", varname = "WWF_BIOME_14")

biomes.list <- list(b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12)

biomes <- raster::stack(b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12)
biomes[biomes == 0] <- NA
biomes <- sum(biomes, na.rm = TRUE)

# Get total area of biome land
areaAllGridcells <- raster::area(biomes)
areaBiomes <- areaAllGridcells * biomes
totalAreaBiomes <- sum(getValues(areaBiomes), na.rm = TRUE)

# Get fractional land area
# FLND <- raster("GC.nc")
FLND <- raster::raster("rsFile.nc", varname = "FLND")

# Specify the minimum fractional area
FLND[FLND < minFracArea] <- NA
FLND <- FLND - FLND + 1

# Exclude certain biomes (deserts and mangroves)
FLND <- FLND * biomes

# png("biomes.png")
# plot(FLND, col = "grey")
# dev.off()

# Initialize and empty vector
xy.allBiomes <- numeric()
xyPerLat <- numeric()

for (i in 1:12) {
  # Get biome and set 0 to NA
  biome <- biomes.list[[i]]
  biome[biome == 0] <- NA
  
  # Calculate the fraction of area covered by biome
  areaEachGridcell <- raster::area(biome)
  areaEachGridcell <- areaEachGridcell * biome
  areaBiome <- sum(getValues(areaEachGridcell), na.rm = TRUE)
  fractionBiome <- areaBiome / totalAreaBiomes
  
  nbiome <- round((fractionBiome * nGridcells))
  
  # Exclude cells based on minimum fractional cover
  biome <- biome * FLND
  # Get coordinates
  xy <- xyFromCell(biome, which(biome[] == 1))
  # Randomly select grid cells for each biome
  
  set.seed(my.seed)
  select <- sample(x = 1:nrow(xy),
                   size = nbiome,
                   replace = FALSE)
  
  # shuffle data
  set.seed(1)
  xy <- xy[sample(nrow(xy), replace = FALSE), ]
  
  # Select points
  xy <- xy[select, ]
  xy <- matrix(xy, ncol = 2)
  
  xy.allBiomes <- rbind(xy.allBiomes, xy)
}

lon <- xy.allBiomes[,1]
lat <- xy.allBiomes[,2]

regular <- "+proj=longlat +ellps=WGS84"
points <- data.frame(lon,lat,1)
coordinates(points) <- ~lon+lat

library(raster)
png("FLND.png")
plot(FLND)
points(points)
dev.off()

# Create a netCDF file

z <- 1

data <- data.frame(xy.allBiomes,z)
data <- rasterFromXYZ(data)
projection(data) <-"+proj=longlat +datum=WGS84"

names(data) <- "FLND"  

data <- extend(data, FLND, NA)

data <- FLND * data

# Flip the raster so that the NH is at the bottom rather than at the top
# This is necessary so that locations are correct when inserting this data
# into the restart file using nco
data <- flip(data, "y")

writeRaster(data, filename="tmp_GC.nc", format = "CDF", varname = "FLND", overwrite=T)

png("tmp_GC.png")
data <- raster("tmp_GC.nc")
plot(data)
dev.off()

# Make a publishable figure that is Atlantic-centered
library(amber)

lon <- xy.allBiomes[,1]
lat <- xy.allBiomes[,2]

lon <- ifelse(lon > 180, lon-360, lon)

regular <- "+proj=longlat +ellps=WGS84"
points <- data.frame(lon,lat,1)
coordinates(points) <- ~lon+lat

my.xlim <- c(-180, 180)
my.ylim <- c(-60, 85)
my.projection <- "+proj=longlat +ellps=WGS84"
shp.filename <- system.file("extdata/ne_110m_land/ne_110m_land.shp", package = "amber")
land <- amber::intFun.coast(my.xlim, my.ylim, my.projection, shp.filename)

# Get high-resolution biome data:
# Create an empty raster:
allBiomes <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_0.5deg.nc", varname = "WWF_BIOME_01")
allBiomes[allBiomes < 10] <- NA

for (i in 1:13) {

ii <- sprintf("%02d", i)

varname <- paste("WWF_BIOME", ii, sep = "_")
biome <- raster("/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/WWF_biomes/WWF_BIOMES_0.5deg.nc", varname = varname)
# biome <- raster::rotate(biome)
biome[biome == 0] <- NA
biome[biome == 1] <- i
allBiomes <- raster::merge(allBiomes, biome)
}


myExtent <- c(-180, 180, -90, 90)

plot.width <- 8
plot.height <- 3.4

png("gridCells.png", width = plot.width, height = plot.height, units = "in", res = 300)
graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(2, 2, 1, 4), lwd = 1, cex = 1)

# colors        
my.breaks <- seq(0.5,13.5,1) # breaks of the colors
my.labels <- seq(1,13,1) # locations where to set the labels
my.col <- terrain.colors(n=length(my.breaks)) # color scale

my.legend <- list(at=my.labels, labels=my.labels, cex.axis=1.0)
my.legend.args <- list(text=expression("Biomes"), side=2, font=1, line=1, cex=1)

plot(allBiomes, col=my.col, breaks=my.breaks, legend=F, main=NA, 
ylim = c(-52,90), xlab="Degrees Longitude", ylab="Degrees Latitude")
plot(allBiomes, legend.only=T, col=my.col, breaks=my.breaks, axis.args = my.legend,
legend.args=my.legend.args, legend.width=2, legend.shrink=1.0, font=1)


raster::plot(land, col = NA, border = 1, add = TRUE)
graphics::points(points, pch = 16, col = "black")

dev.off()

