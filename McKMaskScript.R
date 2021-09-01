###############################
#### Make secr Habitat Mask ###
###############################

#Read in Capture History
McKcapt.comb<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/McKcapt.comb.rds")

#Read in McKenzie boundary and transform to UTMs.
McKHabitat<- readOGR("C:/Users/nelsonj7/Box/secrMcK/McKenzieHabitat_Clipped.shp")
McKHabitat<-spTransform(McKHabitat, CRS("+proj=utm +zone=10 +datum=NAD83")) #transform to UTMs

#scale-down habitat boundary from meters to kilometers
extractCoords <- function(sp.df)
{
  results <- list()
  for(i in 1:length(sp.df@polygons[[1]]@Polygons))
  {
    results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
  }
  results
}

vertices<-extractCoords(McKHabitat)

meters.to.km<-function(meters){
  km<-meters/1000
}

scaled.vertices<-lapply(vertices, meters.to.km)

Polys<-list()
for(i in 1:length(scaled.vertices)){
  Polys[i]<-sp::Polygon(scaled.vertices[[i]])
}
Polys.plural<-sp::Polygons(Polys, ID = "0")
Polys.sp<-sp::SpatialPolygons(list(Polys.plural), proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
McK.spdf<-sp::SpatialPolygonsDataFrame(Polys.sp, data = McKHabitat@data)

##Creating a blank mask with boundary of McK and buffer of 8,200. Then plotting traps (transects)
McKMask.comb<-make.mask(traps(McKcapt.comb), type = "trapbuffer", buffer=8, poly=McK.spdf)
plot(McKMask.comb)
plot(traps(McKcapt.comb), detpar=list(pch=16, cex=0.8), add=TRUE)

##Extracting points from mask
McKMaskcoords<-as.data.frame(McKMask.comb)
McKMaskcoords<-McKMaskcoords[,c(1,2)]*1000
McKMaskcoords<- SpatialPoints(coords=McKMaskcoords, proj4string = CRS("+proj=utm +zone=10"))

###Create data frame for extracted values
McKMaskextracts2018<-McKMaskcoords
McKMaskextracts2019<-McKMaskcoords

##Read in covariate layers and extract values to created data frame.
###DDE
McKDDE<-raster("C:/Users/nelsonj7/Box/secrMcK/MCKENZIE_mean_dde1.tif")
McKDDE<-projectRaster(McKDDE, crs = "+proj=utm +zone=10 +datum=NAD83") #project to UTMs
plot(McKDDE)
DDEex<-extract(McKDDE, McKMaskcoords, method = 'simple')
scaledDDEex<-scale(DDEex) #Z-scaling values to center around 0 with 1 SD

McKMaskextracts2018$DDE<-scaledDDEex
McKMaskextracts2019$DDE<-scaledDDEex

###Percent Slope
McKSlope<-raster("C:/Users/nelsonj7/Box/secrMcK/MCK_mean_slp1.tif")
Slopeex<-extract(McKSlope, McKMaskcoords, method = 'simple')
scaledSlopeex<-scale(Slopeex)

McKMaskextracts2018$PercentSlope<-scaledSlopeex
McKMaskextracts2019$PercentSlope<-scaledSlopeex

###DFE
McKDFE<-raster("C:/Users/nelsonj7/Box/secrMcK/MCKd_edgeint1.tif")
McKDFE<-projectRaster(McKDFE, crs = "+proj=utm +zone=10 +datum=NAD83") #project to UTMs
McKDFEex<-extract(McKDFE, McKMaskcoords, buffer = 1000, fun = mean)
scaledDFEex<-scale(McKDFEex)

McKMaskextracts2018$DFE<-scaledDFEex
McKMaskextracts2019$DFE<-scaledDFEex

###Distance to Roads
#McKRoads<- readOGR("C:/Users/nelsonj7/Box/secrMcK/McK_Roads_retry.shp")
#McKRoads<-spTransform(McKRoads, CRS("+proj=utm +zone=10 +datum=NAD83")) #transform to UTMs
#Template<-McKDFE
#Template[]<-NA
#McKRoadsRaster<-rasterize(McKRoads, Template, field=1)
#saveRDS(McKRoadsRaster, "McKRoadsRaster.rds")
#McKRoads.dist<-distance(McKRoadsRaster)
#saveRDS(McKRoads.dist, file="McKRoadsDist.rds")
McKRoads.dist<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKRoadsDist.rds")
#writeRaster(McKRoads.dist, file = "C:/Users/nelsonj7/Box/secrMcK/McKRoadsDist.tif")
McKRoads.dist.scale<-scale(McKRoads.dist)
scaledRoadsex<-extract(McKRoads.dist.scale, McKMaskcoords, method = 'simple')

McKMaskextracts2018$Roads<-scaledRoadsex
McKMaskextracts2019$Roads<-scaledRoadsex

###Precipitation
McKPrecip.2019<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKenziePrecip_2019_avg.rds")
McKPrecip.2018<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKenziePrecip_2018_avg.rds")

McKPrecip.2018.extract<-extract(McKPrecip.2018, McKMaskcoords, method = 'simple')
McKPrecip.2019.extract<-extract(McKPrecip.2019, McKMaskcoords, method = 'simple')

McKPrecip.2018.scale<-scale(McKPrecip.2018.extract)
McKPrecip.2019.scale<-scale(McKPrecip.2019.extract)

McKMaskextracts2018$Precip<-McKPrecip.2018.scale
McKMaskextracts2019$Precip<-McKPrecip.2019.scale

##Distance to Agriculture
#McKCropsPoly<-readOGR("C:/Users/nelsonj7/Box/secrMcK/McKCrops.shp")
#McKCrops<-spTransform(McKCropsPoly, CRS("+proj=utm +zone=10 +datum=NAD83"))

#Template<-McKDFE
#Template[]<-NA
#McKCropsRaster<-rasterize(McKCrops, Template, field=1)
#saveRDS(McKCropsRaster, "McKCropsRaster.rds")
#McKCrops.dist<-distance(McKCropsRaster)
#saveRDS(McKCrops.dist, file="McKCropsDist.rds")

McKCrops.dist<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKCropsDist.rds")
#writeRaster(McKCrops.dist, file ="C:/Users/nelsonj7/Box/secrMcK/McKCropsDist.tif")
McKCrops.dist.scale<-scale(McKCrops.dist)
scaledCropsex<-extract(McKCrops.dist.scale, McKMaskcoords, method = 'simple')

McKMaskextracts2018$Crops<-scaledCropsex
McKMaskextracts2019$Crops<-scaledCropsex

###All covariates added
head(McKMaskextracts2018)
head(McKMaskextracts2019)

##DDE Values with mean
library(imputeTS)
McKMaskextracts2018$DDE<-na_mean(McKMaskextracts2018$DDE)
McKMaskextracts2019$DDE<-na_mean(McKMaskextracts2019$DDE)

saveRDS(McKMaskextracts2018, file = "McKMaskextracts.2018.comb.km_8000buffer.rds")
saveRDS(McKMaskextracts2019, file = "McKMaskextracts.2019.comb.km_8000buffer.rds")

McKMaskextracts2018<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKMaskextracts.2018.comb_15000buffer4k.rds")
McKMaskextracts2019<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKMaskextracts.2019.comb_15000buffer4k.rds")

###Add extracted values as covariates to mask
covariates(McKMask.comb$McK18)<-McKMaskextracts2018
covariates(McKMask.comb$McK19)<-McKMaskextracts2019
saveRDS(McKMask.comb, file = "McKMask.comb.km_8000buffer.rds")