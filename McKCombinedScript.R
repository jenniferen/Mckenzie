library(secr)
library(rgdal)
library(sp)
library(raster)
library(sf)
library(AICcmodavg)
library(imputeTS)

memory.limit(size = 99999999)
setwd("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie")

captfile2019<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2019McKCH.txt")
trapfile2019<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2019McKTraps.txt")

captfile2018<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018McKCH.txt")
trapfile2018<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018McKTraps.txt")

head(captfile2019)
head(captfile2018)
head(trapfile2018)
head(trapfile2019)
##From 2018 secr script: capt/trap file"McK2018.capthist.4k.txt", "2018McKTraps4k.txt"
#Capture history format is: Session, AnimalID, Occasion,Detector, Sex
#Need to make separate sessions such as: McK18 and McK19
#AnimalIDs need to be separate, add max number of individuals from 18 to every number in 19.
#Occasion stays 1
#X and Y stay the same

#Trap file format is: TRAP ID, X, Y, Character, Julian, effort, AvgPPT
#I think I just need to stack trap files on top of each other and add max # of detectors from 18 to 19
#to correct the TrapID

captfile2019[,2]<-captfile2019[,2]+max(captfile2019[,2])
captfile2019[,4]<-captfile2019[,4]+max(trapfile2018[,1])

trapfile2019[,1]<-trapfile2019[,1] + max(trapfile2018[,1])

combinedcaptfile<-rbind(captfile2018, captfile2019)

trapfile2018[,c(2,3)]<-trapfile2018[,c(2,3)]/1000
trapfile2019[,c(2,3)]<-trapfile2019[,c(2,3)]/1000

write.table(combinedcaptfile, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/combinedcaptfile.txt", row.names = FALSE, quote = FALSE)
write.table(trapfile2018, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018McKtrapfile.txt", row.names = FALSE, quote = FALSE)
write.table(trapfile2019, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2019McKtrapfile.txt", row.names = FALSE, quote = FALSE)
#Add # to first row of both .txt files!

###################################################################################
###################################################################################
####  #  #  #  #SECR#  #  #  #  #####################################################
###################################################################################
###################################################################################
library(secr)
McKcapt.comb<-read.capthist("combinedcaptfile.txt", trapfile = c("2018McKtrapfile.txt", "2019McKtrapfile.txt"), 
                              detector="count", fmt = "trapID", covnames = "Sex", 
                              trapcovnames = c("Julian", "effort", "precip", "Observers"))

summary(McKcapt.comb)
plot(McKcapt.comb, tracks = TRUE)

##Default calculation of sigma
(initialsigma <- RPSV(McKcapt.comb.telem, CC = TRUE))
fit<-secr.fit(McKcapt.comb.telem, buffer = 1.486603, trace=FALSE)
fit
plot(fit, limits = TRUE)


CHt<-read.telemetry(file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018.2019.McKTelemetry.km.txt",
                    covnames = c("Sex"))

CHt
covariates(CHt)
covariates(CHt)$Sex<-as.factor(rep("F", nrow(covariates(CHt))))
McKcapt.comb.telem<-addTelemetry(McKcapt.comb, CHt, type = "independent")

##################
#### Make Mask ###
##################

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
McKMask.comb<-make.mask(traps(McKcapt.comb.telem), type = "trapbuffer", buffer=8, poly=McK.spdf)
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
McKMask.comb<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/McKMask.comb.km_8000buffer.rds")
polyarea(McKMask.comb)

plot(McKMask.comb)
covariates(McKMask.comb)

###Spatial Location Covariate is within the mask as x and y.

###########################
###Adding Telemetry data###
###########################

CHt<-read.telemetry(file = "C:/Users/nelsonj7/Box/secrMcK/2018.2019.McKTelemetry.txt")

McKcapt.comb.telem<-addTelemetry(McKcapt.comb, CHt, type = "independent")
saveRDS(McKcapt.comb.telem, file = "C:/Users/nelsonj7/Box/secrMcK/McKcapt.comb.telem.rds")
McKcapt.comb.telem<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKcapt.comb.telem.rds")

RPSV(CHt, CC = TRUE) #can you get a confidence interval? and model with the mean, lower and upper bounds
#to show if it's super sensitive to this.

######################################
##initial fitting with global models##
######################################

null.model<-secr.fit(McKcapt.comb, model = list(D~1, g0~1, sigma~1), fixed = list(sigma = 1.855),
                     mask = McKMask.comb, binomN = 1, start = c(25, -23),
                     method = "Nelder-Mead")

saveRDS(null.model, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/null.model.rds")
#null.model<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/null.model.rds")

g0.session<-secr.fit(McKcapt.comb, model = list(D~1, g0~Session, sigma~1), fixed = list(sigma = 1.855),
                     mask = McKMask.comb, binomN = 1, start = c(7.5, -3, -2),
                     method = "Nelder-Mead")

saveRDS(g0.session, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.rds")
#g0.session<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.rds")

g0.Julian<-secr.fit(McKcapt.comb, model = list(D~1, g0~Julian, sigma~1), fixed = list(sigma = 1.855),
                    mask = McKMask.comb, binomN = 1, start = c(7.5, -3, -2),
                    method = "Nelder-Mead")

saveRDS(g0.Julian, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.Julian.rds")
#g0.Julian<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.Julian.rds")

g0.effort<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort, sigma~1), fixed = list(sigma = 1.855),
                    mask = McKMask.comb, binomN = 1, start = c(7.5, -3, -2),
                    method = "Nelder-Mead")

saveRDS(g0.effort, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.rds")
#g0.effort<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.rds")

g0.precip<-secr.fit(McKcapt.comb, model = list(D~1, g0~precip, sigma~1), fixed = list(sigma = 1.855),
                    mask = McKMask.comb, binomN = 1, start = c(7.5, -3, -2),
                    method = "Nelder-Mead")

saveRDS(g0.precip, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.precip.rds")
#g0.precip<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.precip.rds")

g0.univariate.models<-secrlist(g0.session, null.model, g0.Julian, g0.effort, g0.precip)

AIC(g0.univariate.models)

g0.effort.Julian<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort+Julian), 
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                           start = c(6.4, -4, 0.35, 0.1), method = "Nelder-Mead")

saveRDS(g0.effort.Julian, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.rds")
#g0.effort.Julian<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.rds")

g0.effort.precip<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort+precip),
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                           start = c(6.4, -4, 0.35, 0.4), method = "Nelder-Mead")

saveRDS(g0.effort.precip, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.precip.rds")
#g0.effort.precip<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.precip.rds")

g0.effort.effort2<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort + effort^2),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, bionomN = 1,
                            start = c(6.4, -4, 0.35, 0.35), method = "Nelder-Mead")
saveRDS(g0.effort.effort2, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.effort2.rds")
#g0.effort.effort2<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.effort2.rds")

g0.Julian.precip<-secr.fit(McKcapt.comb, model = list(D~1, g0~Julian+precip),
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                           start = c(6.4, -4, -0.2, 0.4), method = "Nelder-Mead")
saveRDS(g0.Julian.precip, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.Julian.precip.rds")
#g0.Julian.precip<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.Julian.precip.rds")

g0.effort.Julian.precip<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort + precip + Julian),
                                  fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                                  start = c(6.4, -4, 0.35, 0.4, 0.1), method = "Nelder-Mead")

saveRDS(g0.effort.Julian.precip, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.precip.rds")
#g0.effort.Julian.precip<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.precip.rds")

g0.models<-secrlist(null.model, g0.session, g0.effort, g0.precip, g0.Julian, g0.Julian.precip, g0.effort.Julian, g0.effort.precip,
                    g0.effort.Julian.precip)

AIC(g0.models)

D.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE, g0~effort + precip),
                    fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                    start = c(6.49, 4, -4.3, -0.7, 0.47), method = "Nelder-Mead")

saveRDS(D.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.g0.rds")
#D.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.g0.rds")

D.DFE.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE, g0~effort + precip, sigma~1),
                    fixed = list(sigma = 1.855), mask = McKMask.comb, 
                    binomN = 1,  start = c(5, -1, -2.9, -0.4, 0.3), method = "Nelder-Mead")

saveRDS(D.DFE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.g0.rds")
#D.DFE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.g0.rds")

D.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads, g0~effort + precip, sigma~1),
                      fixed = list(sigma = 1.855), mask = McKMask.comb, 
                      binomN = 1,  start = c(5, -1, -2.9, -0.4, 0.3), method = "Nelder-Mead")

saveRDS(D.Roads.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.g0.rds")
#D.Roads.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.g0.rds")

D.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~PercentSlope, g0~effort + precip, sigma~1),
                      mask = McKMask.comb, binomN = 1)

saveRDS(D.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.g0.rds")
#D.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.g0.rds")

D.precip.g0.<-secr.fit(McKcapt.comb, model = list(D~Precip, g0~effort + precip, sigma~1),
                       mask = McKMask.comb, binomN = 1)

saveRDS(D.precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.g0.rds")
#D.precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.g0.rds")

D.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops, g0~effort + precip, sigma ~ 1),
                      mask = McKMask.comb, binomN = 1)

saveRDS(D.precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.g0.rds")

D.univariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.Crops.g0.)

AIC(D.univariate.models)
D.Crops.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + DDE, g0~effort + precip, sigma ~1),
                          mask = McKMask.comb, binomN = 1)

saveRDS(D.Crops.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.DDE.g0.rds")

D.Crops.DFE.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + DFE, g0~effort + precip, sigma~1),
                        mask = McKMask.comb, binomN = 1)

saveRDS(D.Crops.DFE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.DFE.g0.rds")

D.Crops.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + Roads, g0~effort + precip, sigma~1),
                          mask = McKMask.comb, binomN = 1)

saveRDS(D.Crops.Roads.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.Roads.g0.rds")

D.Crops.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + PercentSlope, g0~effort + precip, sigma~1),
                          mask = McKMask.comb, binomN = 1)

saveRDS(D.Crops.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.slope.g0.rds")

D.Crops.precip.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + Precip, g0~effort + precip, sigma~1),
                          mask = McKMask.comb, binomN = 1)

saveRDS(D.Crops.DFE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.precip.g0.rds")

D.bivariate.models<-secrlist(D.Crops.DDE.g0., D.Crops.DFE.g0., D.Crops.Roads.g0., D.Crops.slope.g0., D.Crops.precip.g0.)
AIC(D.bivariate.models)

D.Roads.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + DFE, g0~effort + precip, sigma~1),
                        mask = McKMask.comb, binomN = 1)

saveRDS(D.Roads.DDE.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.DDE.g0.rds")

D.Roads.DFE.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + DDE, g0~effort + precip, sigma~1),
                          mask = McKMask.comb, binomN = 1)

saveRDS(D.Roads.DFE.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.DFE.g0.rds")

D.Roads.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + PercentSlope, g0~effort + precip, sigma~1),
                            mask = McKMask.comb, binomN = 1)

saveRDS(D.Roads.slope.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.Slope.g0.rds")
#D.Roads.slope.g0.<-readRDS("C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.Slope.g0.rds")

D.Roads.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + Precip, g0~effort + precip, sigma~1),
                            mask = McKMask.comb, binomN = 1)

D.Roads.precip.g0.<-readRDS("C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.precip.g0.rds")

D.Roads.models<-secrlist(D.Roads.DDE.g0., D.Roads.DFE.g0., D.Roads.slope.g0., D.Roads.precip.g0.,D.DDE.g0., D.Roads.g0., D.precip.g0., D.slope.g0.)

AIC(D.Roads.models)

topmodel<-readRDS("C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.DFE.precip.g0.rds")

D.Roads.DFE.precip.g0.<-readRDS("C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.DFE.precip.g0.rds")

D.DFE.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE, g0~effort + precip, sigma~1),
                        mask = McKMask.comb, binomN = 1)

saveRDS(D.DFE.DDE.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.DDE.g0.rds")

D.DFE.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads, g0~effort + precip, sigma~1),
                          mask = McKMask.comb)

saveRDS(D.DFE.Roads.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.Roads.g0.rds")

D.DFE.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + PercentSlope, g0~effort + precip, sigma~1),
                          mask = McKMask.comb)

saveRDS(D.DFE.slope.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.slope.g0.rds")

D.DFE.precip.g0.hn<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Precip, g0~effort + precip, sigma~1),
                             mask = McKMask.comb, binomN = 1)

saveRDS(D.DFE.precip.g0., file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.precip.g0.telem.rds")

D.DFE.precip.g0.effort<-secr.fit(McKcapt, model = list(D~DFE + Precip, g0~effort, sigma~1),
                                 mask = McKMask)

saveRDS(D.DFE.precip.g0.effort, file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.precip.g0.effort.rds")

D.models<-secrlist(D.DDE.g0., D.Roads.g0., D.precip.g0., D.slope.g0., D.DFE.DDE.g0., D.DFE.precip.g0.,
                   D.DFE.Roads.g0., D.DFE.slope.g0.)

D.DFE.precip.g0.sig.Session<-secr.fit(McKcapt.comb, model = list(D~DFE + Precip, g0~effort + precip, sigma ~session),
                                      mask = McKMask.comb, binomN = 1)
saveRDS(D.DFE.precip.g0.sig.Session, file = "C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.precip.g0.sig.Session.rds")

AIC(D.models)

D.models.finalset<-secrlist(topmodel, D.Roads.DDE.g0., D.Roads.DFE.g0., D.Roads.slope.g0., D.Roads.precip.g0.,D.DDE.g0., 
                         D.Roads.g0., D.precip.g0., D.slope.g0., D.Roads.DFE.precip.g0.)

AIC(D.models.finalset)

D.DFE.DDE.g0.
D.DFE.slope.g0.

D.g0.sigma.sex<-secr.fit(McKcapt, model = list(D~DFE, g0~effort + precip, sigma~Sex), CL = TRUE,
                         mask = McKMask)

D.slope.crops.g0.<-secr.fit(McKcapt.comb, model = list(D~PercentSlope + Crops, g0~effort + precip, sigma~1),
                            mask = McKMask.comb, binomN = 1)

D.x.x2.g0.<-secr.fit(McKcapt.comb, model = list(D~x + x2, g0~effort +precip, sigma~1), 
                     mask = McKMask.comb, binomN = 1)

saveRDS(D.slope.crops.g0., "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.crops.g0.rds")

D.slope.crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.crops.g0.rds")

coefficients(D.slope.crops.g0.)

D.slope.crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.crops.g0.rds")

c_hat(D.DFE.g0.)

model<-readRDS("C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/15000buffer_4K/D.DFE.precip.g0.rds")
model$fit
coef(model)

########
##run in parallel
#######
#detection.fits <- par.secr.fit (
#  c('secr.null','secr.t', 'secr.bk','secr.sd', 'secr.bk.sd','secr.t.bk','secr.bait.bk','secr.dna.bk'),
#  ncores = 8) ##755.575 minutes 

# CALCULATE chat

derived(D.Roads.DFE.precip.g0.)
#chat<-NA
#2018
D <- 0.06792655 #density
g0 <- -3.3599118 #
sigma <- 6.9732202
fn <- function (r) r * g0 * exp(-r^2/2/(sigma/100)^2)
expected.nk <- 2 * pi *D * integrate(fn, 0, 5*sigma/100)$value
summary(McKcapt.comb.telem)
nk<-122
#nk <- apply(apply(ch2,c(1,3),sum)>0,2,sum) #realized number of detections
X2 <- sum((nk - expected.nk)^2 / expected.nk)
si <- mean((nk - expected.nk) / expected.nk)
nu <- 121 #degrees of freedom
chat.2018 <- X2/nu / (1 + si) #Fletcher 2012
#out$dispersion <- c(mean.nk = mean(nk), var.nk = sd(nk)^2, chat = chat)
chat.2018

D <- 0.02404881 #density
g0 <- -4.1884166 #
sigma <- 6.6049448
fn <- function (r) r * g0 * exp(-r^2/2/(sigma/100)^2)
expected.nk <- 2 * pi *D * integrate(fn, 0, 5*sigma/100)$value
summary(McKcapt.comb)
nk<-14
#nk <- apply(apply(ch2,c(1,3),sum)>0,2,sum) #realized number of detections
X2 <- sum((nk - expected.nk)^2 / expected.nk)
si <- mean((nk - expected.nk) / expected.nk)
nu <- 14 #degrees of freedom
chat.2019 <- X2/nu / (1 + si) #Fletcher 2012
#out$dispersion <- c(mean.nk = mean(nk), var.nk = sd(nk)^2, chat = chat)
chat.2019


#probability of detection function
d<-seq(0,1,0.01)
lambda.d<-function(d){
  lambda0<-g0*exp(-d^2/2*sigma^2)
  return(lambda0)
}

hazard_halfnormal<-function(lambda.d){
  g0.y<-1-exp(-lambda.d(d))
  return(g0.y)
}

plot(d, hazard_halfnormal(lambda.d(d)))
hazard_halfnormal(d)
lambda.d(d)
hazard_halfnormal(lambda.d(d))
y<-lambda.d(d)
hazard_halfnormal(y)

###Predict Density

topmodelMcK<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.slope.g0.rds")
topmodelDsurface<-predictDsurface(topmodelMcK, mask = McKMask.comb)
plot(topmodelDsurface, col = brewer.pal(8, "YlGn"))
library(RColorBrewer)

McK18.Dsurface<-raster(topmodelDsurface$McK18, covariate = 'D.0', crs = "+proj=utm +zone=10 +datum=NAD83")
McK19.Dsurface<-raster(topmodelDsurface$McK19, covariate = 'D.0', crs = "+proj=utm +zone=10 +datum=NAD83")
writeRaster(McK18.Dsurface, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/McK18.Dsurface.tif", overwrite = TRUE)
writeRaster(McK19.Dsurface, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/McK19.Dsurface.tif", overwrite = TRUE)