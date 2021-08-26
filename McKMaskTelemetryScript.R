library(secr)
library(rgdal)
library(sp)
library(raster)
library(sf)
library(AICcmodavg)


#############################
########ADD??################
##
#D.DDE.g0.
#covariates(TiogaMask)
#covariates(TiogaMask)$DDE
#is.na(covariates(TiogaMask)$DDE)
#covariates(TiogaMask)$DDE[is.na(covariates(TiogaMask$DDE))] <- 0
############################
############################


setwd("C:/Users/nelsonj7/Box/secrMcK")

captfile2019<-read.table("C:/Users/nelsonj7/Box/secrMcK/2019CaptureData.txt")
trapfile2019<-read.table("C:/Users/nelsonj7/Box/secrMcK/2019McKTraps.txt")

captfile2018<-read.table("C:/Users/nelsonj7/Box/secrMcK/2018CaptureData.txt")
trapfile2018<-read.table("C:/Users/nelsonj7/Box/secrMcK/2018McKTraps.txt")

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

captfile2019[,2]<-captfile2019[,2]+max(captfile2018[,2])
captfile2019[,4]<-captfile2019[,4]+max(trapfile2018[,1])

trapfile2019[,1]<-trapfile2019[,1] + max(trapfile2018[,1])

combinedcaptfile<-rbind(captfile2018, captfile2019)

combinedtrapfile<-rbind(trapfile2018, trapfile2019)

write.table(combinedcaptfile, file = "C:/Users/nelsonj7/Box/secrMcK/combinedcaptfile.txt", row.names = FALSE, quote = FALSE)
write.table(combinedtrapfile, file = "C:/Users/nelsonj7/Box/secrMcK/combinedtrapfile.txt", row.names = FALSE, quote = FALSE)

###################################################################################
###################################################################################
####  #  #  #  #SECR#  #  #  #  #####################################################
###################################################################################
###################################################################################
library(secr)
McKcapt.comb<-read.capthist("combinedcaptfile.txt", "combinedtrapfile.txt", 
                            detector="count", fmt = "trapID", covnames = "Sex", 
                            trapcovnames = c("Julian", "effort", "precip"))

summary(McKcapt.comb)
plot(McKcapt.comb, tracks = TRUE)

##Default calculation of sigma
(initialsigma <- RPSV(McKcapt.comb, CC = TRUE))
fit<-secr.fit(McKcapt, buffer = 4*initialsigma, trace=FALSE)
fit
plot(fit, limits = TRUE)


CHt<-read.telemetry(file = "C:/Users/nelsonj7/Box/secrMcK/2018.2019.McKTelemetry.txt",
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

##Creating a blank mask with boundary of McK and buffer of 8,200. Then plotting traps (transects)
McKMask.comb<-make.mask(traps(McKcapt.comb), type = "trapbuffer", buffer=15000, poly=McKHabitat)
plot(McKMask.comb)
plot(traps(McKcapt.comb), detpar=list(pch=16, cex=0.8), add=TRUE)

##Extracting points from mask
McKMaskcoords<-as.data.frame(McKMask.comb)
McKMaskcoords<- SpatialPoints(coords=McKMaskcoords, proj4string = CRS("+proj=utm +zone=10"))

###Create data frame for extracted values
McKMaskextracts2018<-McKMaskcoords
McKMaskextracts2019<-McKMaskcoords

##Read in covariate layers and extract values to created data frame.
###DDE
McKDDE<-raster("MCKENZIE_mean_dde1.tif")
McKDDE<-projectRaster(McKDDE, crs = "+proj=utm +zone=10 +datum=NAD83") #project to UTMs
plot(McKDDE)
DDEex<-extract(McKDDE, McKMaskcoords, method = 'simple')
scaledDDEex<-scale(DDEex) #Z-scaling values to center around 0 with 1 SD

McKMaskextracts2018$DDE<-scaledDDEex
McKMaskextracts2019$DDE<-scaledDDEex

###Percent Slope
McKSlope<-raster("MCK_mean_slp1.tif")
Slopeex<-extract(McKSlope, McKMaskcoords, method = 'simple')
scaledSlopeex<-scale(Slopeex)

McKMaskextracts2018$PercentSlope<-scaledSlopeex
McKMaskextracts2019$PercentSlope<-scaledSlopeex

###DFE
McKDFE<-raster("MCKd_edgeint1.tif")
McKDFE<-projectRaster(McKDFE, crs = "+proj=utm +zone=10 +datum=NAD83") #project to UTMs
McKDFEex<-extract(McKDFE, McKMaskcoords, buffer = 1000, fun = mean)
scaledDFEex<-scale(McKDFEex)

McKMaskextracts2018$DFE<-scaledDFEex
McKMaskextracts2019$DFE<-scaledDFEex

###Distance to Roads
McKRoads<- readOGR("C:/Users/nelsonj7/Box/secrMcK/McK_Roads_retry.shp")
McKRoads<-spTransform(McKRoads, CRS("+proj=utm +zone=10 +datum=NAD83")) #transform to UTMs
Template<-McKDFE
Template[]<-NA
McKRoadsRaster<-rasterize(McKRoads, Template, field=1)
saveRDS(McKRoadsRaster, "McKRoadsRaster.rds")
McKRoads.dist<-distance(McKRoadsRaster)
saveRDS(McKRoads.dist, file="McKRoadsDist.rds")
McKRoads.dist<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKRoadsDist.rds")
writeRaster(McKRoads.dist, file = "C:/Users/nelsonj7/Box/secrMcK/McKRoadsDist.tif")
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
McKCropsPoly<-readOGR("C:/Users/nelsonj7/Box/secrMcK/McKCrops.shp")
McKCrops<-spTransform(McKCropsPoly, CRS("+proj=utm +zone=10 +datum=NAD83"))

Template<-McKDFE
Template[]<-NA
McKCropsRaster<-rasterize(McKCrops, Template, field=1)
saveRDS(McKCropsRaster, "McKCropsRaster.rds")
McKCrops.dist<-distance(McKCropsRaster)
saveRDS(McKCrops.dist, file="McKCropsDist.rds")

McKCrops.dist<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKCropsDist.rds")
writeRaster(McKCrops.dist, file ="C:/Users/nelsonj7/Box/secrMcK/McKCropsDist.tif")
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

saveRDS(McKMaskextracts2018, file = "McKMaskextracts.2018.comb_15000buffer.rds")
saveRDS(McKMaskextracts2019, file = "McKMaskextracts.2019.comb_15000buffer.rds")

McKMaskextracts2018<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKMaskextracts.2018.comb_15000buffer4k.rds")
McKMaskextracts2019<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKMaskextracts.2019.comb_15000buffer4k.rds")

###Add extracted values as covariates to mask
covariates(McKMask.comb$McK18)<-McKMaskextracts2018
covariates(McKMask.comb$McK19)<-McKMaskextracts2019
saveRDS(McKMask.comb, file = "McKMask.comb_15000buffer.rds")
McKMask.comb<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKMask.comb_15000buffer.rds")
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

RPSV(CHt, CC = TRUE)