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
#Capture history format is: Session, AnimalID, Occasion, Detector, Sex

#Trap file format is: TRAP ID, X, Y, Character, Julian, effort, AvgPPT

#Creating unique AnimalID for individuals detected between 2018 and 2019.
captfile2019[,2]<-captfile2019[,2]+max(captfile2019[,2])

#Creating unique detector ID for detectors used between 2018 and 2019
captfile2019[,4]<-captfile2019[,4]+max(trapfile2018[,1])
trapfile2019[,1]<-trapfile2019[,1] + max(trapfile2018[,1])

#Combining 2018 and 2019 capture files into one capture file
combinedcaptfile<-rbind(captfile2018, captfile2019)

#Scaling detector coordinates by 1,000 meters by dividing coordinates by 1,000 m.
trapfile2018[,c(2,3)]<-trapfile2018[,c(2,3)]/1000
trapfile2019[,c(2,3)]<-trapfile2019[,c(2,3)]/1000

#Writing out final .txt files
write.table(combinedcaptfile, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/combinedcaptfile.txt", row.names = FALSE, quote = FALSE)
write.table(trapfile2018, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018McKtrapfile.txt", row.names = FALSE, quote = FALSE)
write.table(trapfile2019, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2019McKtrapfile.txt", row.names = FALSE, quote = FALSE)
#Add # to first row of both .txt files!

###################################
###################################
####            SECR            ###
###################################
###################################
library(secr)
McKcapt.comb<-read.capthist("combinedcaptfile.txt", trapfile = c("2018McKtrapfile.txt", "2019McKtrapfile.txt"), 
                              detector="count", fmt = "trapID", covnames = "Sex", 
                              trapcovnames = c("Julian", "effort", "precip", "Observers"))

summary(McKcapt.comb)
plot(McKcapt.comb, tracks = TRUE)
saveRDS(McKcapt.comb, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/McKcapt.comb.rds")
#McKcapt.comb<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/McKcapt.comb.rds")

##Default calculation of sigma
(initialsigma <- RPSV(McKcapt.comb, CC = TRUE))
fit<-secr.fit(McKcapt.comb, buffer = 1.486603, trace=FALSE)
fit
plot(fit, limits = TRUE)

###########################
###Adding Telemetry data###
###########################

CHt<-read.telemetry(file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018.2019.McKTelemetry.km.txt",
                    covnames = c("Sex"))

CHt
covariates(CHt)
covariates(CHt)$Sex<-as.factor(rep("F", nrow(covariates(CHt))))
McKcapt.comb<-addTelemetry(McKcapt.comb, CHt, type = "independent")

CHt<-read.telemetry(file = "C:/Users/nelsonj7/Box/secrMcK/2018.2019.McKTelemetry.txt")

McKcapt.comb<-addTelemetry(McKcapt.comb, CHt, type = "independent")
saveRDS(McKcapt.comb, file = "C:/Users/nelsonj7/Box/secrMcK/McKcapt.comb.rds")
McKcapt.comb<-readRDS("C:/Users/nelsonj7/Box/secrMcK/McKcapt.comb.rds")

RPSV(CHt, CC = TRUE) #can you get a confidence interval? and model with the mean, lower and upper bounds
#to show if it's super sensitive to this.

##########################
###Read in Habitat Mask###
##########################

McKMask.comb<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/McKMask.comb.km_8000buffer.rds")
polyarea(McKMask.comb)

plot(McKMask.comb)
covariates(McKMask.comb)

##########################
###START RUNNING MODELS###
##########################

##############
##null model##
##############

null.model<-secr.fit(McKcapt.comb, model = list(D~1, g0~1, sigma~1), fixed = list(sigma = 1.855),
                     mask = McKMask.comb, binomN = 1, start = c(25, -23),
                     method = "Nelder-Mead")

saveRDS(null.model, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/null.model.rds")
#null.model<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/null.model.rds")

##########################################
###Probability of detection (g0) models###
##########################################

######################
##Univariate Models###
######################

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

##########################
####Bivariate g0 models###
##########################

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

g0.session.effort<-secr.fit(McKcapt.comb, model = list(D~1, g0~session + effort),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                            start = c(6.4, -4, -0.2, 0.4), method = "Nelder-Mead")
saveRDS(g0.session.effort, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.effort.rds")
#g0.session.effort<-readRDS("C:/Users/nelsonj7/Deskotp/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.effort.rds")

g0.session.Julian<-secr.fit(McKcapt.comb, model = list(D~1, g0~session + Julian),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                            start = c(6.4, -4, -0.2, 0.4), method = "Nelder-Mead")
saveRDS(g0.session.Julian, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.Julian.rds")
#g0.session.Julian<-readRDS("C:/Users/nelsonj7/Deskotp/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.Julian.rds")
  
g0.session.precip<-secr.fit(McKcapt.comb, model = list(D~1, g0~session + precip),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                            start = c(6.4, -4, -0.2, 0.4), method = "Nelder-Mead")
saveRDS(g0.session.precip, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.precip.rds")
#g0.session.precip<-readRDS("C:/Users/nelsonj7/Deskotp/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.session.precip.rds")

##########################
###Trivariate g0 models###
##########################

g0.effort.Julian.precip<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort + precip + Julian),
                                  fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                                  start = c(6.4, -4, 0.35, 0.4, 0.1), method = "Nelder-Mead")

saveRDS(g0.effort.Julian.precip, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.precip.rds")
#g0.effort.Julian.precip<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.precip.rds")

g0.effort.Julian.session<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort + session + Julian),
                                  fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                                  start = c(6.4, -4, 0.35, 0.4, 0.1), method = "Nelder-Mead")

saveRDS(g0.effort.Julian.session, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.session.rds")
#g0.effort.Julian.session<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.session.rds")

g0.effort.precip.session<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort + precip + session),
                                  fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                                  start = c(6.4, -4, 0.35, 0.4, 0.1), method = "Nelder-Mead")

saveRDS(g0.effort.precip.session, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.precip.session.rds")
#g0.effort.precip.session<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.precip.session.rds")

g0.Julian.precip.session<-secr.fit(McKcapt.comb, model = list(D~1, g0~session + precip + Julian),
                                  fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                                  start = c(6.4, -4, 0.35, 0.4, 0.1), method = "Nelder-Mead")

saveRDS(g0.Julian.precip.session, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.Julian.precip.session.rds")
#g0.Julian.precip.session<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.Julian.precip.session.rds")

#####################
###Global g0 model###
#####################

g0.effort.Julian.precip.session<-secr.fit(McKcapt.comb, model = list(D~1, g0~effort + Julian + precip + session),
                                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                          start = c(6.4, -4, 0.35, 0.4, 0.1, -2))

saveRDS(g0.effort.Julian.precip.session, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.precip.session.rds")
#g0.effort.Julian.precip.session<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/g0.effort.Julian.precip.session.rds")

#####################################################
###AIC of all probability of detection (g0) models###
#####################################################

g0.models<-secrlist(null.model, g0.session, g0.effort, g0.precip, g0.Julian, g0.Julian.precip, g0.effort.Julian, g0.effort.precip,
                    g0.session.effort, g0.session.Julian, g0.session.precip, g0.effort.Julian.precip,
                    g0.effort.Julian.session, g0.effort.precip.session, g0.Julian.precip.session, g0.effort.Julian.precip.session)

AIC(g0.models)

####################
###Density models###
####################

###############################
###Univariate Density models###
###############################

D.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE, g0~session + Julian),
                    fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                    start = c(6.49, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.g0.rds")
#D.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.g0.rds")

D.DFE.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE, g0~session + Julian),
                    fixed = list(sigma = 1.855), mask = McKMask.comb, 
                    binomN = 1,  start = c(6.49, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.g0.rds")
#D.DFE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.g0.rds")

D.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads, g0~session + Julian),
                      fixed = list(sigma = 1.855), mask = McKMask.comb, 
                      binomN = 1,  start = c(6.49, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.g0.rds")
#D.Roads.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.g0.rds")

D.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~PercentSlope, g0~session + Julian),
                      fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,  
                      start = c(6.49, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.g0.rds")
#D.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.g0.rds")

D.precip.g0.<-secr.fit(McKcapt.comb, model = list(D~Precip, g0~session + Julian),
                       fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                       start = c(6.49, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.g0.rds")
#D.precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.g0.rds")

D.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops, g0~session + Julian),
                      fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                      start = c(6.49, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.g0.rds")
#D.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.g0.rds")

######################################
###AIC of univariate density models###
######################################

D.univariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.Crops.g0.)

AIC(D.univariate.models)

##############################
###Bivariate Density models###
##############################

D.Crops.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + DDE, g0~session + Julian),
                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                          start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Crops.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.DDE.g0.rds")
#D.Crops.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.DDE.g0.rds")

D.Crops.DFE.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + DFE, g0~session + Julian),
                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                          start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Crops.DFE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.DFE.g0.rds")
#D.Crops.DFE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.DFE.g0.rds")

D.Crops.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + Roads, g0~session + Julian),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                            start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Crops.Roads.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.Roads.g0.rds")
#D.Crops.Roads.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.Roads.g0.rds")

D.Crops.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + PercentSlope, g0~session + Julian),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                            start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Crops.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.slope.g0.rds")
#D.Crops.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.slope.g0.rds")

D.Crops.precip.g0.<-secr.fit(McKcapt.comb, model = list(D~Crops + Precip, g0~session + Julian),
                             fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                             start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Crops.precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.precip.g0.rds")
#D.Crops.precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Crops.precip.g0.rds")

D.Roads.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + DFE, g0~session + Julian),
                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                          start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.DDE.g0.rds")
#D.Roads.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.DDE.g0.rds")

D.Roads.DFE.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + DDE, g0~session + Julian),
                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                          start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.DFE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.DFE.g0.rds")
#D.Roads.DFE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.DFE.g0.rds")

D.Roads.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + PercentSlope, g0~session + Julian),
                            fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                            start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.Slope.g0.rds")
#D.Roads.slope.g0.<-readRDS("C:/Users/nelsonj7/Box/secr/CombinedYears/Results/August2021/Telemetry/D.Roads.Slope.g0.rds")

D.Roads.precip.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + Precip, g0~session + Julian),
                             fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                             start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.precip.g0., "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.precip.g0.rds")
#D.Roads.precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.precip.g0.rds")

D.DFE.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + DDE, g0~session + Julian),
                        fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                        start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.DDE.g0.rds")
#D.DFE.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.DDE.g0.rds")

D.DFE.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads, g0~session + Julian),
                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                          start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.g0.rds")
#D.DFE.Roads.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.g0.rds")

D.DFE.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + PercentSlope, g0~session + Julian),
                          fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                          start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.slope.g0.rds")
#D.DFE.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.slope.g0.rds")

D.DFE.precip.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Precip, g0~session + Julian),
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                           start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.precip.g0.rds")
#D.DFE.precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.precip.g0.rds")

D.precip.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~Precip + PercentSlope, g0~session + Julian),
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                           start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.precip.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.slope.g0.rds")
#D.precip.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.slope.g0.rds")


D.precip.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~Precip + DDE, g0~session + Julian),
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                           start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.precip.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.DDE.g0.rds")
#D.precip.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.precip.DDE.g0.rds")


D.slope.DDE.g0.<-secr.fit(McKcapt.comb, model = list(D~PercentSlope + DDE, g0~session + Julian),
                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                           start = c(6.49, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.slope.DDE.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.DDE.g0.rds")
#D.slope.DDE.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.slope.DDE.g0.rds")

###############################
###Trivariate Density models###
###############################

D.DDE.DFE.Roads.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads, g0~session + Julian),
                                      fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                      start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.g0.rds")
#D.DDE.DFE.Roads.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.g0.rds")

D.DDE.DFE.Precip.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Precip, g0~session + Julian),
                               fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                               start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Precip.g0.rds")
#D.DDE.DFE.Precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Precip.g0.rds")

D.DDE.DFE.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + PercentSlope, g0~session + Julian),
                              fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                              start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.slope.g0.rds")
#D.DDE.DFE.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.slope.g0.rds")

D.DDE.DFE.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Crops, g0~session + Julian),
                              fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                              start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Crops.g0.rds")
#D.DDE.DFE.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Crops.g0.rds")

D.DFE.Roads.Precip.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + Precip, g0~session + Julian),
                                 fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                 start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.Precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.g0.rds")
#D.DFE.Roads.Precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.g0.rds")

D.DFE.Roads.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + PercentSlope, g0~session + Julian),
                                fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.slope.g0.rds")
#D.DFE.Roads.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.slope.g0.rds")

D.DFE.Roads.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + Crops, g0~session + Julian),
                                fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Crops.g0.rds")
#D.DFE.Roads.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Crops.g0.rds")

D.DFE.Precip.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Precip + Percent Slope, g0~session + Julian),
                                 fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                 start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")
###
saveRDS(D.DFE.Precip.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.slope.g0.rds")
#D.DFE.Precip.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.slope.g0.rds")

D.DFE.Precip.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Precip + Crops, g0~session + Julian),
                                 fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                 start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Precip.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.Crops.g0.rds")
#D.DFE.Precip.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.Crops.g0.rds")

D.DFE.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + PercentSlope + Crops, g0~session + Julian),
                                fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.slope.Crops.g0.rds")
#D.DFE.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.slope.Crops.g0.rds")

D.Roads.Precip.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + Precip + PercentSlope, g0~session + Julian),
                               fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                               start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Precip.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.slope.g0.rds")
#D.DFE.Precip.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.slope.g0.rds")

D.Roads.Precip.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + Precip + Crops, g0~session + Julian),
                                   fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                   start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.Precip.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.Precip.Crops.g0.rds")
#D.Roads.Precip.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.Precip.Crops.g0.rds")

D.Roads.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + PercentSlope + Crops, g0~session + Julian),
                                  fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                  start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.Precip.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.Precip.Crops.g0.rds")
#D.Roads.Precip.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Precip.Crops.g0.rds")

D.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~Precip + PercentSlope + Crops, g0~session + Julian),
                                   fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                   start = c(6.49, 2.3, 1.2, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Precip.slope.Crops.g0.rds")
#D.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Precip.slope.Crops.g0.rds")

###############################
###4-variable Density models###
###############################

D.DDE.DFE.Roads.Precip.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads + Precip, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.Precip.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.g0.rds")
#D.DDE.DFE.Roads.Precip.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.g0.rds")

D.DDE.DFE.Roads.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads + PercentSlope, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.slope.g0.rds")
#D.DDE.DFE.Roads.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.slope.g0.rds")

D.DDE.DFE.Roads.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads + Crops, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Crops.g0.rds")
#D.DDE.DFE.Roads.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Crops.g0.rds")

D.DDE.Roads.Precip.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + Roads + Precip + PercentSlope, g0~session + Julian),
                                       fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                       start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.Roads.Precip.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Roads.Precip.slope.g0.rds")
#D.DDE.Roads.Precip.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Roads.Precip.slope.g0.rds")

D.DDE.Roads.Precip.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + Roads + Precip + Crops, g0~session + Julian),
                                       fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                       start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.Roads.Precip.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Roads.Precip.Crops.g0.rds")
#D.DDE.Roads.Precip.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Roads.Precip.Crops.g0.rds")

D.DDE.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + Roads + Precip + PercentSlope, g0~session + Julian),
                              fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                              start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Precip.slope.Crops.g0.rds")
#D.DDE.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Precip.slope.Crops.g0.rds")

D.DFE.Roads.Precip.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + Precip + PercentSlope, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.Precip.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.slope.g0.rds")
#D.DFE.Roads.Precip.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.slope.g0.rds")

D.DFE.Roads.Precip.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + Precip + Crops, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.Precip.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.Crops.g0.rds")
#D.DFE.Roads.Precip.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.Crops.g0.rds")

D.DFE.Roads.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + PercentSlope + Crops, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.slope.Crops.g0.rds")
#D.DFE.Roads.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.slope.Crops.g0.rds")

D.DFE.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Precip + PercentSlope + Crops, g0~session + Julian),
                                      fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                      start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.slope.Crops.g0.rds")
#D.DFE.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Precip.slope.Crops.g0.rds")

D.Roads.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~Roads + Precip + PercentSlope + Crops, g0~session + Julian),
                                     fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                     start = c(6.49, 2.3, 1.2, -1, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.Roads.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.Precip.slope.Crops.g0.rds")
#D.Roads.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.Roads.Precip.slope.Crops.g0.rds")

###############################
###5-variable Density models###
###############################

D.DDE.DFE.Roads.Precip.slope.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads + Precip + PercentSlope, g0~session + Julian),
                                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                           start = c(6.49, 2.3, 1.2, -1, -3, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.Precip.slope.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.slope.g0.rds")
#D.DDE.DFE.Roads.Precip.slope.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.slope.g0.rds")

D.DDE.DFE.Roads.Precip.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads + Precip + Crops, g0~session + Julian),
                                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                           start = c(6.49, 2.3, 1.2, -1, -3, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.Precip.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.Crops.g0.rds")
#D.DDE.DFE.Roads.Precip.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.Crops.g0.rds")

D.DDE.Roads.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + Roads + Precip + PercentSlope + Crops, g0~session + Julian),
                                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                           start = c(6.49, 2.3, 1.2, -1, -3, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.Roads.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Roads.Precip.slope.Crops.g0.rds")
#D.DDE.Roads.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.Roads.Precip.slope.Crops.g0.rds")

D.DFE.Roads.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DFE + Roads + Precip + PercentSlope + Crops, g0~session + Julian),
                                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                           start = c(6.49, 2.3, 1.2, -1, -3, 4, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DFE.Roads.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.slope.Crops.g0.rds")
#D.DFE.Roads.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DFE.Roads.Precip.slope.Crops.g0.rds")

##########################
###Global Density Model###
##########################

D.DDE.DFE.Roads.Precip.slope.Crops.g0.<-secr.fit(McKcapt.comb, model = list(D~DDE + DFE + Roads + Precip + PercentSlope + Crops, g0~session + Julian),
                                           fixed = list(sigma = 1.855), mask = McKMask.comb, binomN = 1,
                                           start = c(6.49, 2.3, 1.2, -1, -3, 4, 2, -3.7, 0.9, -0.55), method = "Nelder-Mead")

saveRDS(D.DDE.DFE.Roads.Precip.slope.Crops.g0., file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.slope.Crops.g0.rds")
#D.DDE.DFE.Roads.Precip.slope.Crops.g0.<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/Results/2021August/D.DDE.DFE.Roads.Precip.slope.Crops.g0.rds")



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
summary(McKcapt.comb)
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