library(dplyr)
library(raster)
library(sp)
library(sf)

setwd("~/OSU/ElkSCR/McKenzie/")
data<-read.csv("2019McKenzieTracksCombined.csv")
head(data)

###What do I want to do? Group each transect, calculate distances between points,
##choose distance to separate by and extract those points. Let's start with 200 m.
data.sp<-SpatialPointsDataFrame(coords = data[,c(4:5)], 
                                data = data,
                                proj4string = CRS(as.character("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")))

data.sp$dist<-NA

#For each group of transects, Calculate distance between consecutive points, place 0 if first point.
for(i in 1:nrow(data.sp)) {
  if(i > 1 && data.sp$Transect[i] == data.sp$Transect[i-1]) {
    data.sp$dist[i] <- sqrt(((data.sp$UTM_X[i] - data.sp$UTM_X[i-1]) ^ 2) + ((data.sp$UTM_Y[i] - data.sp$UTM_Y[i-1]) ^ 2))
  } else {
    data.sp$dist[i] <- NA
  }
}

data.sp$dist.sum<-NA
for(i in 1:nrow(data.sp)) {
  if(i > 1 && data.sp$Transect[i] == data.sp$Transect[i-1]) {
    data.sp$dist.sum[i] <- data.sp$dist.sum[i-1] + data.sp$dist[i]
  } else {
    data.sp$dist.sum[i] <- 0
  }
}

data.sp$dist.sum
max(data.sp$dist.sum)

#Probelm transects that need subtraction: C31T2, C6T1
#C6T1: 1002:1020 - 1104
#C31T2: 4954:4981 -1500

data.sp$dist.sum[1002:1020]<-data.sp$dist.sum[1002:1020]-1104
data.sp$dist.sum[4954:4981]<-data.sp$dist.sum[4954:4981]-1500

#Now, I want to extract points that occur in increments of x-specified meters

seg.200<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 200))), ]))
seg.400<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 400))), ]))
seg.500<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 500))), ]))
seg.600<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 600))), ]))
seg.800<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 800))), ]))
seg.1000<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1000))), ]))
seg.1200<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1200))), ]))
seg.1400<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1400))), ]))
seg.1500<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1500))), ]))
seg.1600<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1600))), ]))
seg.1800<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1800))), ]))
seg.2000<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2000))), ]))
seg.2200<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2200))), ]))
seg.2400<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2400))), ]))
seg.2500<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2500))), ]))
seg.2600<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2600))), ]))
seg.2800<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2800))), ]))
seg.3000<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3000))), ]))
seg.3200<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3200))), ]))
seg.3400<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3400))), ]))
seg.3500<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3500))), ]))
seg.3600<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3600))), ]))
seg.3800<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3800))), ]))
seg.4000<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4000))), ]))
seg.4200<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4200))), ]))
seg.4400<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4400))), ]))
seg.4500<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4500))), ]))
seg.4600<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4600))), ]))
seg.4800<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4800))), ]))
seg.5000<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5000))), ]))
seg.5200<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5200))), ]))
seg.5400<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5400))), ]))
seg.5500<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5500))), ]))
seg.5600<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5600))), ]))
seg.5800<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5800))), ]))
seg.6000<-do.call('rbind', by(data.sp, data.sp$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6000))), ]))

#200 meters
#seg.comb<-rbind(seg.200, seg.400, seg.600, seg.800, seg.1000, seg.1200, seg.1400, seg.1600, seg.1800, seg.2000, seg.2200,
#                seg.2400, seg.2600, seg.2800, seg.3000, seg.3200, seg.3400, seg.3600, seg.3800, seg.4000, seg.4200, seg.4400,
#                seg.4600, seg.4800, seg.5000, seg.5200, seg.5400, seg.5600, seg.5800)

#400 meters
#seg.comb<-rbind(seg.400, seg.800, seg.1200, seg.1600, seg.1800, seg.2200,
#                seg.2600, seg.3000, seg.3400, seg.3800, seg.4200,
#                seg.4600, seg.5000, seg.5400, seg.5800)

#500 meters
#seg.comb<-rbind(seg.500, seg.1000, seg.1500, seg.2000, seg.2500, seg.3000, seg.3500, seg.4000, seg.4500, seg.5000,
#                seg.5500, seg.6000)

#600 meters
seg.comb<-rbind(seg.600, seg.1200, seg.1800, seg.2400, seg.3000, seg.3600, seg.4200, seg.4800, seg.5600)

#800 meters



#1000 meters
seg.comb<-
  
seg.comb.trim<-seg.comb[!duplicated(seg.comb$OBJECTID),]

#Create trap file
traps<-seg.comb.trim[,-c(1,3,7:15)]
head(traps)

#Need TrapID, X, Y, character, effort, Julian, precipitation, and # of observers.
TrapID<-seq(1:nrow(traps)) #Trap ID
traps$TrapID<-TrapID
traps$char<-rep("/", nrow(traps))
head(traps)
names(traps)<-c("Transect", "X", "Y", "Date", "TrapID", "char")

#effort covariate: want to take the max dist of each transect, scale and add to each point.
Effort<-NULL
Transects<-as.list(unique(seg.comb.trim$Transect))
for(i in Transects){
  sub<-subset(seg.comb.trim, seg.comb.trim$Transect == i)
  Effort$Transect[i]<- i
  Effort$dist[i]<-max(sub$dist.sum)
}

head(Effort)
Effort$dist.scale<-scale(Effort$dist)

traps$effort <- Effort$dist.scale[match(traps$Transect, Effort$Transect)]
head(traps)

#Julian Day
Jday<-as.POSIXlt(traps$Date, format = "%m/%d/%y")
traps$JulianDay<-Jday$yday-131
head(traps)

#Precipitation
Precip.julian.end<-traps$JulianDay+13
Precip.julian.start<-traps$JulianDay
Precip.df<-as.data.frame(cbind(traps$X, traps$Y, Precip.julian.end, Precip.julian.start))
#
coordinates(Precip.df)<-Precip.df[,c(1,2)]
precip.cov<-rep(NA, nrow(Precip.df))
test.precip<-NULL

for(i in 1:nrow(Precip.df)){
  test.precip<-calc(x2_clip.proj.2019[[Precip.df$Precip.julian.start[i]:Precip.df$Precip.julian.end[i]]], fun = mean, na.rm = T)
  precip.cov[i]<-extract(test.precip, Precip.df[i,], method = 'simple')
}

traps$precip<-scale(precip.cov)

#Number of observers
McK2019obs<-c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2,2,2,1,2,3,2,2,2,2,2,2,3,3,3,2,1,2,2,2,3,3,3,2,2,
                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,1,1,1,2,
                2,2,1,1,1,1,1,1,2,2,2,3,3,3)

obs.df<-as.data.frame(cbind(unlist(Transects), McK2019obs))
names(obs.df)<-c("Transect", "McK2019obs")
traps$Observers<-obs.df$McK2019obs[match(traps$Transect, obs.df$Transect)]
traps$Observers<-as.numeric(traps$Observers)
traps$Observers<-scale(traps$Observers)
head(traps)

#Organize and export traps
traps<-traps[,c(5,2,3,6,7:10)]
head(traps)

###Scale Julian Day
traps$JulianDay<-scale(traps$JulianDay)

write.table(traps, file = "2019McKTraps.txt", quote = FALSE, row.names = FALSE)
#Make sure to add # in front of column names in .txt file
