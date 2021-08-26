library(dplyr)
library(raster)
library(sp)
library(sf)

setwd("~/OSU/ElkSCR/McKenzie")
data<-read.csv("McKTracksOnepersonpertrack.csv")
head(data)

###What do I want to do? Group each transect, calculate distances between points,
##choose distance to separate by and extract those points. Let's start with 200 m.
data.sp<-SpatialPointsDataFrame(coords = data[,c(9,8)], 
                                data = data,
                                proj4string = CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))

data.sp.utm<- spTransform(data.sp, CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"))

data.sp.utm$UTM<-coordinates(data.sp.utm)
data.sp.utm$X<-data.sp.utm$UTM[,1]
data.sp.utm$Y<-data.sp.utm$UTM[,2]
data.sp.utm$dist<-NA

#For each group of transects, Calculate distance between consecutive points, place 0 if first point.
for(i in 1:nrow(data.sp.utm)) {
  if(i > 1 && data.sp.utm$Transect[i] == data.sp.utm$Transect[i-1]) {
    data.sp.utm$dist[i] <- sqrt(((data.sp.utm$X[i] - data.sp.utm$X[i-1]) ^ 2) + ((data.sp.utm$Y[i] - data.sp.utm$Y[i-1]) ^ 2))
  } else {
    data.sp.utm$dist[i] <- NA
  }
}

data.sp.utm$dist.sum<-NA
for(i in 1:nrow(data.sp.utm)) {
  if(i > 1 && data.sp.utm$Transect[i] == data.sp.utm$Transect[i-1]) {
    data.sp.utm$dist.sum[i] <- data.sp.utm$dist.sum[i-1] + data.sp.utm$dist[i]
  } else {
    data.sp.utm$dist.sum[i] <- 0
  }
}

data.sp.utm$dist.sum
max(data.sp.utm$dist.sum)

#Now, I want to extract points that occur in increments of x-specified meters
x<-200

a <- data.frame(id = c("A","A","A","B","B","B"),
                b = c(1.2, 1.5, 1.8, 1.1, 1.6, 1.4))

seg.200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 200))), ]))
seg.400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 400))), ]))
seg.500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 500))), ]))
seg.600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 600))), ]))
seg.800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 800))), ]))
seg.1000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1000))), ]))
seg.1200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1200))), ]))
seg.1400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1400))), ]))
seg.1500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1500))), ]))
seg.1600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1600))), ]))
seg.1800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1800))), ]))
seg.2000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2000))), ]))
seg.2200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2200))), ]))
seg.2400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2400))), ]))
seg.2500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2500))), ]))
seg.2600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2600))), ]))
seg.2800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2800))), ]))
seg.3000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3000))), ]))
seg.3200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3200))), ]))
seg.3400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3400))), ]))
seg.3500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3500))), ]))
seg.3600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3600))), ]))
seg.3800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3800))), ]))
seg.4000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4000))), ]))
seg.4200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4200))), ]))
seg.4400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4400))), ]))
seg.4500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4500))), ]))
seg.4600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4600))), ]))
seg.4800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4800))), ]))
seg.5000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5000))), ]))
seg.5200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5200))), ]))
seg.5400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5400))), ]))
seg.5500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5500))), ]))
seg.5600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5600))), ]))
seg.5800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5800))), ]))
seg.6000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6000))), ]))

#200 meters
#seg.comb<-rbind(seg.200, seg.400, seg.600, seg.800, seg.1000, seg.1200, seg.1400, seg.1600, seg.1800, seg.2000, seg.2200,
#                seg.2400, seg.2600, seg.2800, seg.3000, seg.3200, seg.3400, seg.3600, seg.3800, seg.4000, seg.4200, seg.4400,
#                seg.4600, seg.4800, seg.5000, seg.5200, seg.5400, seg.5600, seg.5800, seg.6000)

#400 meters
#seg.comb<-rbind(seg.400, seg.800, seg.1200, seg.1600, seg.1800, seg.2200,
#                seg.2600, seg.3000, seg.3400, seg.3800, seg.4200,
#                seg.4600, seg.5000, seg.5400, seg.5800)

#500 meters
#seg.comb<-rbind(seg.500, seg.1000, seg.1500, seg.2000, seg.2500, seg.3000, seg.3500, seg.4000, seg.4500, seg.5000,
#                seg.5500, seg.6000)

#600 neters
seg.comb<-rbind(seg.600, seg.1200, seg.1800, seg.2400, seg.3000, seg.3600, seg.4200, seg.4800, seg.5600)

seg.comb.trim<-seg.comb[!duplicated(seg.comb$OBJECTID),]

#Create trap file
traps<-seg.comb.trim[,-c(2:9,11:12,15:19)]
traps<-traps[,c(1,3:4,2)]
head(traps)

#Need TrapID, X, Y, character, effort, Julian, precipitation, and # of observers.
TrapID<-seq(1:nrow(traps)) #Trap ID
traps$TrapID<-TrapID
traps$char<-rep("/", nrow(traps))

#effort covariate: want to take the max dist of each transect, scale and add to each point.
Effort<-NULL
Transects<-as.list(unique(seg.comb.trim$Transect))
for(i in Transects){
  sub<-subset(seg.comb.trim, seg.comb.trim$Transect == i)
  Effort$Transect[i]<- i
  Effort$dist[i]<-max(sub$dist.sum)
}

Effort$dist.scale<-scale(Effort$dist)

traps$effort <- Effort$dist.scale[match(traps$Transect, Effort$Transect)]

#Julian Day
Jday<-as.POSIXlt(traps$Date, format = "%m/%d/%y")
traps$JulianDay<-Jday$yday-136

#Precipitation
Precip.julian.end<-traps$JulianDay+15
Precip.julian.start<-traps$JulianDay+1
Precip.df<-as.data.frame(cbind(traps$X, traps$Y, Precip.julian.end, Precip.julian.start))
#
coordinates(Precip.df)<-Precip.df[,c(1,2)]
precip.cov<-rep(NA, nrow(Precip.df))
test.precip<-NULL
x2_proj.2018<-"C:/Users/jnelson/Documents/OSU/ElkSCR/PRISM/McKPrecip_2018_daily.tif"

for(i in 1:nrow(Precip.df)){
  test.precip<-calc(x2_proj.2018.crop[[Precip.df$Precip.julian.start[i]:Precip.df$Precip.julian.end[i]]], fun = mean, na.rm = T)
  precip.cov[i]<-extract(test.precip, Precip.df[i,], method = 'simple')
}

summary(precip.cov)
Precip.df
traps$precip<-scale(precip.cov)

#Number of observers
McK2018obs<-c(2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,
              1,1,1,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,
              1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,2,2,2,
              1,1,1,2,2,2,2,2,2,2,2,2,1,1,1)

obs.df<-as.data.frame(cbind(unlist(Transects), McK2018obs))
names(obs.df)<-c("Transect", "McK2018obs")
traps$Observers<-obs.df$McK2018obs[match(traps$Transect, obs.df$Transect)]
traps$Observers<-as.numeric(traps$Observers)
traps$Observers<-scale(traps$Observers)
head(traps)

#Organize and export traps
traps<-traps[,c(5,2,3,6,7:10)]
head(traps)

###Scale Julian Day
traps$JulianDay<-scale(traps$JulianDay)

write.table(traps, file = "2018McKTraps.txt", quote = FALSE, row.names = FALSE)
#Make sure to add # in front of column names in .txt file
