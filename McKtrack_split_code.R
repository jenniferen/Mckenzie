library(dplyr)

setwd("~/OSU/Thesis/Analysis/McKenzie")
data<-as_tibble(read.csv("McKTracksOnepersonpertrack.csv"))

data<-data %>% arrange(Order)

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

k<-7
i=1
transects<-unique(data$Transect)
out<-NULL

for(i in 1:length(transects)){

a<-subset(data, Transect==transects[i])
b<-chunk(a$Order, k)


c<-data.frame(Transect=rep(transects[i],k), 
              section.ID=c(1:k),
              start.Order=unlist(lapply(b,first)),
              end.Order=unlist(lapply(b,last)),
              med.Order=floor(unlist(lapply(b,median))),
              stringsAsFactors = FALSE
              )
out<-rbind(out, c)
              }

out$med.Lat<-data$Latitude[match(out$med.Order,data$Order)]
out$med.Lon<-data$Longitude[match(out$med.Order,data$Order)]
out$Date<-data$Date[match(out$med.Order,data$Order)]
head(out)
#-----------------------------------------------------------------------------
#Creating Trap matrix by assigning each trap a number and converting coordinates
#to UTM (zone 10)
#Adding Trap ID to dataframe and including 
out$TrapID<-c(1:length(out$Transect))
head(out)
#reordering
dim(out)
Traps<-out[,c(9,1,2,7,6,8)]
head(Traps)
library(sf)

df.SP<-st_as_sf(Traps, coords=c("med.Lon", "med.Lat"), crs=4326) #makes points spatial
df.SP<-st_transform(x=df.SP, crs=32610) #transforming from WGS84 to UTM zone 10, 
#can look up other CRS for other coordinate systems
df.SP$utm_E<-st_coordinates(df.SP)[,1] #get coordinates
df.SP$utm_N<-st_coordinates(df.SP)[,2] #get coordinates

#back to lat-long
df.SP<-st_transform(x=df.SP, crs=4326)
df.SP$lon<-st_coordinates(df.SP)[,1]#adding to dataset
df.SP$lat<-st_coordinates(df.SP)[,2]#adding to dataset

#add Julian Day
library(date)
df.SP$Julian<-julian.Date(df.SP$Date, origin=as.Date("2018-01-01"))+17531
df.SP$Julian<-scale(df.SP$Julian)

##add Effort
head(df.SP)
head(effortcov)


df.SP$effort <- effortcov$effort[match(df.SP$Transect, effortcov$transectsunique)]
head(df.SP)

#MCKENZIE DDE/Habitat Use TO TRAP COVARIATE
library(raster)
head(df.SP)
#DDEprojected is McKenzie DDE in UTMs

DDEex<-extract(DDEprojected, df.SP, buffer = 1000, fun=mean)
scaledDDEex<-scale(DDEex)
df.SP$DDE<-scaledDDEex

#TAPS is percent slope in UTMS
slopeex<-extract(TAPS, df.SP, method = 'simple')
scaledslopeex<-scale(slopeex) #Z-scaling values to center around 0 with 1 SD
df.SP$PercentSlope<-scaledslopeex

#Roads is distance to raods in meters
Roadsex<-extract(McKRoads, df.SP, method = 'simple')
scaledRoadsex<-scale(Roadsex) #Z-scaling values to center around 0 with 1 SD
df.SP$DistRoads<-scaledRoadsex

#DForestEdge is distance to forest edge in meters
DFEdgeex<-extract(DForestEdge, df.SP, method='simple')
scaledDFEdgeex<-scale(DFEdgeex)
df.SP$DFEdge<-scaledDFEdgeex

#Habitat Use
HabitatUseex<-extract(HabitatUse, df.SP, buffer = 1000, fun = mean)
scaledHabitatUseex<-scale(HabitatUseex)
df.SP$HabitatUse<-scaledHabitatUseex

df.SP<-st_set_geometry(df.SP, NULL) #makes non-spatial


traps<-df.SP
head(traps)
traps<-traps[,c(1,2,3,5,6,9,10,11,12,13,14,15)]
head(traps)
setwd("~/OSU/Thesis/Analysis/2018McKenzieGenotypes")
write.csv(traps, file="McKTraps.csv")
write.table(traps, file="McKTraps.txt")

#----------------------------------------------------------------
#map
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldTopoMap", group = "Topo") %>%
  addProviderTiles("Esri.WorldImagery", group = "ESRI Aerial") %>%
  addCircleMarkers(data=df.SP, group="Transects", radius = 4, opacity=1, fill = "darkblue",stroke=TRUE,
                   fillOpacity = 0.75, weight=2, fillColor = "yellow",
                   popup = paste0("Transect Name: ", df.SP$Transect)) %>%
  addLayersControl(
    baseGroups = c("Topo","ESRI Aerial"),
    overlayGroups = c("Hot SPrings"),
    options = layersControlOptions(collapsed = T))
