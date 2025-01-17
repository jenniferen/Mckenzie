##########################################
####Build Capture History File############
##########################################

setwd("~/OSU/ElkSCR/McKenzie")

#Read in trap file
traps<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018McKTraps.txt")
head(traps)
colnames(traps)<-c("TrapID", "X", "Y", "sep", "Effort", "Julian", "Precip", "Observers")
head(traps)

#Read in coordinates and sex of identified elk
McKElk2018<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2018McKElkLocations.csv", header=TRUE)
head(McKElk2018)

#Changing collection date to Julian day to gather precipitation data from PRISM script for Habitat Mask
McKElk2018$CollectionDate<-as.Date(McKElk2018$Collection.Date, "%m/%d/%y")

require(lubridate)
McKElk2018$CollectionDate = yday(McKElk2018$CollectionDate)
head(McKElk2018)

McKElk2018$CollectionDate<-McKElk2018$CollectionDate - 133
head(McKElk2018)

####Putting in Session, Animal ID, occassion, X, Y format
McKElk2018$Session<-rep("McK18", nrow(McKElk2018))
McKElk2018$Occasion<-rep(1, nrow(McKElk2018))
head(McKElk2018)
McKElk2018<-McKElk2018[,c(8,2,9,4,3,7)]
head(McKElk2018)
colnames(McKElk2018)<-c("Session", "AnimalID", "Occasion", "X", "Y", "Collection Date")

McKElk2018.NoDate<-McKElk2018[,c(1:5)]

write.csv(McKElk2018, "2018McKElk.csv", row.names=FALSE, quote = FALSE)


#Calculate distances between individuals and nearest trap
library(sf)
projcrs <- "+proj=utm +zone=10 +datum=NAD83"

a<-st_as_sf(McKElk2018.NoDate, coords=c("X", "Y"), crs=st_crs(4326))
a<-st_transform(a, crs="+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ") #WGS 84 / UTM zone 10N
b<-st_as_sf(traps, coords=c("X", "Y"), crs="+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
closest<-list() #creates a blank list to store data
for(i in seq_len(nrow(a))){ #for each Elk in a
  closest[[i]]<-b[which.min( #write in each row of closest, b (Trap ID) which is the 
    st_distance(b, a[i,]) #minimum distance between all trap IDs and the current Elk, i.
  ),]
}

##Calculate distances between individuals
huzzah<-st_distance(a[-1,],a[-nrow(a),],by_element=TRUE)
write.csv(huzzah, file = "McK2018DistancesbetweenIndividualsmeters.csv", row.names = FALSE)
###

closestTrapID<-sapply(closest, "[[", 1) #Extracting TrapIDs from list
closestTrapID #viewing TrapIDs

#####
c<-c(0, huzzh
d<-cbind(c, huzzah)
test<-cbind(huzzah, closestTrapID)
###

head(McKElk2018)
McKElk2018$TrapID<-closestTrapID #creating column in McK Elk
head(McKElk2018)

McK2018.capthist<-McKElk2018[,c(1,2,3,7)]
names(McK2018.capthist)<-c("Session", "ID", "Occasion", "Detector")
head(McK2018.capthist)
write.csv(McK2018.capthist, file = "2018McKCH.csv", row.names = FALSE, quote = FALSE)

#Find this file, examine individuals and delete
#whenever an individual was detected at a trap more than once. Delete first column (doesn't have a header)
#Then save as a tab-deliminated txt. Open txt and add # and one space before Session.