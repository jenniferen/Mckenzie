##########################################
####Build Capture History File############
##########################################

setwd("~/OSU/ElkSCR/McKenzie")

#Read in trap file
traps<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2019McKTraps.txt")
head(traps)
colnames(traps)<-c("TransectID", "X", "Y", "sep", "Effort", "Julian", "Precip", "Observers")
head(traps)

#Read in coordinates and sex of identified elk
McKElk2019<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/McKenzie/2019McKenzieSampleCoordinatesUTMIndividualsAdded.csv", header=TRUE)
head(McKElk2019)

#Adding Unknown to Sex column for individuals with unknown Sex
levels<-levels(factor(McKElk2019$Sex)) ###Making NAs Unknowns
McKElk2019$Sex<-factor(McKElk2019$Sex, levels = levels)
head(McKElk2019)

#Changing collection date to Julian day to gather precipitation data from PRISM script for Habitat Mask
McKElk2019$CollectionDate<-as.Date(McKElk2019$Collection.Date, "%m/%d/%y")

require(lubridate)
McKElk2019$CollectionDate = yday(McKElk2019$CollectionDate)
head(McKElk2019)

McKElk2019$CollectionDate<-McKElk2019$CollectionDate - 134
head(McKElk2019)

####Putting in Session, Animal ID, occassion, X, Y format
McKElk2019$Session<-rep("McK19", nrow(McKElk2019))
McKElk2019$Occasion<-rep(1, nrow(McKElk2019))
head(McKElk2019)
McKElk2019<-McKElk2019[,c(9,7,10,4,5,8)]
head(McKElk2019)
colnames(McKElk2019)<-c("Session", "AnimalID", "Occasion", "X", "Y", "Collection Date")

McKElk2019.NoDate<-McKElk2019[,c(1:5)]

##Convert to UTM
#coordinates(McKElk2019)<-McKElk2019[,c(5,4)]
#crs(McKElk2019)<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
#McKElk2019<-spTransform(McKElk2019, CRS("+proj=utm +zone=10 +datum=NAD83"))

#McKElk2019.UTM<-as.data.frame(McKElk2019)

#McKElk2019.UTM<-McKElk2019.UTM[,c(1,2,3,7,8,6)]
#colnames(McKElk2019.UTM)<-c("Session", "AnimalID", "Occasion", "X", "Y", "Collection Date")

write.csv(McKElk2019, "2019McKElk.csv", row.names=FALSE, quote = FALSE)

#Calculate distances between individuals and nearest trap
library(sf)
a<-st_as_sf(McKElk2019, coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83 
            +units=m +no_defs")
b<-st_as_sf(traps, coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83
            +units=m +no_defs")
closest<-list() #creates a blank list to store data
for(i in seq_len(nrow(a))){ #for each Elk in a
  closest[[i]]<-b[which.min( #write in each row of closest, b (Trap ID) which is the 
    st_distance(b, a[i,]) #minimum distance between all trap IDs and the current Elk, i.
  ),]
}

##Calculate distances between individuals
huzzah<-st_distance(a[-1,],a[-nrow(a),],by_element=TRUE)
write.csv(huzzah, file = "McK2019DistancesbetweenIndividualsmeters.csv", row.names = FALSE)
###

closestTrapID<-sapply(closest, "[[", 1) #Extracting TrapIDs from list
closestTrapID #viewing TrapIDs

head(McKElk2019)
McKElk2019$TrapID<-closestTrapID #creating column in McK Elk
head(McKElk2019)

McK2019.capthist<-McKElk2019[,c(1,2,3,7)]
names(McK2019.capthist)<-c("Session", "ID", "Occasion", "Detector")
head(McK2019.capthist)
write.csv(McK2019.capthist, file = "2019McKCH.csv", row.names = FALSE, quote = FALSE)

#Find this file, examine individuals and delete
#whenever an individual was detected at a trap more than once. Delete first column (doesn't have a header)
#Then save as a tab-deliminated txt. Open txt and add # and one space before Session.