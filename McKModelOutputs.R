######################################
##Top Model Calculations and Outputs##
######################################

# CALCULATE chat

derived(top.model) #Use these outputs below

#2018 c-hat
D <- 0.06792655 #density
g0 <- -3.3599118 #
sigma <- 6.9732202
fn <- function (r) r * g0 * exp(-r^2/2/(sigma/100)^2)
expected.nk <- 2 * pi *D * integrate(fn, 0, 5*sigma/100)$value
summary(Tiogacapt.comb.telem)
nk<-122 #realized number of detections

X2 <- sum((nk - expected.nk)^2 / expected.nk)
si <- mean((nk - expected.nk) / expected.nk)
nu <- 121 #degrees of freedom
chat.2018 <- X2/nu / (1 + si) #Fletcher 2012
#out$dispersion <- c(mean.nk = mean(nk), var.nk = sd(nk)^2, chat = chat)
chat.2018

#2019 c-hat
D <- 0.09298142 #density
g0 <- -3.3599118 #probability of detection
sigma <- 6.9732202
fn <- function (r) r * g0 * exp(-r^2/2/(sigma/100)^2)
expected.nk <- 2 * pi *D * integrate(fn, 0, 5*sigma/100)$value
summary(Tiogacapt.comb.telem)
nk<-167
#nk <- apply(apply(ch2,c(1,3),sum)>0,2,sum) #realized number of detections
X2 <- sum((nk - expected.nk)^2 / expected.nk)
si <- mean((nk - expected.nk) / expected.nk)
nu <- 166 #degrees of freedom
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

topmodelMcK<-readRDS("C:/Users/nelsonj7/Box/secrMcK/CombinedYears/D.Crops.slope.g0.rds")
topmodelDsurface<-predictDsurface(topmodelMcK, mask = McKMask.comb)
plot(topmodelDsurface, col = brewer.pal(8, "YlGn"))
library(RColorBrewer)

McK18.Dsurface<-raster(topmodelDsurface$McK18, covariate = 'D.0', crs = "+proj=utm +zone=10 +datum=NAD83")
McK19.Dsurface<-raster(topmodelDsurface$McK19, covariate = 'D.0', crs = "+proj=utm +zone=10 +datum=NAD83")
writeRaster(McK18.Dsurface, file = "C:/Users/nelsonj7/Box/secrMcK/CombinedYears/McK18.Dsurface.tif", overwrite = TRUE)
writeRaster(McK19.Dsurface, file = "C:/Users/nelsonj7/Box/secrMcK/CombinedYears/McK19.Dsurface.tif", overwrite = TRUE)

########
##run in parallel
#######
#detection.fits <- par.secr.fit (
#  c('secr.null','secr.t', 'secr.bk','secr.sd', 'secr.bk.sd','secr.t.bk','secr.bait.bk','secr.dna.bk'),
#  ncores = 8) ##755.575 minutes 