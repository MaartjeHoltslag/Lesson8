#Load required libraries
library(raster)
#Set working directory
setwd("/home/user/git/geoScripting/Lesson8")
dir.create("./data")
#Download required data
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/GewataB1.rda", destfile= "data/GewataB1.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/GewataB2.rda", destfile= "data/GewataB2.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/GewataB3.rda", destfile= "data/GewataB3.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/GewataB4.rda", destfile= "data/GewataB4.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/GewataB5.rda", destfile= "data/GewataB5.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/GewataB7.rda", destfile= "data/GewataB7.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/vcfGewata.rda", destfile= "data/vcfGewata.rda", method = "wget")
download.file(url = "https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/trainingPoly.rda", destfile= "data/trainingPoly.rda", method = "wget")
#Load required data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")
load("data/trainingPoly.rda")
#Remove cloud and water errors in the data
vcfGewata[vcfGewata > 100] <- NA
#Brick raster data
Gewata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
GewataBands <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7)
#Create data frame
GewataValues <- as.data.frame(getValues(Gewata))
#Plot Gewata pairs
pairs(Gewata)
#Create correlation model
correlationModel <- lm(vcf2000Gewata ~ gewataB1+ gewataB2 + gewataB3 + gewataB4 + gewataB5 + gewataB7, data = GewataValues) 
summary(correlationModel)
#Create predicted tree cover raster based on the correlation model
predictTreeC <- predict(GewataBands, model=correlationModel, na.rm=TRUE)
predictTreeC[predictTreeC < 0] <- NA
predictTreeCdf <- as.data.frame(getValues(predictTreeC))
plot(predictTreeC, main="Predicted tree cover", zlim=c(0,100))
#Compare predicted tree cover raster to vcf raster
pairs(brick(predictTreeC, vcfGewata))
#Calculate RMSE to compare predicted tree cover raster to vcf raster
RMSE <- (sqrt(mean((GewataValues$vcf2000Gewata-predictTreeCdf)^2, na.rm = T)))
#Calculate NDVI
ndvi <- overlay(GewataB3, GewataB4, fun=function(x,y){(y-x)/(x+y)})
#Make trainingPoly classes integers
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
#Create classes raster
classes <- rasterize(trainingPoly, ndvi, field='Code')
#Calculate RMSE for each class
difference <- (vcfGewata-predictTreeC)^2
TreeCovDif <- zonal(difference, classes)
RMSEclasses <- (sqrt(TreeCovDif[,2]))
RMSEclasses
