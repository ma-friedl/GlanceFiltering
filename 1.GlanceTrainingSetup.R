######################################################
# GLanCE Training Data Analysis and Filtering Script #
######################################################

# First step in training data filtering.
# Ingest data, generate some basic plots;
# Screen out some basic problematic cases;
# Set up data for more involved filtering.

# Set up working dir
setwd("/projectnb/measures/users/rkstan/glance/")

# load required libraries
library(randomForest)
library(ggmap)
library(maptools)
library(rgdal)
library(corrplot)
library(caret)
library(ranger)
library(mlr)
library(MLmetrics)
library(e1071)
library(doMC)

# source script with function definitions
source('GlanceFiltering/0.GlanceFunctionDefs.R')

glancedat <- read.csv('training_data/NA_Training_Master_V1_2021_11_23_predictors_final.csv',header=TRUE)
#glancedat <- read.csv('training_data/NA_Training_Master_V1_NO_LCMAP_2021_07_09_predictors_final.csv', header=TRUE)

# select continent
cont <- 'NA'      # NA = North America, AF = Africa.....

# read training data from csv file
#if (cont=='AF') glancedat <- read.csv('data/AF_Training_Master_V1_2020_08_24_Coef_float.csv',header=TRUE)

#if (cont=='NA') {  
  # glancedat <- read.csv('data/NA_Training_Master_V1_NO_LCMAP_2021_03_17_predictors.csv',header=TRUE)[,2:142] # old code
  # horrible kludge to make this work with new colnames format.
#  c.names <- colnames(read.csv('training_data/NA_Training_Master_V1_NO_LCMAP_2021_07_09_predictors_final.csv',header=TRUE))
#  glancedat <- read.csv('training_data/NA_Training_Master_V1_NO_LCMAP_2021_07_09_predictors_final.csv',header=TRUE)[,c(c.names,"LC_Confidence")][,-1]
#}

## get rid of cases with 'unfilled (4) or snow/ice' (1) as LC class in AF
#if (cont == 'AF') {
#  glancedat <- subset(glancedat,glancedat[,'LC_Class']!='Unfilled' & glancedat[,'LC_Class']!='Snow/Ice')
#}

# replace missing confidence values with -99 and convert to factor
#glancedat[is.na(glancedat[,'LC_Confidence']),'LC_Confidence']=-999
#glancedat[,'LC_Confidence']=factor(glancedat[,'LC_Confidence'])

# map features to columns - varies by continent
#if (cont == 'AF') {
#  blue <- 3:17
#  green <- 25:39
#  red <- 65:79
 # nir <- 48:62
 # swir1 <- 80:94
 # swir2 <- 95:109
 # lst <- 112:126
 # topo <- c(2,20,23)
#  categorical <- c(18,19,46,47,127,128)
 # years <- c(110,24,130)
 # orphans <- c(21,22,24,40:43,63)
#  lc.class <- c(44:45)
 # climate <- c(64,111)
#}


blue <- grep("BLUE_*", colnames(glancedat))
green <- grep("GREEN_*", colnames(glancedat))
red <- grep("RED_*", colnames(glancedat))
nir <- grep("NIR_*", colnames(glancedat))
swir1 <- grep("SWIR1_*", colnames(glancedat))
swir2 <- grep("SWIR2_*", colnames(glancedat))
lst <- grep("TEMP_.", colnames(glancedat))
topo <- c(grep("ASPECT", colnames(glancedat)), grep("DEM_SLOPE", colnames(glancedat)), 
          grep("ELEVATION", colnames(glancedat)),
          grep("MIN_LSZA", colnames(glancedat)), grep("MAX_LSZA", colnames(glancedat)))
#categorical <- c(19,20,46,47)
years <- grep("*Year*", colnames(glancedat))
#datasetID <- 22
#orphans <- c(23,42,43,49,66,134:137)
lc.class <- grep("*Class*", colnames(glancedat))
#site.id <- c(44,45)
#ecoreg <- grep("*Ecoregion", colnames(glancedat))
#climate <- c(grep("TEMPERATURE", colnames(glancedat)), grep("RAINFALL", colnames(glancedat)))
climate <- c(grep("*Temp", colnames(glancedat)), grep("precip", colnames(glancedat)), grep("waterDef", colnames(glancedat)))
lat.lon <- c(grep("lat", colnames(glancedat)), grep("lon", colnames(glancedat)))
aux_vars <- c(grep("DEVELOPED", colnames(glancedat)), grep("WATER_OCCURRENCE", colnames(glancedat)),
              grep("recentMag", colnames(glancedat)))

###### NOTES: (1) MAPPING OF FEATURES TO COLS NEEDS TO BE UPDATED FOR AFRICA; 
#             (2) NUMERICAL SITE ID (COL 45) IS NOT UNIQUE TO EACH SITE
#                 ASSUME SAME SITE, DIFFERENT YEARS

# extract lat/long for each sample - format varies by continent
#if (cont == 'AF') {
#  geo <- glancedat[,'.geo']
#  nc <- nchar(geo)
#  lat.lon <- strsplit(substr(geo,32,(nc-2)),',')
#  lat.lon <- matrix(as.numeric(unlist(lat.lon)),ncol=2,byrow=T)[,c(2,1)]
#}

#if (cont == 'NA') {
#  lat.lon <- glancedat[,lat.lon]
#}

lat.lon <- glancedat[,lat.lon]

# assign row and col names to lat and lon data
rownames(lat.lon) <- rownames(glancedat)
colnames(lat.lon) <- c('Lat','Lon')

# Convert LC character string classes to integers for convenience
# LC1: 1=Bare; 2=Developed; 3=Forest; 4=Herbaceous; 
#      5=Shrub; 6=Snow/Ice; 7 = unfilled; 8 = Water
lc1.f <- factor(glancedat[,'LC_Class'])
table(lc1.f)  
lc1 <- as.integer(factor(glancedat[,'LC_Class']))
plot(lc1.f,main="Land Cover Level 1",col='lightgreen')

# LC2:  1 = Beach/Sand; 2=Deciduous; 3=Developed 4=Evergreen
#       5 = Grassland; 6 = Mixed;  7=Moss/Lichen; 8 = Rock
#       9 = Row Crops; 10=Shrub; 11= Snow/Ice; 12= Soil
#       13 = Unfilled; 14 = Water
#lc2.f <- factor(glancedat[,'LC_Class2'])
#table(lc2.f)
#lc2 <- as.integer(factor(glancedat[,'LC_Class2']))
#plot(lc2.f,main="Land Cover Level 2",col='lightblue',cex.lab=0.5)

# extract surface reflectance data 
sr.data <- glancedat[,c(blue,green,red,nir,swir1,swir2)]

# identify and get rid of cases witn NA's in sr features
sr.data <- na.omit(sr.data)

# and remove these cases from glancedat
glancedat <- glancedat[rownames(sr.data),]

# compute reflectances and NDVI from CCDC coefs 
blue.sr <- t(apply(glancedat[,blue],1,doPixSr))
green.sr <- t(apply(glancedat[,green],1,doPixSr))
red.sr <- t(apply(glancedat[,red],1,doPixSr))
nir.sr <- t(apply(glancedat[,nir],1,doPixSr))
swir1.sr <- t(apply(glancedat[,swir1],1,doPixSr))
swir2.sr <- t(apply(glancedat[,swir2],1,doPixSr))
lst.K <- t(apply(glancedat[,lst],1,doPixSr))*100
ndvi <- (nir.sr-red.sr)/(nir.sr+red.sr)

# now plot max values across all sites for each DOY - looking for out of range!
plot(apply(nir.sr,2,max),type='l',col='green',lwd=2,ylab='Max NIR')
plot(apply(red.sr,2,min),type='l',col='green',lwd=2,ylab='Min Red')
plot(apply(lst.K,2,max),type='l',col='green',lwd=2,ylab='Max LST')
plot(apply(lst.K,2,min),type='l',col='green',lwd=2,ylab='Min LST')

# filter cases where fitted SR data are <0 or >1, suggesting bad ccdc results
blue.oor <- getOor(blue.sr)
green.oor <- getOor(green.sr)
red.oor <- getOor(red.sr)
nir.oor <- getOor(nir.sr)
swir1.oor <- getOor(swir1.sr)
swir2.oor <- getOor(swir2.sr)

# put results all together in a single data frame, flag rows with no OOR values
oor.df <- data.frame(blue.oor,green.oor,red.oor,nir.oor,swir1.oor,swir2.oor)
n.oor.bypix <- rowSums(oor.df)
good.pix <- n.oor.bypix==0

# barplot class distribution for OOR cases
plot(lc1.f,main="OOR Land Cover",col='lightgreen')
plot(lc1.f[!good.pix],col='red',add=T)

# remove oor cases from glancedat and surface reflectances 
# (~472 cases in North America)
glancedat <- glancedat[good.pix,]
sr.data <- sr.data[good.pix,]
blue.sr <- blue.sr[good.pix,]
green.sr <- green.sr[good.pix,]
red.sr <- red.sr[good.pix,]
nir.sr <- nir.sr[good.pix,]
swir1.sr <- swir1.sr[good.pix,]
swir2.sr <- swir2.sr[good.pix,]
ndvi <- ndvi[good.pix,]
lst.K <- lst.K[good.pix,]

# compute mean time series for ndvi, red, nir, & swir, by class
ndvi.mn.forest <- getMeanSpec(LC.class='Forest',sr.dat=ndvi)
ndvi.mn.shrub <- getMeanSpec(LC.class='Shrub',sr.dat=ndvi)
ndvi.mn.herb <- getMeanSpec(LC.class='Herbaceous',sr.dat=ndvi)
ndvi.mn.bare <- getMeanSpec(LC.class='Bare',sr.dat=ndvi)
ndvi.mn.dev <- getMeanSpec(LC.class='Developed',sr.dat=ndvi)
ndvi.mn.water <- getMeanSpec(LC.class='Water',sr.dat=ndvi)

red.mn.forest <- getMeanSpec(LC.class='Forest',sr.dat=red.sr)
red.mn.shrub <- getMeanSpec(LC.class='Shrub',sr.dat=red.sr)
red.mn.herb <- getMeanSpec(LC.class='Herbaceous',sr.dat=red.sr)
red.mn.bare <- getMeanSpec(LC.class='Bare',sr.dat=red.sr)
red.mn.dev <- getMeanSpec(LC.class='Developed',sr.dat=red.sr)
red.mn.water <- getMeanSpec(LC.class='Water',sr.dat=red.sr)

nir.mn.forest <- getMeanSpec(LC.class='Forest',sr.dat=nir.sr)
nir.mn.shrub <- getMeanSpec(LC.class='Shrub',sr.dat=nir.sr)
nir.mn.herb <- getMeanSpec(LC.class='Herbaceous',sr.dat=nir.sr)
nir.mn.bare <- getMeanSpec(LC.class='Bare',sr.dat=nir.sr)
nir.mn.dev <- getMeanSpec(LC.class='Developed',sr.dat=nir.sr)
nir.mn.water <- getMeanSpec(LC.class='Water',sr.dat=nir.sr)

swir1.mn.forest <- getMeanSpec(LC.class='Forest',sr.dat=swir1.sr)
swir1.mn.shrub <- getMeanSpec(LC.class='Shrub',sr.dat=swir1.sr)
swir1.mn.herb <- getMeanSpec(LC.class='Herbaceous',sr.dat=swir1.sr)
swir1.mn.bare <- getMeanSpec(LC.class='Bare',sr.dat=swir1.sr)
swir1.mn.dev <- getMeanSpec(LC.class='Developed',sr.dat=swir1.sr)
swir1.mn.water <- getMeanSpec(LC.class='Water',sr.dat=swir1.sr)

# plot results for NDVI
plot(1:365,ndvi.mn.forest,
     type='l',
     col='darkgreen',
     lwd=2,
     ylim=c(-0.1,0.7),
     xlab='DOY',
     ylab='NDVI')
lines(1:365,ndvi.mn.shrub,col='green',lwd=2)
lines(1:365,ndvi.mn.herb,col='lightgreen',lwd=2)
lines(1:365,ndvi.mn.bare,col='brown',lwd=2)
lines(1:365,ndvi.mn.dev,col='red',lwd=2)
lines(1:365,ndvi.mn.water,col='blue',lwd=2)
title('Mean Modeled NDVI by LC')

# Generate EDA plots for each class, write as pdf to file.
class.list=unique(glancedat[,'LC_Class'])
mask <- class.list!="Snow/Ice"
class.list=class.list[mask]

pdf(file=paste(cont,'GlanceEDA.pdf',sep=''))
for (cls in class.list) {
  edaPlotsByClass(data=sr.data,
                  lcDat=glancedat[,'LC_Class'],
                  lcClass=cls)
}
dev.off()

# Now look at/flag extreme cases in features 
# First, set threshold for 'extreme' 
sd.thresh <- 3     # 3 SD from mean in each class = 0.27% of  data for N(0,1)

# now loop through each class in data set
filterData <- NULL
for (cls in class.list) {
  # apply 'filterFeatures' to data stratified by LC class
  clean.cls <- filterFeatures(data=sr.data,
                              lcDat=glancedat[,'LC_Class'],
                              lcClass=cls,
                              n.sd=sd.thresh)
  # filterData is matrix where all entries more than sd.thresh standard
  # deviations from mean replaced with NAs
  filterData <- rbind(filterData,clean.cls)
}

# Inspect results: 1. n extreme features in each training pixel? 
n.extreme.vals.byrow <- rowSums(is.na(filterData))

# plot number of extreme features in each training sample - slow to plot
barplot(n.extreme.vals.byrow,
        col='red',
        main="n Features > 3 SD from mean, each sample (by class)")

# histogram of # extreme features in each pixel
n.ex.row.hist <- hist(n.extreme.vals.byrow,breaks=0:max(n.extreme.vals.byrow),plot=F)
plot(n.ex.row.hist, 
     col='red',
     xlab='N Extreme Values',
     main="Distribution of Extreme Feature Values/Pixel in Training Data")

# cumulative distribution
plot(n.ex.row.hist$breaks[-1]-1,cumsum(n.ex.row.hist$density),
     type='l',
     col='red',
     xlab='N Extreme Values',
     ylab='CDF Extreme Values',
     main="Distribution of Extreme Feature Values/Pixel in Training Data")

# Inspect results: 2. how many extreme pixels in each feature? 
n.extreme.vals.bycol <- colSums(is.na(filterData))

barplot(n.extreme.vals.bycol,
        col='red',
        main="n Features > 3 SD from mean, each feature (by class)")

# Histogram of extreme cases for each feature
n.ex.col.hist <- hist(n.extreme.vals.bycol,nclass=6,plot=F)
plot(n.ex.col.hist, 
     col='red',
     xlab='N Extreme Values',
     main="Distribution of Extreme Values/Feature in Training Data")

# 33% worst offending features?
sort(n.extreme.vals.bycol)[64:96]
mean(n.extreme.vals.bycol[64:96])/dim(glancedat)[1]  # worst features have ~1.5% extreme vals

# now do cumulative distribution
n.ex.col.hist <- hist(n.extreme.vals.bycol,breaks=0:max(n.extreme.vals.bycol),plot=F)
plot(n.ex.col.hist$breaks[-1]-1,cumsum(n.ex.col.hist$density),
     type='l',
     col='red',
     xlab='N Extreme Values',
     ylab='CDF Extreme Values',
     ylim=c(0,1),
     main="Distribution of Extreme Feature Values/Feature in Training Data")

# Compute percentage of rows with n features more than sd.thresh SDs from mean
sum(n.extreme.vals.byrow>1)/dim(sr.data)[1]  # 22% of training cases have more than 1 feature > sd,thresh from mean
sum(n.extreme.vals.byrow>2)/dim(sr.data)[1]  # 16 of training cases have more than 2 features > sd.thresh from mean
sum(n.extreme.vals.byrow>3)/dim(sr.data)[1]  # 13% training cases have more than 3 features > sd.thresh from mean
sum(n.extreme.vals.byrow>4)/dim(sr.data)[1]  # 10% of training cases have more than 4 features > sd.thresh from mean
sum(n.extreme.vals.byrow>5)/dim(sr.data)[1]  # 8% of training cases have more than 5 features > sd.thresh from mean
sum(n.extreme.vals.byrow<1)/dim(sr.data)[1]  # 67% of training cases have no features > sd.thresh from mean

# let's look at classes
clean.cases <- n.extreme.vals.byrow<1        # look at pixels w/no extreme values
lc.clean.cases <- glancedat[clean.cases,'LC_Class']

# compare with original
plot(lc1.f,main="Land Cover Level 1, Original",col='lightgreen')
barplot(table(factor(lc.clean.cases)),
        ylab='Frequency',
        col='green',
        add=TRUE,
        main='After Removing Pixels with 1 Feature > Threshold SD')

# NOTE, AS OF RIGHT NOW, INFORMATION/RESULTS FROM THIS LAST SECTION NOT USED - PURELY EDA

