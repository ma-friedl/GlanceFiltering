######################################################
# GLanCE Training Data Filtering Script              #
######################################################

###################  PART 1  ######################### 
# Ingest data, screen some basic problematic cases   #
# Set up data for more involved filtering            #
######################################################

# Set up working dir and functions
setwd("~/Dropbox (BOSTON UNIVERSITY)/Main/Rwork/GLanCETraining")

# Source script with function definitions
source('Rscripts/0.GlanceFunctionDefs.R')

# load required libraries
library(randomForest)
library(ranger)

# Parameters for script
cont <- 'NA'      # NA = North America, AF = Africa.....
inputFile <- 'data/NA_Training_Master_V1_2021_11_23_predictors_final.csv'
ntrees=1000
low.margin.thresh <- 0.05
quartile.thresh <- 2      # 2=exclude upper 2 quartiles; 3 = exclude uppermost quartile
n.regions <- 9            # Number of clusters/regions

# Read training data
glancedat <- read.csv(inputFile,header=TRUE)

# Extract bands and features
blue <- grep("BLUE_*", colnames(glancedat))
green <- grep("GREEN_*", colnames(glancedat))
red <- grep("RED_*", colnames(glancedat))
nir <- grep("NIR_*", colnames(glancedat))
swir1 <- grep("SWIR1_*", colnames(glancedat))
swir2 <- grep("SWIR2_*", colnames(glancedat))
lst <- grep("TEMP_.", colnames(glancedat))
topo <- c(grep("ASPECT", colnames(glancedat)), 
          grep("DEM_SLOPE", colnames(glancedat)), 
          grep("ELEVATION", colnames(glancedat)),
          grep("MIN_LSZA", colnames(glancedat)), 
          grep("MAX_LSZA", colnames(glancedat)))
years <- grep("*Year*", colnames(glancedat))
lc.class <- grep("*Class*", colnames(glancedat))
Level1.class <- grep("*LC_Class*", colnames(glancedat))
climate <- c(grep("*Temp", colnames(glancedat)), 
             grep("precip", colnames(glancedat)), 
             grep("waterDef", colnames(glancedat)))
Lat.Lon <- c(grep("lat", colnames(glancedat)), 
             grep("lon", colnames(glancedat)))
aux_vars <- c(grep("DEVELOPED", colnames(glancedat)), 
              grep("WATER_OCCURRENCE", colnames(glancedat)),
              grep("recentMag", colnames(glancedat)))

# assign row and col names to lat and lon data
lat.lon <- glancedat[,Lat.Lon]
rownames(lat.lon) <- rownames(glancedat)
colnames(lat.lon) <- c('Lat','Lon')

# extract surface reflectance data 
sr.data <- glancedat[,c(blue,green,red,nir,swir1,swir2)]

# identify and get rid of cases with NA's in sr features
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

# remove oor cases from glancedat and surface reflectances 
# (~a few hundred cases in North America)
glancedat <- glancedat[good.pix,]

###################  PART 2  ######################### 
# Filter training data based on Random Forest        #
######################################################

# First build climate/geographic regions via clustering
# select climate data, scale th  data, and add lat/long
clust.dat <- scale(glancedat[,c('maxTemp','minTemp','precip')])
clust.dat <- data.frame(glancedat[,'lon'],glancedat[,'lat'],clust.dat)

# estimate using hierarchical clustering 
# d.mat <- dist(clust.dat)    # compute distance matrix
# clust.full <- hclust(d.mat,method='ward.D2')
# region.labels <- cutree(clust2.full, k=n.regions)

# estimate using kmeans
region.labels <- kmeans(clust.dat,centers=n.regions,iter.max=100)$cluster

# subset LC class and features
featureDat <- glancedat[,c(Level1.class,blue,green,red,nir,swir1,swir2,lst,topo,climate)]

# convert LC class to a factor for classification
featureDat[,1] <- factor(featureDat[,'LC_Class'])

# remove sin/cos coefs from features
cos.sins <- c(grep("COS", colnames(glancedat)),
              grep("SIN", colnames(glancedat)))
good.cols <- !((1:dim(featureDat)[2]) %in% cos.sins)
featureDat <- featureDat[,good.cols]
nfeat <- dim(featureDat)[2]-1

# create empty object to store filtered rows
clean.rows <- NULL

# Filter data in each region
for (region in 1:n.regions) {
  # subset data for region
  featureDat.subReg <- subset(featureDat,region.labels==region)
  print(paste("Filtering Region",region, ":", dim(featureDat.subReg)[1],"Cases",sep=" "))
  
  # estimate and plot variable importance
  ranger.perm <- ranger(LC_Class~.,data=featureDat.subReg,num.trees=ntrees,importance='permutation')
  
  # select 50 most important features - based on trial and error - not super sensitive
  best.features <- names(sort(ranger.perm$variable.importance))[(nfeat-49):nfeat]
  featureDat.subReg <- featureDat.subReg[,c('LC_Class',best.features)]
  
  # Estimate base model 
  ranger.base <- ranger(LC_Class~.,
                        data=featureDat.subReg,
                        num.trees=ntrees,
                        importance='none',
                        probability=FALSE)
  names(ranger.base$predictions) <- rownames(featureDat.subReg)
  
  ########## Filtering Step 1: get rid of bat-shit crazy confusion cases  ###########
  bad.pix <-(featureDat.subReg[,'LC_Class'] == 'Bare' & ranger.base$predictions == 'Forest') |
    (featureDat.subReg[,'LC_Class'] == 'Bare' & ranger.base$predictions == 'Water') |
    (featureDat.subReg[,'LC_Class'] == 'Developed' & ranger.base$predictions == 'Water') |
    (featureDat.subReg[,'LC_Class'] == 'Developed' & ranger.base$predictions == 'Snow/Ice') |
    (featureDat.subReg[,'LC_Class'] == 'Forest' & ranger.base$predictions == 'Bare') |
    (featureDat.subReg[,'LC_Class'] == 'Forest' & ranger.base$predictions == 'Water') |
    (featureDat.subReg[,'LC_Class'] == 'Forest' & ranger.base$predictions == 'Herbaceous') |
    (featureDat.subReg[,'LC_Class'] == 'Herbaceous' & ranger.base$predictions == 'Forest') |
    (featureDat.subReg[,'LC_Class'] == 'Herbaceous' & ranger.base$predictions == 'Water') |
    (featureDat.subReg[,'LC_Class'] == 'Herbaceous' & ranger.base$predictions == 'Snow/Ice') |
    (featureDat.subReg[,'LC_Class'] == 'Shrub' & ranger.base$predictions == 'Water') |
    (featureDat.subReg[,'LC_Class'] == 'Shrub' & ranger.base$predictions == 'Snow/Ice') |
    (featureDat.subReg[,'LC_Class'] == 'Snow/Ice' & ranger.base$predictions != 'Snow/Ice') |
    (featureDat.subReg[,'LC_Class'] == 'Water' & ranger.base$predictions != 'Water') 
  
  # remove BSC cases from training data
  featureDat.subReg <- featureDat.subReg[!bad.pix,]
  
  # and re-estimate model with new subsetted training data
  ranger.sub1 <- ranger(LC_Class~.,
                        data=featureDat.subReg,
                        num.trees=ntrees,
                        importance='none',
                        probability=FALSE)
  names(ranger.sub1$predictions) <- rownames(featureDat.subReg)
  
  # Estimate ranger in probability mode to so we can get margins
  ranger.sub1.probs <- ranger(LC_Class~.,
                              data=featureDat.subReg, 
                              num.trees=ntrees,
                              importance='none',
                              probability=TRUE)$predictions
  ranger.probs.preds <- colnames(ranger.sub1.probs)[apply(ranger.sub1.probs,1,which.max)]
  
  # compute margin on first, second most likely class, predicted class, true class
  margins.all <- getRangerMarg(trueClass=featureDat.subReg[,1],
                               predClass=ranger.sub1$predictions,
                               probs=ranger.sub1.probs)
  
  # get upper value on upper 1/2 of margins for misclassified!
  egreg.thresh.high <- boxplot(margins.all[,'Margin']~margins.all[,'missClass'],plot=F)$stats[(quartile.thresh+1),2] # upper quartiles 
  
  ##### Filtering step 2: remove based on egregious & low margin cases  ####
  # split training data into egregious/low margin versus non-egregious sufficient margin
  bad.pix <- (margins.all[,'missClass'] & (margins.all[,'Margin'] > egreg.thresh.high)) | # egregiously misclassified
    (margins.all[,'Margin'] < low.margin.thresh)                                          # indistinguishable
  
  # Remove "bad pixels"
  goodTraining <- featureDat.subReg[!bad.pix,]
  clean.rows <- c(clean.rows,rownames(goodTraining))
  
}

# Subset data set based on filtering
clean.dat.allfeatures <- glancedat[clean.rows,]

# write to file
outfile <- paste(cont,"_Filtered.csv")
write.csv(clean.dat.allfeatures,file=outfile,row.names=FALSE)


