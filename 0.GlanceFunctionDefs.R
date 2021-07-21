############################  function defs  ###########################

edaPlotsByClass <- function(data=sr.data,lcDat=lc1,lcClass='Forest') {
  # generate EDA plots for surface reflectances, for requested class
  lc.sub <- subset(data,lcDat==lcClass)
  
  # generate boxplots
  boxplot(scale(lc.sub),cex=0.2,col='red',main=lcClass)
  for (i in 1:dim(data)[2]){
    hist(data[,i],nclass=100,
         main=paste(lcClass,colnames(data)[i]),
         col='green',
         xlab='Feature Data')
  }
  
  # plot PCA
  sr.pca <- princomp(data,cor=T)
  barplot(sr.pca$sdev,col='blue')
  plot(sr.pca$scores[,1],sr.pca$scores[,2],
       xlab='PC1',
       ylab='PC2',
       main=lcClass)
  lc.pc.sub=subset(sr.pca$scores,lcDat==lcClass)
  points(lc.pc.sub[,1],lc.pc.sub[,2],pch=16,col='red')
}

# doPixSr <- function(pix.coefs,bnd) {
#   int <- paste(bnd,'INTP',sep='_')
#   amp1 <- paste(bnd,'AMPLITUDE',sep='_')
#   amp2 <- paste(bnd,'AMPLITUDE2',sep='_')
#   amp3 <- paste(bnd,'AMPLITUDE3',sep='_')
#   phs1 <- paste(bnd,'PHASE',sep='_')
#   phs2 <- paste(bnd,'PHASE2',sep='_')
#   phs3 <- paste(bnd,'PHASE3',sep='_')
#   cs1 <- paste(bnd,'COS',sep='_')
#   cs2 <- paste(bnd,'COS2',sep='_')
#   cs3 <- paste(bnd,'COS3',sep='_')
#   sn1 <- paste(bnd,'SIN',sep='_')
#   sn2 <- paste(bnd,'SIN2',sep='_')
#   sn3 <- paste(bnd,'SIN3',sep='_')
#   slp <- paste(bnd,'SLP',sep='_')
#   doys <- 1:365
#   
#   sr <- pix.coefs[,int]+
#     pix.coefs[,cs1]*cos((doys+pix.coefs[,phs1])*2*pi/365)+
#     pix.coefs[,cs2]*cos((doys+pix.coefs[,phs2])*4*pi/365)+
#     pix.coefs[,cs3]*cos((doys+pix.coefs[,phs3])*6*pi/365)+
#     pix.coefs[,sn1]*sin((doys+pix.coefs[,phs1])*2*pi/365)+
#     pix.coefs[,sn2]*sin((doys+pix.coefs[,phs2])*4*pi/365)+
#     pix.coefs[,sn3]*sin((doys+pix.coefs[,phs3])*6*pi/365)
#   return(sr)
# }

doPixSr <- function(pix.coefs=glancedat[,blue]) {
  # compute reflectance for given band from CCDC coefs
  # do for every day of year
  doys <- 1:365
  sr <- pix.coefs[7]+
    pix.coefs[4]*cos((doys+pix.coefs[8])*2*pi/365)+
    pix.coefs[5]*cos((doys+pix.coefs[9])*4*pi/365)+
    pix.coefs[6]*cos((doys+pix.coefs[10])*6*pi/365)+
    pix.coefs[12]*sin((doys+pix.coefs[8])*2*pi/365)+
    pix.coefs[13]*sin((doys+pix.coefs[9])*4*pi/365)+
    pix.coefs[14]*sin((doys+pix.coefs[10])*6*pi/365)
  return(sr)
}

filterFeatures <- function(data=sr.data,lcDat=glancedat[,lc.class][,1],lcClass='Forest',n.sd=1.96){
  # replace cases in each band outside of n.sd standard deviations in given LC class 
  lc.sub <- subset(data,lcDat==lcClass)
  out.vals <- abs(scale(lc.sub)) > n.sd 
  barplot(colSums(out.vals),main=paste("N outliers, Each Feature",lcClass))
  lc.sub[out.vals] <- NA
  return(lc.sub)
}

plotSpecs <- function(line.n) {
  # annual time series of each band for requested line number
  plot(blue.sr[line.n,],type='l',col='blue',ylim=c(0,0.8),
       ylab='Reflectance',xlab='DOY',main=glancedat[line.n,'LC_Class'])
  lines(1:365,green.sr[line.n,],col='green')
  lines(1:365,red.sr[line.n,],col='red')
  lines(1:365,nir.sr[line.n,],col='black')
  lines(1:365,swir1.sr[line.n,],col='magenta')
  lines(1:365,swir2.sr[line.n,],col='purple')
  lines(1:365,ndvi[line.n,],col='red',lty=2)
  
}

getMeanSpec <- function(LC.class='Forest',sr.dat=blue.sr,lc.vec=glancedat[,'LC_Class']) {
  # compute mean annual reflectance time series for band, and class 
  mn.spec <- colMeans(subset(sr.dat,lc.vec==LC.class))
  return(mn.spec)
  
}

getrfMarg <- function(trueClass=featureDat[,1],mod=rfmodel.all) {
  # compute margins from estimated Random Forest model
  predClass <- predict(mod)
  missClass <- trueClass!=predClass
  probs <- predict(mod,type='prob')
  nclass <- dim(probs)[2]
  sort.probs <- t(apply(probs,1,sort))
  Margin <- sort.probs[,nclass]-sort.probs[,(nclass-1)]
  return(data.frame(trueClass,predClass,missClass,Margin))
}

getRangerMarg <- function(trueClass=featureDat[,1],
                          predClass=ranger.all$predictions,
                          probs=ranger.all.probs) {
  # estimate RF model and return margins
  missClass <- trueClass!=predClass
  nclass <- dim(probs)[2]
  sort.probs <- t(apply(probs,1,sort))
  Margin <- sort.probs[,nclass]-sort.probs[,(nclass-1)]
  return(data.frame(trueClass,predClass,missClass,Margin))
}

getOor <- function(dat=res.sr) {
  # flag cases where reflectance is out of range
  rw.mx <- apply(dat,1,max)
  rw.mn <- apply(dat,1,min)
  oor <- (rw.mx > 1 | rw.mn < 0)
  return(oor)
}

writePointsKmz <- function(lc.points=lat.lon,
                           Train.class=glancedat[,'LC_Class'],
                           Pred.class=glancedat[,'LC_Class'],
                           ID=glancedat[,'system.index'],
                           year=glancedat[,years],
                           outfil='LCpoints.kml') {
  # write kmz file with misclassified locations
  lc.points <- data.frame(lc.points,Train.class,Pred.class,year,ID)
  
  #Build a SpatialPointsData Frame
  coordinates(lc.points)<-c("Lon","Lat") 
  #names(lc.points) <- ID
  proj4string(lc.points)<-CRS("+proj=longlat +datum=WGS84")
  
  ##write SpatialPointsDataFrame to a KML File
  writeOGR(lc.points, dsn=outfil, layer= "lc.points", driver="KML",overwrite_layer = TRUE)
}

assignP <- function(clss){
  # assign prescribed prior likelihood for training site with land cover label = 'clss'
  if (clss == 'Bare') {
    pval <- 0.17
  } else if ( clss == 'Developed') {
    pval <- 0.15
  } else if (clss == 'Forest') {
    pval <- 0.08
  } else if (clss == 'Herbaceous') {
    pval <- 0.05
  } else if (clss == 'Shrub') {
    pval <- 0.10
  } else if (clss == 'Snow/Ice') {
    pval <- 0.15
  } else {
    pval <- 0.15
  }
  return(pval)
}

makeGrid <- function(xmin=-170,xmax=-50, ymin=15,ymax=80,grd.sz=10) {
  # function to estimate center point of regular grid in lat-lon
  # xmin, xmax, ymin, ymax = bounding box in degrees; 
  # increment determine size of grids (incr=5 c)
  
  # set up grid
  x.incr <- y.incr <- grd.sz/2
  if (xmin>xmax) x.incr <- -grd.sz/2
  if (ymin>ymax) y.incr <- -grd.sz/2
  x <- seq(xmin,xmax,x.incr)
  y <- seq(ymin,ymax,y.incr)
  
  # compute number of cells in grid
  ncells <- (length(x)-1)*(length(y)-1)
  
  # create array to store results
  grid <- array(data=NA,dim=c(ncells,1,2)) #,dimnames = c("Grid","Lon","Lat"))
  
  # initialize cell counter
  cell <- 0
  
  # loop through and assign grid center point
  for (i in 1:(length(x)-1)){
    for (j in 1:(length(y)-1)) {
      cell <- cell+1
      x.center <- (x[i]+x[i+1])/2
      y.center <- (y[j]+y[j+1])/2
      grid[cell,,] <- c(x.center,y.center) 
    }
  }
  return(grid) 
}

getSampsInWindow <- function(center.point,rng,latlon) {
  # get rownames for training samples in prescribed radius from center point
  grid.samps <- subset(latlon,latlon[,'Lat'] > (center.point[2] - rng) & 
                latlon[,'Lat'] < (center.point[2] + rng) &
                latlon[,'Lon'] > (center.point[1] - rng) &
                latlon[,'Lon'] < (center.point[1] + rng))
  # return rownames
  return(row.names(grid.samps))
    
}


filterECdat <- function(f.dat = goodTraining,
                        low.marg.thresh = 0.05,
                        pc.test=0.1,
                        dontSplit=FALSE) {
  # function to clean data set, by ecoregion, and return train and test set
  
  # convert LC label into factor
  f.dat[,1]=factor(f.dat[,1])
  
  # estimate probabilities for each case from so we can get margins
  ranger.all.probs <- ranger(LC_Class~.,
                             data=f.dat, 
                             num.trees=ntrees,
                             importance='none',
                             probability=TRUE)$predictions
  
  # predicted class based on probabilities
  ranger.probs.preds <- colnames(ranger.all.probs)[apply(ranger.all.probs,1,which.max)]
  
  # get margins: trueClass, predClass, missClass, Margin
  margins.all <- getRangerMarg(trueClass=f.dat[,1],
                               predClass=ranger.probs.preds,
                               probs=ranger.all.probs)
  
  # identify cases that are egregiously mis-classified
  egreg.high.thresh <- boxplot(margins.all[,'Margin']~margins.all[,'missClass'],plot=F)$stats[3,2] # upper 2 quartiles 
  egreg.high <- margins.all[,'missClass']==TRUE & margins.all[,'Margin'] > egreg.high.thresh
  
  # now subset cases that (1) are not egregious and (2) do not have excessively low margin
  f.dat <- subset(f.dat,!egreg.high & margins.all[,'Margin'] > low.marg.thresh)
  
  if (dontSplit==TRUE) {
    return(f.dat) 
  } else {
    # split into train and test
    test.smpl <- sample(1:dim(f.dat)[1],size=round(pc.test*dim(f.dat)[1]))
    train.indx <- !((1:dim(f.dat)[1]) %in% test.smpl)
    f.test.dat <- f.dat[test.smpl,]
    f.dat <- f.dat[train.indx,]
    
    # create list and return
    train.test <- NULL
    train.test$train <- f.dat
    train.test$test <- f.test.dat
    return(train.test)
  }
}

doEcClass <- function(ec.dat=featureDat) {
  ranger.ec <- ranger(LC_Class~.,
                      data=ec.dat,
                      num.trees=500,
                      importance='none',
                      probability=FALSE)
  return(ranger.ec)
}

plotOmCom <- function(cf.mat) {
  com.errs <- rowSums(cf.mat)-diag(cf.mat)
  om.errs <- colSums(cf.mat)-diag(cf.mat)
  plot(com.errs,om.errs,pch=16,col='red',xlab='Commission',ylab='Omission')
  abline(0,1)
  text(com.errs+1,om.errs+1,names(om.errs))
}


##########################  end function defs  ###########################
