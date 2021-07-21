##########################################
# Script to filter GLanCE training data  #
##########################################

# run set up in read and pre-process GLanCE training data
source('Rscripts/1.GlanCETrainingSetup.R')

# subset LC class and features
featureDat <- glancedat[,c(lc.class[1],blue,green,red,nir,swir1,swir2,lst,topo,climate)]

# convert LC class to a factor for classification
featureDat[,1] <- factor(featureDat[,'LC_Class'])

# remove sin/cos coefs from features
cos.sins <- c(5:7,13:15,
              21:23,29:31,
              37:39,45:47,
              53:55,61:63,
              69:71,77:79,
              85:87,93:95,
              101:103,109:111)
good.cols <- !((1:dim(featureDat)[2]) %in% cos.sins)
featureDat <- featureDat[,good.cols]
nfeat <- dim(featureDat)[2]-1

# save a version of original - need to be able to restore later
featureDat.Base <- featureDat

# set number of trees for Random Forest
ntrees=500

# set number of classes
nclass <- length(unique(featureDat[,1]))

# estimate and plot variable importance
ranger.perm <- ranger(LC_Class~.,data=featureDat,num.trees=ntrees,importance='permutation')
barplot(sort(ranger.perm$variable.importance),horiz = TRUE,las=1,col='red',cex.names=0.5)
title('Variable Importance from Ranger')

# select 50 most important features - based on trial and error - not super sensitive
best.features <- names(sort(ranger.perm$variable.importance))[(nfeat-49):nfeat]
featureDat <- featureDat[,c('LC_Class',best.features)]

# Estimate base model 
ranger.base <- ranger(LC_Class~.,
                     data=featureDat,
                     num.trees=ntrees,
                     importance='none',
                     probability=FALSE)
names(ranger.base$predictions) <- rownames(featureDat)

# confusion matrix and accuracy
ranger.base$confusion.matrix
1-ranger.base$prediction.error

########## Filtering Step 1: get rid of bat-shit crazy confusion cases  ###########
bad.pix <- (featureDat[,'LC_Class'] == 'Bare' & ranger.base$predictions == 'Forest') |
  (featureDat[,'LC_Class'] == 'Bare' & ranger.base$predictions == 'Water') |
  (featureDat[,'LC_Class'] == 'Developed' & ranger.base$predictions == 'Water') |
  (featureDat[,'LC_Class'] == 'Forest' & ranger.base$predictions == 'Bare') |
  (featureDat[,'LC_Class'] == 'Forest' & ranger.base$predictions == 'Water') |
  (featureDat[,'LC_Class'] == 'Herbaceous' & ranger.base$predictions == 'Forest') |
  (featureDat[,'LC_Class'] == 'Herbaceous' & ranger.base$predictions == 'Water') |
  (featureDat[,'LC_Class'] == 'Shrub' & ranger.base$predictions == 'Water') |
  (featureDat[,'LC_Class'] == 'Snow/Ice' & ranger.base$predictions == 'Forest') |
  (featureDat[,'LC_Class'] == 'Water' & ranger.base$predictions == 'Forest') 

# look at confidence of bad.dat pixels
plot(glancedat[rownames(featureDat[bad.pix,]),'LC_Confidence'],col='green')
title('Distribution of Confidence in "Bad Pixels" ')

# Remove bad pixels 
featureDat <- featureDat[!bad.pix,]
print(paste('Deleted',dim(featureDat.Base)[1]-dim(featureDat)[1],'Training Cases'))

# look at change in distribution
barplot(table(featureDat.Base[,'LC_Class']),col='Green', main="Training Dist Before vs After Removing BS Cases")      
barplot(table(featureDat[,'LC_Class']),col='blue',add=T)

# and re-estimate model with new subsetted training data
ranger.sub1 <- ranger(LC_Class~.,
                     data=featureDat,
                     num.trees=ntrees,
                     importance='none',
                     probability=FALSE)
names(ranger.sub1$predictions) <- rownames(featureDat)

# look at and save confusion matrix and accuracy
ranger.sub1$confusion.matrix
1-ranger.sub1$prediction.error
ranger.sub1.cf <- ranger.sub1$confusion.matrix
ranger.sub1.acc <- 1-ranger.sub1$prediction.error

# Estimate ranger in probability mode to estimate margins
ranger.sub1.probs <- ranger(LC_Class~.,
                           data=featureDat, 
                           num.trees=ntrees,
                           importance='none',
                           probability=TRUE)$predictions
ranger.probs.preds <- colnames(ranger.sub1.probs)[apply(ranger.sub1.probs,1,which.max)]

# now look at accuracy and confusion matrix from RF in probability mode
ranger.probs.cf <- table(featureDat[,'LC_Class'],ranger.probs.preds)
ranger.probs.acc <- sum(diag(ranger.probs.cf))/sum(ranger.probs.cf)
ranger.probs.cf
ranger.probs.acc

# confusion matrix for classification vs probability model
class.prob.cf <- table(ranger.sub1$predictions,ranger.probs.preds)
class.prob.cf
sum(diag(class.prob.cf))/sum(class.prob.cf) # 90% agreement (10% disagreement!)

# compute margin on first, second most likely class, predicted class, true class
margins.all <- getRangerMarg(trueClass=featureDat[,1],
                             predClass=ranger.sub1$predictions,
                             probs=ranger.sub1.probs)

# look at result
boxplot(margins.all[,'Margin']~margins.all[,'missClass'],
        xlab='Misclassified',
        ylab='Margin',
        col='red')

# let's look at mis-classified
misclassified.pix <- subset(margins.all,margins.all[,'missClass'])  # ~3700 pixels
hist(misclassified.pix[,'Margin'],
     col='red',
     xlab='Margin',
     main='Freq Distribution for Margin on Misclassified Cases')

# look at margins across classes
boxplot(misclassified.pix[,'Margin']~misclassified.pix[,'trueClass'],
        ylab='Margin',
        xlab='',
        col='red',
        main='Margins on Misclassified: Training Label Class')

boxplot(misclassified.pix[,'Margin']~misclassified.pix[,'predClass'],
        ylab='Margin',
        xlab='',
        col='red',
        main='Margins on Misclassified: Predicted Label Class')

# get upper value on upper 1/2 of margins for misclassified!
egreg.thresh.high <- boxplot(margins.all[,'Margin']~margins.all[,'missClass'],plot=F)$stats[3,2] # upper 2 quartiles 
egreg <- subset(misclassified.pix,misclassified.pix[,'Margin']>egreg.thresh.high)                # egregious = margin > 0.19

# now lets look at how different subsets look
egreg.rows <- rownames(egreg)
egreg.pix <- data.frame(lat.lon[egreg.rows,],
                        glancedat[egreg.rows,'LC_Class'],
                        ranger.sub1$predictions[egreg.rows],
                        glancedat[egreg.rows,'Dataset'])
colnames(egreg.pix) <- c('Lat','Lon','Train.Class','Pred.Class','Source')

barplot(table(as.factor(glancedat[,'Dataset'])),
        col='green',
        main='Distribution of Egregious Errors wrt Source')
barplot(table(as.factor(egreg.pix[,'Source'])),add=T,col='red')

# look at distribution of egregious errors across confidence
all.cf <- table(glancedat[,'LC_Confidence'])
eg.cf <- table(glancedat[egreg.rows,'LC_Confidence'])
eg.cf/all.cf

plot(glancedat[,'LC_Confidence'],col='green')
plot(glancedat[egreg.rows,'LC_Confidence'],add=T,col='red')
title('Confidence: Full Training Data vs Egregious Mislabels')
barplot(eg.cf/all.cf,col='red',main='Egregious Mislabels As Proportions of Each Confidence Level')


# now look at patterns in each source of training data

# ABoVE
abv.sub <- subset(egreg.pix,egreg.pix[,'Source']=='ABoVE')
table(abv.sub[,3:4])

# Clustering
clust.sub <- subset(egreg.pix,egreg.pix[,'Source']=='CLUSTERING')
table(clust.sub[,3:4])

# STEP
step.sub <- subset(egreg.pix,egreg.pix[,'Source']=='STEP')
table(step.sub[,3:4])

# AUGMENT
aug.sub <- subset(egreg.pix,egreg.pix[,'Source']=='Training_augment')
table(aug.sub[,3:4])

# Generate margins, but include class predicted by RF in probability mode
margins.all.p <- getRangerMarg(trueClass=featureDat[,1],
                               predClass=ranger.probs.preds,
                               probs=ranger.sub1.probs)

# combine results for margins for both classifiers into a single data frame
margins.comp <- data.frame(margins.all,margins.all.p)[,c(1,2,6,3,7,8)]
colnames(margins.comp) <- c('True.Class','RF.Class.Pred','RF.Prob.Pred','RF.Class.MissClass', 'RF.Prob.MissClass','Margin')

# look at margins for cases where standard RF vs RF in probability mode disagree
preds.probs.delt=subset(margins.comp,margins.comp[,'RF.Class.Pred']!=margins.comp[,'RF.Prob.Pred'])
hist(preds.probs.delt[,'Margin'],
     breaks=seq(0,0.75,0.01),
     xlab='Margin',
     col='red',
     main='Histogram of Margins Where RF.p vs RF.c Disagree n~ 800')

# overall confusion matrix
egreg.cf <- table(true=egreg[,1],pred=egreg[,2])         # confusion matrix for egregious errors;
n <- sum(egreg.cf)
egreg.cf <- rbind(egreg.cf,colSums(table(egreg[,1:2])))  
egreg.cf <- cbind(egreg.cf,c(rowSums(table(egreg[,1:2])),n))  # rows = labels; cols = RF preds

# plot frequency of commission error classes with freq in underlying class data
com.errs <- colSums(table(egreg[,1:2]))
om.errs <- rowSums(table(egreg[,1:2]))
barplot(table(lc1.f),col='green',main='Commission Errors - Egregious Cases')
barplot(com.errs,add=T,col='red')
barplot(table(lc1.f),col='green',main='Omission Errors - Egregious Cases')
barplot(om.errs,add=T,col='blue')

# examine balance in commission vs omission
plot(com.errs,om.errs,pch=16,col='red',
     xlab='Commission',ylab='Omission', 
     main='Omission vs Commission (Egregious) Errors')
abline(0,1)
text(com.errs+10,om.errs+25,names(om.errs))

# create kmz with egregious issues
colnames(egreg.pix)=c('Lat','Lon','Train','Pred','Source')
writePointsKmz(lc.points=egreg.pix[,c('Lat','Lon')], 
               Train.class = egreg.pix[,'Train'], 
               Pred.class = egreg.pix[,'Pred'],
               ID = glancedat[rownames(egreg.pix),'system.index'],
               year = glancedat[rownames(egreg.pix),years],
               outfil = paste('kmls/',cont,'badTrain.kml',sep=''))


##### Filtering step 2: remove based on egregious & low margin cases  ####
low.margin.thresh <- 0.05

# split training data into egregious/low margin versus non-egregious sufficient margin
bad.pix <- (margins.all[,'missClass'] & (margins.all[,'Margin'] > egreg.thresh.high)) | # egregiously misclassified
  (margins.all[,'Margin'] < low.margin.thresh)                                          # indistinguishable

# subset data into 'good' and 'bad' training
badTraining <- featureDat[bad.pix,]
goodTraining <- featureDat[!bad.pix,]

# look at distribution of screened data
plot(as.factor(featureDat[,'LC_Class']),col='green',main = 'Full Training vs "Good" Training')
plot(as.factor(goodTraining[,'LC_Class']),col='blue', add=T)
plot(as.factor(featureDat[,'LC_Class']),col='green',main = 'Full Training vs "Bad" Training')
plot(as.factor(badTraining[,'LC_Class']),col='blue', add=T)

# re-estimate model
ranger.goodT <- ranger(LC_Class~.,
                       data=goodTraining,
                       num.trees=ntrees,
                       importance='none',
                       probability=FALSE)

# Generate CF and accuracy statistics
ranger.goodT.cf <- ranger.goodT$confusion.matrix
ranger.goodT.acc <- 1-ranger.goodT$prediction.error
ranger.goodT.cf
ranger.goodT.acc
com.errs <- rowSums(ranger.goodT.cf)-diag(ranger.goodT.cf)
om.errs <- colSums(ranger.goodT.cf)-diag(ranger.goodT.cf)
plot(com.errs,om.errs,pch=16,col='red',xlab='Commission',ylab='Omission')
abline(0,1)
text(com.errs+10,om.errs+5,names(om.errs))

# Assign likelihoods to each training case 
pvec <- matrix(NA,dim(goodTraining)[1])
for (i in 1:dim(goodTraining)[1]) {
  pvec[i] <- assignP(goodTraining[i,'LC_Class'])
}

###  At this point have removed: BS-crazy, egregious, low margin  ###

# generate new sample from original data rebalanced based on likelihoods
nsamps <- 10000
smpl <- sample(1:dim(goodTraining)[1],size=nsamps,replace=FALSE,prob=pvec)

# Plot resampled vs original training set
plot(as.factor(featureDat[,'LC_Class']),col='green',main = 'Full Training Set vs Screened & Rebalanced Training')
plot(as.factor(goodTraining[smpl,'LC_Class']),col='blue',add=TRUE)

# Plot resampled vs training set filtered for egregious errors
plot(as.factor(goodTraining[,'LC_Class']),col='green',main = 'Filtered Training Set vs Screened & Rebalanced Training')
plot(as.factor(goodTraining[smpl,'LC_Class']),col='blue',add=TRUE)

# percentage of each class in re-sampled training data
table(goodTraining[smpl,'LC_Class'])/sum(table(goodTraining[smpl,'LC_Class']))

# re-estimate model, with re-sampled data, estimate confusion matrix and accuracy
ranger.resamp <- ranger(LC_Class~.,
                      data=goodTraining[smpl,],
                      num.trees=ntrees,
                      importance='none',
                      probability=FALSE)
ranger.resamp$confusion.matrix
1-ranger.resamp$prediction.error

# and plot result
om.errs <- rowSums(ranger.resamp$confusion.matrix)-diag(ranger.resamp$confusion.matrix)
com.errs <- colSums(ranger.resamp$confusion.matrix)-diag(ranger.resamp$confusion.matrix)
plot(com.errs,om.errs,pch=16,col='red',xlab='Commission',ylab='Omission',ylim=c(0,400),xlim=c(0,400))
abline(0,1)
text(com.errs+10,om.errs+5,names(om.errs))

# cross-validation
clean.dat <- goodTraining[smpl,]
n.sampl <- dim(clean.dat)[1]
n.cv <- 20

# create sample and data frame for results; set up DF as list of factors
smpl.cv <- sample(1:n.cv,size=n.sampl,replace=T)
preds.cv <- data.frame(matrix('Bare',n.sampl))
true.cv <- data.frame(matrix('Bare',n.sampl))
preds.cv[,1] <- factor(preds.cv[,1],levels=c("Bare","Developed","Forest","Herbaceous","Shrub","Snow/Ice","Water"))
true.cv[,1] <- factor(true.cv[,1],levels=c("Bare","Developed","Forest","Herbaceous","Shrub","Snow/Ice","Water"))

# and do the cross validation
for (i in 1:n.cv) {
  test <- smpl.cv==i
  train.dat <- clean.dat[!test & smpl,]
  test.dat <- clean.dat[test,]
  ranger.cv <- ranger(LC_Class~.,
                      data=train.dat,
                      num.trees=ntrees,
                      importance='none',
                      probability=FALSE)
  preds.cv[test,] <- predict(ranger.cv,data=test.dat)$predictions
  true.cv[test,] <- clean.dat[test,'LC_Class']
}

# summarize and plot results.
cv.cf <- table(true.cv[,1],preds.cv[,1])
cv.acc <- sum(diag(cv.cf))/sum(cv.cf)
cv.cf
cv.acc
om.errs <- rowSums(cv.cf)-diag(cv.cf)
com.errs <- colSums(cv.cf)-diag(cv.cf)
plot(com.errs,om.errs,pch=16,col='red',xlab='Commission',ylab='Omission',ylim=c(0,400),xlim=c(0,400))
abline(0,1)
text(com.errs+10,om.errs+5,names(om.errs))

# write data sets to file.
clean.rows <- rownames(clean.dat)
clean.dat.allfeatures <- glancedat[clean.rows,]
write.csv(clean.dat[,c(1,51:2)],file='data/NA_Filtered_70Features.csv',row.names=FALSE)
write.csv(clean.dat.allfeatures,file='data/NA_Filtered_AllFeatures.csv',row.names=FALSE)

