
##########################

# save, then delete cases w/ecoregion = 0
EC.0.dat <- subset(glancedat,glancedat[,'Level1_Ecoregion']==0)
glancedat <- subset(glancedat,glancedat[,'Level1_Ecoregion']!=0)

# Define aggregated ecoregions
ec.grp1 <- glancedat[,'Level1_Ecoregion'] == 2 | glancedat[,'Level1_Ecoregion'] == 4 # high latitude, low stature, wet
ec.grp2 <- glancedat[,'Level1_Ecoregion'] == 3 | glancedat[,'Level1_Ecoregion'] == 5 # high latitude forests
ec.grp3 <- glancedat[,'Level1_Ecoregion'] == 6 | glancedat[,'Level1_Ecoregion'] == 7 # western forests
ec.grp4 <- glancedat[,'Level1_Ecoregion'] == 8                                       # eastern forest
ec.grp5 <- glancedat[,'Level1_Ecoregion'] == 9 | glancedat[,'Level1_Ecoregion'] == 10 | glancedat[,'Level1_Ecoregion'] == 11 # arid west
ec.grp6 <- glancedat[,'Level1_Ecoregion'] == 12 | glancedat[,'Level1_Ecoregion'] == 13 # temperate + arid highlands/Mex
ec.grp7 <- glancedat[,'Level1_Ecoregion'] == 14 | glancedat[,'Level1_Ecoregion'] == 15 # tropical/sub-tropical forests

# extract rows from original data corresponding to each group
ec.grp1.rws <- rownames(glancedat[ec.grp1,])
ec.grp2.rws <- rownames(glancedat[ec.grp2,])
ec.grp3.rws <- rownames(glancedat[ec.grp3,])
ec.grp4.rws <- rownames(glancedat[ec.grp4,])
ec.grp5.rws <- rownames(glancedat[ec.grp5,])
ec.grp6.rws <- rownames(glancedat[ec.grp6,])
ec.grp7.rws <- rownames(glancedat[ec.grp7,])

# look at distribution of classes in each aggregated ecoregion
barplot(table(glancedat[ec.grp1.rws,'LC_Class']),col='blue', main='High Latitude/Tundra/Wetlands')
barplot(table(glancedat[ec.grp2.rws,'LC_Class']),col='green', main='Boreal + Taiga Forests')
barplot(table(glancedat[ec.grp3.rws,'LC_Class']),col='brown', main='Western Forests')
barplot(table(glancedat[ec.grp4.rws,'LC_Class']),col='orange', main='Eastern Forests')
barplot(table(glancedat[ec.grp5.rws,'LC_Class']),col='red', main='Arid West')
barplot(table(glancedat[ec.grp6.rws,'LC_Class']),col='black', main='Sub-Tropical/Semi Arid Northern Mexico')
barplot(table(glancedat[ec.grp7.rws,'LC_Class']),col='purple', main='Tropical/Sub-Tropical Forests')

barplot(c(sum(ec.grp1),sum(ec.grp2),sum(ec.grp3),sum(ec.grp4),sum(ec.grp5),sum(ec.grp6),sum(ec.grp7)),
          names.arg=c('EC1','EC2','EC3','EC4','EC5','EC6','EC7'),
        col='green',
        main='Number of Samples in Each Aggregated Ecoregion')
  
#  break out into super-ecoregions
ec1.gdat <- glancedat[ec.grp1.rws,]
ec2.gdat <- glancedat[ec.grp2.rws,]
ec3.gdat <- glancedat[ec.grp3.rws,]
ec4.gdat <- glancedat[ec.grp4.rws,]
ec5.gdat <- glancedat[ec.grp5.rws,]
ec6.gdat <- glancedat[ec.grp6.rws,]
ec7.gdat <- glancedat[ec.grp7.rws,]

# filter each super-region
splt.flag=FALSE                # dontSplit= FALSE -> split data into train & test 
ec1.dat <- filterECdat(f.dat=ec1.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag) 
ec2.dat <- filterECdat(f.dat=ec2.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag)
ec3.dat <- filterECdat(f.dat=ec3.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag)
ec4.dat <- filterECdat(f.dat=ec4.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag)
ec5.dat <- filterECdat(f.dat=ec5.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag)
ec6.dat <- filterECdat(f.dat=ec6.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag)
ec7.dat <- filterECdat(f.dat=ec7.gdat[,c('LC_Class',best.features)],
                       low.marg.thresh = 0.05,
                       pc.test=0.1,
                       dontSplit=splt.flag)

# create a single data frame with all the data
if (splt.flag) {
  filtered.data.final <- rbind(ec1.dat,
                               ec2.dat,
                               ec3.dat,
                               ec4.dat,
                               ec5.dat,
                               ec6.dat,
                               ec7.dat)
} else {
  filtered.data.final <- rbind(ec1.dat$train,ec1.dat$test,
                               ec2.dat$train,ec2.dat$test,
                               ec3.dat$train,ec3.dat$test,
                               ec4.dat$train,ec4.dat$test,
                               ec5.dat$train,ec5.dat$test,
                               ec6.dat$train,ec6.dat$test,
                               ec7.dat$train,ec7.dat$test)
}

# look at class dist in each super-region
barplot(table(glancedat[ec.grp1.rws,'LC_Class']),col='lightblue', main='High Latitude/Tundra/Wetlands')
barplot(table(ec1.dat$train[,1]),col='blue',add=T)

barplot(table(glancedat[ec.grp2.rws,'LC_Class']),col='lightgreen', main='Boreal + Taiga Forests')
barplot(table(ec2.dat$train[,1]),col='green',add=T)

barplot(table(glancedat[ec.grp3.rws,'LC_Class']),col='brown', main='Western Forests')
barplot(table(ec3.dat$train[,1]),col='brown1',add=T)

barplot(table(glancedat[ec.grp4.rws,'LC_Class']),col='orange', main='Eastern Forests')
barplot(table(ec4.dat$train[,1]),col='orange3',add=T)

barplot(table(glancedat[ec.grp5.rws,'LC_Class']),col='darkred', main='Arid West')
barplot(table(ec5.dat$train[,1]),col='red',add=T)

barplot(table(glancedat[ec.grp6.rws,'LC_Class']),col='gray', main='Sub-Tropical/Semi Arid Northern Mexico')
barplot(table(ec6.dat$train[,1]),col='black',add=T)

barplot(table(glancedat[ec.grp7.rws,'LC_Class']),col='purple', main='Tropical/Sub-Tropical Forests')
barplot(table(ec7.dat$train[,1]),col='purple3',add=T)

# now do classifications
ec1.rf <- doEcClass(ec.dat=ec1.dat$train)
ec2.rf <- doEcClass(ec.dat=ec2.dat$train)
ec3.rf <- doEcClass(ec.dat=ec3.dat$train)
ec4.rf <- doEcClass(ec.dat=ec4.dat$train)
ec5.rf <- doEcClass(ec.dat=ec5.dat$train)
ec6.rf <- doEcClass(ec.dat=ec6.dat$train)
ec7.rf <- doEcClass(ec.dat=ec7.dat$train)

# compile predictions from training data into single data frames
all.preds <- unlist(list(ec1.rf$predictions,
                         ec2.rf$predictions,
                         ec3.rf$predictions,
                         ec4.rf$predictions,
                         ec5.rf$predictions,
                         ec6.rf$predictions,
                         ec7.rf$predictions))

all.true <- unlist(list(ec1.dat$train[,'LC_Class'],
                        ec2.dat$train[,'LC_Class'],
                        ec3.dat$train[,'LC_Class'],
                        ec4.dat$train[,'LC_Class'],
                        ec5.dat$train[,'LC_Class'],
                        ec6.dat$train[,'LC_Class'],
                        ec7.dat$train[,'LC_Class']))

all.cf.mat <- table(all.true,all.preds)
all.cf.mat
sum(diag(all.cf.mat))/sum(all.cf.mat)

# plot errors of omission/commission
om.errs <- rowSums(all.cf.mat)-diag(all.cf.mat)
com.errs <- colSums(all.cf.mat)-diag(all.cf.mat)
plot(com.errs,om.errs,pch=16,col='red',xlab='Commission',ylab='Omission',ylim=c(0,550),xlim=c(0,550))
abline(0,1)
text(com.errs+10,om.errs+5,names(om.errs))

# now look at test data
ec1.test.preds <- predict(ec1.rf,data=ec1.dat$test)
ec2.test.preds <- predict(ec2.rf,data=ec2.dat$test)
ec3.test.preds <- predict(ec3.rf,data=ec3.dat$test)
ec4.test.preds <- predict(ec4.rf,data=ec4.dat$test)
ec5.test.preds <- predict(ec5.rf,data=ec5.dat$test)
ec6.test.preds <- predict(ec6.rf,data=ec6.dat$test)
ec7.test.preds <- predict(ec7.rf,data=ec7.dat$test)

# put predictions from test data in a single vector
all.test.preds <- unlist(list(ec1.test.preds$predictions,
                              ec2.test.preds$predictions,
                              ec3.test.preds$predictions,
                              ec4.test.preds$predictions,
                              ec5.test.preds$predictions,
                              ec6.test.preds$predictions,
                              ec7.test.preds$predictions))

# put training labels from test data in a single vector
all.test.true <- unlist(list(ec1.dat$test[,'LC_Class'],
                             ec2.dat$test[,'LC_Class'],
                             ec3.dat$test[,'LC_Class'],
                             ec4.dat$test[,'LC_Class'],
                             ec5.dat$test[,'LC_Class'],
                             ec6.dat$test[,'LC_Class'],
                             ec7.dat$test[,'LC_Class']))

# estimate accuracy
all.test.cf <- table(all.test.true,all.test.preds)
all.test.cf
sum(diag(all.test.cf))/sum(all.test.cf)
plotOmCom(all.test.cf)

# now, put it all back in a single data set, estimate a single model, and see what changes?
all.dat <- rbind(ec1.dat$train,ec2.dat$train,ec3.dat$train,ec4.dat$train,ec5.dat$train)
dim(all.dat)
all.rf <- doEcClass(all.dat)
1-all.rf$prediction.error
all.rf$confusion.matrix
plotOmCom(all.rf$confusion.matrix)

# predict on test data using model based on 'all the data' (ie., does not include EcoRegions)
all.dat.test <- rbind(ec1.dat$test,ec2.dat$test,ec3.dat$test,ec4.dat$test,ec5.dat$test)
test.preds <- predict(all.rf,data=all.dat.test)
test.cf <- table(all.dat.test[,'LC_Class'],test.preds$predictions)
sum(diag(test.cf))/sum(test.cf)

# write csv files
# Top 50 features
write.csv(filtered.data.final,
          file='data/NA_Filtered_ByEcoRegionsTop50Features.csv',
          row.names=FALSE)

# All features
write.csv(glancedat[rownames(filtered.data.final),],
          file='data/NA_Filtered_ByEcoRegionsAllFeatures.csv',
          row.names=FALSE)

