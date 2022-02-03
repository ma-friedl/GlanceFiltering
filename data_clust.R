# Example script to cluster Data into climate/geographic regions

# Set up working dir
setwd("/projectnb/measures/users/rkstan/glance/")

# First read data
cont <- 'EU'      # NA = North America, AF = Africa.....
inputFile <- 'training_data/EU_Training_Master_V1_2021_12_22_predictors_BU_sub.csv'
glancedat <- read.csv(inputFile,header=TRUE)
n.regions <- 9

# select climate data attributes, scale, and add lat/long
clust.dat <- scale(glancedat[,c('maxTemp','minTemp','precip')])
clust.dat <- data.frame(glancedat[,'lon'],glancedat[,'lat'],clust.dat)
pdf(paste("training_data/results/",cont,"_kmclusters.pdf", sep=""))
# first do kmeans

for (i in 3:n.regions) {
  clust.km <- kmeans(clust.dat,centers=i,iter.max=100)$cluster
  plot(clust.dat[,1],clust.dat[,2],col=clust.km,xlab='Long',ylab='Lat', main=paste(i,"Clusters"))
  hist(clust.km,plot=F,breaks=seq(0.5,i+0.5,1))$counts
  plot(as.factor(clust.km),col=1:i)
}
dev.off()




## select climate data attributes, scale, and add lat/long
#clust.dat <- scale(glancedat[,c('maxTemp','minTemp','precip')])
#clust.dat <- data.frame(glancedat[,'lon'],glancedat[,'lat'],clust.dat)
## first do kmeans
#clust.km <- kmeans(clust.dat,centers=n.regions,iter.max=100)$cluster
#plot(clust.dat[,1],clust.dat[,2],col=clust.km,xlab='Long',ylab='Lat')
#hist(clust.km,plot=F,breaks=seq(0.5,n.regions+0.5,1))$counts
#plot(as.factor(clust.km),col=1:n.regions)

## now do hierarchical - first compute distance matrix
#d.mat <- dist(clust.dat)
## estimate using hierarchical clustering
#clust.hr <- hclust(d.mat,method='ward.D2')
#clust.labels <- cutree(clust.hr, k=n.regions)
#plot(clust.dat[,1],clust.dat[,2],col=clust.labels,xlab='Long',ylab='Lat')
##hist(clust.hr,plot=F)
#plot(as.factor(clust.labels),col=1:n.regions+0.5)
