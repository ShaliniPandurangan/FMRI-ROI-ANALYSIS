setwd("C:/Users/ernst/Dropbox/ernst/Genetics/TalksMeetings/2017/11-Startup-Research-Siena/04-Data/data")
library(igraph)
library(huge)
library(scatterplot3d)
library(SparseTSCGM)
library(ergm)
library(statnet)

# D: DTI Structural network data (anatomical interactions) 
# 24 individuals
# 2 measurements per individual
# Measurements are number of anatomical connections between 70 ROIs in Desikan atlas (note: 2 are unknown)
# Note: the block structure corresponds with the left and right region of brain
load("C:/Users/ernst/Dropbox/ernst/Genetics/TalksMeetings/2017/11-Startup-Research-Siena/04-Data/data/DTI-connectome.RData")
dim(D)
head(D[,,1,1])
image(D[,,1,1])

# Y: fMRI connectome data (dynamic functional interactions)
# 24 individuals
# 2 measurements per individual
# 70 times series of changes in blood-oxygen levels for each ROI of length 404 timepoints 
load("C:/Users/ernst/Dropbox/ernst/Genetics/TalksMeetings/2017/11-Startup-Research-Siena/04-Data/data/fMRI-ROI-time-series.RData")
dim(Y)
plot(Y[1,,2,2])
lines(smooth.spline(x=1:403,y=Y[1,-404,2,2], df=80))


# W: fMRI connectome data (dynamic functional interactions)
# 24 individuals
# 2 measurements per individual
# Measurements is synchronization of temporal fMRI profiles between pairs of ROI regions. 
#            in Desikan atlas (note: 2 are unknown)
# synchronization = correlation
# Note: secundary diagonal, which corresponds with the LEFT and RIGHT region of brain
load("C:/Users/ernst/Dropbox/ernst/Genetics/TalksMeetings/2017/11-Startup-Research-Siena/04-Data/data/fMRI-connectome.RData")
dim(W)
head(W[,,2,2])
image(W[,,2,2])
range(W[,,2,2])
# Note: synchronization = correlation
cor(Y[1,,2,2],Y[2,,2,2],use="complete.obs")
W[1,2,2,2]


# ROI Covariates
# 3 covariates for the 70 ROIs in the study
# - which hemisphere (left, right)
# - which lobe (temporal, limbic, frontal, inter-hemispheric occipital, parietal)
# - 3D locations
roi.cov<-read.table("ROI-covariates.txt",header=T)
head(roi.cov)
tabulate(roi.cov$lobes)
unique(roi.cov$lobes)
scatterplot3d(x=roi.cov$centroid_x,y=roi.cov$centroid_y,z=roi.cov$centroid_z)

# SUBJECT covariates
# 6 covariates for the 24 subjects in the study
# - age
# - handedness (Left or Right)
# - current diagnosis (YES, NO, NA)
# - type current diagnosis (4 people with some CURRENT abuse/depression)
# - lifetime diagnosis (YES, NO, NA)
# - type lifetime diagnois (7 people that have had some form of abuse/depression)

sub.cov<-read.table("SUBJ-covariates.txt",header=T)
head(sub.cov)
dim(sub.cov)


# --------------------------------------------------------------
# Temporal slices 
tmp22<-apply(Y[,-404,2,2],1,function(y,x,df){smooth.spline(x,y,df=df)$y},x=1:403,df=80)


# ----------------
# SPARSE TSCGM
Y.long<-list()
tscgm<-list()
tms<-1:403
for (i in 1:24){
  if (is.na(Y[1,1,i,1])){
    Y.long[[i]]<-NA
  } else {
    if (is.na(Y[1,1,i,2])){
      hlp<-t(Y[,-404,i,1])
      rp<-1
    } else {
      hlp<-t(matrix(c(rbind(Y[,-404,i,1],Y[,-404,i,2])),nrow=70))
      rp<-2
    }
    print(i)
    Y.long[[i]]<-as.longitudinal(hlp,repeats=rep(rp,403))
    tscgm[[i]]<-sparse.tscgm(data=Y.long[[i]],optimality="bic")
  }
}
plot(Y.long[[2]])
tscgm2<-sparse.tscgm(data=Y.long[[2]],lam1=.8,lam2=.5,control=list(silent=F))
plot.tscgm(tscgm2,mat="precision")
plot.tscgm(tscgm2,mat="autoregression")

tscgm2<-sparse.tscgm(data=Y.long[[2]],control=list(silent=F),optimality = "bic")
plot.tscgm(tscgm2,mat="precision")
plot.tscgm(tscgm2,mat="autoregression")

plot.tscgm(tscgm[[3]],mat="precision")
plot.tscgm(tscgm[[3]],mat="autoregression")
