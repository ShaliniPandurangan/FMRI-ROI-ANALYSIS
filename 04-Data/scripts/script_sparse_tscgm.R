library(imputeTS)
load("/Users/usi/Desktop/04-Data/fMRI-ROI-time-series/fMRI-ROI-time-series.RData")
library(SparseTSCGM)
library(longitudinal)
library(MVN)

n=24
T=404
p=70
sum(is.na(t(Y[,,2,1])))


Ycomplete <-  Y[,,which(apply(Y[,,,1],3,function(x) sum(is.na(x)))!=28280),1]
Ycomplete <-  apply(Ycomplete,c(3,1),na_ma)

Ytrans<-qnorm(rank((Ycomplete)+.5)/(sum(!is.na(Ycomplete))+1))
dim(Ytrans)<-dim(Ycomplete)

list_timep = do.call(rbind,lapply(1:T, function(i) Ytrans[i,,]))

# list_timep = do.call(rbind,lapply(1:T, function(i) Ycomplete[i,,]))

sprseTSCGMdata = as.longitudinal(list_timep, repeats = 22)
sprseTSCGMdata = scale(sprseTSCGMdata, scale=FALSE, center=TRUE)
sparse.tscgm(data = sprseTSCGMdata,model="ar1", penalty="lasso",optimality="bic",
             control=list(maxit.out = 5, maxit.in = 50, tol.out = 1e-03))
