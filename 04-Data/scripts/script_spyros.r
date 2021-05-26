library(imputeTS)

n=24
T=404
p=70
sum(is.na(t(Y[,,2,1])))


Ycomplete <-  Y[,,which(apply(Y[,,,1],3,function(x) sum(is.na(x)))!=28280),1]
Ycomplete <-  apply(Ycomplete,c(3,1),na_ma)
list_timep = do.call(rbind,lapply(1:T, function(i) Ycomplete[i,,]))
sprseTSCGMdata = as.longitudinal(list_timep, repeats = 22)
sparse.tscgm(data = sprseTSCGMdata,model="ar1", penalty="lasso",optimality="bic",
             control=list(maxit.out = 10, maxit.in = 100))
