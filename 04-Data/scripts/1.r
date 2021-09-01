sMS_AIC <- summary(ms_AIC)
sMS_BIC <- summary(ms_BIC)
##sBIC <- summary(BIC)
sarmod <- summary(armod)
##resultAIC <- list()
##resultBIC <- list()
resultMS_AIC <- list()
resultMS_BIC <- list()
##resultAIC <- sAIC
resultMS_AIC <- sMS_AIC
resultMS_BIC <- sMS_BIC
##resultBIC <- sBIC
sarmod <- summary(armod)
ARcoef <- sarmod$coefficients[-1,1]
MScoef_AIC <- sMS_AIC$coefficients[-1,1]
MScoef_BIC <- sMS_BIC$coefficients[-1,1]
ARcoefPval <- sarmod$coefficients[-1,4]
MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
Interc <- ARcoef[1]
IntercMS_AIC <- MScoef_AIC[1]
IntercMS_A=BIC <- MScoef_BIC[1]
resultAR <- list()
resultMS_AIC <- list()
resultMS_BIC <- list()
resultAR$COEF <- ARcoef
resultMS_AIC $COEF <- MScoef_AIC
resultMS_BIC $COEF <- MScoef_BIC
resultAR$PVAL <- ARcoefPval
resultMS_AIC $PVAL <- MScoefPval_AIC
resultMS_BIC $PVAL <- MScoefPval_BIC
resultAR$INTERC <- Interc
resultMS_AIC $INTERC <- IntercMS_AIC
resultMS_BIC $INTERC <- IntercMS_BIC
##resultAR$sAIC <- resultAIC
##resultAR$sMS_AIC <- resultMS_AIC
##resultAR$sMS_BIC <- resultMS_BIC
##resultAR$sBIC <- resultBIC
resultAR$RESID <- armod$residuals
resultAR$FITTED <- armod$fitted.values
resultAR$STD.DEV <- sarmod$sigma
resultAR$R.SQUARED <- sarmod$r.squared
resultMS_AIC$RESID <- ms_AIC$residuals
resultMS_AIC$FITTED <- ms_AIC$fitted.values
resultMS_AIC$STD.DEV <- sMS_AIC $sigma
resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
resultMS_BIC$RESID <- ms_BIC$residuals
resultMS_BIC$FITTED <- ms_BIC$fitted.values
resultMS_BIC$STD.DEV <- sMS_BIC$sigma
resultMS_BIC$R.SQUARED <- sMS_BIC$r.squared
})
model_results <-lapply(1:nROI, function(i){
cat(i)
y <- Yfit[,i]
y <- unname(y)
x <- Xfit
complete_obs <- complete.cases(y) & complete.cases(x)
y <- y[complete_obs]
x <- x[complete_obs,]
data1 <-data.frame(cbind(y,x))
##y_name <- gsub( "_.*$", "", names(data1)[i] )
colnames(x) <- gsub("-",".",colnames(x))
##names(data1)[1] <- y_name
model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))
## armod <- lm(y ~ x)
lw <- formula(paste("~",paste(colnames(x)[1:10],collapse = "+")))
up <- formula(paste("~",paste(colnames(x)[1:70],collapse = "+")))
armod <- lm(model_formulae, data=data1)
ms_AIC <- mStep(armod, scope = list(upper=up),direction="both", k=2)
##ms_BIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
##AIC <- step(armod,direction = "both", k = 2)
##LARS <- lars(model_formulae, type="lasso")
##BIC <- step(armod,direction = "both", k = log(70))
##AIC <- stepAIC(armod, direction=c("both"))
## BIC <- BIC(armod, direction=c("both"))
##sAIC <- summary(AIC)
sMS_AIC <- summary(ms_AIC)
##sMS_BIC <- summary(ms_BIC)
##sBIC <- summary(BIC)
sarmod <- summary(armod)
##resultAIC <- list()
##resultBIC <- list()
resultMS_AIC <- list()
##resultMS_BIC <- list()
##resultAIC <- sAIC
resultMS_AIC <- sMS_AIC
##resultMS_BIC <- sMS_BIC
##resultBIC <- sBIC
sarmod <- summary(armod)
ARcoef <- sarmod$coefficients[-1,1]
MScoef_AIC <- sMS_AIC$coefficients[-1,1]
##MScoef_BIC <- sMS_BIC$coefficients[-1,1]
ARcoefPval <- sarmod$coefficients[-1,4]
MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
##MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
Interc <- ARcoef[1]
IntercMS_AIC <- MScoef_AIC[1]
##IntercMS_A=BIC <- MScoef_BIC[1]
resultAR <- list()
resultMS_AIC <- list()
##resultMS_BIC <- list()
resultAR$COEF <- ARcoef
resultMS_AIC $COEF <- MScoef_AIC
##resultMS_BIC $COEF <- MScoef_BIC
resultAR$PVAL <- ARcoefPval
resultMS_AIC $PVAL <- MScoefPval_AIC
##resultMS_BIC $PVAL <- MScoefPval_BIC
resultAR$INTERC <- Interc
resultMS_AIC $INTERC <- IntercMS_AIC
##resultMS_BIC $INTERC <- IntercMS_BIC
##resultAR$sAIC <- resultAIC
##resultAR$sMS_AIC <- resultMS_AIC
##resultAR$sMS_BIC <- resultMS_BIC
##resultAR$sBIC <- resultBIC
resultAR$RESID <- armod$residuals
resultAR$FITTED <- armod$fitted.values
resultAR$STD.DEV <- sarmod$sigma
resultAR$R.SQUARED <- sarmod$r.squared
resultMS_AIC$RESID <- ms_AIC$residuals
resultMS_AIC$FITTED <- ms_AIC$fitted.values
resultMS_AIC$STD.DEV <- sMS_AIC $sigma
resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
return(resultMS_AIC)
})
VAR_mat_AIC <- do.call(rbind, lapply(model_results, function(x) x$COEF))
View(VAR_mat_AIC)
VAR_mat_AIC <-  t(VAR_mat_AIC)
library(qgraph)
qgraph_net <- qgraph(VAR_mat_AIC,
posCol="blue",
negCol="purple")
View(VAR_mat)
library(qgraph)
qgraph_net <- qgraph(VAR_mat_AIC,
posCol="blue",
negCol="purple")
Yfit <- do.call(rbind, lapply(Data_YX, function(x) x$Y))
library(SparseTSCGM)
library(mgm)
library(igraph)
library(huge)
library(scatterplot3d)
library(SparseTSCGM)
library(ergm)
library(MASS)
library(fbasics)
library(dplyr)
library(tidypredict)
library(statnet)
library(mgm)
## Get the data
Data_Y <- Y
Data_Y <- Y[,,,1]
## Get the data dimensions
dat_dims <- dim(Data_Y)
nTime <- dat_dims[2]
nPersons <- dat_dims[3]
nROI <- dat_dims[1]
## Specify the lag of the proces
lag <- 1
Data_YX <- lapply(1:nPersons, function(i){
result <- list()
result$Y <- NA
result$X <- NA
get_lags <- mgm:::lagData(data = t(Data_Y[,,i]), lags = 1:lag)
result$Y <- get_lags$data_response
result$X <- do.call(cbind, get_lags$l_data_lags)
return(result)
})
Yfit <- do.call(rbind, lapply(Data_YX, function(x) x$Y))
dimnames(Yfit)[[1]] <- NULL
Xfit <- do.call(rbind, lapply(Data_YX, function(x) x$X))
dimnames(Xfit) <- list(NULL, paste(paste(colnames(Yfit),"_lag_",sep = ""), 1:lag, sep = ""))
## Fit linear models and get useful output
model_results <-lapply(1:nROI, function(i){
cat(i)
y <- Yfit[,i]
y <- unname(y)
x <- Xfit
complete_obs <- complete.cases(y) & complete.cases(x)
y <- y[complete_obs]
x <- x[complete_obs,]
data1 <-data.frame(cbind(y,x))
##y_name <- gsub( "_.*$", "", names(data1)[i] )
colnames(x) <- gsub("-",".",colnames(x))
##names(data1)[1] <- y_name
model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))
## armod <- lm(y ~ x)
lw <- formula(paste("~",paste(colnames(x)[1:10],collapse = "+")))
up <- formula(paste("~",paste(colnames(x)[1:70],collapse = "+")))
armod <- lm(model_formulae, data=data1)
ms_AIC <- mStep(armod, scope = list(upper=up),direction="both", k=2)
ms_BIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
##AIC <- step(armod,direction = "both", k = 2)
##LARS <- lars(model_formulae, type="lasso")
##BIC <- step(armod,direction = "both", k = log(70))
##AIC <- stepAIC(armod, direction=c("both"))
##BIC <- BIC(armod, direction=c("both"))
##sAIC <- summary(AIC)
sMS_AIC <- summary(ms_AIC)
sMS_BIC <- summary(ms_BIC)
##sBIC <- summary(BIC)
sarmod <- summary(armod)
##resultAIC <- list()
##resultBIC <- list()
resultMS_AIC <- list()
resultMS_BIC <- list()
##resultAIC <- sAIC
resultMS_AIC <- sMS_AIC
resultMS_BIC <- sMS_BIC
##resultBIC <- sBIC
sarmod <- summary(armod)
ARcoef <- sarmod$coefficients[-1,1]
MScoef_AIC <- sMS_AIC$coefficients[-1,1]
MScoef_BIC <- sMS_BIC$coefficients[-1,1]
ARcoefPval <- sarmod$coefficients[-1,4]
MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
Interc <- ARcoef[1]
IntercMS_AIC <- MScoef_AIC[1]
IntercMS_A=BIC <- MScoef_BIC[1]
resultAR <- list()
resultMS_AIC <- list()
resultMS_BIC <- list()
resultAR$COEF <- ARcoef
resultMS_AIC $COEF <- MScoef_AIC
resultMS_BIC $COEF <- MScoef_BIC
resultAR$PVAL <- ARcoefPval
resultMS_AIC $PVAL <- MScoefPval_AIC
resultMS_BIC $PVAL <- MScoefPval_BIC
resultAR$INTERC <- Interc
resultMS_AIC $INTERC <- IntercMS_AIC
resultMS_BIC $INTERC <- IntercMS_BIC
##resultAR$sAIC <- resultAIC
##resultAR$sMS_AIC <- resultMS_AIC
##resultAR$sMS_BIC <- resultMS_BIC
##resultAR$sBIC <- resultBIC
resultAR$RESID <- armod$residuals
resultAR$FITTED <- armod$fitted.values
resultAR$STD.DEV <- sarmod$sigma
resultAR$R.SQUARED <- sarmod$r.squared
resultMS_AIC$RESID <- ms_AIC$residuals
resultMS_AIC$FITTED <- ms_AIC$fitted.values
resultMS_AIC$STD.DEV <- sMS_AIC $sigma
resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
resultMS_BIC$RESID <- ms_BIC$residuals
resultMS_BIC$FITTED <- ms_BIC$fitted.values
resultMS_BIC$STD.DEV <- sMS_BIC$sigma
resultMS_BIC$R.SQUARED <- sMS_BIC$r.squared
})
ls()
View(AIC_res)
View(resultAIC)
VAR_mat_AIC <- do.call(rbind, lapply(model_results, function(x) x$COEF))
View(VAR_mat_AIC)
library(mst)
library(mStep)
library(qtlmt)
data(etrait)
mdf<- data.frame(traits,markers)
## Not run:
mlm<- lm(cbind(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16) ~
m1 + m2 + m3 + m4 + m5, data=mdf)
loww<- formula(paste("~ ", paste("m",1:3,collapse=" + ",sep="")))
upp<- formula(paste("~", paste("m",1:15,collapse=" + ",sep="")))
ob<- mStep(mlm, scope=list(lower=loww), k=99, direction="backward", data=mdf)
of<- mStep(mlm, scope=list(upper=upp), k=5, direction="forward", data=mdf)
o1<- mStep(mlm, scope=list(upper=upp), k=5, direction="both", data=mdf)
o2<- mStep(o1, scope=list(upper=upp), k=2, direction="forward", data=mdf)
View(o2)
o2[["terms"]]
dimnames(VAR_mat_AIC) <-  list(colnames(Yfit),colnames(Yfit))
View(VAR_mat_AIC)
View(VAR_mat)
model_results <-lapply(1:nROI, function(i){
cat(i)
y <- Yfit[,i]
y <- unname(y)
x <- Xfit
complete_obs <- complete.cases(y) & complete.cases(x)
y <- y[complete_obs]
x <- x[complete_obs,]
data1 <-data.frame(cbind(y,x))
##y_name <- gsub( "_.*$", "", names(data1)[i] )
colnames(x) <- gsub("-",".",colnames(x))
##names(data1)[1] <- y_name
model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))
## armod <- lm(y ~ x)
lw <- formula(paste("~",paste(colnames(x)[1:10],collapse = "+")))
up <- formula(paste("~",paste(colnames(x)[1:70],collapse = "+")))
armod <- lm(model_formulae, data=data1)
##ms_AIC <- mStep(armod, scope = list(upper=up),direction="both", k=2)
ms_BIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
##AIC <- step(armod,direction = "both", k = 2)
##LARS <- lars(model_formulae, type="lasso")
##BIC <- step(armod,direction = "both", k = log(70))
##AIC <- stepAIC(armod, direction=c("both"))
##BIC <- BIC(armod, direction=c("both"))
##sAIC <- summary(AIC)
##sMS_AIC <- summary(ms_AIC)
sMS_BIC <- summary(ms_BIC)
##sBIC <- summary(BIC)
##sarmod <- summary(armod)
##resultAIC <- list()
##resultBIC <- list()
##resultMS_AIC <- list()
resultMS_BIC <- list()
##resultAIC <- sAIC
##resultMS_AIC <- sMS_AIC
resultMS_BIC <- sMS_BIC
##resultBIC <- sBIC
sarmod <- summary(armod)
ARcoef <- sarmod$coefficients[-1,1]
##MScoef_AIC <- sMS_AIC$coefficients[-1,1]
MScoef_BIC <- sMS_BIC$coefficients[-1,1]
ARcoefPval <- sarmod$coefficients[-1,4]
##MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
Interc <- ARcoef[1]
##IntercMS_AIC <- MScoef_AIC[1]
IntercMS_A=BIC <- MScoef_BIC[1]
resultAR <- list()
##resultMS_AIC <- list()
resultMS_BIC <- list()
resultAR$COEF <- ARcoef
##resultMS_AIC $COEF <- MScoef_AIC
resultMS_BIC $COEF <- MScoef_BIC
resultAR$PVAL <- ARcoefPval
##resultMS_AIC $PVAL <- MScoefPval_AIC
resultMS_BIC $PVAL <- MScoefPval_BIC
resultAR$INTERC <- Interc
##resultMS_AIC $INTERC <- IntercMS_AIC
resultMS_BIC$INTERC <- IntercMS_BIC
##resultAR$sAIC <- resultAIC
##resultAR$sMS_AIC <- resultMS_AIC
resultAR$sMS_BIC <- resultMS_BIC
##resultAR$sBIC <- resultBIC
resultAR$RESID <- armod$residuals
resultAR$FITTED <- armod$fitted.values
resultAR$STD.DEV <- sarmod$sigma
resultAR$R.SQUARED <- sarmod$r.squared
##resultMS_AIC$RESID <- ms_AIC$residuals
##resultMS_AIC$FITTED <- ms_AIC$fitted.values
##resultMS_AIC$STD.DEV <- sMS_AIC $sigma
##resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
resultMS_BIC$RESID <- ms_BIC$residuals
resultMS_BIC$FITTED <- ms_BIC$fitted.values
resultMS_BIC$STD.DEV <- sMS_BIC $sigma
resultMS_BIC$R.SQUARED <- sMS_BIC $r.squared
return(resultMS_BIC)
})
model_results <-lapply(1:nROI, function(i){
cat(i)
y <- Yfit[,i]
y <- unname(y)
x <- Xfit
complete_obs <- complete.cases(y) & complete.cases(x)
y <- y[complete_obs]
x <- x[complete_obs,]
data1 <-data.frame(cbind(y,x))
##y_name <- gsub( "_.*$", "", names(data1)[i] )
colnames(x) <- gsub("-",".",colnames(x))
##names(data1)[1] <- y_name
model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))
## armod <- lm(y ~ x)
lw <- formula(paste("~",paste(colnames(x)[1:10],collapse = "+")))
up <- formula(paste("~",paste(colnames(x)[1:70],collapse = "+")))
armod <- lm(model_formulae, data=data1)
##ms_AIC <- mStep(armod, scope = list(upper=up),direction="both", k=2)
ms_BIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
##AIC <- step(armod,direction = "both", k = 2)
##LARS <- lars(model_formulae, type="lasso")
##BIC <- step(armod,direction = "both", k = log(70))
##AIC <- stepAIC(armod, direction=c("both"))
##BIC <- BIC(armod, direction=c("both"))
##sAIC <- summary(AIC)
##sMS_AIC <- summary(ms_AIC)
sMS_BIC <- summary(ms_BIC)
##sBIC <- summary(BIC)
##sarmod <- summary(armod)
##resultAIC <- list()
##resultBIC <- list()
##resultMS_AIC <- list()
resultMS_BIC <- list()
##resultAIC <- sAIC
##resultMS_AIC <- sMS_AIC
resultMS_BIC <- sMS_BIC
##resultBIC <- sBIC
sarmod <- summary(armod)
ARcoef <- sarmod$coefficients[-1,1]
##MScoef_AIC <- sMS_AIC$coefficients[-1,1]
MScoef_BIC <- sMS_BIC$coefficients[-1,1]
ARcoefPval <- sarmod$coefficients[-1,4]
##MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
Interc <- ARcoef[1]
##IntercMS_AIC <- MScoef_AIC[1]
IntercMS_BIC <- MScoef_BIC[1]
resultAR <- list()
##resultMS_AIC <- list()
resultMS_BIC <- list()
resultAR$COEF <- ARcoef
##resultMS_AIC $COEF <- MScoef_AIC
resultMS_BIC $COEF <- MScoef_BIC
resultAR$PVAL <- ARcoefPval
##resultMS_AIC $PVAL <- MScoefPval_AIC
resultMS_BIC $PVAL <- MScoefPval_BIC
resultAR$INTERC <- Interc
##resultMS_AIC $INTERC <- IntercMS_AIC
resultMS_BIC$INTERC <- IntercMS_BIC
##resultAR$sAIC <- resultAIC
##resultAR$sMS_AIC <- resultMS_AIC
resultAR$sMS_BIC <- resultMS_BIC
##resultAR$sBIC <- resultBIC
resultAR$RESID <- armod$residuals
resultAR$FITTED <- armod$fitted.values
resultAR$STD.DEV <- sarmod$sigma
resultAR$R.SQUARED <- sarmod$r.squared
##resultMS_AIC$RESID <- ms_AIC$residuals
##resultMS_AIC$FITTED <- ms_AIC$fitted.values
##resultMS_AIC$STD.DEV <- sMS_AIC $sigma
##resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
resultMS_BIC$RESID <- ms_BIC$residuals
resultMS_BIC$FITTED <- ms_BIC$fitted.values
resultMS_BIC$STD.DEV <- sMS_BIC $sigma
resultMS_BIC$R.SQUARED <- sMS_BIC $r.squared
return(resultMS_BIC)
})
VAR_mat_BIC <- do.call(rbind, lapply(model_results, function(x) x$COEF))
dimnames(VAR_mat_BIC) <-  list(colnames(Yfit),colnames(Yfit))
View(VAR_mat_BIC)
return(resultMS_BIC)
VAR_mat_BIC <- do.call(rbind, lapply(model_results, function(x) x$COEF))
View(VAR_mat_BIC)
model_results <-lapply(1:nROI, function(i){
cat(i)
y <- Yfit[,i]
y <- unname(y)
x <- Xfit
complete_obs <- complete.cases(y) & complete.cases(x)
y <- y[complete_obs]
x <- x[complete_obs,]
data1 <-data.frame(cbind(y,x))
##y_name <- gsub( "_.*$", "", names(data1)[i] )
colnames(x) <- gsub("-",".",colnames(x))
##names(data1)[1] <- y_name
model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))
## armod <- lm(y ~ x)
lw <- formula(paste("~",paste(colnames(x)[1:10],collapse = "+")))
up <- formula(paste("~",paste(colnames(x)[1:70],collapse = "+")))
armod <- lm(model_formulae, data=data1)
##ms_AIC <- mStep(armod, scope = list(upper=up),direction="both", k=2)
ms_BIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
##AIC <- step(armod,direction = "both", k = 2)
##LARS <- lars(model_formulae, type="lasso")
##BIC <- step(armod,direction = "both", k = log(70))
##AIC <- stepAIC(armod, direction=c("both"))
##BIC <- BIC(armod, direction=c("both"))
##sAIC <- summary(AIC)
##sMS_AIC <- summary(ms_AIC)
sMS_BIC <- summary(ms_BIC)
##sBIC <- summary(BIC)
##sarmod <- summary(armod)
##resultAIC <- list()
##resultBIC <- list()
##resultMS_AIC <- list()
resultMS_BIC <- list()
##resultAIC <- sAIC
##resultMS_AIC <- sMS_AIC
##resultMS_BIC <- sMS_BIC
##resultBIC <- sBIC
sarmod <- summary(armod)
ARcoef <- sarmod$coefficients[-1,1]
##MScoef_AIC <- sMS_AIC$coefficients[-1,1]
MScoef_BIC <- sMS_BIC$coefficients[-1,1]
ARcoefPval <- sarmod$coefficients[-1,4]
##MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
Interc <- ARcoef[1]
##IntercMS_AIC <- MScoef_AIC[1]
IntercMS_BIC <- MScoef_BIC[1]
resultAR <- list()
##resultMS_AIC <- list()
##resultMS_BIC <- list()
resultAR$COEF <- ARcoef
##resultMS_AIC $COEF <- MScoef_AIC
resultMS_BIC$COEF <- MScoef_BIC
resultAR$PVAL <- ARcoefPval
##resultMS_AIC $PVAL <- MScoefPval_AIC
resultMS_BIC$PVAL <- MScoefPval_BIC
resultAR$INTERC <- Interc
##resultMS_AIC $INTERC <- IntercMS_AIC
resultMS_BIC$INTERC <- IntercMS_BIC
##resultAR$sAIC <- resultAIC
##resultAR$sMS_AIC <- resultMS_AIC
##resultAR$sMS_BIC <- resultMS_BIC
##resultAR$sBIC <- resultBIC
resultAR$RESID <- armod$residuals
resultAR$FITTED <- armod$fitted.values
resultAR$STD.DEV <- sarmod$sigma
resultAR$R.SQUARED <- sarmod$r.squared
##resultMS_AIC$RESID <- ms_AIC$residuals
##resultMS_AIC$FITTED <- ms_AIC$fitted.values
##resultMS_AIC$STD.DEV <- sMS_AIC $sigma
##resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
resultMS_BIC$RESID <- ms_BIC$residuals
resultMS_BIC$FITTED <- ms_BIC$fitted.values
resultMS_BIC$STD.DEV <- sMS_BIC $sigma
resultMS_BIC$R.SQUARED <- sMS_BIC $r.squared
return(resultMS_BIC)
})
savehistory("~/Desktop/1.r")
