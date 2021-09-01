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
library(qtlmt)

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

## Get the effective observations :
## T_effective = nTime - lag

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
model_results <-lapply(1, function(i){
  cat(i)
  y <- Yfit[,i]
  y <- unname(y)
  
  x <- Xfit
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  data1 <-data.frame(cbind(y,x))
  
  colnames(x) <- gsub("-",".",colnames(x))
  model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))
  
  ##lw <- paste("y","~",paste(colnames(x)[1:10],collapse = "+"))
  up <- paste("y","~",paste(colnames(x),collapse = "+"))
  
  armod <- lm(model_formulae, data=data1)
  
  ## AIC using qtlmt package
  ms_AIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
  sMS_AIC <- summary(ms_AIC)
  
  sarmod <- summary(armod)
  resultMS_AIC <- list()
  View(sMS_AIC)
  ##ARcoef <- sarmod$coefficients[-1,1]
  MScoef_AIC <- sMS_AIC$coefficients[-1,1]
  
  ##ARcoefPval <- sarmod$coefficients[-1,4]
  MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
  
  ##resultAR <- list()

  ##resultAR$COEF <- ARcoef
  resultMS_AIC$COEF <- MScoef_AIC
  
  ##resultAR$PVAL <- ARcoefPval
  resultMS_AIC$PVAL <- MScoefPval_AIC
  
  return(resultMS_AIC)
})