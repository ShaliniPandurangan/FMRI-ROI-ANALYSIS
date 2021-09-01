#load("/Users/shalinipandurangan/Desktop/04-Data/fMRI-ROI-time-series/fMRI-ROI-time-series.RData")
#load("D:/WORK/SUPERVISION/Shalini/04-Data/fMRI-ROI-time-series/fMRI-ROI-time-series/fMRI-ROI-time-series.RData")
# Y is the data frame of fMRI-ROI-time-series data
# y values of all the subjects of all 70 regions of interest
# Sub_covariates <- read_delim("/Users/usi/Desktop/04-Data/covariates/SUBJ-covariates.txt",delim = " ")
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
library(qgraph)
library(HDeconometrics)
library(abind)
library(fastDummies)
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
Yfit_SC <- do.call(rbind, lapply(Data_YX, function(x) x$Y))
dimnames(Yfit)[[1]] <- NULL
Xfit_SC <- do.call(rbind, lapply(Data_YX, function(x) x$X))
dimnames(Xfit) <- list(NULL, paste(paste(colnames(Yfit),"_lag_",sep = ""), 1:lag, sep = ""))

## Linear Regression, Linear regression model selection with AIC and BIC
model_results <-lapply(1:nROI, function(i){
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
  armod <- lm(model_formulae, data=data1)
  ms_AIC <- stepAIC(armod,direction="both")
  ms_BIC <- stepAIC(armod,direction="both", k= log(70))
  sMS_AIC <- summary(ms_AIC)
  sMS_BIC <- summary(ms_BIC)
  sarmod <- summary(armod)
  ARcoef <- sarmod$coefficients[-1,1]
  ARcoefPval <- sarmod$coefficients[-1,4]
  Interc <- ARcoef[1]
  MScoef_AIC <- sMS_AIC$coefficients[-1,1]
  MScoefPval_AIC <- sMS_AIC$coefficients[-1,4]
  IntercMS_AIC <- MScoef_AIC[1]
  MScoef_BIC <- sMS_BIC$coefficients[-1,1]
  MScoefPval_BIC <- sMS_BIC$coefficients[-1,4]
  IntercMS_BIC <- MScoef_BIC[1]
  resultAR <- list()
  resultMS_AIC <- list()
  resultMS_BIC <- list()
  resultAR$COEF <- ARcoef
  resultAR$PVAL <- ARcoefPval
  resultAR$INTERC <- Interc
  resultAR$RESID <- armod$residuals
  resultAR$FITTED <- armod$fitted.values
  resultAR$STD.DEV <- sarmod$sigma
  resultAR$R.SQUARED <- sarmod$r.squared
  resultMS_AIC $COEF <- MScoef_AIC
  resultMS_AIC $PVAL <- MScoefPval_AIC
  resultMS_AIC $INTERC <- IntercMS_AIC
  resultMS_AIC$RESID <- ms_AIC$residuals
  resultMS_AIC$FITTED <- ms_AIC$fitted.values
  resultMS_AIC$STD.DEV <- sMS_AIC $sigma
  resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
  resultMS_BIC $COEF <- MScoef_BIC
  resultMS_BIC $PVAL <- MScoefPval_BIC
  resultMS_BIC $INTERC <- IntercMS_AIC
  resultMS_BIC$RESID <- ms_BIC$residuals
  resultMS_BIC$FITTED <- ms_BIC$fitted.values
  resultMS_BIC$STD.DEV <- sMS_BIC $sigma
  resultMS_BIC$R.SQUARED <- sMS_BIC $r.squared
  return(list(resultMS_AIC,resultMS_BIC,resultAR))
})

## VAR matrix after the model selection methods using AIC and BIC
VAR_mat_AIC <- do.call(rbind, lapply(model_results[[1]], function(x) x["COEF"][[1]]))

VAR_mat_BIC <- do.call(rbind, lapply(model_results[[2]], function(x) x["COEF"][[1]]))

VAR_mat <- do.call(rbind, lapply(model_results[[3]], function(x) x["COEF"][[1]]))


dimnames(VAR_mat_AIC) <-  list(colnames(Yfit),colnames(Yfit))

dimnames(VAR_mat_BIC) <-  list(colnames(Yfit),colnames(Yfit))

dimnames(VAR_mat) <-  list(colnames(Yfit),colnames(Yfit))


## Improving sparcity using P Values. 
VAR_mat_Pval_AIC <- do.call(rbind, lapply(model_results[[1]], function(x) x["PVAL"][[1]]))
VAR_mat_AIC[!VAR_mat_Pval_AIC < 0.05] <- 0

VAR_mat_Pval_BIC <- do.call(rbind, lapply(model_results[[2]], function(x) x["PVAL"][[1]]))
VAR_mat_BIC[!VAR_mat_Pval_BIC < 0.05] <- 0

VAR_mat_Pval <- do.call(rbind, lapply(model_results[[3]], function(x) x["PVAL"][[1]]))
VAR_mat[!VAR_mat_Pval < 0.05] <- 0

## Transpose to get the adjacency 
VAR_mat_AIC <-  t(VAR_mat_AIC)

VAR_mat_BIC <-  t(VAR_mat_BIC)

VAR_mat <-  t(VAR_mat)

## Get the covariance matrix using the model selection methods AIC and BIC 
Sigma_diag_AIC <- diag(sapply(model_results[[1]], function(x) x["STD.DEV"][[1]])^2)

Sigma_diag_BIC <- diag(sapply(model_results[[2]], function(x) x["STD.DEV"][[1]])^2)

Sigma_diag <- diag(sapply(model_results[[3]], function(x) x["STD.DEV"][[1]])^2)


## Get the residuals using the model selection methods AIC and BIC 
residuals_VAR_AIC <-  do.call(cbind,lapply(model_results[[1]], function(x) x["RESID"][[1]]))

residuals_VAR_BIC <-  do.call(cbind,lapply(model_results[[2]], function(x) x["RESID"][[1]]))

residuals_VAR <-  do.call(cbind,lapply(model_results[[3]], function(x) x["RESID"][[1]]))


## Get the sample covariance matrix of the residuals using the model selection methods AIC and BIC 
Sigma_residuals_AIC <-  cov(do.call(cbind,lapply(model_results[[1]], function(x) x["RESID"][[1]])))

Sigma_residuals_BIC <-  cov(do.call(cbind,lapply(model_results[[2]], function(x) x["RESID"][[1]])))

Sigma_residuals <-  cov(do.call(cbind,lapply(model_results[[3]], function(x) x["RESID"][[1]])))


model_results_ridge <-lapply(1:nROI, function(i){
  cat(i)
  y <- Yfit[,i]
  y <- unname(y)  
  x <- Xfit
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  ## Ridge regression cross validation 
  cv_ridge <- cv.glmnet(x, y, alpha = 0)
  print(cv_ridge$lambda.min)
  model_ridge <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)
  print(coef(model_ridge))
  coef_ridge <- coef(model_ridge)[-1]
  coef_ridge_reg <- list()
  coef_ridge_reg $COEF <- coef_ridge
  return(coef_ridge_reg)
})


model_results_lasso <-lapply(1:nROI, function(i){
  cat(i)
  y <- Yfit[,i]
  y <- unname(y) 
  x <- Xfit
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  ## lasso regression cross validation 
  cv_lasso <- cv.glmnet(x, y, alpha = 1)
  print(cv_lasso$lambda.min)
  model_lasso <- glmnet(x, y, alpha = 1, lambda = cv_lasso$lambda.min)
  print(coef(model_lasso))
  coef_lasso <- coef(model_lasso)[-1]
  coef_lasso_reg <- list()
  coef_lasso_reg $COEF <- coef_lasso
  return(coef_lasso_reg)
})

model_results_lasso_BIC <-lapply(1:nROI, function(i){
  cat(i)
  y <- Yfit[,i]
  y <- unname(y) 
  x <- Xfit
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  ## lasso regression BIC
  lasso_bic=ic.glmnet(x,y)
  print(lasso_bic$lambda[length(lasso_bic$lambda)])
  model_lasso_bic <- ic.glmnet(x, y, alpha = 1, lambda = lasso_bic$lambda[length(lasso_bic$lambda)])
  print(model_lasso_bic$nvar)
  coef_lasso_BIC <- coef(model_lasso_bic)[-1]
  coef_lasso_BIC_reg <- list()
  n_vars_lasso_BIC <- list()
  n_vars_lasso_BIC$val <- model_lasso_bic$nvar
  ##print(model_lasso_bic$nvar)
  coef_lasso_BIC_reg $COEF <- coef_lasso_BIC
  return(coef_lasso_BIC_reg)
})

## plotting random plot of mean squared error for a range of lambdas with alpha. 
for (i in 1:5){
    j <- sample.int(nROI, 1)
    y <- Yfit[,j]
    y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    cv1=cv.glmnet(x,y,nfolds=10,alpha=1)
    plot(cv1)
    cv.5=cv.glmnet(x,y,nfolds=10,alpha=.5)
    cv0=cv.glmnet(x,y,nfolds=10,alpha=0)
    plot(cv1);plot(cv.5);plot(cv0)
    plot(log(cv1$lambda),cv1$cvm,pch=19,col="green",xlab="log(Lambda)",ylab=cv1$name)
    points(log(cv.5$lambda),cv.5$cvm,pch=19,col="pink")
    points(log(cv0$lambda),cv0$cvm,pch=19,col="yellow")
    legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("green","pink","yellow"))
}


## plotting network grqphs 

qgraph_net <- qgraph(VAR_mat,
                     posCol="blue",
                     negCol="purple")

qgraph_net <- qgraph(VAR_mat_AIC,
                     posCol="blue",
                     negCol="purple")

qgraph_net <- qgraph(VAR_mat_BIC,
                     posCol="blue",
                     negCol="purple")



# model_results_cov <-lapply(1:nROI, function(i){
#     cat(i)
#     y <- Yfit[,i]
#     y <- unname(y) 
#     print(y)
#     x <- Xfit
#     complete_obs <- complete.cases(y) & complete.cases(x)
#     y <- y[complete_obs]
#     x <- x[complete_obs,]
#     data1 <-data.frame(cbind(y,x))
#     ##model_formulae <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",paste("factor(",(dimnames(x)[[2]][72]),")"),"+",paste("factor(",(dimnames(x)[[2]][73]),")")))
#     armod <- lm(model_formulae, data=data1)
#     ##up <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",paste("factor(",(dimnames(x)[[2]][72]),")"),"+",paste("factor(",(dimnames(x)[[2]][73]),")")))
#     ##armod_cov <- mStep(armod, scope = list(upper=up),direction="both", k=2)
#     sarmod <- summary(armod)
#     ARcoef <- sarmod$coefficients[-1,1]
#     ARcoefPval <- sarmod$coefficients[-1,4]
#     Interc <- ARcoef[1]
#     resultAR <- list()
#     resultAR$COEF <- ARcoef
#     resultAR$PVAL <- ARcoefPval
#     resultAR$INTERC <- Interc
#     resultAR$RESID <- armod$residuals
#     resultAR$FITTED <- armod$fitted.values
#     resultAR$STD.DEV <- sarmod$sigma
#     resultAR$R.SQUARED <- sarmod$r.squared
#     return(resultAR)                
# })

## model results with covariates entire population
model_results_cov <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit_SC)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit_SC)[i]
    y <- Yfit_SC[,i]
    #y <- unname(y) 
    x <- Xfit_SC
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit_SC[complete_obs,i],1,x))
    # weights = rep(7/17.nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    data_fit$handedness <- factor(data_fit$handedness)
    model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
    model_lower <- lm(model_fl_lower, data_fit)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit)
    
    ms_BIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower), 
                      k= log(70))
    sbic<- summary(ms_BIC)
    resultBIC <- list()
    resultBIC$COEF <- sbic$coefficients[-1,1]
    resultBIC$PVAL <- sbic$coefficients[-1,4]
    resultBIC$INTERC <- sbic$coefficients[1,1]
    resultBIC$RESID <- ms_BIC$residuals
    resultBIC$FITTED <- ms_BIC$fitted.values
    resultBIC$STD.DEV <- sbic$sigma
    resultBIC$R.SQUARED <- sbic$r.squared
    return(resultBIC)                
})

# ## diagnosed group 
# model_results_diag <-lapply(1:nROI, function(i){
#     cat(i)
#     y <- Yfit_diagnosed[,i]
#     y <- unname(y) 
#     x <- Xfit_diagnosed
#     complete_obs <- complete.cases(y) & complete.cases(x)
#     y <- y[complete_obs]
#     x <- x[complete_obs,]
#     data1 <-data.frame(cbind(y,x))
#     model_formulae <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",dimnames(x)[[2]][72],"+",dimnames(x)[[2]][73]))
#     armod <- lm(model_formulae, data=data1)
#     up <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",dimnames(x)[[2]][72],"+",dimnames(x)[[2]][73]))
#     ms_BIC <- mStep(armod, scope = list(upper=up),direction="both", k= log(70))
#     ##up <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",paste("factor(",(dimnames(x)[[2]][72]),")"),"+",paste("factor(",(dimnames(x)[[2]][73]),")")))
#     ##armod_cov <- mStep(armod, scope = list(upper=up),direction="both", k=2)
#     sbic<- summary(ms_BIC)
#     ARcoef <- sbic$coefficients[-1,1]
#     ARcoefPval <- sbic$coefficients[-1,4]
#     Interc <- ARcoef[1]
#     resultAR <- list()
#     resultAR$COEF <- ARcoef
#     resultAR$PVAL <- ARcoefPval
#     resultAR$INTERC <- Interc
#     resultAR$RESID <- ms_BIC$residuals
#     resultAR$FITTED <- ms_BIC$fitted.values
#     resultAR$STD.DEV <- sbic$sigma
#     resultAR$R.SQUARED <- sbic$r.squared
#     return(resultAR)                
# })

## non_diagnosed group 
model_results_non_diag <-lapply(1:nROI, function(i){
    cat(i)
    y <- Yfit_non_diagnosed[,i]
    y <- unname(y) 
    x <- Xfit_non_diagnosed
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    data1 <-data.frame(cbind(y,x))
    model_formulae <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",dimnames(x)[[2]][72],"+",dimnames(x)[[2]][73]))
    armod <- lm(model_formulae, data=data1)
    up <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",dimnames(x)[[2]][72],"+",dimnames(x)[[2]][73]))
    armod_cov <- mStep(armod, scope = list(upper=up),direction="both", k=2)
    sarmod <- summary(armod_cov)
    ARcoef <- sarmod$coefficients[-1,1]
    ARcoefPval <- sarmod$coefficients[-1,4]
    Interc <- ARcoef[1]
    resultAR <- list()
    resultAR$COEF <- ARcoef
    resultAR$PVAL <- ARcoefPval
    resultAR$INTERC <- Interc
    resultAR$RESID <- armod_cov$residuals
    resultAR$FITTED <- armod_cov$fitted.values
    resultAR$STD.DEV <- sarmod$sigma
    resultAR$R.SQUARED <- sarmod$r.squared
    return(resultAR)                
})

VAR_mat_NON_DIA <- do.call(rbind, lapply(model_results_non_diag, function(x) x[["COEF"]]))

VAR_mat_BIC <- do.call(rbind, lapply(model_results[[2]], function(x) x["COEF"][[1]]))

VAR_mat <- do.call(rbind, lapply(model_results[[3]], function(x) x["COEF"][[1]]))



## Linear Regression, Linear regression model selection with AIC and BIC
model_results <-lapply(1:nROI, function(i){
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
  armod <- lm(model_formulae, data=data1)
  ms_AIC <- stepAIC(armod,direction="both")
  ms_BIC <- stepAIC(armod,direction="both", k= log(70))
  sMS_AIC <- summary(ms_AIC)
  sMS_BIC <- summary(ms_BIC)
  sarmod <- summary(armod)
  resultAR <- list()
  resultMS_AIC <- list()
  resultMS_BIC <- list()
  resultAR$COEF <- sarmod$coefficients[-1,1]
  resultAR$PVAL <- sarmod$coefficients[-1,4]
  resultAR$INTERC <- ARcoef[1]
  resultAR$RESID <- armod$residuals
  resultAR$FITTED <- armod$fitted.values
  resultAR$STD.DEV <- sarmod$sigma
  resultAR$R.SQUARED <- sarmod$r.squared
  resultMS_AIC $COEF <- sMS_AIC$coefficients[-1,1]
  resultMS_AIC $PVAL <- sMS_AIC$coefficients[-1,4]
  resultMS_AIC $INTERC <- MScoef_AIC[1]
  resultMS_AIC$RESID <- ms_AIC$residuals
  resultMS_AIC$FITTED <- ms_AIC$fitted.values
  resultMS_AIC$STD.DEV <- sMS_AIC $sigma
  resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
  resultMS_BIC $COEF <- sMS_BIC$coefficients[-1,1]
  resultMS_BIC $PVAL <- sMS_BIC$coefficients[-1,4]
  resultMS_BIC $INTERC <- MScoef_BIC[1]
  resultMS_BIC$RESID <- ms_BIC$residuals
  resultMS_BIC$FITTED <- ms_BIC$fitted.values
  resultMS_BIC$STD.DEV <- sMS_BIC $sigma
  resultMS_BIC$R.SQUARED <- sMS_BIC $r.squared
  return(list(resultMS_AIC,resultMS_BIC,resultAR))
})



## diagnosed group 
model_results_diag <-lapply(1:nROI, function(i){
    cat(i)
    y <- Yfit_diagnosed[,i]
    y <- unname(y) 
    x <- Xfit_diagnosed
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    data1 <-data.frame(cbind(y,x))
    model_formulae <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",dimnames(x)[[2]][72],"+",dimnames(x)[[2]][73]))
    armod <- lm(model_formulae, data1)
    ms_AIC <- stepAIC(armod,direction="both")    
    ms_BIC <- stepAIC(armod,direction="both", k= log(70))
    ##up <- paste("y","~",paste(colnames(x[,1:70]),collapse = "+"),"+",paste("factor(",(dimnames(x)[[2]][71]),")","+",paste("factor(",(dimnames(x)[[2]][72]),")"),"+",paste("factor(",(dimnames(x)[[2]][73]),")")))
    ##armod_cov <- mStep(armod, scope = list(upper=up),direction="both", k=2)
    sbic<- summary(ms_BIC)
    resultAR <- list()
    resultAR$COEF <- sbic$coefficients[-1,1]
    resultAR$PVAL <- sbic$coefficients[-1,4]
    resultAR$INTERC <- sbic$coefficients[-1,1][1]
    resultAR$RESID <- ms_BIC$residuals
    resultAR$FITTED <- ms_BIC$fitted.values
    resultAR$STD.DEV <- sbic$sigma
    resultAR$R.SQUARED <- sbic$r.squared
    return(resultAR)                
})

## diagnosed group(FINAL)
model_results_diag <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit_diagnosed)[i]
    y <- Yfit_diagnosed[,i]
    #y <- unname(y) 
    x <- Xfit_diagnosed
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit_diagnosed[complete_obs,i],1,x))
    # weights = rep(7/17.nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    data_fit$handedness <- factor(data_fit$handedness)
    model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
    model_lower <- lm(model_fl_lower, data_fit)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit)
    
    ms_BIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower), 
                      k= log(70))
    sbic<- summary(ms_BIC)
    resultAR <- list()
    resultAR$COEF <- sbic$coefficients[-1,1]
    resultAR$PVAL <- sbic$coefficients[-1,4]
    resultAR$INTERC <- sbic$coefficients[-1,1][1]
    resultAR$RESID <- ms_BIC$residuals
    resultAR$FITTED <- ms_BIC$fitted.values
    resultAR$STD.DEV <- sbic$sigma
    resultAR$R.SQUARED <- sbic$r.squared
    return(resultAR)                
})


## diagnosed group(FINAL)
model_results_non_diag <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit_diagnosed)[i]
    y <- Yfit_diagnosed[,i]
    #y <- unname(y) 
    x <- Xfit_diagnosed
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit_diagnosed[complete_obs,i],1,x))
    # weights = rep(7/17.nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    data_fit$handedness <- factor(data_fit$handedness)
    model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
    model_lower <- lm(model_fl_lower, data_fit)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit)
    
    ms_BIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower), 
                      k= log(70))
    sbic<- summary(ms_BIC)
    resultAR <- list()
    resultAR$COEF <- sbic$coefficients[-1,1]
    resultAR$PVAL <- sbic$coefficients[-1,4]
    resultAR$INTERC <- sbic$coefficients[-1,1][1]
    resultAR$RESID <- ms_BIC$residuals
    resultAR$FITTED <- ms_BIC$fitted.values
    resultAR$STD.DEV <- sbic$sigma
    resultAR$R.SQUARED <- sbic$r.squared
    return(resultAR)                
})


## Linear Regression, penalised Linear regression lasso without covariates
model_results_lasso <-lapply(1:nROI, function(i){
    cat(i)
    y <- Yfit[,i]
    y <- unname(y)  
    x <- Xfit[,1:73]
    weig <- Xfit[,74]
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    x_weig <- weig[complete_obs]
    fit = cv.glmnet(x, y, alpha = 1, weights = x_weig , nlambda = 20,type.measure = "mse", nfolds = 20)
    sarmod <- fit
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    return(list(resultAR))
})

NETWORK_LASSO <- sapply(model_results_lasso, function(y){
  coef_result <- rep(0,length(colnames(Xfit)))
  names(coef_result) <- model_results_diag_lasso[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_diag_lasso[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
  })

## Linear Regression, Linear regression model selection with AIC and BIC(FINAL)
model_results <-lapply(1:nROI, function(i){
    cat(i)
    y <- Yfit[,i]
    y <- unname(y)  
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    armod <- lm(y~., as.data.frame(cbind(y,1,x)))
    ms_AIC <- stepAIC(armod,direction="both")
    ms_BIC <- stepAIC(armod,direction="both", k= log(70))
    sMS_AIC <- summary(ms_AIC)
    sMS_BIC <- summary(ms_BIC)
    sarmod <- summary(armod)
    resultAR <- list()
    resultMS_AIC <- list()
    resultMS_BIC <- list()
    resultAR$COEF <- sarmod$coefficients[-1,1]
    resultAR$PVAL <- sarmod$coefficients[-1,4]
    resultAR$INTERC <- ARcoef[1]
    resultAR$RESID <- armod$residuals
    resultAR$FITTED <- armod$fitted.values
    resultAR$STD.DEV <- sarmod$sigma
    resultAR$R.SQUARED <- sarmod$r.squared
    resultMS_AIC $COEF <- sMS_AIC$coefficients[-1,1]
    resultMS_AIC $PVAL <- sMS_AIC$coefficients[-1,4]
    resultMS_AIC $INTERC <- MScoef_AIC[1]
    resultMS_AIC$RESID <- ms_AIC$residuals
    resultMS_AIC$FITTED <- ms_AIC$fitted.values
    resultMS_AIC$STD.DEV <- sMS_AIC $sigma
    resultMS_AIC$R.SQUARED <- sMS_AIC $r.squared
    resultMS_BIC $COEF <- sMS_BIC$coefficients[-1,1]
    resultMS_BIC $PVAL <- sMS_BIC$coefficients[-1,4]
    resultMS_BIC $INTERC <- MScoef_BIC[1]
    resultMS_BIC$RESID <- ms_BIC$residuals
    resultMS_BIC$FITTED <- ms_BIC$fitted.values
    resultMS_BIC$STD.DEV <- sMS_BIC $sigma
    resultMS_BIC$R.SQUARED <- sMS_BIC $r.squared
    return(list(resultMS_AIC,resultMS_BIC,resultAR))
})