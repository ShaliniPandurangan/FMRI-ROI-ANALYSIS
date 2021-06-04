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

dum <- array(NA, c(3,404,24))

Data_Y_SC <- abind(Data_Y,dum,along = 1)
dimnames(Data_Y_SC)[[1]][71] <- "handedness"
dimnames(Data_Y_SC)[[1]][72] <- "current_diagnosis"
dimnames(Data_Y_SC)[[1]][73] <- "lifetime_diagnosis"

Sub_covariates$handedness.factor <- factor(Sub_covariates$handedness)
Sub_covariates$current_diagnosis.factor <- factor(Sub_covariates$current_diagnosis)
Sub_covariates$lifetime_diagnosis.factor <- factor(Sub_covariates$lifetime_diagnosis)

for (i in 1:nPersons){
Data_Y_SC[71,,i] <- as.numeric(Sub_covariates$handedness.factor[i])
Data_Y_SC[72,,i] <- as.numeric(Sub_covariates$current_diagnosis.factor[i])
Data_Y_SC[73,,i] <- as.numeric(Sub_covariates$lifetime_diagnosis.factor[i])
}

Data_diagnosed <- Data_Y[,,Sub_covariates$sub_id[!is.na(Sub_covariates$type_lifetime_diagnosis)]]

Data_non_diagnosed <- Data_Y[,,Sub_covariates$sub_id[is.na(Sub_covariates$type_lifetime_diagnosis)]]

Data_diagnosed_sc <- Data_Y_SC[,,Sub_covariates$sub_id[!is.na(Sub_covariates$type_lifetime_diagnosis)]]

Data_non_diagnosed_sc <- Data_Y_SC[,,Sub_covariates$sub_id[is.na(Sub_covariates$type_lifetime_diagnosis)]]


Data_YX_diagnosed <- lapply(1:7, function(i){
  result <- list()
  result$Y <- NA
  result$X <- NA
  get_lags <- mgm:::lagData(data = t(Data_diagnosed[,,i]), lags = 1:lag) 
  result$Y <- get_lags$data_response
  result$X <- do.call(cbind, get_lags$l_data_lags) 
  return(result)
}) 

Yfit_diagnosed <- do.call(rbind, lapply(Data_YX_diagnosed, function(x) x$Y))

Data_YX_non_diagnosed <- lapply(1:17, function(i){
  result <- list()
  result$Y <- NA
  result$X <- NA
  get_lags <- mgm:::lagData(data = t(Data_non_diagnosed[,,i]), lags = 1:lag) 
  result$Y <- get_lags$data_response
  result$X <- do.call(cbind, get_lags$l_data_lags) 
  return(result)
}) 

Yfit_non_diagnosed <- do.call(rbind, lapply(Data_YX_non_diagnosed, function(x) x$Y))

Data_YX_non_diagnosed_X <- lapply(1:17, function(i){
  result <- list()
  result$Y <- NA
  result$X <- NA
  get_lags <- mgm:::lagData(data = t(Data_non_diagnosed_sc[,,i]), lags = 1:lag) 
  result$Y <- get_lags$data_response
  result$X <- do.call(cbind, get_lags$l_data_lags) 
  return(result)
}) 

Xfit_non_diagnosed <- do.call(rbind, lapply(Data_YX_non_diagnosed_X, function(x) x$X))

Data_YX_diagnosed_X <- lapply(1:7, function(i){
  result <- list()
  result$Y <- NA
  result$X <- NA
  get_lags <- mgm:::lagData(data = t(Data_diagnosed_sc[,,i]), lags = 1:lag) 
  result$Y <- get_lags$data_response
  result$X <- do.call(cbind, get_lags$l_data_lags) 
  return(result)
}) 
Xfit_diagnosed <- do.call(rbind, lapply(Data_YX_diagnosed_X, function(x) x$X))


