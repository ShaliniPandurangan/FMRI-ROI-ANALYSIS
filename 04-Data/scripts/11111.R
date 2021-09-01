library(sna)
library(igraph)
library(networktools)
library(network)
library(intergraph)

VAR_MAT_NET <- graph_from_adjacency_matrix(VAR_mat,mode="directed",weighted = TRUE)

centralization.degree(VAR_MAT_NET)

centralization.degree(VAR_MAT_NET)

centralization.degree(VAR_MAT_NET,loops = FALSE, mode = "all")


dum <- array(NA, c(1,404,24))

Data_Y_SC <- abind(Data_Y,dum,along = 1)
dimnames(Data_Y_SC)[[1]][71] <- "sub_id"
colnames(x)

BIC_NETWORK <- sapply(model_results, function(y){
  coef_result <- rep(0,length(colnames(x)))
  names(coef_result) <- colnames(x)
  res <- y[[2]]$COEF
  indexes <- match(names(res),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_PVAL <- sapply(model_results, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed)[1:70]))
  names(coef_result) <- colnames(Xfit_diagnosed)[1:70]
  res <- y[[2]]$PVAL
  names(res) <- stri_sub(names(res), 2, -2)
  indexes <- match(names(res),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_1 <- BIC_NETWORK

BIC_NETWORK_1[!BIC_NETWORK_PVAL < 0.05] <- 0

match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR))[1]
Match <- match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR))

length(match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR)))

le <- length(match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR)))

## for key_value matching

for(i in 2:70){
  for(j in 1:length(match(names(model_results[[i]][[1]][["COEF"]]),colnames(VAR_MAR)))){
    Match <- match(names(model_results[[i]][[1]][["COEF"]]),colnames(VAR_MAR))
    VAR_MAR[i, colnames(VAR_MAR)[Match[j]]] <- model_results[[i]][[1]][["COEF"]][[j]]}}


## Diagnosed
model_results_diag <-lapply(1:2, function(i){
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



BIC_NETWORK_DIAG <- sapply(model_results_diag, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_DIAG_PVAL <- sapply(model_results_diag, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

## diagnosed cross validation lasso 
dum <- array(NA, c(2828,1))

Xfit_diagnosed <- abind(Xfit_diagnosed,dum,along = 2)
dimnames(Xfit_diagnosed)[[2]][74] <- "weights"
Xfit_diagnosed[,74] <- rep(7/17, length(Xfit_diagnosed[,74]))
fit = glmnet(x, y, alpha = 0.2, weights = c(rep(1,50),rep(2,50)), nlambda = 20)

## diagnosed cross validation lasso 
model_results_diag_lasso <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_diagnosed)[i]
  y <- Yfit_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_diagnosed[,1:73]
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  fit = cv.glmnet(x, y, alpha = 1, weights = rep(7/17,dim(x)[1]), nlambda = 20,type.measure = "mse", nfolds = 20)
  resultAR <- list()
  resultAR$COEF <- coef(fit, s = "lambda.min")
  resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
  resultAR$fit <- fit
  return(resultAR)                
})


NETWORK_DIAG_LASSO <- sapply(model_results_diag_lasso, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed)))
  names(coef_result) <- model_results_diag_lasso[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_diag_lasso[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

# non diagnosed cross validation lasso 
model_results_non_diag_lasso <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_non_diagnosed)[i]
  y <- Yfit_non_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_non_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  fit = cv.glmnet(x, y, alpha = 1, nlambda = 20,type.measure = "mse", nfolds = 20)
  resultAR <- list()
  resultAR$COEF <- coef(fit, s = "lambda.min")
  resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
  resultAR$fit <- fit
  return(resultAR)                
})

NETWORK_NON_DIAG_LASSO <- sapply(model_results_non_diag_lasso, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed)))
  names(coef_result) <- model_results_non_diag_lasso[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_non_diag_lasso[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

## non diagnosed 
model_results_non_diag <-lapply(1:2, function(i){
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

BIC_NETWORK_NON_DIAG <- sapply(model_results_non_diag, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+2)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_NON_DIAG_PVAL <- sapply(model_results_non_diag, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+2)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})


