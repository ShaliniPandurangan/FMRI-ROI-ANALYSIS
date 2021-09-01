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

BIC_NETWORK_DIAG <- sapply(model_results, function(y){
  coef_result <- rep(0,length(colnames(x)))
  names(coef_result) <- colnames(x)
  res <- y[[2]]$COEF
  indexes <- match(names(res),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

VAR_MAT_COV_BIC <- sapply(model_results_cov, function(y){
  coef_result <- rep(0,length(colnames(Xfit_SC))+2)
  names(coef_result) <- paste( "`",colnames(Xfit_SC),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result  
})

VAR_MAT_COV_BIC_PVAL <- sapply(model_results_cov, function(y){
  coef_result <- rep(0,length(colnames(Xfit_SC))+2)
  names(coef_result) <- paste( "`",colnames(Xfit_SC),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

names(model_results_cov[[1]][["COEF"]])


match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR))[1]
Match <- match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR))

length( match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR)))

le <- length( match(names(model_results[[5]][[1]][["COEF"]]),colnames(VAR_MAR)))

## for key_value matching

for(i in 2:70){
  for(j in 1:length(match(names(model_results[[i]][[1]][["COEF"]]),colnames(VAR_MAR)))){
    Match <- match(names(model_results[[i]][[1]][["COEF"]]),colnames(VAR_MAR))
    VAR_MAR[i, colnames(VAR_MAR)[Match[j]]] <- model_results[[i]][[1]][["COEF"]][[j]]}}


## diagnosed AIC
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


## non diagnosed BIC
model_results_non_diag <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_non_diagnosed)[i]
  y <- Yfit_non_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_non_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_non_diagnosed[complete_obs,i],1,x))
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