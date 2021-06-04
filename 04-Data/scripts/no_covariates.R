## WITHOUT COVARIATES

load("Yfit_NC.Rdata")
load("Xfit_NC.Rdata")
load("Xfit_NC_ND.Rdata")
load("Xfit_NC_ND.Rdata")
load("Yfit.Rdata")
load("Xfit.Rdata")


## DIAGNOSED_GROUP WITHOUT WEIGHTS BIC

model_results_diag_wo_co_nw <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_NC)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC)[i]
  y <- Yfit_NC[,i]
  #y <- unname(y) 
  x <- Xfit_NC
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  data_fit <- as.data.frame(cbind(Yfit_NC[complete_obs,i],1,x))
  
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  model_fl_lower <- as.formula(paste(yname,"~1"))
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



BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_diag_wo_co_nw, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC)))
    names(coef_result) <- paste( "`",colnames(Xfit_NC),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL <- sapply(model_results_diag_wo_co_nw, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC)))
    names(coef_result) <- paste( "`",colnames(Xfit_NC),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[!BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL < 0.05] <- 0

BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_DIAG_WEIGHTS_WITHOUT_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

Residuals_BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_diag_wo_co_nw, function(y){
    res <- y$RESID
    })
# BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- {qgraph(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[1:70,],
#                                posCol="blue",
#                                negCol="purple")
#   title("BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES")}

## graph density without weights and no covariates_BIC_diagnosed
graph.density(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)

## betweenness centrality max
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_betw(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_betw(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_degree(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_degree(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]


## DIAGNOSED_GROUP WITHOUT WEIGHTS AIC

model_results_diag_wo_co_nw_aic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_NC)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC)[i]
  y <- Yfit_NC[,i]
  #y <- unname(y) 
  x <- Xfit_NC
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  data_fit <- as.data.frame(cbind(Yfit_NC[complete_obs,i],1,x))
  
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  model_fl_lower <- as.formula(paste(yname,"~1"))
  model_lower <- lm(model_fl_lower, data_fit)
  model_fl_upper <- as.formula(paste(yname,"~."))
  model_upper <- lm(model_fl_upper, data_fit)
  
  ms_AIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower))

  saic<- summary(ms_AIC)
  resultAIC <- list()
  resultAIC$COEF <- saic$coefficients[-1,1]
  resultAIC$PVAL <- saic$coefficients[-1,4]
  resultAIC$INTERC <- saic$coefficients[-1,1][1]
  resultAIC$RESID <- ms_AIC$residuals
  resultAIC$FITTED <- ms_AIC$fitted.values
  resultAIC$STD.DEV <- saic$sigma
  resultAIC$R.SQUARED <- saic$r.squared
  return(resultAIC)                
})

AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_diag_wo_co_nw_aic, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC)))
    names(coef_result) <- paste( "`",colnames(Xfit_NC),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL <- sapply(model_results_diag_wo_co_nw_aic, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC)))
    names(coef_result) <- paste( "`",colnames(Xfit_NC),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[!AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL < 0.05] <- 0

AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

# AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- {qgraph(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[1:70,],
#                                posCol="blue",
#                                negCol="purple")
#   title("AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES")}

Residuals_AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_diag_wo_co_nw_aic, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_betw(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_betw(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_degree(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_degree(AIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]

  


## diagnosed cross validation lasso without weights. 
model_results_diag_lasso_ww_nc <-lapply(1:nROI, function(i){
  cat(i)
#   colnames(Yfit_NC)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC)[i]
  y <- Yfit_NC[,i]
  #y <- unname(y) 
  x <- Xfit_NC
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  fit = cv.glmnet(x, y, alpha = 1, nlambda = 20,type.measure = "mse", nfolds = 20)
  fitted_values <- predict(fit, x, s = 'lambda.min')
  resultAR <- list()
  resultAR$COEF <- coef(fit, s = "lambda.min")
  resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
  resultAR$fit <- fit
  resultAR$residuals <-  y - fitted_values
  return(resultAR)                
})

MATRIX_DIAG_LASSO_WW_NC <- sapply(model_results_diag_lasso_ww_nc, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC))+1)
    names(coef_result) <- model_results_diag_lasso_ww_nc[[1]][["COEF"]]@Dimnames[[1]]
    res <- y$COEF
    indexes <- match(dimnames(model_results_diag_lasso_ww_nc[[1]][["COEF"]])[[1]],names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

Residuals_DIAG_LASSO_WW_NC <- sapply(model_results_diag_lasso_ww_nc, function(y){
    res <- y$residuals
    })

NETWORK_DIAG_LASSO_WW_NC_igraph <- graph_from_adjacency_matrix(MATRIX_DIAG_LASSO_WW_NC[2:71,],mode = "directed",weighted = TRUE)

save(Residuals_DIAG_LASSO_WW_NC,file = "Residuals_DIAG_LASSO_WITHOUT_WEIGHTS_NC.Rdata")
# NETWORK_DIAG_LASSO_WW_NC <- {qgraph(MATRIX_DIAG_LASSO_WW_NC[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_DIAG_LASSO_WW")}

## graph density without weights and no covariates_lasso_diagnosed
graph.density(NETWORK_DIAG_LASSO_WW_NC_igraph)

## betweenness centrality max
dimnames(MATRIX_DIAG_LASSO_WW_NC)[[1]][which.max(centr_betw(NETWORK_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_DIAG_LASSO_WW_NC)[[1]][which.min(centr_betw(NETWORK_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_DIAG_LASSO_WW_NC)[[1]][which.max(centr_degree(NETWORK_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_DIAG_LASSO_WW_NC)[[1]][which.min(centr_degree(NETWORK_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]

## non diagnosed cross validation lasso with weights
model_results_non_diag_lasso_nc <-lapply(1:nROI, function(i){
    cat(i)
    y <- Yfit_NC_ND[,i]
    #y <- unname(y) 
    x <- Xfit_NC_ND
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    fit = cv.glmnet(x, y, alpha = 1, nlambda = 20, weights = rep(7/17,dim(x)[1]), type.measure = "mse", nfolds = 20)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

save(model_results_non_diag_lasso_nc,file="model_results_non_diag_lasso_nc.Rdata")

MATRIX_NON_DIAG_LASSO_NO_COVARIATES <- sapply(model_results_non_diag_lasso_nc, function(y){
  coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
  names(coef_result) <- model_results_non_diag_lasso_nc[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_non_diag_lasso_nc[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

NETWORK_NON_DIAG_LASSO_NC_igraph <- graph_from_adjacency_matrix(MATRIX_NON_DIAG_LASSO_NO_COVARIATES[2:71,],mode = "directed",weighted = TRUE)


Residuals_DIAG_LASSO_NC <- sapply(model_results_non_diag_lasso_nc, function(y){
    res <- y$residuals
    })

# NETWORK_NON_DIAG_LASSO_NC <- {qgraph(MATRIX_NON_DIAG_LASSO_NO_COVARIATES[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_NON_DIAG_LASSO_NC")}

## graph density with weights and no covariates_lasso_non_diagnosed
graph.density(NETWORK_NON_DIAG_LASSO_NC_igraph)


## betweenness centrality max
dimnames(MATRIX_NON_DIAG_LASSO_NO_COVARIATES)[[1]][which.max(centr_betw(NETWORK_NON_DIAG_LASSO_NC_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_NON_DIAG_LASSO_NO_COVARIATES)[[1]][which.min(centr_betw(NETWORK_NON_DIAG_LASSO_NC_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_NON_DIAG_LASSO_NO_COVARIATES)[[1]][which.max(centr_degree(NETWORK_NON_DIAG_LASSO_NC_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_NON_DIAG_LASSO_NO_COVARIATES)[[1]][which.min(centr_degree(NETWORK_NON_DIAG_LASSO_NC_igraph)[["res"]])+1]

## non diagnosed cross validation lasso without weights
model_results_non_diag_lasso_ww_nc<-lapply(1:nROI, function(i){
    cat(i)
    y <- Yfit_NC_ND[,i]
    #y <- unname(y) 
    x <- Xfit_NC_ND
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    fit = cv.glmnet(x, y, alpha = 1, nlambda = 20, type.measure = "mse", nfolds = 20)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES <- sapply(model_results_non_diag_lasso_ww_nc, function(y){
  coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
  names(coef_result) <- model_results_non_diag_lasso_ww_nc[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_non_diag_lasso_ww_nc[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

NETWORK_NON_DIAG_LASSO_WW_NC_igraph <- graph_from_adjacency_matrix(MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

Residuals_NON_DIAG_LASSO_NC_WW <- sapply(model_results_non_diag_lasso_ww_nc, function(y){
    res <- y$residuals
    })
# NETWORK_NON_DIAG_LASSO_WW_NC <- {qgraph(MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_NON_DIAG_LASSO_WW_NC")}

## graph density without weights and no covariates_lasso_non_diagnosed
graph.density(NETWORK_NON_DIAG_LASSO_WW_NC_igraph)


## betweenness centrality max
dimnames(MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES)[[1]][which.max(centr_betw(NETWORK_NON_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES)[[1]][which.min(centr_betw(NETWORK_NON_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES)[[1]][which.max(centr_degree(NETWORK_NON_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_NON_DIAG_LASSO_WW_NO_COVARIATES)[[1]][which.min(centr_degree(NETWORK_NON_DIAG_LASSO_WW_NC_igraph)[["res"]])+1]

## NON_DIAGNOSED_GROUP WITHOUT WEIGHTS BIC

model_results_non_diag_wo_co_nw <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_NC_ND)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC_ND)[i]
  y <- Yfit_NC_ND[,i]
  #y <- unname(y) 
  x <- Xfit_NC_ND
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  data_fit <- as.data.frame(cbind(Yfit_NC_ND[complete_obs,i],1,x))
  
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  model_fl_lower <- as.formula(paste(yname,"~1"))
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



BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co_nw, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL <- sapply(model_results_diag_wo_co_nw, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[!BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL < 0.05] <- 0

BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

# BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- {qgraph(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph[1:70,],
#                                posCol="blue",
#                                negCol="purple")
#   title("BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES")}
Residuals_BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co_nw, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_betw(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_betw(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_degree(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_degree(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]

## graph density without weights and no covariates_BIC_diagnosed
graph.density(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)

## NON_DIAGNOSED_GROUP WITHOUT WEIGHTS AIC

model_results_non_diag_wo_co_nw_aic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_NC_ND)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC_ND)[i]
  y <- Yfit_NC_ND[,i]
  #y <- unname(y) 
  x <- Xfit_NC_ND
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  data_fit <- as.data.frame(cbind(Yfit_NC_ND[complete_obs,i],1,x))
  
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  model_fl_lower <- as.formula(paste(yname,"~1"))
  model_lower <- lm(model_fl_lower, data_fit)
  model_fl_upper <- as.formula(paste(yname,"~."))
  model_upper <- lm(model_fl_upper, data_fit)
  
  ms_AIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower))

  saic<- summary(ms_AIC)
  resultAIC <- list()
  resultAIC$COEF <- saic$coefficients[-1,1]
  resultAIC$PVAL <- saic$coefficients[-1,4]
  resultAIC$INTERC <- saic$coefficients[-1,1][1]
  resultAIC$RESID <- ms_AIC$residuals
  resultAIC$FITTED <- ms_AIC$fitted.values
  resultAIC$STD.DEV <- saic$sigma
  resultAIC$R.SQUARED <- saic$r.squared
  return(resultAIC)                
})



AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co_nw_aic, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL <- sapply(model_results_diag_wo_co_nw_aic, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[!AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_PVAL < 0.05] <- 0

AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

Residuals_AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co_nw_aic, function(y){
    res <- y$RESID
    })
# AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES <- {qgraph(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph[1:70,],
#                                posCol="blue",
#                                negCol="purple")
#   title("AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES")}

## betweenness centrality max
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_betw(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_betw(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.max(centr_degree(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES)[[1]][which.min(centr_degree(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)[["res"]])]

## NON_DIAGNOSED_GROUP WEIGHTS BIC

model_results_non_diag_wo_co <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_NC_ND)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC_ND)[i]
  y <- Yfit_NC_ND[,i]
  #y <- unname(y) 
  x <- Xfit_NC_ND
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  
  data_fit <- as.data.frame(cbind(Yfit_NC_ND[complete_obs,i],1,x))
  weights_1 = rep(7/17,nrow(data_fit))
  names(data_fit)[1] <- paste("Y",i,sep = "")
  model_fl_lower <- as.formula(paste(yname,"~1"))
  model_lower <- lm(model_fl_lower, data_fit,weights =weights_1)
  model_fl_upper <- as.formula(paste(yname,"~."))
  model_upper <- lm(model_fl_upper, data_fit,weights =weights_1)
  
  ms_BIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower), 
                    k= log(70))

  sbic<- summary(ms_BIC)
  resultBIC <- list()
  resultBIC$COEF <- sbic$coefficients[-1,1]
  resultBIC$PVAL <- sbic$coefficients[-1,4]
  resultBIC$INTERC <- sbic$coefficients[-1,1][1]
  resultBIC$RESID <- ms_BIC$residuals
  resultBIC$FITTED <- ms_BIC$fitted.values
  resultBIC$STD.DEV <- sbic$sigma
  resultBIC$R.SQUARED <- sbic$r.squared
  return(resultBIC)                
})



BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_PVAL <- sapply(model_results_non_diag_wo_co, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES[!BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_PVAL < 0.05] <- 0

BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

# BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES <- {qgraph(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph[1:70,],
#                                posCol="blue",
#                                negCol="purple")
#   title("BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES")}

## graph density without weights and no covariates_BIC_diagnosed
graph.density(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)

Residuals_BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.max(centr_betw(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.min(centr_betw(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.max(centr_degree(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.min(centr_degree(BIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]

## NON_DIAGNOSED_GROUP WEIGHTS AIC

model_results_non_diag_wo_co_aic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_NC_ND)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_NC_ND)[i]
  y <- Yfit_NC_ND[,i]
  #y <- unname(y) 
  x <- Xfit_NC_ND
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  data_fit <- as.data.frame(cbind(Yfit_NC_ND[complete_obs,i],1,x))
  weights_1 = rep(7/17,nrow(data_fit))
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  model_fl_lower <- as.formula(paste(yname,"~1"))
  model_lower <- lm(model_fl_lower, data_fit,weights =weights_1)
  model_fl_upper <- as.formula(paste(yname,"~."))
  model_upper <- lm(model_fl_upper, data_fit,weights =weights_1)
  
  ms_AIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower))

  saic<- summary(ms_AIC)
  resultAIC <- list()
  resultAIC$COEF <- saic$coefficients[-1,1]
  resultAIC$PVAL <- saic$coefficients[-1,4]
  resultAIC$INTERC <- saic$coefficients[-1,1][1]
  resultAIC$RESID <- ms_AIC$residuals
  resultAIC$FITTED <- ms_AIC$fitted.values
  resultAIC$STD.DEV <- saic$sigma
  resultAIC$R.SQUARED <- saic$r.squared
  return(resultAIC)                
})



AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co_aic, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_PVAL <- sapply(model_results_non_diag_wo_co_aic, function(y){
    coef_result <- rep(0,length(colnames(Xfit_NC_ND))+1)
    names(coef_result) <- paste( "`",colnames(Xfit_NC_ND),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES[!AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_PVAL < 0.05] <- 0

AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES[1:70,],mode = "directed",weighted = TRUE)

# AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES <- {qgraph(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph[1:70,],
#                                posCol="blue",
#                                negCol="purple")
#   title("AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES")}

Residuals_AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES <- sapply(model_results_non_diag_wo_co_aic, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.max(centr_betw(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.min(centr_betw(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.max(centr_degree(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES)[[1]][which.min(centr_degree(AIC_NETWORK_NON_DIAG_WITHOUT_COVARIATES_igraph)[["res"]])]

## total population. 
## total population without weights lasso

model_results_without_weights_no_covariates <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit[,1:70]
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    fit = cv.glmnet(x, y, alpha = 1, nlambda = 20,type.measure = "mse", nfolds = 20)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})



MATRIX_LASSO_WW_NC <- sapply(model_results_without_weights_no_covariates, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+4)
    names(coef_result) <- model_results_without_weights_no_covariates[[1]][["COEF"]]@Dimnames[[1]]
    res <- y$COEF
    indexes <- match(dimnames(model_results_without_weights_no_covariates[[1]][["COEF"]])[[1]],names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

NETWORK_MATRIX_LASSO_WW_NC_igraph <- graph_from_adjacency_matrix(MATRIX_LASSO_WW_NC[1:70,],mode = "directed",weighted = TRUE)


Residuals_TOTAL_POPU_LASSO_NC_WW <- sapply(model_results_without_weights_no_covariates, function(y){
    res <- y$residuals
    })
## graph density without weights and no covariates_lasso_total_population
graph.density(NETWORK_MATRIX_LASSO_WW_NC_igraph)


## betweenness centrality max
dimnames(MATRIX_LASSO_WW_NC)[[1]][which.max(centr_betw(NETWORK_MATRIX_LASSO_WW_NC_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_LASSO_WW_NC)[[1]][which.min(centr_betw(NETWORK_MATRIX_LASSO_WW_NC_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_LASSO_WW_NC)[[1]][which.max(centr_degree(NETWORK_MATRIX_LASSO_WW_NC_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_LASSO_WW_NC)[[1]][which.min(centr_degree(NETWORK_MATRIX_LASSO_WW_NC_igraph)[["res"]])+1]


## total population weights lasso
model_results_weights_nc <-lapply(1:nROI, function(i){
    cat(i)

    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    weight <- x[,74]
    x <- x[,-(71:74)]
    fit = cv.glmnet(x, y, alpha = 1, weights = weight, nlambda = 20,type.measure = "mse", nfolds = 20)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

MATRIX_LASSO_WEIGHT_NC <- sapply(model_results_weights_nc, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- model_results_weights_nc[[1]][["COEF"]]@Dimnames[[1]]
    res <- y$COEF
    indexes <- match(dimnames(model_results_weights_nc[[1]][["COEF"]])[[1]],names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

NETWORK_LASSO_WEIGHT_NC_igraph <- graph_from_adjacency_matrix(MATRIX_LASSO_WEIGHT_NC[1:70,],mode = "directed",weighted = TRUE)


Residuals_TOTAL_POPU_LASSO_NC <- sapply(model_results_weights_nc, function(y){
    res <- y$residuals
    })

## graph density weights and no covariates_lasso_total_population
graph.density(NETWORK_LASSO_WEIGHT_NC_igraph)



## betweenness centrality max
dimnames(MATRIX_LASSO_WEIGHT_NC)[[1]][which.max(centr_betw(NETWORK_LASSO_WEIGHT_NC_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_LASSO_WEIGHT_NC)[[1]][which.min(centr_betw(NETWORK_LASSO_WEIGHT_NC_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_LASSO_WEIGHT_NC)[[1]][which.max(centr_degree(NETWORK_LASSO_WEIGHT_NC_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_LASSO_WEIGHT_NC)[[1]][which.min(centr_degree(NETWORK_LASSO_WEIGHT_NC_igraph)[["res"]])+1]


##total population without weights BIC
model_results_BIC_NC <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    x <- x[,-(71:74)]
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    ##weights_1 = rep(7/17,nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    model_fl_lower <- as.formula(paste(yname,"~1"))
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

BIC_NETWORK_POPU_NC <- sapply(model_results_BIC_NC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_POPU_NC_PVAL <- sapply(model_results_BIC_NC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


BIC_NETWORK_POPU_NC[!BIC_NETWORK_POPU_NC_PVAL < 0.05] <- 0

BIC_NETWORK_POPU_NC_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_POPU_NC[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_POPU_NC <- {qgraph(BIC_NETWORK_POPU_NC[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_POPU_NC")}

Residuals_BIC_NETWORK_POPU_NC <- sapply(model_results_BIC_NC, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(BIC_NETWORK_POPU_NC)[[1]][which.max(centr_betw(BIC_NETWORK_POPU_NC_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_POPU_NC)[[1]][which.min(centr_betw(BIC_NETWORK_POPU_NC_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_POPU_NC)[[1]][which.max(centr_degree(BIC_NETWORK_POPU_NC_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_POPU_NC)[[1]][which.min(centr_degree(BIC_NETWORK_POPU_NC_igraph)[["res"]])]

## graph density without weights and no covariates_bic_total_population
graph.density(BIC_NETWORK_POPU_NC_igraph)

##total population weights BIC
model_results_BIC_NC_WEIGHTS <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    weight <- x[,74]
    x <- x[,-(71:74)]
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    model_fl_lower <- as.formula(paste(yname,"~1"))
    model_lower <- lm(model_fl_lower, data_fit, weights = weight)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit, weights = weight)
    
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

BIC_NETWORK_POPU_WEIGHTS_NC <- sapply(model_results_BIC_NC_WEIGHTS, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_POPU_WEIGHTS_NC_PVAL <- sapply(model_results_BIC_NC_WEIGHTS, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


BIC_NETWORK_POPU_WEIGHTS_NC[!BIC_NETWORK_POPU_WEIGHTS_NC_PVAL < 0.05] <- 0

BIC_NETWORK_POPU_WEIGHTS_NC_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_POPU_WEIGHTS_NC[1:70,],mode = "directed",weighted = TRUE)

Residuals_BIC_NETWORK_POPU_WEIGHTS_NC <- sapply(model_results_BIC_NC_WEIGHTS, function(y){
    res <- y$RESID
    })
    
# NETWORK_POPU_NC_WITH_WEIGHTS <- {qgraph(BIC_NETWORK_POPU_WEIGHTS_NC[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_POPU_WEIGHTS")}

## graph density weights and no covariates_bic_total_population
graph.density(BIC_NETWORK_POPU_WEIGHTS_NC_igraph)



## betweenness centrality max
dimnames(BIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.max(centr_betw(BIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.min(centr_betw(BIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.max(centr_degree(BIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.min(centr_degree(BIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]

##total population without weights AIC
model_results_AIC_NC <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    x <- x[,-(71:74)]
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    ##weights_1 = rep(7/17,nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    model_fl_lower <- as.formula(paste(yname,"~1"))
    model_lower <- lm(model_fl_lower, data_fit)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit)
    
    ms_AIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower))
    saic<- summary(ms_AIC)
    resultAIC <- list()
    resultAIC$COEF <- saic$coefficients[-1,1]
    resultAIC$PVAL <- saic$coefficients[-1,4]
    resultAIC$INTERC <- saic$coefficients[-1,1][1]
    resultAIC$RESID <- ms_AIC$residuals
    resultAIC$FITTED <- ms_AIC$fitted.values
    resultAIC$STD.DEV <- saic$sigma
    resultAIC$R.SQUARED <- saic$r.squared
    return(resultAIC)                
})

AIC_NETWORK_POPU_NC <- sapply(model_results_AIC_NC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_POPU_NC_PVAL <- sapply(model_results_AIC_NC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


AIC_NETWORK_POPU_NC[!AIC_NETWORK_POPU_NC_PVAL < 0.05] <- 0

AIC_NETWORK_POPU_NC_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_POPU_NC[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_POPU_NC_AIC <- {qgraph(AIC_NETWORK_POPU_NC[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_POPU_NC")}

## graph density without weights and no covariates_bic_total_population
graph.density(BIC_NETWORK_POPU_NC_igraph)

Residuals_AIC_NETWORK_POPU_NC <- sapply(model_results_AIC_NC, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(AIC_NETWORK_POPU_NC)[[1]][which.max(centr_betw(AIC_NETWORK_POPU_NC_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_POPU_NC)[[1]][which.min(centr_betw(AIC_NETWORK_POPU_NC_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_POPU_NC)[[1]][which.max(centr_degree(AIC_NETWORK_POPU_NC_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_POPU_NC)[[1]][which.min(centr_degree(AIC_NETWORK_POPU_NC_igraph)[["res"]])]

##total population weights AIC
model_results_AIC_NC_WEIGHTS <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    weight <- x[,74]
    x <- x[,-(71:74)]
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    model_fl_lower <- as.formula(paste(yname,"~1"))
    model_lower <- lm(model_fl_lower, data_fit, weights = weight)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit, weights = weight)
    
    ms_AIC <- stepAIC(model_lower,direction="both",scope=list(upper=model_upper,lower=model_lower), 
                      k= log(70))
    saic<- summary(ms_AIC)
    resultAIC <- list()
    resultAIC$COEF <- saic$coefficients[-1,1]
    resultAIC$PVAL <- saic$coefficients[-1,4]
    resultAIC$INTERC <- saic$coefficients[-1,1][1]
    resultAIC$RESID <- ms_AIC$residuals
    resultAIC$FITTED <- ms_AIC$fitted.values
    resultAIC$STD.DEV <- saic$sigma
    resultAIC$R.SQUARED <- saic$r.squared
    return(resultAIC)                
})

AIC_NETWORK_POPU_WEIGHTS_NC <- sapply(model_results_AIC_NC_WEIGHTS, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_POPU_WEIGHTS_NC_PVAL <- sapply(model_results_AIC_NC_WEIGHTS, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


AIC_NETWORK_POPU_WEIGHTS_NC[!AIC_NETWORK_POPU_WEIGHTS_NC_PVAL < 0.05] <- 0

AIC_NETWORK_POPU_WEIGHTS_NC_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_POPU_WEIGHTS_NC[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_POPU_NC_WITH_WEIGHTS_AIC <- {qgraph(AIC_NETWORK_POPU_WEIGHTS_NC[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_POPU_NC_WITH_WEIGHTS_AIC")}

Residuals_AIC_NETWORK_POPU_NC_WEIGHTS <- sapply(model_results_AIC_NC_WEIGHTS, function(y){
    res <- y$RESID
    })

## betweenness centrality max
dimnames(AIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.max(centr_betw(AIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.min(centr_betw(AIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.max(centr_degree(AIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_POPU_WEIGHTS_NC)[[1]][which.min(centr_degree(AIC_NETWORK_POPU_WEIGHTS_NC_igraph)[["res"]])]



graph.density(AIC_NETWORK_POPU_WEIGHTS_NC_igraph)
## graph density without weights and no covariates_BIC_diagnosed
graph.density(BIC_NETWORK_DIAG_NO_WEIGHTS_WITHOUT_COVARIATES_igraph)

## graph density without weights and no covariates_lasso_diagnosed
graph.density(NETWORK_DIAG_LASSO_WW_NC_igraph)

## graph density with weights and no covariates_lasso_non_diagnosed
graph.density(NETWORK_NON_DIAG_LASSO_NC_igraph)

## graph density without weights and no covariates_lasso_non_diagnosed
graph.density(NETWORK_NON_DIAG_LASSO_WW_NC_igraph)

## graph density without weights and no covariates_bic_total_population
graph.density(BIC_NETWORK_POPU_NC_igraph)

## graph density weights and no covariates_bic_total_population
graph.density(BIC_NETWORK_POPU_WEIGHTS_NC_igraph)

