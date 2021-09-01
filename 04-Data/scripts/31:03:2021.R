## WITH COVARIATES

## DIAGNOSED_GROUP WITHOUT WEIGHTS BIC
model_results_diag_nw_bic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_diagnosed)[i]
  y <- Yfit_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  x <- x[,-(72:74)]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_diagnosed[complete_obs,i],1,x))
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
  resultBIC$INTERC <- sbic$coefficients[-1,1][1]
  resultBIC$RESID <- ms_BIC$residuals
  resultBIC$FITTED <- ms_BIC$fitted.values
  resultBIC$STD.DEV <- sbic$sigma
  resultBIC$R.SQUARED <- sbic$r.squared
  return(resultBIC)                
})


BIC_NETWORK_DIAG_NO_WEIGHTS <- sapply(model_results_diag_nw_bic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_DIAG_NO_WEIGHTS_PVAL <- sapply(model_results_diag_nw_bic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_DIAG_NO_WEIGHTS[!BIC_NETWORK_DIAG_NO_WEIGHTS_PVAL < 0.05] <- 0

BIC_NETWORK_DIAG_NO_WEIGHTS_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_DIAG_NO_WEIGHTS[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_DIAG_NO_WEIGHTS_BIC <- {qgraph(BIC_NETWORK_DIAG_NO_WEIGHTS[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_DIAG_NO_WEIGHTS")}

graph.density(BIC_NETWORK_DIAG_NO_WEIGHTS_igraph)

## betweenness centrality max
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.max(centr_betw(BIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.min(centr_betw(BIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.max(centr_degree(BIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.min(centr_degree(BIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]



## DIAGNOSED_GROUP WITHOUT WEIGHTS AIC
model_results_diag_nw_aic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_diagnosed)[i]
  y <- Yfit_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  x <- x[,-(72:74)]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_diagnosed[complete_obs,i],1,x))
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  data_fit$handedness <- factor(data_fit$handedness)
  model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
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


AIC_NETWORK_DIAG_NO_WEIGHTS <- sapply(model_results_diag_nw_aic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

AIC_NETWORK_DIAG_NO_WEIGHTS_PVAL <- sapply(model_results_diag_nw_aic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

AIC_NETWORK_DIAG_NO_WEIGHTS[!AIC_NETWORK_DIAG_NO_WEIGHTS_PVAL < 0.05] <- 0

AIC_NETWORK_DIAG_NO_WEIGHTS_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_DIAG_NO_WEIGHTS[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_DIAG_NO_WEIGHTS_AIC <- {qgraph(AIC_NETWORK_DIAG_NO_WEIGHTS[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("AIC_NETWORK_DIAG_NO_WEIGHTS")}

graph.density(AIC_NETWORK_DIAG_NO_WEIGHTS_igraph)

## betweenness centrality max
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.max(centr_betw(AIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.min(centr_betw(AIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.max(centr_degree(AIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_DIAG_NO_WEIGHTS)[[1]][which.min(centr_degree(AIC_NETWORK_DIAG_NO_WEIGHTS_igraph)[["res"]])]



## diagnosed cross validation lasso without weights. 
model_results_diag_lasso_ww <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_diagnosed)[i]
  y <- Yfit_diagnosed[,i]
  #y <- unname(y) 
  x <- as.matrix(dummy_cols(Xfit_diagnosed[,1:71],select_columns = dimnames(Xfit_diagnosed)[[2]][71],ignore_na = T))
  p.fac = rep(1, 20)
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  p.fac = rep(1, dim(x)[2])
  p.fac[c(71:74)] =0
  fit = cv.glmnet(x, y, alpha = 1, nlambda = 20,type.measure = "mse", nfolds = 20, penalty.factor = p.fac)
  fitted_values <- predict(fit, x, s = 'lambda.min')
  resultAR <- list()
  resultAR$COEF <- coef(fit, s = "lambda.min")
  resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
  resultAR$fit <- fit
  resultAR$residuals <-  y - fitted_values
  return(resultAR)                
})

MATRIX_DIAG_LASSO_WW <- sapply(model_results_diag_lasso_ww, function(y){
  coef_result <- rep(0,length(colnames(Xfit_diagnosed))+3)
  names(coef_result) <- model_results_diag_lasso_ww[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_diag_lasso_ww[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

NETWORK_DIAG_LASSO_WW_igraph <- graph_from_adjacency_matrix(MATRIX_DIAG_LASSO_WW[2:71,],mode = "directed",weighted = TRUE)

Residuals_DIAG_LASSO_WW <- sapply(model_results_diag_lasso_ww, function(y){
    res <- y$residuals
    })

# NETWORK_DIAG_LASSO_WW <- {qgraph(MATRIX_DIAG_LASSO_WW[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_DIAG_LASSO_WW")}


## betweenness centrality max
dimnames(MATRIX_DIAG_LASSO_WW)[[1]][which.max(centr_betw(NETWORK_DIAG_LASSO_WW_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_DIAG_LASSO_WW)[[1]][which.min(centr_betw(NETWORK_DIAG_LASSO_WW_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_DIAG_LASSO_WW)[[1]][which.max(centr_degree(NETWORK_DIAG_LASSO_WW_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_DIAG_LASSO_WW)[[1]][which.min(centr_degree(NETWORK_DIAG_LASSO_WW_igraph)[["res"]])+1]

graph.density(NETWORK_DIAG_LASSO_WW_igraph)

graph.density(NETWORK_DIAG_LASSO_WITH_WEIGHTS)

## NON DIAGNOSED_GROUP WITHOUT WEIGHTS BIC
model_results_non_diag_nw_bic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_non_diagnosed)[i]
  y <- Yfit_non_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_non_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  x <- x[,-(72:74)]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_non_diagnosed[complete_obs,i],1,x))
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
  resultBIC$INTERC <- sbic$coefficients[-1,1][1]
  resultBIC$RESID <- ms_BIC$residuals
  resultBIC$FITTED <- ms_BIC$fitted.values
  resultBIC$STD.DEV <- sbic$sigma
  resultBIC$R.SQUARED <- sbic$r.squared
  return(resultBIC)                
})


BIC_NETWORK_NON_DIAG_NO_WEIGHTS <- sapply(model_results_non_diag_nw_bic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_NON_DIAG_NO_WEIGHTS_PVAL <- sapply(model_results_non_diag_nw_bic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_NON_DIAG_NO_WEIGHTS[!BIC_NETWORK_NON_DIAG_NO_WEIGHTS_PVAL < 0.05] <- 0

BIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_NON_DIAG_NO_WEIGHTS[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_DIAG_NO_WEIGHTS_BIC <- {qgraph(BIC_NETWORK_NON_DIAG_NO_WEIGHTS[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_NON_DIAG_NO_WEIGHTS")}


## betweenness centrality max
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.max(centr_betw(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.min(centr_betw(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.max(centr_degree(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.min(centr_degree(BIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]

## NON DIAGNOSED_GROUP WITHOUT WEIGHTS AIC
model_results_non_diag_nw_aic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_non_diagnosed)[i]
  y <- Yfit_non_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_non_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  x <- x[,-(72:74)]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_non_diagnosed[complete_obs,i],1,x))
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  data_fit$handedness <- factor(data_fit$handedness)
  model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
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


AIC_NETWORK_NON_DIAG_NO_WEIGHTS <- sapply(model_results_non_diag_nw_aic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

AIC_NETWORK_NON_DIAG_NO_WEIGHTS_PVAL <- sapply(model_results_non_diag_nw_aic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

AIC_NETWORK_NON_DIAG_NO_WEIGHTS[!AIC_NETWORK_NON_DIAG_NO_WEIGHTS_PVAL < 0.05] <- 0

AIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_NON_DIAG_NO_WEIGHTS[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_DIAG_NO_WEIGHTS_AIC <- {qgraph(AIC_NETWORK_NON_DIAG_NO_WEIGHTS[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("AIC_NETWORK_NON_DIAG_NO_WEIGHTS")}



## betweenness centrality max
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.max(centr_betw(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.min(centr_betw(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.max(centr_degree(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_NON_DIAG_NO_WEIGHTS)[[1]][which.min(centr_degree(AIC_NETWORK_NON_DIAG_NO_WEIGHTS_igraph)[["res"]])]


## NON DIAGNOSED_GROUP WEIGHTS AIC
model_results_non_diag_aic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_non_diagnosed)[i]
  y <- Yfit_non_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_non_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  weight <- x[,74]
  x <- x[,-(72:74)]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_non_diagnosed[complete_obs,i],1,x))
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  data_fit$handedness <- factor(data_fit$handedness)
  model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
  model_lower <- lm(model_fl_lower, data_fit,weights =weight)
  model_fl_upper <- as.formula(paste(yname,"~."))
  model_upper <- lm(model_fl_upper, data_fit,weights =weight)
  
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


AIC_NETWORK_NON_DIAG <- sapply(model_results_non_diag_aic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

AIC_NETWORK_NON_DIAG_PVAL <- sapply(model_results_non_diag_aic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

AIC_NETWORK_NON_DIAG[!AIC_NETWORK_NON_DIAG_PVAL < 0.05] <- 0

AIC_NETWORK_NON_DIAG_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_NON_DIAG[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_DIAG_AIC <- {qgraph(AIC_NETWORK_NON_DIAG[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("AIC_NETWORK_NON_DIAG")}


## betweenness centrality max
dimnames(AIC_NETWORK_NON_DIAG)[[1]][which.max(centr_betw(AIC_NETWORK_NON_DIAG_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_NON_DIAG)[[1]][which.min(centr_betw(AIC_NETWORK_NON_DIAG_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_NON_DIAG)[[1]][which.max(centr_degree(AIC_NETWORK_NON_DIAG_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_NON_DIAG)[[1]][which.min(centr_degree(AIC_NETWORK_NON_DIAG_igraph)[["res"]])]


## NON DIAGNOSED_GROUP WITH WEIGHTS BIC
model_results_non_diag_bic <-lapply(1:nROI, function(i){
  cat(i)
  colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
  
  yname <- colnames(Yfit_non_diagnosed)[i]
  y <- Yfit_non_diagnosed[,i]
  #y <- unname(y) 
  x <- Xfit_non_diagnosed
  complete_obs <- complete.cases(y) & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  weight <- <- x[,74]
  x <- x[,-(72:74)]
  colnames(x)[71] <- "handedness"
  data_fit <- as.data.frame(cbind(Yfit_non_diagnosed[complete_obs,i],1,x))
  names(data_fit)[1] <- paste("Y",i,sep = "")
  
  data_fit$handedness <- factor(data_fit$handedness)
  model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
  model_lower <- lm(model_fl_lower, data_fit,weights =weight)
  model_fl_upper <- as.formula(paste(yname,"~."))
  model_upper <- lm(model_fl_upper, data_fit,weights =weight)
  
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


BIC_NETWORK_NON_DIAG <- sapply(model_results_non_diag_bic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$COEF
  indexes <- match(names(y[["COEF"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_NON_DIAG_PVAL <- sapply(model_model_results_non_diag_bicresults_non_diag_nw_bic, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+1)
  names(coef_result) <- paste( "`",colnames(Xfit_non_diagnosed),"`",sep = "")
  names(coef_result)[74] <- "handedness2"
  names(coef_result)[75] <- "handedness3"
  res <- y$PVAL
  indexes <- match(names(y[["PVAL"]]),names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

BIC_NETWORK_NON_DIAG[!BIC_NETWORK_NON_DIAG_PVAL < 0.05] <- 0

BIC_NETWORK_NON_DIAG_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_NON_DIAG[1:70,],mode = "directed",weighted = TRUE)

# NETWORK_DIAG_BIC <- {qgraph(BIC_NETWORK_NON_DIAG[1:70,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_NON_DIAG")}

## betweenness centrality max
dimnames(BIC_NETWORK_NON_DIAG)[[1]][which.max(centr_betw(BIC_NETWORK_NON_DIAG_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_NON_DIAG)[[1]][which.min(centr_betw(BIC_NETWORK_NON_DIAG_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_NON_DIAG)[[1]][which.max(centr_degree(BIC_NETWORK_NON_DIAG_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_NON_DIAG)[[1]][which.min(centr_degree(BIC_NETWORK_NON_DIAG_igraph)[["res"]])]


## non diagnosed cross validation lasso with weights
model_results_non_diag_lasso <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit_non_diagnosed)[i]
    y <- Yfit_non_diagnosed[,i]
    #y <- unname(y) 
    x <- as.matrix(dummy_cols(Xfit_non_diagnosed[,1:71],select_columns = dimnames(Xfit_non_diagnosed)[[2]][71],ignore_na = T))
    p.fac = rep(1, 20)
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    p.fac = rep(1, dim(x)[2])
    p.fac[c(71:74)] =0
    fit = cv.glmnet(x, y, alpha = 1, nlambda = 20, weights = rep(7/17,dim(x)[1]), type.measure = "mse", nfolds = 20, penalty.factor = p.fac)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

MATRIX_NON_DIAG_LASSO <- sapply(model_results_non_diag_lasso, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed))+4)
  names(coef_result) <- model_results_non_diag_lasso[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_non_diag_lasso[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

NETWORK_NON_DIAG_LASSO_igraph <- graph_from_adjacency_matrix(MATRIX_NON_DIAG_LASSO[2:71,],mode = "directed",weighted = TRUE)

# NETWORK_NON_DIAG_LASSO <- {qgraph(MATRIX_NON_DIAG_LASSO[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_NON_DIAG_LASSO")}

Residuals_NON_DIAG_LASSO <- sapply(model_results_non_diag_lasso, function(y){
    res <- y$residuals
    })



## betweenness centrality max
dimnames(MATRIX_NON_DIAG_LASSO)[[1]][which.max(centr_betw(NETWORK_NON_DIAG_LASSO_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_NON_DIAG_LASSO)[[1]][which.min(centr_betw(NETWORK_NON_DIAG_LASSO_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_NON_DIAG_LASSO)[[1]][which.max(centr_degree(NETWORK_NON_DIAG_LASSO_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_NON_DIAG_LASSO)[[1]][which.min(centr_degree(NETWORK_NON_DIAG_LASSO_igraph)[["res"]])+1]


## non diagnosed cross validation lasso with out weights
model_results_non_diag_lasso_ww <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit_non_diagnosed)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit_non_diagnosed)[i]
    y <- Yfit_non_diagnosed[,i]
    #y <- unname(y) 
    x <- as.matrix(dummy_cols(Xfit_non_diagnosed[,1:71],select_columns = dimnames(Xfit_non_diagnosed)[[2]][71],ignore_na = T))
    p.fac = rep(1, 20)
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    p.fac = rep(1, dim(x)[2])
    p.fac[c(71:74)] =0
    fit = cv.glmnet(x, y, alpha = 1, nlambda = 20, type.measure = "mse", nfolds = 20, penalty.factor = p.fac)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

MATRIX_NON_DIAG_LASSO_WW <- sapply(model_results_non_diag_lasso_ww, function(y){
  coef_result <- rep(0,length(colnames(Xfit_non_diagnosed)))
  names(coef_result) <- model_results_non_diag_lasso[[1]][["COEF"]]@Dimnames[[1]]
  res <- y$COEF
  indexes <- match(dimnames(model_results_non_diag_lasso[[1]][["COEF"]])[[1]],names(coef_result))
  coef_result[indexes] <- res
  coef_result
})

NETWORK_NON_DIAG_LASSO_WW_igraph <- graph_from_adjacency_matrix(MATRIX_NON_DIAG_LASSO_WW[2:71,],mode = "directed",weighted = TRUE)

# NETWORK_NON_DIAG_LASSO_WW <- {qgraph(MATRIX_NON_DIAG_LASSO_WW[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("NETWORK_NON_DIAG_LASSO_WW")}

Residuals_NON_DIAG_LASSO_WW <- sapply(model_results_non_diag_lasso_ww, function(y){
    res <- y$residuals
    })

## betweenness centrality max
dimnames(MATRIX_NON_DIAG_LASSO_WW)[[1]][which.max(centr_betw(NETWORK_NON_DIAG_LASSO_WW_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_NON_DIAG_LASSO_WW)[[1]][which.min(centr_betw(NETWORK_NON_DIAG_LASSO_WW_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_NON_DIAG_LASSO_WW)[[1]][which.max(centr_degree(NETWORK_NON_DIAG_LASSO_WW_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_NON_DIAG_LASSO_WW)[[1]][which.min(centr_degree(NETWORK_NON_DIAG_LASSO_WW_igraph)[["res"]])+1]

## total population. 
## total population without weights.

model_results_without_weights <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- as.matrix(dummy_cols(Xfit,select_columns = dimnames(Xfit)[[2]][71],ignore_na = T))
    p.fac = rep(1, 20)
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    x <- x[,-(71:74)]
    p.fac = rep(1, dim(x)[2])
    p.fac[c(71:74)] =0
    fit = cv.glmnet(x, y, alpha = 1, nlambda = 20,type.measure = "mse", nfolds = 20, penalty.factor = p.fac)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

save(model_results_without_weights,file = "model_results_without_weights.Rdata")


MATRIX_LASSO_WW <- sapply(model_results_without_weights, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+4)
    names(coef_result) <- model_results_without_weights[[1]][["COEF"]]@Dimnames[[1]]
    res <- y$COEF
    indexes <- match(dimnames(model_results_without_weights[[1]][["COEF"]])[[1]],names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

NETWORK_MATRIX_LASSO_WW_igraph <- graph_from_adjacency_matrix(MATRIX_LASSO_WW[2:71,],mode = "directed",weighted = TRUE)


Residuals_POPU_LASSO_WW <- sapply(model_results_without_weights, function(y){
    res <- y$residuals
    })



## betweenness centrality max
dimnames(MATRIX_LASSO_WW)[[1]][which.max(centr_betw(NETWORK_MATRIX_LASSO_WW_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_LASSO_WW)[[1]][which.min(centr_betw(NETWORK_MATRIX_LASSO_WW_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_LASSO_WW)[[1]][which.max(centr_degree(NETWORK_MATRIX_LASSO_WW_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_LASSO_WW)[[1]][which.min(centr_degree(NETWORK_MATRIX_LASSO_WW_igraph)[["res"]])+1]


## total population weights.
model_results_weights <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- as.matrix(dummy_cols(Xfit,select_columns = dimnames(Xfit)[[2]][71],ignore_na = T))
    p.fac = rep(1, 20)
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    weight <- x[,74]
    x <- x[,-(71:74)]
    p.fac = rep(1, dim(x)[2])
    p.fac[c(71:74)] =0
    fit = cv.glmnet(x, y, alpha = 1, weights = weight, nlambda = 20,type.measure = "mse", nfolds = 20, penalty.factor = p.fac)
    fitted_values <- predict(fit, x, s = 'lambda.min')
    resultAR <- list()
    resultAR$COEF <- coef(fit, s = "lambda.min")
    resultAR$COEF_1se <- coef(fit, s = "lambda.1se")
    resultAR$fit <- fit
    resultAR$residuals <-  y - fitted_values
    return(resultAR)                
})

MATRIX_LASSO_WEIGHT <- sapply(model_results_weights, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+1)
    names(coef_result) <- model_results_weights[[1]][["COEF"]]@Dimnames[[1]]
    res <- y$COEF
    indexes <- match(dimnames(model_results_weights[[1]][["COEF"]])[[1]],names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

NETWORK_LASSO_WEIGHT_igraph <- graph_from_adjacency_matrix(MATRIX_LASSO_WEIGHT[2:71,],mode = "directed",weighted = TRUE)


Residuals_POPU_LASSO <- sapply(model_results_weights, function(y){
    res <- y$residuals
    })

## betweenness centrality max
dimnames(MATRIX_LASSO_WEIGHT)[[1]][which.max(centr_betw(NETWORK_LASSO_WEIGHT_igraph)[["res"]])+1]
## betweenness centrality min
dimnames(MATRIX_LASSO_WEIGHT)[[1]][which.min(centr_betw(NETWORK_LASSO_WEIGHT_igraph)[["res"]])+1]
## degree centrality max
dimnames(MATRIX_LASSO_WEIGHT)[[1]][which.max(centr_degree(NETWORK_LASSO_WEIGHT_igraph)[["res"]])+1]
## degree centrality min
dimnames(MATRIX_LASSO_WEIGHT)[[1]][which.min(centr_degree(NETWORK_LASSO_WEIGHT_igraph)[["res"]])+1]


##total population without weights BIC
model_results_BIC <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    x <- x[,-(72:74)]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    ##weights_1 = rep(7/17,nrow(data_fit))
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

BIC_NETWORK_POPU <- sapply(model_results_BIC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_POPU_PVAL <- sapply(model_results_BIC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


BIC_NETWORK_POPU[!BIC_NETWORK_POPU_PVAL < 0.05] <- 0

BIC_NETWORK_POPU_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_POPU[2:71,],mode = "directed",weighted = TRUE)

# NETWORK_POPU <- {qgraph(BIC_NETWORK_POPU[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_DIAG_NO_WEIGHTS")}

## betweenness centrality max
dimnames(BIC_NETWORK_POPU)[[1]][which.max(centr_betw(BIC_NETWORK_POPU_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_POPU)[[1]][which.min(centr_betw(BIC_NETWORK_POPU_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_POPU)[[1]][which.max(centr_degree(BIC_NETWORK_POPU_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_POPU)[[1]][which.min(centr_degree(BIC_NETWORK_POPU_igraph)[["res"]])]

##total population weights BIC
model_results_BIC_WW <-lapply(1:nROI, function(i){
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
    x <- x[,-(72:74)]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    ##weights_1 = rep(7/17,nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    data_fit$handedness <- factor(data_fit$handedness)
    model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
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

BIC_NETWORK_POPU_WEIGHTS <- sapply(model_results_BIC_WW, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

BIC_NETWORK_POPU_WEIGHTS_PVAL <- sapply(model_results_BIC_WW, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


BIC_NETWORK_POPU_WEIGHTS[!BIC_NETWORK_POPU_WEIGHTS_PVAL < 0.05] <- 0

BIC_NETWORK_POPU_WEIGHTS_igraph <- graph_from_adjacency_matrix(BIC_NETWORK_POPU_WEIGHTS[2:71,],mode = "directed",weighted = TRUE)

# NETWORK_POPU <- {qgraph(BIC_NETWORK_POPU_WEIGHTS[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("BIC_NETWORK_POPU_WEIGHTS")}


## betweenness centrality max
dimnames(BIC_NETWORK_POPU_WEIGHTS)[[1]][which.max(centr_betw(BIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]
## betweenness centrality min
dimnames(BIC_NETWORK_POPU_WEIGHTS)[[1]][which.min(centr_betw(BIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]
## degree centrality max
dimnames(BIC_NETWORK_POPU_WEIGHTS)[[1]][which.max(centr_degree(BIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]
## degree centrality min
dimnames(BIC_NETWORK_POPU_WEIGHTS)[[1]][which.min(centr_degree(BIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]

##total population without weights AIC
model_results_AIC <-lapply(1:nROI, function(i){
    cat(i)
    colnames(Yfit)[1:70] <- paste("Y",1:70,sep = "")
    
    yname <- colnames(Yfit)[i]
    y <- Yfit[,i]
    #y <- unname(y) 
    x <- Xfit
    complete_obs <- complete.cases(y) & complete.cases(x)
    y <- y[complete_obs]
    x <- x[complete_obs,]
    x <- x[,-(72:74)]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    ##weights_1 = rep(7/17,nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    data_fit$handedness <- factor(data_fit$handedness)
    model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
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

AIC_NETWORK_POPU <- sapply(model_results_BIC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_POPU_PVAL <- sapply(model_results_BIC, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


AIC_NETWORK_POPU[!AIC_NETWORK_POPU_PVAL < 0.05] <- 0

AIC_NETWORK_POPU_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_POPU[2:71,],mode = "directed",weighted = TRUE)

# NETWORK_POPU_AIC <- {qgraph(AIC_NETWORK_POPU[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("AIC_NETWORK_POPU")}


## betweenness centrality max
dimnames(AIC_NETWORK_POPU)[[1]][which.max(centr_betw(AIC_NETWORK_POPU_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_POPU)[[1]][which.min(centr_betw(AIC_NETWORK_POPU_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_POPU)[[1]][which.max(centr_degree(AIC_NETWORK_POPU_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_POPU)[[1]][which.min(centr_degree(AIC_NETWORK_POPU_igraph)[["res"]])]

##total population weights AIC
model_results_AIC_WW <-lapply(1:nROI, function(i){
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
    x <- x[,-(72:74)]
    colnames(x)[71] <- "handedness"
    data_fit <- as.data.frame(cbind(Yfit[complete_obs,i],1,x))
    ##weights_1 = rep(7/17,nrow(data_fit))
    names(data_fit)[1] <- paste("Y",i,sep = "")
    
    data_fit$handedness <- factor(data_fit$handedness)
    model_fl_lower <- as.formula(paste(yname,"~1+handedness"))
    model_lower <- lm(model_fl_lower, data_fit, weights = weight)
    model_fl_upper <- as.formula(paste(yname,"~."))
    model_upper <- lm(model_fl_upper, data_fit, weights = weight)
    
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

AIC_NETWORK_POPU_WEIGHTS <- sapply(model_results_AIC_WW, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$COEF
    indexes <- match(names(y[["COEF"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})

AIC_NETWORK_POPU_WEIGHTS_PVAL <- sapply(model_results_AIC_WW, function(y){
    coef_result <- rep(0,length(colnames(Xfit))+3)
    names(coef_result) <- paste( "`",colnames(Xfit),"`",sep = "")
    names(coef_result)[75] <- "handedness2"
    names(coef_result)[76] <- "handedness3"
    res <- y$PVAL
    indexes <- match(names(y[["PVAL"]]),names(coef_result))
    coef_result[indexes] <- res
    coef_result
})


AIC_NETWORK_POPU_WEIGHTS[!AIC_NETWORK_POPU_WEIGHTS_PVAL < 0.05] <- 0

AIC_NETWORK_POPU_WEIGHTS_igraph <- graph_from_adjacency_matrix(AIC_NETWORK_POPU_WEIGHTS[2:71,],mode = "directed",weighted = TRUE)

# NETWORK_POPU_AIC <- {qgraph(AIC_NETWORK_POPU_WEIGHTS[2:71,],
#                                    posCol="blue",
#                                    negCol="purple")
#   title("AIC_NETWORK_POPU_WEIGHTS")}


## betweenness centrality max
dimnames(AIC_NETWORK_POPU_WEIGHTS)[[1]][which.max(centr_betw(AIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]
## betweenness centrality min
dimnames(AIC_NETWORK_POPU_WEIGHTS)[[1]][which.min(centr_betw(AIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]
## degree centrality max
dimnames(AIC_NETWORK_POPU_WEIGHTS)[[1]][which.max(centr_degree(AIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]
## degree centrality min
dimnames(AIC_NETWORK_POPU_WEIGHTS)[[1]][which.min(centr_degree(AIC_NETWORK_POPU_WEIGHTS_igraph)[["res"]])]

## graph density DIAGNOSED_GROUP WITH WEIGHTS_BIC
graph.density(BIC_NETWORK_DIAG_WEIGHTS_igraph)

## graph density DIAGNOSED_GROUP WITHOUT WEIGHTS_BIC
graph.density(BIC_NETWORK_DIAG_NO_WEIGHTS_igraph)

## graph density DIAGNOSED_GROUP cross validation WITH_WEIGHTS
graph.density(NETWORK_DIAG_LASSO_WITH_WEIGHTS)

## graph density DIAGNOSED_GROUP cross validation WITH_WEIGHTS
graph.density(NETWORK_DIAG_LASSO_WW_igraph)

## graph density NON_DIAGNOSED_GROUP cross validation WITH WEIGHTS
graph.density(NETWORK_NON_DIAG_LASSO_igraph)

## graph density NON_DIAGNOSED_GROUP cross validation WITHOUT WEIGHTS
graph.density(NETWORK_NON_DIAG_LASSO_WW_igraph)


## graph density TOTAL GROUP cross validation  without weights
graph.density(NETWORK_MATRIX_LASSO_WW_igraph)


## graph density TOTAL GROUP cross validation  with weights
graph.density(NETWORK_LASSO_WEIGHT_igraph)

## graph density TOTAL GROUP BIC  without weights
graph.density(NETWORK_MATRIX_LASSO_WW_igraph)

## graph density TOTAL GROUP BIC  with weights
graph.density(BIC_NETWORK_POPU_WEIGHTS_igraph)