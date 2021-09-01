#load("/Users/shalinipandurangan/Desktop/04-Data/fMRI-ROI-time-series/fMRI-ROI-time-series.RData")
load("D:/WORK/SUPERVISION/Shalini/04-Data/fMRI-ROI-time-series/fMRI-ROI-time-series/fMRI-ROI-time-series.RData")
# Y is the data frame of fMRI-ROI-time-series data
# y values of all the subjects of all 70 regions of interest

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

model_results <-lapply(1:nROI, function(i){    cat(i)    y <- Yfit[,i]    y <- unname(y)        x <- Xfit    complete_obs <- complete.cases(y) & complete.cases(x)    y <- y[complete_obs]    x <- x[complete_obs,]    data1 <-data.frame(cbind(y,x))    ##y_name <- gsub( "_.*$", "", names(data1)[i] )    colnames(x) <- gsub("-",".",colnames(x))    ##names(data1)[1] <- y_name    model_formulae <- paste("y","~",paste(colnames(x),collapse = "+"))    ##lambda_seq <- 10^seq(2, -2, by = -.1)    ##cv_output <- cv.glmnet(x, y,                           ##alpha = 1, lambda = lambda_seq,                            ##nfolds = 5)
    lamdaseq(x,y,weight = NA, alpha = 1, standardis = TRUE, lamdaRatio = 1e-04, nLamda = 100, addZeroLamda = FALSE)                       
    best_lam <- list() 
    lamda <- lamdaseq(x,y)$lamda
    ##best_lam <- cv_output$lambda.min                                                                       
    print(lamda)})

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
  up <- paste("y","~",paste(colnames(x)[1:70],collapse = "+"))
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

## Get the Autoregressive matrix 
VAR_mat <- do.call(rbind, lapply(model_results, function(x) x$COEF))
dimnames(VAR_mat) <-  list(colnames(Yfit),colnames(Yfit))

## Lets set some elements to zero based on the signifance of the p values
VAR_mat_Pval <- do.call(rbind, lapply(model_results, function(x) x$PVAL))
VAR_mat[!VAR_mat_Pval < 0.05] <- 0

## Transpose to get the adjacency
VAR_mat <-  t(VAR_mat)

## Get the covariance matrix
Sigma_diag <- diag(sapply(model_results, function(x) x$STD.DEV)^2)

## Get the residuals
residuals_VAR <-  do.call(cbind,lapply(model_results, function(x) x$RESID))

## Get the sample covariance matrix of the residuals
Sigma_residuals <-  cov(do.call(cbind,lapply(model_results, function(x) x$RESID)))


##############################################
############## PLOT the results ##############
##############################################

################################## FIRST WAY OF PLOTTING NETWORKS : IGRAPH 

library(igraph)
igraph_net <- graph_from_adjacency_matrix(VAR_mat, mode = "directed", weighted = TRUE )

## Basic plot
plot(igraph_net)
### Change some graphical parameters if you like

## For example you may want to change the size of the nodes according to 
## their degree centrality
deg_centr <- centr_degree(igraph_net)


### Layout 1 
plot(igraph_net, 
     edge.arrow.width = 0.1,
     edge.arrow.size = 0.5,
     size= deg_centr$res/10, 
     layout = layout.fruchterman.reingold(igraph_net), ## Change to a circle 
     edge.color= ifelse(E(igraph_net)$weight > 0, "green", "red" ) ## change edge color (green=positive, red=negative)
)

### Layout 1 
plot(igraph_net, 
     edge.arrow.width = 0.1,
     edge.arrow.size = 0.5,
     size= deg_centr$res/10, 
     layout = layout.circle(igraph_net), ## Change to a circle 
     edge.color= ifelse(E(igraph_net)$weight > 0, "green", "red" ) ## change edge color (green=positive, red=negative)
)



################################## SECOND WAY OF PLOTTING NETWORKS : NETWORK

library(network)

network_net <- network(VAR_mat, directed = TRUE )

## add attributed to vertices and edges 
network_net %v% "size" <-  deg_centr$res
network_net %e% "weight" <-  c(VAR_mat)
network_net %e% "color" <-  ifelse(network_net %e% "weight" >0 , "green", "tomato")

## Plot 
plot(network_net, 
     label  = network.vertex.names(network_net), # labels
     pad=0,
     label.pad=0,
     vertex.cex = get.vertex.attribute(network_net, "size")/20, # size of the nodes
     edge.col = network::get.edge.value(network_net, "color") # color of the edges
)



################################## Third WAY OF PLOTTING NETWORKS : QGRAPH


library(qgraph)

qgraph_net <- qgraph(VAR_mat,
                     posCol="blue",
                     negCol="purple")



################################## Fourth way of plotting graphs  : GGNET
library(GGally)
library()
library(network)
#library(ggnet2)

network_net <- network(VAR_mat[1,], directed = TRUE )

## add attributed to vertices and edges 
network_net %v% "size" <-  deg_centr$res
network_net %e% "weight" <-  c(VAR_mat)
network_net %e% "color" <-  ifelse(network_net %e% "weight" >0 , "green", "tomato")



ggnet2(network_net, node.size = "degree", 
       arrow.size = 8, arrow.gap = 0.025,
       label = TRUE, color = "darkgreen", label.color = "black",
       edge.color = "color")+
  guides( size = FALSE, color=FALSE)