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
  complete_obs <- complete.cases(y)
 & complete.cases(x)
  y <- y[complete_obs]
  x <- x[complete_obs,]
  armod <- lm(y ~ x)
  sarmod <- summary(armod)
  ARcoef <- sarmod$coefficients[-1,1]
  ARcoefPval <- sarmod$coefficients[-1,4]
  Interc <- ARcoef[1]
  resultAR <- list()
  resultAR$COEF <- ARcoef
  resultAR$PVAL <- ARcoefPval
  resultAR$INTERC <- Interc
  resultAR$RESID <- armod$residuals
  resultAR$FITTED <- armod$fitted.values
  resultAR$STD.DEV <- sarmod$sigma
  resultAR$R.SQUARED <- sarmod$r.squared
  return(resultAR)
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
?ggnet2

network_net <- network(VAR_mat, directed = TRUE )

## add attributed to vertices and edges 
network_net %v% "size" <-  deg_centr$res
network_net %e% "weight" <-  c(VAR_mat)
network_net %e% "color" <-  ifelse(network_net %e% "weight" >0 , "green", "tomato")



ggnet2(network_net, node.size = "degree", 
       arrow.size = 8, arrow.gap = 0.025,
       label = TRUE, color = "darkgreen", label.color = "black",
       edge.color = "color")+
  guides( size = FALSE, color=FALSE)