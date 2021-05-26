library(ezsim)
library(MVN)
library(vars)
library(sarima)
library(gtools)
library(car)
library(ggplot2)
library(grid)
library(gridExtra)
library(plot.matrix)

diagnostics <-  function(x, data.fit, AR.MAT, lm.fit,  plot, uni.acf.test="LjungBox"){
  K=70
  p=1
  if (missing(x)) stop("x must be provided")
  if (missing(data.fit)) stop("data.fit must be provided")
  if (missing(AR.MAT)) stop("AR.MAT must be provided")
  if (missing(lm.fit)) stop("lm.fit must be provided")
  
  
  if (missing(plot)) plot <- TRUE
  if (missing(uni.acf.test)) uni.acf.test <- "LjungBox"
  if (!uni.acf.test %in% c("LjungBox", "LiMcLeod","Breusch-Godfrey")) 
  stop("The acf test can be one of: 'LjungBox','LiMcLeod' or 'Breusch-Godfrey'")
  #if (x$call$model.select == "none") stop("No model has been selected for diagnostic checks. Use the function in order to select an shVAR model.")
  
  ID <- NULL
  out <- vector("list", 3L)
  names(out) <- c("stationarity","multivariateNormality", "autocorrelation")
  out$stationarity <- list()
  out$stationarity$roots <- {
    A <- t(AR.MAT)[1:70,1:70]
    companion <- matrix(0, nrow = K * p, ncol = K * p)
    companion[1:K, 1:(K * p)] <- A
    if (p > 1) {
          j <- 0
        for (i in (K + 1):(K * p)) {
            j <- j + 1
            companion[i, j] <- 1
          }
      }
    roots <- eigen(companion)$values
    Mod(roots)
  }

  if (all(out$stationarity$roots < 1)){
    out$stationarity$is.stationary <- TRUE
    cat("The moduli of the roots of the characteristic polynomial is less than 1. The process is stationary.\n") 
  }else{
    out$stationarity$is.stationary <- FALSE
    cat("The moduli for some of the roots of the characteristic polynomial is greater than or equal to 1. The process is non stationary.\n")
  }
  
  ## Multivariate normality
  normality_tests <- data.frame(matrix(NA, ncol = 4, nrow = 7))
  dimnames(normality_tests)[[2]] <- c("Test", "Statistic", "p value", "package")
  normality_tests[1:5,4] <-  rep("mvn",5)
  normality_tests[6:7,4] <-  rep("vars",2)
  test_mard <- mvn(x, mvnTest = "mardia")
  test_hz <- mvn(x, mvnTest = "hz")
  test_dh <- mvn(x, mvnTest = "dh")
  test_JB <- jbtest(x = x, obs = nrow(x), K=70,obj.name ="testJB")
  normality_tests[1:2,2:3] <- sapply(test_mard$multivariateNormality[1:2,2:3],function(y) as.numeric(as.character(y)))
  normality_tests[1:2,1] <- as.character(test_mard$multivariateNormality[1:2,1]) 
  normality_tests[3,2:3] <- as.numeric(test_hz$multivariateNormality[2:3])
  normality_tests[3,1] <- as.character(test_hz$multivariateNormality[1,1])
  normality_tests[4,2:3] <- as.numeric(test_dh$multivariateNormality[c(2,4)])
  normality_tests[4,1] <- as.character(test_dh$multivariateNormality[1,1])
  normality_tests[5,2] <- as.numeric(test_JB$JB$statistic)
  normality_tests[5,3] <- as.numeric(test_JB$JB$p.value)
  normality_tests[5,1] <- "Jarque-Bera"
  normality_tests[6,1] <- "Jarque-Bera Skewness"
  normality_tests[7,1] <- "Jarque-Bera Kurtosis"
  normality_tests[6,2] <-  test_JB$Skewness$statistic
  normality_tests[7,2] <-  test_JB$Kurtosis$statistic
  normality_tests[6,3] <-  test_JB$Skewness$p.value
  normality_tests[7,3] <-  test_JB$Kurtosis$p.value
  
  out$multivariateNormality <- normality_tests
  
  #l.matrix[upper.tri(l.matrix)] <- (70+2):(70+1+sum(upper.tri(l.matrix)))
  
  resid_mod <- x
  nresid <-  nrow(resid_mod)
  effective_timep <- nresid / length(unique(data.fit$ID))
  dtacf <- data.frame(resid_mod,ID=rep(1:7, each=effective_timep))
  
  autoc <- do.call(rbind,lapply(1:7, function(i){
    dti <- subset(dtacf, ID==i)[,1:70]
    pttest <- serial_test(resi =as.matrix(dti) ,p=1, K=70, type = "PT.adjusted")
    VEOUT <- c(pttest$statistic,pttest$p.value)
    names(VEOUT) <- c("statistic","p.value")
    VEOUT
    
  }))
  rownames(autoc) <-  paste("ID",1:7)
  out$autocorrelation <- autoc
  
  library("tseries")    
  comb_ts <- cbind(Y[1,,2,1], Y[1,,2,2]) # please make sure the length of both your timeseries
  color<-rainbow(ncol(comb_ts))
  ts.plot(comb_ts, plot.type = "single",col=color)
  
  if (plot){
    
    distances <- mahalanobis(dtacf[complete.cases(dtacf),1:70], colMeans(dtacf[,1:70], na.rm = T), cov(dtacf[complete.cases(dtacf),1:70],))
    qqPlot(distances, distribution = "chisq", df = mean(distances,na.rm = T),  lwd = 1, grid = FALSE,cex=2, 
           col=rainbow(length(unique(dtacf[,-seq(1,70)])))[as.numeric( as.factor(dtacf[,-seq(1,70)]))],
           main = "Multivariate normality Q-Q plot", xlab = expression(chi^2 * " quantiles"), ylab = expression("Mahalanobis distances "^2))
    l.matrix <- matrix(0, nrow = 70, ncol = 70)
    l.matrix[lower.tri(l.matrix, diag=FALSE)] <- 1:sum(lower.tri(l.matrix, diag=FALSE))
    l.matrix <- t(l.matrix)
    diag(l.matrix) <- max(l.matrix[upper.tri(l.matrix)]) + 1:70
    l.matrix[lower.tri(l.matrix, diag=FALSE)] <- max(diag(l.matrix)) + 1
    order_mat <- which(lower.tri(l.matrix, diag = FALSE), arr.ind = TRUE)
    maxdiag <- max(diag(l.matrix)) + 1
    l.matrix[order_mat[order(order_mat[, 1]), ]]  <- seq(maxdiag,maxdiag-1+sum(l.matrix==maxdiag))
    layout(l.matrix)
    cmbcols <- combinations(70, 2, v=paste("X",1:70,sep=""))
    cmbcols2 <- combinations(70, 2, v=paste("X",rev(1:70),sep="") ,set = FALSE)
    
    # #Multivariate normality test
    # apply(cmbcols,1, function(y){
    #   d = data.matrix(x[,y])
    #   mvn(data = d[complete.cases(d),], multivariatePlot = "contour", desc = FALSE)
    # })
    
    #
    sapply(1:10, function(i){
      qqnorm(x[,i], pch = 1, frame = TRUE, ylim=c(-4,4),
             main =paste("Normal Q-Q Plot for variable:", paste("X",1:70,sep="")[i]))
      qqline(x[,i], col = "steelblue", lwd = 2)
      tst1 <- shapiro.test(x[,i])
      tst2 <- wilcox.test(x[,i])
      text(2,-3,paste(paste("Shapiro normality test:", round(tst1$p.value,3)),"\n",
                      paste("Wilcoxon normality test:", round(tst2$p.value,3))))})
    
    apply(cmbcols2[rev(1:nrow(cmbcols2)),][1:10,],1, function(y){
      plot(x[,y[1]],x[,y[2]],scales = "free",
           xlab = paste("Residuals for",paste("X",1:70,sep="")[y[1]]),
           ylab = paste("Residuals for",paste("X",1:70,sep="")[y[2]]))
      abline(lm(x[,y[2]]~x[,y[1]]),col="red")
      text(2,-3,paste("Correlation:",round(cor(x[,y[1]],x[,y[2]]),3)))
      
    })
    
    ## fitted values vs the residuals plots
    par(mfrow=c(1,1))
    gglist_resvsfitted <- lapply(1:(70), function(i){
      dtplot <- data.frame(cbind(x[complete.cases(x),i],lm.fit[[i]]$FITTED ))
      colnames(dtplot) <- c("residuals", "fitted")
      p1<-ggplot(dtplot, aes(fitted, residuals))+geom_point()
      p1<-p1+stat_smooth(method="loess", formula = 'y ~ x')+geom_hline(yintercept=0, col="red", linetype="dashed")
      p1<-p1+xlab("Fitted values")+ylab("Residuals")
      p1<-p1+ggtitle(paste("Residual vs Fitted Plot for variable:",paste("X",1:70,sep="")[i]))+theme_bw()
      invisible(p1) 
    })

    #par(mfrow=c(ceiling(70 / 4) ,4))
    # nCol <- floor(sqrt(70))
    # do.call("grid.arrange", c(gglist_resvsfitted, ncol=nCol))
    
    ## auto correlation of the residuals
    ## The Breusch–Godfrey test is a test for autocorrelation in the errors in a regression model. 
    ## It makes use of the residuals from the model being considered in a regression analysis, and a test statistic is derived from these. The null hypothesis is that there is no serial correlation of any order up to p
    
    if (uni.acf.test=="Breusch-Godfrey"){
      cat("When uni.acf.test='Breusch-Godfrey' no plot for the pacf is returned ")
      acf_test_plot <- t(do.call("rbind",lapply(1:7, function(i){
        dti <- subset(dtacf, ID==i)[,1:70]
        par(mfrow=c(1,1))
        tests <- sapply(1:70,function(j){
          y <- dti[,j]
          Mlags <- cbind(
            filter(y, c(0,1), method= "conv", sides=1),
            filter(y, c(0,0,1), method= "conv", sides=1),
            filter(y, c(0,0,0,1), method= "conv", sides=1))
          colnames(Mlags) <- paste("lag", seq_len(ncol(Mlags)))
          # store p-value of the Breusch-Godfrey test
          if (any(c(bgtest(y ~ 1+Mlags, order=1, type="F", fill=NA)$p.value,
                    bgtest(y ~ 1+Mlags, order=2, type="F", fill=NA)$p.value,
                    bgtest(y ~ 1+Mlags, order=3, type="F", fill=NA)$p.value,
                    bgtest(y ~ 1+Mlags, order=4, type="F", fill=NA)$p.value) < 0.05)) return(1) else return(0)
          
        })
        return(tests)
      })))
      colnames(acf_test_plot) <- paste("ID",1:7)
      rownames(acf_test_plot) <- paste("X",1:70,sep="")
      plot(acf_test_plot==0,col=c('red', 'green'),xlab="ID", ylab="Variable",key = NULL,
           main="Violations of univariate autocorrelation for each ID ")
    }else{
      dotests <- lapply(1:7, function(i){
        #dti <- subset(dtacf[complete.cases(dtacf),1:70], ID==i)[,1:70]
        dti <- subset(dtacf[complete.cases(dtacf),], ID==i)[,1:70]
        checkpacf <- apply(dti,2,partialAutocorrelations,max.lag=10)
        checkacf <- apply(dti,2,autocorrelations,max.lag=10)
        
        testacf <- sapply(1:70,function(r){
          whiteNoiseTest(checkacf[[r]],h0="iid",nlags=5, x=dti[,r],n=length(dti[,r]),
                         method=uni.acf.test)$test[3]
        })
        testpacf <- sapply(1:70,function(r){
          whiteNoiseTest(checkpacf[[r]],h0="iid",nlags=5, x=dti[,r],n=length(dti[,r]),
                         method=uni.acf.test)$test[3]
        })
        return(list(acf=ifelse(testacf <0.05,1,0),pacf=ifelse(testpacf<0.05,1,0)))
      })
      acf_test_plot <- t(do.call("rbind",lapply(dotests, function(t) t$acf)))
      pacf_test_plot <- t(do.call("rbind",lapply(dotests, function(t) t$pacf)))
      colnames(acf_test_plot) <- paste("ID",1:7)
      rownames(acf_test_plot) <- paste("X",1:70,sep="")
      colnames(pacf_test_plot) <- paste("ID",1:7)
      rownames(pacf_test_plot) <- paste("X",1:70,sep="")
      plot(acf_test_plot==0,col=c('red', 'green'),xlab="ID", ylab="Variable",key = NULL,
                  main="Violations of univariate autocorrelation for each ID ")
      plot(pacf_test_plot==0,col=c('red', 'green'),xlab="ID", ylab="Variable",key = NULL,
                  main="Violations of univariate partial autocorrelation for each ID ")
      # plot_matrix(acf_test_plot==0,col=c('red', 'green'),xlab="ID", ylab="Variable",key = NULL,
      #             main="Violations of univariate autocorrelation for each ID ")
      # plot_matrix(pacf_test_plot==0,col=c('red', 'green'),xlab="ID", ylab="Variable",key = NULL,
      #             main="Violations of univariate partial autocorrelation for each ID ")
    }
    
    
    
  }
  return(out)
  
}
