JarqueBera <- function (x) 
{
  obj.name <- deparse(substitute(x))
  K <- x$call$descriptives$p
  resid <- x$estimates$residuals
  obs <- nrow(resid)
  resids <- scale(resid, scale = FALSE)
  jbm.resids <- jbtest(resids, obs = obs, K = K, obj.name = obj.name)
  return(jbm.resids)
}

jbtest <- function (x, obs, K, obj.name) 
{
  P <- chol(crossprod(data.matrix(x[complete.cases(x),]))/obs)
  resids.std <- data.matrix(x[complete.cases(x),]) %*% solve(P)
  b1 <- apply(resids.std, 2, function(x) sum(x^3)/obs)
  b2 <- apply(resids.std, 2, function(x) sum(x^4)/obs)
  s3 <- obs * t(b1) %*% b1/6
  s4 <- obs * t(b2 - rep(3, K)) %*% (b2 - rep(3, K))/24
  STATISTIC <- s3 + s4
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- 2 * K
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "JB-Test (multivariate)"
  result1 <- list(statistic = STATISTIC, parameter = PARAMETER, 
                  p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", 
                                                                     obj.name))
  class(result1) <- "htest"
  STATISTIC <- s3
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- K
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Skewness only (multivariate)"
  result2 <- list(statistic = STATISTIC, parameter = PARAMETER, 
                  p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", 
                                                                     obj.name))
  class(result2) <- "htest"
  STATISTIC <- s4
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- K
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Kurtosis only (multivariate)"
  result3 <- list(statistic = STATISTIC, parameter = PARAMETER, 
                  p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", 
                                                                     obj.name))
  class(result3) <- "htest"
  result <- list(JB = result1, Skewness = result2, Kurtosis = result3)
  return(result)
}

serial_test <- function (resi, p, K, lags.pt = 16, lags.bg = 5, type = c("PT.asymptotic", 
                                                                         "PT.adjusted", "BG", "ES")) 
{
  #obj.name <- deparse(substitute(x))
  type <- match.arg(type)
  K <- K
  resids <- resi
  obs <- nrow(resids)
  if ((type == "PT.asymptotic") || (type == "PT.adjusted")) {
    lags.pt <- abs(as.integer(lags.pt))
    ptm <- pt_multi(p=p, K = K, obs = obs, lags.pt = lags.pt, resids = resids)
    ifelse(type == "PT.asymptotic", test <- ptm[[1]], test <- ptm[[2]])
  }
  else {
    lags.bg <- abs(as.integer(lags.bg))
    bgm <- bg_serial(x, K = K, obs = obs, lags.bg = lags.bg, 
                     obj.name = obj.name, resids = resids)
    ifelse(type == "BG", test <- bgm[[1]], test <- bgm[[2]])
  }
  return(test)
}

pt_multi <- function (p, K, obs, lags.pt, resids) 
{
  C0 <- crossprod(resids)/obs
  C0inv <- solve(C0)
  tracesum <- rep(NA, lags.pt)
  for (i in 1:lags.pt) {
    Ut.minus.i <- matlag1(resids, lag = i)[-c(1:i), ]
    Ut <- resids[-c(1:i), ]
    Ci <- crossprod(Ut, Ut.minus.i)/obs
    tracesum[i] <- sum(diag(t(Ci) %*% C0inv %*% Ci %*% C0inv))
  }
  vec.adj <- obs - (1:lags.pt)
  Qh <- obs * sum(tracesum)
  Qh.star <- obs^2 * sum(tracesum/vec.adj)
  nstar <- K^2 * p
  STATISTIC <- Qh
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- (K^2 * lags.pt - nstar)
  
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Portmanteau Test (asymptotic)"
  PT1 <- list(statistic = STATISTIC, parameter = PARAMETER, 
              p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object"))
  class(PT1) <- "htest"
  STATISTIC <- Qh.star
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- (K^2 * lags.pt - nstar)
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Portmanteau Test (adjusted)"
  PT2 <- list(statistic = STATISTIC, parameter = PARAMETER, 
              p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object"))
  class(PT2) <- "htest"
  result <- list(PT1 = PT1, PT2 = PT2)
  return(result)
}


matlag1 <- function (x, lag = 1) 
{
  totcols <- ncol(x)
  nas <- matrix(NA, nrow = lag, ncol = totcols)
  x <- rbind(nas, x)
  totrows <- nrow(x)
  x <- x[-c((totrows - lag + 1):totrows), ]
  return(x)
}
qqPlot<-function(x, ...) {
  UseMethod("qqPlot")
}

qqPlot.default <- function(x, distribution="norm", groups, layout,
                           ylim=range(x, na.rm=TRUE), ylab=deparse(substitute(x)),
                           xlab=paste(distribution, "quantiles"), glab=deparse(substitute(groups)),
                           main=NULL, las=par("las"),
                           envelope=.95,
                           col=carPalette()[1], col.lines=carPalette()[2], lwd=2, pch=1, cex=par("cex"),
                           line=c("quartiles", "robust", "none"), id=TRUE, grid=TRUE, ...){
  if (!missing(groups)){
    if (isTRUE(id)) id <- list(n=2)
    if (is.null(id$labels)) id$labels <- seq(along=x)
    grps <- levels(as.factor(groups))
    if (missing(layout)) layout <- mfrow(length(grps), max.plots=12)
    if (prod(layout) < length(grps)) stop("layout cannot accomodate ", length(grps), " plots")
    oldpar <- par(mfrow=layout)
    on.exit(par(oldpar))
    for (group in grps){

      id.gr <- id
      if (!isFALSE(id)) id.gr$labels <- id$labels[groups == group]
      qqPlot.default(x[groups == group], distribution=distribution, ylim=ylim, ylab=ylab,
                     xlab=xlab, main=paste(glab, "=", group), las=las, envelope=envelope, col=col,
                     col.lines=col.lines, pch=pch, cex=cex, line=line, id=id.gr, grid=grid, ...)
    }
    return(invisible(NULL))
  }
  id <- applyDefaults(id, defaults=list(method="y", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
  if (isFALSE(id)){
    id.n <- 0
    id.method <- "none"
    labels <- id.cex <- id.col <- id.location <- NULL
  }
  else{
    labels <- id$labels
    if (is.null(labels)) labels <- if(!is.null(names(x))) names(x) else seq(along=x)
    id.method <- id$method
    id.n <- if ("identify" %in% id.method) Inf else id$n
    id.cex <- id$cex
    id.col <- id$col
    id.location <- id$location
  }
  line <- match.arg(line)
  index <- seq(along=x)
  good <- !is.na(x)
  ord <- order(x[good])
  if (length(col) == length(x)) col <- col[good][ord]
  if (length(pch) == length(x)) pch <- pch[good][ord]
  if (length(cex) == length(x)) cex <- cex[good][ord]
  ord.x <- x[good][ord]
  ord.lab <- labels[good][ord]
  q.function <- eval(parse(text=paste("q", distribution, sep="")))
  d.function <- eval(parse(text=paste("d", distribution, sep="")))
  n <- length(ord.x)
  P <- ppoints(n)
  z <- q.function(P, ...)
  plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las, ylim=ylim)
  if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
  points(z, ord.x, col=col, pch=pch, cex=cex)
  if (line == "quartiles" || line == "none"){
    Q.x <- quantile(ord.x, c(.25,.75))
    Q.z <- q.function(c(.25,.75), ...)
    b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
    a <- Q.x[1] - b*Q.z[1]
    if (line == "quartiles") abline(a, b, col=col.lines, lwd=lwd)
  }
  if (line=="robust") {
    coef <- coef(rlm(ord.x ~ z))
    a <- coef[1]
    b <- coef[2]
    abline(a, b, col=col.lines, lwd=lwd)
  }
  conf <- if (envelope == FALSE) .95 else envelope
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/d.function(z, ...))*sqrt(P*(1 - P)/n)
  fit.value <- a + b*z
  upper <- fit.value + zz*SE
  lower <- fit.value - zz*SE
  if (envelope != FALSE) {
    lines(z, upper, lty=2, lwd=lwd, col=col.lines)
    lines(z, lower, lty=2, lwd=lwd, col=col.lines)
  }
  extreme <- showLabels(z, ord.x, labels=ord.lab,
                        method = id.method, n = id.n, cex=id.cex, col=id.col, location=id.location)
  if (is.numeric(extreme)){
    nms <- names(extreme)
    extreme <- index[good][ord][extreme]
    if (!all(as.character(extreme) == nms)) names(extreme) <- nms
  }
  if (length(extreme) > 0) extreme else invisible(NULL)
}

qqPlot.formula <- function (formula, data, subset, id=TRUE, ylab, glab, ...) {
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
    m$data <- as.data.frame(data)
  m$formula <- m$... <- m$id <- m$ylab <- NULL
  m[[1]] <- as.name("model.frame")
  if (missing(ylab)) ylab <- as.character(formula[[2]])
  if (length(formula) == 2){
    groups <- FALSE
  }
  else if (length(formula) == 3){
    groups <- TRUE
    if(missing(glab)) glab <- as.character(formula[[3]])
  }
  m$formula <-formula
  if (missing(data)){
    x <- na.omit(eval(m, parent.frame()))
  }
  else{
    x <- eval(m, parent.frame())
  }
  if (!isFALSE(id)){
    if (isTRUE(id)){
      id <- list(n=2)
    }
    if (is.null(id$labels)){
      id$labels <- rownames(x)
    }
  }
  if (!groups && ncol(x) > 1) stop("more than one variable specified")
  if (groups && ncol(x) != 2) stop("formula must be of the form variable ~ groups")
  if (!groups) qqPlot(x[, 1], id=id, ylab=ylab, ...)
  else qqPlot(x[, 1], groups=x[, 2], id=id, ylab=ylab, glab=glab, ...)
}

applyDefaults <- function(args, defaults, type=""){
  if (isFALSE(args)) return(FALSE)
  names <- names(args)
  names <- names[names != ""]
  if (!isTRUE(args) && !is.null(args) && length(names) != length(args)) warning("unnamed ", type, " arguments, will be ignored")
  if (isTRUE(args) || is.null(names)) defaults
  else defaults[names] <- args[names]
  as.list(defaults)
}
