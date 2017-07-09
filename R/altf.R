
altf <- function (y,x,window=NULL,initial.period=NULL,d=NULL)
  {

### computes some forecast quality measures for some alternative forecasts
 
### requires "forecast", "MSwM", "stats" and "xts" packages

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different independent variables

### window - a size of a rolling regression window, a number of observations,
###          if not specified 10% of all observations are taken

### initial.period - a number of observation since which forecast quality measures are computed,
###                  by default the whole sample is used, i.e., initial.period = 1

### d - logical, used for hit.ratio calculation,
###     d = FALSE for level time-series,
###     d = TRUE if time-series represent changes,
###     by default d = FALSE


### checking initial data 

if (missing(y)) { stop("please, specify y") }
if (missing(x)) { stop("please, specify x") }
if (! (is.numeric(y) || is.matrix(y))) { stop("y must be numeric or matrix") }
if (is.matrix(y) && ! (ncol(y) == 1)) { stop("y must be a one column matrix") }
if (! is.matrix(x)) { stop("x must be a matrix") }
if (is.null(colnames(x))) 
  { 
    colnames(x) <- colnames(x, do.NULL = FALSE, prefix = "X")
    warning('column names of x were automatically created') 
  }
if (anyNA(colnames(x))) { stop("x must have column names") }
if (is.matrix(y) && is.null(colnames(y))) 
  { 
    warning('column name of y was automatically created') 
    colnames(y) <- colnames(y, do.NULL = FALSE, prefix = "Y") 
  }
if (is.matrix(y) && anyNA(colnames(y))) 
  { 
    warning('column name of y was automatically created') 
    colnames(y) <- "Y1" 
  }
if (! length(y) == nrow(x)) { stop("y and x must have the same number of observations") }
if (anyNA(y)) { stop("missing values in y") }
if (anyNA(x)) { stop("missing values in x") }
if (is.null(window)) { window <- floor(length(y)/10) }
if (window < 1) { window <- 1 }
if (!(is.numeric(window))) { stop("window must be numeric") }
if ((window < 0) || (window > length(y))) { stop("window must be a positive number less then the total number of observations") }
if (length(y) < 3) { stop("time-series too short: there have to be more than 3 observations") }
if (is.null(initial.period)) { initial.period <- 1 }
if (! is.numeric(initial.period)) { stop("initial.period must be a number") }
if ((initial.period <= 0) || (initial.period > length(y))) { stop("initial.period must be greater than or equal to 1, and less than the number of observations") }
if (is.null(d)) { d <- FALSE }
if (! is.logical(d)) { stop("d must be logical, i.e., TRUE or FALSE") }
if (requireNamespace('forecast')) { } else { stop('package >>forecast<< is required') }
if (requireNamespace('stats')) { } else { stop('package >>stats<< is required') }
if (requireNamespace('MSwM')) { } else { stop('package >>MSwM<< is required') }
requireNamespace('xts')

y <- as.matrix(y)
x <- as.matrix(x)

######################### naive
##################################################

y.naive <- (c(NA,as.vector(y)))[-(1+length(as.vector(y)))]

######################### OLS
##################################################

y.ols <- lm(y ~ x)$fitted.values

######################### recursive OLS
##################################################

y.rec.ols <- vector()

y.rec.ols[1] <- lm(y[1] ~ t(x[1,]))$fitted.values[1]

for (i in 2:nrow(x))
  {
    y.rec.ols[i] <- lm(y[1:i] ~ x[1:i,])$fitted.values[i]
  }

######################### rolling OLS
##################################################

y.roll.ols <- roll.reg(y=y,x=x,window=window)

######################### TVP
##################################################

y.tvp <- tvp(y=y,x=x,V=1,lambda=0.99)

######################### AR(1)
##################################################

x1 <- (as.vector(y))[1:(length(as.vector(y))-1)]

yy <- (as.vector(y))[-1]

y.ar1 <- lm(yy ~ x1)$fitted.values
y.ar1 <- c(NA,y.ar1)

rm(x1,yy)

######################### AR(2)
##################################################

x1 <- (as.vector(y))[2:(length(as.vector(y))-1)]
x2 <- (as.vector(y))[1:(length(as.vector(y))-2)]
x2 <- cbind(x1,x2)

yy <- (as.vector(y))[-c(1,2)]

y.ar2 <- lm(yy ~ x2)$fitted.values
y.ar2 <- c(NA,NA,y.ar2)

rm(x1,yy,x2)

######################### auto ARIMA
##################################################

y.auto.arima <- y - auto.arima(as.vector(y))$residuals 

######################### TVP-AR(1)
##################################################

x1 <- (as.vector(y))[1:(length(as.vector(y))-1)]
yy <- (as.vector(y))[-1]

y.tvp.ar1 <- tvp(y=yy,x=x1,V=1,lambda=0.99)
y.tvp.ar1 <- c(NA,y.tvp.ar1)

rm(x1,yy)

######################### TVP-AR(2)
##################################################

x1 <- (as.vector(y))[2:(length(as.vector(y))-1)]
x2 <- (as.vector(y))[1:(length(as.vector(y))-2)]
x2 <- cbind(x1,x2)

yy <- (as.vector(y))[-c(1,2)]

y.tvp.ar2 <- tvp(y=yy,x=x1,V=1,lambda=0.99)
y.tvp.ar2 <- c(NA,NA,y.tvp.ar2)

rm(x1,yy,x2)

######################### Markov Switching Models
##################################################

ob <- lm(y ~ x)
msm <- msmFit(object=ob,k=2,sw=rep(TRUE,(ncol(x)+2)))

y.ms <- msm@model$fitted.values

rm(ob,msm)

##################################################

fq <- rbind(
accuracy(f=(as.vector(y.naive))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.ols))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.rec.ols))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.roll.ols))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.tvp))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.ar1))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
forecast::accuracy(f=(as.vector(y.ar2))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.auto.arima))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.tvp.ar1))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.tvp.ar2))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))]),
accuracy(f=(as.vector(y.ms))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))])
)

hr <- rbind(
hit.ratio(y=as.vector(y),y.hat=as.vector(y.naive),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.ols),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.rec.ols),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.roll.ols),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.tvp),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.ar1),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.ar2),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.auto.arima),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.tvp.ar1),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.tvp.ar2),d=d),
hit.ratio(y=as.vector(y),y.hat=as.vector(y.ms),d=d)
)

fq <- cbind(fq,hr)

rownames(fq) <- c("naive","OLS","rec. OLS","roll. OLS","TVP","AR(1)","AR(2)","auto ARIMA","TVP-AR(1)","TVP-AR(2)","MS")
colnames(fq)[6] <- c("HR")

return(round(fq,digits=4))
  
  }
  