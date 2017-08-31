\name{altf2}
\alias{altf2}

\title{Computes a Few Alternative Forecasts Based on Model Averaging.}

\description{
It is necessary to compare a given forecast method with some alternative ones. This function computes selected forecast quality measures for a few selected forecast methods (which might be treated as alternative ones to Dynamic Model Averaging, Dynamic Model Selection, etc.). 

ME (Mean Error), RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), MPE (Mean Percentage Errror) and MAPE (Mean Absolute Percentage Error) are computed as \code{\link[forecast]{accuracy}}. HR (Hit Ratio) is computed as \code{\link{hit.ratio}}.
}

\usage{
altf2(y,x,mods.incl=NULL,av=NULL,window=NULL,initial.period=NULL,d=NULL,
f=NULL,fmod=NULL,parallel=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns correspond to different independent variables}
\item{mods.incl}{optional, \code{\link[base]{matrix}} indicating which models will be used in averaging, if not specified all possible models will be used, see \code{\link{fDMA}}}
\item{av}{optional, a method for model averaging, \code{av="ord"} corresponds to equal weights for each model, \code{av="aic"} corresponds to information theoretic model averaging based on Akaike Information Criterion, \code{av="bic"} corresponds to information theoretic model averaging based on Bayesian Information Criterion, \code{av="mse"} corresponds to setting weights proportional to the inverse of the models Mean Squared Error, if not specified \code{av="ord"} is used}
\item{window}{optional, \code{\link[base]{numeric}}, a size of a rolling regression window (a number of observations), if not specified 10\% of all observations are taken}
\item{initial.period}{optional, \code{\link[base]{numeric}}, a number of observation since which forecast quality measures are computed, if not specified the whole sample is used, i.e., \code{initial.period=1}}
\item{d}{optional, \code{\link[base]{logical}}, a parameter used for HR (Hit Ratio) calculation, should be \code{d=FALSE} for level time-series and \code{d=TRUE} if time-series represent changes, if not specified \code{d=FALSE}}
\item{f}{optional, \code{\link[base]{logical}} vector, indicating which of alternative forecast -- av. OLS, av. rec. OLS, av. roll. OLS and av. TVP -- should be averaged, if not specified \code{f=c(\link[base]{rep}(TRUE,4)}, i.e., all alternative forecast are computed}
\item{fmod}{optional, class \code{dma} object, a model to be compared with alternative forecast}
\item{parallel}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
}

\details{For each \code{av} method, in the initial period equal weights for each model are taken, and then successively updated based on the chosen criterion. For OLS models weights are not updated. The same weight for each model (estimated from the whole sample) is taken for each period.}

\value{class \code{altf2} object, \code{\link[base]{list}} of
\item{$summary}{\code{\link[base]{matrix}} of forecast quality measures ordered by columns, forecast methods are ordered by rows}
\item{$y.hat}{\code{\link[base]{list}} of predicted values from all forecasting methods which were applied}
\item{$y}{\code{y}, forecasted time-series}
\item{$coeff.}{\code{\link[base]{list}} of coefficients from all forecasting methods which were applied}
\item{$weights}{\code{\link[base]{list}} of weights of models used in averaging for all forecasting methods which were applied}
\item{$p.val.}{\code{\link[base]{list}} of p-values (averaged with respect to suitable weights) for t-test of statistical significance for coefficients from all forecasting methods which were applied (for TVP they are not computed)}
}

\examples{
\dontrun{
# models for untransformed data
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a1 <- altf2(y=wti,x=drivers)

# do not include first 12 observations for forecast quality measures,
# i.e., treat first 12 observations as a ''training set''
a2 <- altf2(y=wti,x=drivers,initial.period=12)

# models for log-differenced data
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
a3 <- altf2(y=ld.wti,x=ld.drivers,d=TRUE)

fcomp <- c(TRUE,TRUE,TRUE,FALSE)
a4 <- altf2(y=ld.wti,x=ld.drivers,d=TRUE,f=fcomp)
a5 <- altf2(y=ld.wti,x=ld.drivers,d=TRUE,f=fcomp,av="aic")

m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
a6 <- altf2(y=ld.wti,x=ld.drivers,d=TRUE,f=fcomp,fmod=m1)

# models just with one independent variable and a constant will be averaged
mds <- diag(1,ncol(ld.drivers),ncol(ld.drivers))
mds <- cbind(rep(1,ncol(ld.drivers)),mds)
a7 <- altf2(y=ld.wti,x=ld.drivers,d=TRUE,mods.incl=mds)
}
}

\references{
Burnham, K. P., Anderson, D. R., 2004. Multimodel inference: Understanding AIC and BIC in model selection. \emph{Sociological Methods & Research} \bold{33}, 261--304.

Gelman, A., Hwang, J., Vehtari, A., 2014. Understanding predictive information criteria for Bayesian models. \emph{Statistics and Computing} \bold{24}, 997--1016.

Kapetanios, G., Labhard, V., Price, S., 2008. Forecasting using Bayesian and information-theoretic model averaging. \emph{Journal of Business & Economic Statistics} \bold{26}, 33--41. 

Timmermann, A., 2006. Forecast combinations. In: Elliott, G., et al. (eds.), \emph{Handbook of Economic Forecasting}, Elsevier.
}

\seealso{\code{\link{plot.altf2}}, \code{\link{print.altf2}}, \code{\link{summary.altf2}}, \code{\link{rec.reg}}, \code{\link{roll.reg}}, \code{\link{tvp}}, \code{\link{altf}}, \code{\link{altf3}}, \code{\link{altf4}}.}
