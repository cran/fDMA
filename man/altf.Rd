\name{altf}
\alias{altf}

\title{Computes a Few Alternative Forecasts.}

\description{
It is necessary to compare a given forecast method with some alternative ones. This function computes selected forecast quality measures for a few selected forecast methods (which might be treated as alternative ones to Dynamic Model Averaging, Dynamic Model Selection, etc.). 

Naive forecast (naive) is computed in a way that all forecasts are set to be the value of the last observation.

For rolling OLS forecast (roll. OLS) for the first periods (until the size of a \code{window} is obtained) are estimated through recursive OLS (rec. OLS). 

Autoregressive models (AR(1) and AR(2)) are computed by ordinary least squares method. 

Time-varying parameters models (TVP, TVP-AR(1) and TVP-AR(2)) are computed as \code{\link{tvp}} with \code{V = 1} and \code{lambda = 0.99}.

Auto ARIMA (auto ARIMA) is computed as \code{\link[forecast]{auto.arima}}.

Markov Switching Models (MS) are computed as \code{\link[MSwM]{msmFit}}. The object is a simple linear regression model. 2 regimes are assumed. All coefficients have switching. 

ME (Mean Error), RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), MPE (Mean Percentage Errror) and MAPE (Mean Absolute Percentage Error) are computed as \code{\link[forecast]{accuracy}}. HR (Hit Ratio) is computed as \code{\link{hit.ratio}}.
}

\usage{
altf(y,x,window=NULL,initial.period=NULL,d=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns correspond to different independent variables}
\item{window}{optional, \code{\link[base]{numeric}}, a size of a rolling regression window (a number of observations), if not specified 10\% of all observations are taken}
\item{initial.period}{optional, \code{\link[base]{numeric}}, a number of observation since which forecast quality measures are computed, if not specified the whole sample is used, i.e., \code{initial.period=1}}
\item{d}{optional, \code{\link[base]{logical}}, a parameter used for HR (Hit Ratio) calculation, should be \code{d=FALSE} for level time-series and \code{d=TRUE} if time-series represent changes, if not specified \code{d=FALSE}}
}

\value{\code{\link[base]{matrix}}, 
forecast quality measures are ordered by columns, 
forecast methods are ordered by rows, i.e., alternative forecasts are ordered by rows}

\examples{
\dontrun{
# models for untransformed data
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a1 <- altf(y=wti,x=drivers)

# do not include first 12 observations for forecast quality measures,
# i.e., treat first 12 observations as a ''training set''
a2 <- altf(y=wti,x=drivers,initial.period=12)

# models for log-differenced data
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
a3 <- altf(y=ld.wti,x=ld.drivers,d=TRUE)
}
}
