\name{fDMA}
\alias{fDMA}

\title{Computes Dynamic Model Averaging.}

\description{
The function estimates Dynamic Model Averaging (and some of its variations). The method is described in \cite{Raftery et al. (2010)}. 
}

\usage{
fDMA(y,x,alpha,lambda,initvar,W=NULL,initial.period=NULL,V.meth=NULL,kappa=NULL, 
gprob=NULL,omega=NULL,model=NULL,parallel=NULL,m.prior=NULL,mods.incl=NULL, 
DOW=NULL,DOW.nmods=NULL,DOW.type=NULL,bm=NULL)
}

\arguments{
\item{y}{a numeric or a column matrix of a dependent variable, if \code{y} is a xts object, then plots will have time index on the x axis}
\item{x}{a matrix of independent variables, different columns correspond to different variables}
\item{alpha}{a forgetting factor \eqn{\alpha} between 0 and 1 used in probabilities estimations}
\item{lambda}{a forgetting factor \eqn{\lambda} between 0 and 1 used in variance approximations}
\item{initvar}{initial variance in the state space equation, i.e., the number by which the unit matrix is multiplied}
\item{W}{optional, a method for setting the initial values of variance for the models equations, \code{W="reg"} corresponds to the method based on the linear regression as in the paper by \cite{Raftery et al. (2010)}, alternatively an arbitrary positive number can be specified, by default the method of \cite{Raftery et al. (2010)} is used}
\item{initial.period}{optional, a number of observation since which MSE (Mean Squared Error) and MAE (Mean Absolute Error) are computed, by default the whole sample is used, i.e., \code{initial.period=1}}
\item{V.meth}{optional, a method for the state space equation variance updating, \code{V.meth="rec"} corresponds to the recursive moment estimator, as in the paper by \cite{Raftery et al. (2010)}, \code{V.meth = "ewma"} corresponds to the exponentially weighted moving average as in, for example, \cite{Koop and Korobilis (2012)}, by default \code{V.meth = "rec"} is used}
\item{kappa}{optional, a parameter in the exponentially weighted moving average, between 0 and 1, used if \code{V.meth = "ewma"}}
\item{gprob}{optional, a matrix of Google probabilities as in \cite{Koop and Onorante (2014)}, columns should correspond to columns of x}
\item{omega}{optional, a parameter between 0 and 1 used in probabilities estimations, used if \code{gprob} is specified}
\item{model}{optional, \code{model="dma"} for Dynamic Model Averaging, \code{model="dms"} for Dynamic Model Selection, or \code{model="med"} for Median Probability Model as in \cite{Barbieri and Berger (2004)}, by default \code{model="dma"} is used}
\item{parallel}{optional, a logical, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{m.prior}{optional, a parameter for general model prior \cite{(Mitchell and Beauchamp, 1988)}, by default \code{m.prior=0.5}, which corresponds to the uniform distribution, i.e., non-informative priors, see also \cite{Eicher et al. (2011)}}
\item{mods.incl}{optional, a matrix indicating which models should be used for estimation, the first column indicates inclusion of a constant, by default all possible models with a constant are used, inclusion of a variable is indicated by 1, omitting by 0}
\item{DOW}{optional, a threshold for Dynamic Occam's Window \cite{(Onorante and Raftery, 2016)}, should be a number between 0 and 1, if \code{DOW=0}, then no Dynamic Occam's Window is applied, by default \code{DOW=0}, Dynamic Occam's Window can be applied only if all possible models with a constant are used}
\item{DOW.nmods}{optional, initial number of models for Dynamic Occam's Window, should be less than the number of all possible models and larger than or equal to 2, they are randomly chosen, if \code{DOW.nmods=0}, then initially models with exactly one variable are taken, by default \code{DOW.nmods=0}}
\item{DOW.type}{optional, \code{DOW.type="r"} corresponds to DMA-R from \cite{Onorante and Raftery (2016)}, \code{DOW.type="e"} to DMA-E, by default \code{DOW.type="r"}}
\item{bm}{optional, indicate whether benchmark forecast should be computed, these benchmarks are naive forecast (all forecasts are set to be the value of the last observation) and Auto Arima \code{\link[forecast]{auto.arima}}, by default \code{bm=FALSE}}
}

\value{class \code{dma} object
\item{$y.hat}{forecasted values}
\item{$post.incl}{posteriori inclusion probabilities for independent variables}
\item{$MSE}{Mean Squared Error of forecast}
\item{$MAE}{Mean Absolute Error of forecast} 
\item{$models}{models included in estimations}
\item{$post.mod}{posteriori probabilities of all used models}
\item{$exp.var}{expected number of variables (incl. constant)}
\item{$exp.coef.}{expected values of regression coefficients}
\item{$parameters}{parameters of the estimated model} 
\item{$yhat.all.mods}{predictions from every sub-model used in estimations} 
\item{$y}{independent variable}
\item{$benchmarks}{Mean Squared Error and Mean Absolute Error of naive and auto ARIMA forecast}
\item{$DOW.init.mods}{models initially selected to Dynamic Occam's Window, if this method has been selected}
}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,model="dms")
fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,V.meth="ewma",kappa=0.9)
fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,DOW=0.7)
}
\dontshow{
s <- ld.drivers[1:20,1:4]
fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10)
fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,model="dms")
fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,V.meth="ewma",kappa=0.9)
fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,DOW=0.7)
}
# Google trends are available since 2004
data(trends)
\dontrun{
gp <- trends/100
s <- ld.drivers['2004-01-01/']
fDMA(y=ld.wti['2004-01-01/'],x=s,alpha=0.99,lambda=0.90,initvar=10,gprob=gp,omega=0.5)
}
\dontshow{
s <- ld.drivers['2004-01-01/'][1:20,1:4]
gp <- trends/100
gp <- gp[1:20,1:4]
fDMA(y=ld.wti['2004-01-01/'][1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,gprob=gp,omega=0.5)
}
# models just with one independent variable and a constant will be averaged
mds <- diag(1,ncol(ld.drivers),ncol(ld.drivers))
mds <- cbind(rep(1,ncol(ld.drivers)),mds)
\dontrun{fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,mods.incl=mds)}
\dontshow{
s <- ld.drivers[1:20,1:4]
fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,mods.incl=mds[1:5,1:5])
}
}

\source{
Raftery, A. E., Karny, M., Ettler, P., 2010. Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.
}

\references{
Barbieri, M. M., Berger, J. O., 2004. Optimal predictive model selection. \emph{The Annals of Statistics} \bold{32}, 870--897.

Eicher, T. S., Papageorgiou, C., Raftery, A. E., 2011. Default priors and predictive performance in Bayesian Model Averaging, with application to growth determinants. \emph{Journal of Applied Econometrics} \bold{26}, 30--55.

Koop, G., Korobilis, D., 2012. Forecasting inflation using Dynamic Model Averaging. \emph{International Economic Review} \bold{53}, 867--886.

Koop, G., Onorante, L., 2014. Macroeconomic nowcasting using Google probabilities. \url{http://www.ecb.europa.eu/events/pdf/conferences/140407/OnoranteKoop_MacroeconomicNowcastingUsingGoogleProbabilities.pdf}.

Mitchell, T. J., Beauchamp, J. J., 1988. Bayesian variable selection in linear regression (with discussion). \emph{Journal of the American Statistical Association} \bold{83}, 1023--1036.

Onorante, L., Raftery, A. E., 2016. Dynamic model averaging in large model spaces using
dynamic Occam's window. \emph{European Economic Review} \bold{81}, 2--14.
}
