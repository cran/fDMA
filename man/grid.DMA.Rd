\name{grid.DMA}
\alias{grid.DMA}

\title{Computes \link{fDMA} Function for Multiple Values of \code{alpha} and \code{lambda}.}

\description{
Sometimes it is necessary to consider various values of parameters \eqn{\alpha} and \eqn{\lambda} in Dynamic Model Averaging (or Dynamic Model Selection, etc.). This function computes fDMA function for all combinations of \code{alpha} and \code{lambda} for given grids.

This function is a wrapper of \link{fDMA}.
}

\usage{
grid.DMA(y,x,grid.alpha,grid.lambda,initvar,W=NULL,initial.period=NULL,V.meth=NULL,
kappa=NULL,gprob=NULL,omega=NULL,model=NULL,parallel.grid=NULL,m.prior=NULL,
mods.incl=NULL,DOW=NULL,DOW.nmods=NULL,DOW.type=NULL)
}

\arguments{
\item{y}{see \link{fDMA}}
\item{x}{see \link{fDMA}}
\item{grid.alpha}{a numeric vector}
\item{grid.lambda}{a numeric vector}
\item{initvar}{see \link{fDMA}}
\item{W}{see \link{fDMA}}
\item{initial.period}{see \link{fDMA}}
\item{V.meth}{see \link{fDMA}}
\item{kappa}{see \link{fDMA}}
\item{gprob}{see \link{fDMA}}
\item{omega}{see \link{fDMA}}
\item{model}{see \link{fDMA}}
\item{parallel.grid}{optional, a logical, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{m.prior}{see \link{fDMA}}
\item{mods.incl}{see \link{fDMA}}
\item{DOW}{see \link{fDMA}}
\item{DOW.nmods}{see \link{fDMA}}
\item{DOW.type}{see \link{fDMA}}
}

\value{an object of class \code{grid.dma}

\item{$models}{a list of list of models}
\item{$MSE}{a matrix with Mean Squared Error (MSE) for all estimated models}
\item{$MAE}{a matrix with Mean Absolute Error (MAE) for all estimated models}
}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)
g <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=10)
}
\dontshow{
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)
g <- grid.DMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],grid.alpha=gra,grid.lambda=grl,initvar=10)
}
# extract model with alpha=0.97 and lambda=0.95
g$models[[3]][[2]]
}
