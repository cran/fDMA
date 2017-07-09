\name{tvp}
\alias{tvp}

\title{Computes Time-Varying Parameters Regression.}

\description{
This function computes Time-varying Parameters Regression (TVP) with the updating procedure as in \cite{Raftery et. al (2010)}. 
}

\usage{
tvp(y,x,V,lambda,W=NULL)
}

\arguments{
\item{y}{a numeric or a column matrix of a dependent variable}
\item{x}{a matrix of independent variables, different columns should correspond to different variables}
\item{V}{initial variance in the state space equation for the recursive moment estimator updating method, as in \cite{Raftery et al. (2010)}}
\item{lambda}{a forgetting factor between 0 and 1 used in variance approximations}
\item{W}{optional, initial value of variance for the model equations, if not specified the method based on the linear regression, as in \cite{Raftery et al. (2010)} is used}
}

\value{a vector of fitted (forecasted) values}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99)}
\dontshow{tvp(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],V=1,lambda=0.99)}
\dontrun{tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99,W=1)}
\dontshow{tvp(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],V=1,lambda=0.99,W=1)}
}

\references{
Raftery, A. E., Karny, M., Ettler, P., 2010. Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.
}
