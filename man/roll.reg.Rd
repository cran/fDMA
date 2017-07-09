\name{roll.reg}
\alias{roll.reg}

\title{Computes Rolling Regression.}

\description{
This function computes Rolling Regression. For the first \code{window}-1 observations Recursive Regression is computed. Since \code{window}-th observation the rolling is performed. 
}

\usage{
roll.reg(y,x,window)
}

\arguments{
\item{y}{a numeric or a column matrix of a dependent variable}
\item{x}{a matrix of independent variables, different columns should correspond to different variables}
\item{window}{a numeric, a size of a window for rolling}
}

\value{a vector of fitted (forecasted) values}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
roll <- roll.reg(y=ld.wti,x=ld.drivers,window=100)
}
