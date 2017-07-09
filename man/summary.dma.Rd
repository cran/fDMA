\name{summary.dma}
\alias{summary.dma}

\title{Summarizes Outcomes from \code{dma} Object.}

\description{The function summarizes outcomes obtained from \link{fDMA}.}

\usage{\method{summary}{dma}(object, ...)}

\arguments{
\item{object}{an object of \code{dma} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{print}.

Additionally:

If \code{object} comes from Dynamic Model Averaging (DMA), it shows how often (in comparision to the whole analyzed period) a posteriori inclusion probability for a given variable exceeds 1/2. It also shows minimum, maximum and mean posteriori inclusion probability for every variable throught the analyzed period. 

If \code{object} comes from Dynamic Model Selection (DMS) or Median Probability Model (MED), it shows how often (in comparision to the whole analyzed period) a given variable is present in the selected model. 
}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dma")
m2 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dms")
}
\dontshow{
m1 <- fDMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],alpha=0.99,lambda=0.99,initvar=1,model="dma")
m2 <- fDMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],alpha=0.99,lambda=0.99,initvar=1,model="dms")
}
summary(m1)
summary(m2)
}
