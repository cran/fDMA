\name{crudeoil}
\alias{crudeoil}

\title{Crude Oil Data.}

\description{Selected data from oil market}

\usage{data(crudeoil)}

\format{\code{crudeoil} is a \code{xts} objects such that
\itemize{
\item \code{crudeoil$WTI} -- WTI spot price in USD per barrel
\item \code{crudeoil$MSCI} -- MSCI World Index
\item \code{crudeoil$TB3MS} -- U.S. 3-month treasury bill secondary market rate in \%
\item \code{crudeoil$CSP} -- Crude steel production in thousand tonnes
\item \code{crudeoil$TWEXM} -- Trade weighted U.S. dollar index (Mar, 1973 = 100)
\item \code{crudeoil$PROD} -- U.S. product supplied for crude oil and petroleum products in thousands of barrels
\item \code{crudeoil$CONS} -- Total consumption of petroleum products in OECD in quad BTU
\item \code{crudeoil$VXO} -- Implied volatility of S&P 100
}
}

\details{The data are in monthly frequency. They cover the period between Jan, 1990 and Dec, 2016.
}

\source{The data are provided by CBOE, Federal Reserve Bank of St. Louis, MSCI, U.S. Energy Information Administration and World Steel Association. 

\url{http://www.cboe.com/products/vix-index-volatility/vix-options-and-futures/vix-index/vix-historical-
data}

http://www.eia.gov

\url{http://fred.stlouisfed.org}

\url{http://www.msci.com/end-of-day-data-search}

\url{http://www.worldsteel.org/steel-by-topic/statistics/Statistics-monthly-crude-steel-and-iron-data-
/steel-archive.html}
}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{m <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dma")}
}
