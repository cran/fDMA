\name{normalize}
\alias{normalize}

\title{Normalizes a Numeric Matrix by Columns.}

\description{
For a variable considered to be used in Dynamic Model Averaging (or Dynamic Model Selection, etc.), sometimes it is desirable to have all its values between 0 and 1. This function rescales the values to fit between 0 and 1.

If the argument is not a matrix, the function tries to convert the object into a matrix. For example, it works smoothly for a \code{xts} objects.
}

\usage{
normalize(data)
}

\arguments{
\item{data}{a matrix, observations are put in rows, and variables are grouped by columns}
}

\value{a matrix}

\examples{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
normalize(wti)
normalize(drivers)
normalize(cbind(c(0,1,2),c(1,2,3),c(0,1,3)))
}
