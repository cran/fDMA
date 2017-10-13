\name{plot.altf4}
\alias{plot.altf4}

\title{Plots Selected Outcomes from \code{altf4} Object.}

\description{The function plots selected outcomes from \code{altf4} object.}

\usage{\method{plot}{altf4}(x, ...)}

\arguments{
\item{x}{an object of \code{altf4} class}
\item{...}{not used}
}

\details{
After executing the command, the User is asked to choose

1 - for plotting expected coefficients in separate png files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a png file in the current working directory),

2 - for plotting weights of all models used in averaging,

3 - for plotting expected window size.


Chosing 0 exits the \code{plot} command.
}

\note{
It is suggested to execute \code{\link[grDevices]{graphics.off}} before exectuing \code{plot} command for \code{altf4} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link[grDevices]{graphics.off}} is not executed before plotting \code{altf4} object, sometimes a legend might cover the important parts of the plot.
}

\examples{
\dontrun{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf4(y=wti,x=drivers,windows=c(36,100,150))

plot(a)
}
}
