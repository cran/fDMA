\name{plot.dma}
\alias{plot.dma}

\title{Plots Selected Outcomes from \link{fDMA} Function.}

\description{The function plots selected outcomes from \link{fDMA}.}

\usage{\method{plot}{dma}(x, ...)}

\arguments{
\item{x}{an object of \code{dma} class}
\item{...}{not used}
}

\details{
If \code{x} comes from estimation of Dynamic Model Averaging (DMA), after executing the command, the User is asked to choose

1 - for plotting actual vs. predicted values,

2 - for plotting residuals,

3 - for plotting the expected number of variables (including constant),

4 - for plotting posteriori inclusion probabilities (including constant) on one plot,

5 - for plotting posteriori inclusion probabilities (including constant) in separate \code{png} files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a \code{png} file in the current working directory),

6 - for plotting expected coefficients (including constant) on one plot,

7 - for plotting expected coefficients (including constant) in separate \code{png} files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a \code{png} file in the current working directory).

Chosing 0 exits the \code{plot} command.

---------------------------------------

If {x} comes from estimation of Dynamic Model Selection (DMS) or Median Probability Model (MED), after executing \code{plot} the User is asked to choose

1 - for plotting actual vs. predicted values,

2 - for plotting residuals,

3 - for plotting the expected number of variables (including constant),

4 - for producing a plot showing which variables (including constant) are included in the DMS or MED model in each time,

5 - for plotting expected coefficients (including constant) on one plot,

6 - for plotting expected coefficients (including constant) in separate \code{png} files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a \code{png} file in the current working directory).

Chosing 0 exits the \code{plot} command.
}

\note{
It is suggested to execute \code{\link{graphics.off}} before exectuing \code{plot} command for \code{dma} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link{graphics.off}} is not executed before plotting \code{dma} object, sometimes a legend might cover the important parts of the plot.
}

\examples{
\dontrun{
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dma")
m2 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dms")
# graphics.off()
plot(m1)
# graphics.off()
plot(m2)
}
}
