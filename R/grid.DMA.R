
grid.DMA <- function (y,x,grid.alpha,grid.lambda,initvar,
                      W=NULL,initial.period=NULL,V.meth=NULL,kappa=NULL,gprob=NULL,omega=NULL,model=NULL,parallel.grid=NULL,
                      m.prior=NULL,mods.incl=NULL,DOW=NULL,DOW.nmods=NULL,DOW.type=NULL)
{

### this is a wrapper of fDMA() 

### "foreach", "doParallel", "stats", "forecast" and "xts" packages are required

### grid.alpha and grid.lambda must be numeric vectors, 
### fDMA for all combinations of alpha and lambda 
### possible to be obtained from the values in grid.alpha and grid.lambda respectively are computed

#############################################

if (missing(grid.alpha)) { stop("please, specify grid.alpha") }
if (! missing(grid.alpha) && ! is.numeric(grid.alpha)) { stop("grid.alpha must be a collection of numbers") }
if (anyNA(grid.alpha)) { stop("missing values in grid") }
for (i in 1:length(grid.alpha))
  {
    if ((grid.alpha[i] <= 0) || (grid.alpha[i] > 1)) { stop("grid.alpha must be a collection of numbers greater than 0, and less than or equal to 1") }
  }
rm(i)
if (missing(grid.lambda)) { stop("please, specify grid.lambda") }
if (! missing(grid.lambda) && ! is.numeric(grid.lambda)) { stop("grid.lambda must be a collection of numbers") }
if (anyNA(grid.lambda)) { stop("missing values in grid.lambda") }
for (i in 1:length(grid.lambda))
  {
    if ((grid.lambda[i] <= 0) || (grid.lambda[i] > 1)) { stop("grid.lambda must be a collection of numbers greater than 0, and less than or equal to 1") }
  }
rm(i)
if (missing(initvar)) { stop("please, specify initvar (i.e., initial variance)") }
if (! missing(initvar) && ! is.numeric(initvar)) { stop("initvar must be a number") }
if (initvar <= 0) { stop("variance (initvar) must be positive") }
if (is.null(parallel.grid)) { parallel.grid <- FALSE }
if (! is.logical(parallel.grid)) { stop("parallel.grid must be logical, i.e., TRUE or FALSE") }
if (requireNamespace('foreach')) 
  {
  } 
else 
  {
    stop("package >>foreach<< is required")
  }
if (parallel.grid == TRUE)
  {
    if (requireNamespace('doParallel')) 
      {
      } 
    else 
      {
        stop("for parallel computations package >>doParallel<< is required")
      }
  }

#############################################

grid.alpha <- unique(grid.alpha)
grid.lambda <- unique(grid.lambda)
grid.alpha <- sort(grid.alpha, decreasing=TRUE)
grid.lambda <- sort(grid.lambda, decreasing=TRUE)

i.alpha <- NULL
j.lambda <- NULL

if (parallel.grid == TRUE)
  {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    gDMA <- foreach(i.alpha=grid.alpha, .packages=c("xts","fDMA")) %:% 
        foreach(j.lambda=grid.lambda, .packages=c("xts","fDMA")) %dopar%
          {
            fDMA(y=y,x=x,alpha=i.alpha,lambda=j.lambda,initvar=initvar,W=W,initial.period=initial.period,V.meth=V.meth,kappa=kappa,
                gprob=gprob,omega=omega,model=model,parallel=FALSE,m.prior=m.prior,mods.incl=mods.incl,DOW=DOW,DOW.nmods=DOW.nmods,DOW.type=DOW.type)
          }
          
    stopCluster(cl)
  }
else
  {
    gDMA <- foreach(i.alpha=grid.alpha, .packages=c("xts","fDMA")) %:% 
        foreach(j.lambda=grid.lambda, .packages=c("xts","fDMA")) %do%
          {
            fDMA(y=y,x=x,alpha=i.alpha,lambda=j.lambda,initvar=initvar,W=W,initial.period=initial.period,V.meth=V.meth,kappa=kappa,
                gprob=gprob,omega=omega,model=model,parallel=FALSE,m.prior=m.prior,mods.incl=mods.incl,DOW=DOW,DOW.nmods=DOW.nmods,DOW.type=DOW.type)
          }
  }

mse <- matrix(NA,nrow=length(grid.lambda),ncol=length(grid.alpha))
mae <- matrix(NA,nrow=length(grid.lambda),ncol=length(grid.alpha))
rownames(mse) <- grid.lambda
colnames(mse) <- grid.alpha
rownames(mae) <- grid.lambda
colnames(mae) <- grid.alpha

for (i in 1:length(grid.alpha))
  {
    for (j in 1:length(grid.lambda))
      {
        mse[j,i] <- gDMA[[i]][[j]]$RMSE
        mae[j,i] <- gDMA[[i]][[j]]$MAE
      }
  }

temp <- list(gDMA,mse,mae)
names(temp) <- c("models","RMSE","MAE")
class(temp) <- "grid.dma"
return(temp)

}
