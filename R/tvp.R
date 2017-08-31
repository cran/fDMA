
tvp <- function(y,x,V,lambda,W=NULL,kappa=NULL,c=NULL)
{

### estimates time-varying parameters regression model

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different variables

### V - initial variance in the state space equation for the recursive moment estimator updating method, 
###     as in the paper by Raftery et al. (2010),

### lambda - a forgetting factor between 0 and 1 used in variance approximations

### W - a method for setting the initial values of variance for the models equations,
###     by default (if W is not specified) the method based on the linear regression,
###     as in the paper by Raftery et al. (2010) is used,
###     alternatively an arbitrary positive number can be specified

### kappa - a parameter in the exponentially weighted moving average, between 0 and 1,
###         if not specified, the method from Raftery et al. (2010) is applied

### c - a parameter indicating whether constant is included,
###     by default c=TRUE (constant is included),
###     it is not possible to set c=FALSE if ncol(x)=0 


requireNamespace('xts')
requireNamespace('stats')

if ( is.null(c) ) { c <- TRUE }

x <- as.matrix(x)
if ( ncol(x) == 0 ) { c <- TRUE }

xe <- cbind(rep.int(1,nrow(x)),x)


if ( c == TRUE )
  {
    theta <- matrix(0,ncol=1, nrow=ncol(xe))
  }
else
  {
    theta <- matrix(0,ncol=1, nrow=ncol(xe)-1)
  }
thetas <- theta

if (is.null(W))
{
  if (ncol(x)>0)
    {
      beta <- lm(y~x)$coefficients[1]
    }
  else
    {
      beta <- mean(y)
    }
  if (length(y)>1)
    {
      beta <- as.numeric(beta * beta + var(as.numeric(y)))
    }
  else
    {
      beta <- as.numeric(beta * beta)
    }

  E <- diag(beta,ncol=1+ncol(x), nrow=1+ncol(x))

  v <- vector()
  if (ncol(x)>0 && length(y)>1)
    {
      for (i in 1:ncol(x))
        {
          v[i] <- var(x[,i])
        }
    }
  if (length(y)>1)
    {
      vary <- as.numeric(var(as.numeric(y)))
    }
  
  if (ncol(x)>0 && length(y)>1) 
    {
      for (j in 1:length(v))
        {
          E[j+1,j+1] <- vary/v[j]
        }
    }
  else
    {
      if (ncol(x)>0)
        {
          for (j in 1:ncol(x))
            {
              E[j+1,j+1] <- 0
            }
        }
    }
}

if (!is.null(W))
{
  E <- diag(W,ncol=1+ncol(x), nrow=1+ncol(x))
}

y.tvp <- vector()
pdensi <- vector()

#########################
#########################

if (c == FALSE) 
  { 
    xe <- xe[,-1,drop=FALSE] 
    E <- E[-1,-1]
  }

for (t in 1:nrow(x))
{
  xx <- xe[t,]
  yhat <- as.numeric(crossprod(xx,theta))
  ei <- as.numeric(as.vector(y)[t] - yhat)
  R <- E / lambda
  tv <- as.numeric(crossprod(xx,R) %*% xx)
  Vu <- as.numeric(V + tv)
  E <- R - (R %*% xx) %*%  crossprod(xx,R) / Vu  
  pdensi[t] <- exp(-0.5 * ei * ei / Vu ) / sqrt(2*pi*Vu)
  if (is.null(kappa))
    {
      temp <- ( (t-1) * V + (ei * ei - tv) ) / t 
      if (temp>0) 
        { 
          theta <- theta + ( R %*% xx ) * ei / (temp + tv) 
          V <- temp
        }
      else 
        { 
          theta <- theta + ( R %*% xx ) * ei / Vu 
        }
     }
  else
     {
        temp <- V * kappa + (1-kappa) * ei * ei 
        theta <- theta + ( R %*% xx ) * ei / (temp + tv) 
     }
  y.tvp[t] <- yhat
  thetas <- cbind(thetas,theta)
}

rm(xx,xe,yhat,ei,R,E,tv,Vu,temp,t,V,theta)

thetas <- t(thetas[,-1])
if (!is.null(colnames(x))) 
  { 
    if ( c == TRUE ) 
      {
        colnames(thetas) <- c("const",colnames(x)) 
      }
    else
      {
        if (ncol(x)==1) { thetas <- t(thetas) }
        colnames(thetas) <- colnames(x)
      }
  }
if (ncol(x)==0)
  { 
    thetas <- t(thetas) 
    colnames(thetas) <- "const"
  }

ret <- list(y.tvp,thetas,pdensi,as.matrix(y))
names(ret) <- c("y.hat","thetas","pred.dens.","y")
class(ret) <- "tvp"

return(ret)

}
