
tvp <- function(y,x,V,lambda,W=NULL)
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

requireNamespace('xts')
requireNamespace('stats')

x <- as.matrix(x)

xe <- cbind(rep.int(1,nrow(x)),x)

theta <- matrix(0,ncol=1, nrow=ncol(xe))

if (is.null(W))
{
  beta <- lm(y~xe)$coefficients[1]
  beta <- as.numeric(beta * beta + var(as.numeric(y)))

  v <- vector()
  for (i in 1:ncol(x))
    {
      v[i] <- var(x[,i])
    }
  vary <- var(as.numeric(y))
  E <- diag(beta,ncol=1+length(v), nrow=1+length(v))
  if (length(v)>1) 
    {
      for (j in 2:length(v))
        {
          E[j,j] <- vary/v[j-1]
        }
    }

  rm(v,beta,vary)
}

if (!is.null(W))
{
  E <- diag(W,ncol=1+ncol(x), nrow=1+ncol(x))
}

y.tvp <- vector()
  
#########################
#########################

for (t in 1:nrow(x))
{
  xx <- xe[t,]
  yhat <- as.numeric(crossprod(xx,theta))
  ei <- as.numeric(as.vector(y)[t] - yhat)
  R <- E / lambda
  tv <- as.numeric(crossprod(xx,R) %*% xx)
  Vu <- as.numeric(V + tv)
  E <- R - (R %*% xx) %*%  crossprod(xx,R) / Vu  
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
  y.tvp[t] <- yhat
}

rm(xx,xe,yhat,ei,R,E,tv,Vu,temp,t,V,theta)

return(y.tvp)

}
