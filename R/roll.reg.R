
roll.reg <- function(y,x,window)
{

### estimates rolling regression model

### y - a numeric or a column matrix of a dependent variable

### x - a matrix of independent variables (drivers), different columns correspond to different variables

### window - a numeric, a size of a window for rolling

requireNamespace('xts')
requireNamespace('stats')

x <- as.matrix(x)

y.roll.ols <- vector()

for (i in 1:(window-1))
  {
    if (i==1) 
      { 
        y.roll.ols[1] <- lm(y[1] ~ t(x[1,]))$fitted.values[1]
      }
    else
      {
        y.roll.ols[i] <- lm(y[1:i] ~ x[1:i,])$fitted.values[i]
      }
  }

for (i in window:nrow(x))
  {
    y.roll.ols[i] <- lm(y[(i-window+1):i] ~ x[(i-window+1):i,])$fitted.values[window]
  }

return(y.roll.ols)

}
