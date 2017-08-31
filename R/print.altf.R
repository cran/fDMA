
print.altf <- function(x, ...)
  {
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary)
  }
