
print.altf4 <- function(x, ...)
  {
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary)
  }
