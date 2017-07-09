
print.grid.dma <- function(x, ...)
  {
   cat("MSE: ")
   cat("\n")
   print(round(x$MSE,digits=4))
   cat("\n")
   cat("\n")
   cat("MAE: ")
   cat("\n")
   print(round(x$MAE,digits=4))
   cat("\n")
   cat("* alphas by columns, lambdas by rows")
  }
