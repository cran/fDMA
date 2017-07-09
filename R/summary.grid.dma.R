
summary.grid.dma <- function(object, ...)
  {
   x <- object
   
   cat("MSE: ")
   cat("\n")
   print(round(x$MSE,digits=4))
   cat("\n")
   cat("Indices of the model minimising MSE:")
   temp <- which(x$MSE == min(x$MSE), arr.ind = TRUE)[1,]
   names(temp) <- c("","")
   print(temp)
   cat("\n")
   cat("MAE: ")
   cat("\n")
   print(round(x$MAE,digits=4))
   cat("\n")
   cat("Indices of the model minimising MAE:")
   temp <- which(x$MAE == min(x$MAE), arr.ind = TRUE)[1,]
   names(temp) <- c("","")
   print(temp)   
   cat("\n")
   cat("* alphas by columns, lambdas by rows")
  }
