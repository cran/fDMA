
summary.dma <- function(object, ...)
  {
   
   x <- object
   
   E.coef <- round(colMeans(x$exp.coef.),digits=4)
   probs <- round(colMeans(x$post.incl),digits=2)
   min.probs <- round(apply(x$post.incl,2,min),digits=2)
   max.probs <- round(apply(x$post.incl,2,max),digits=2)
   
   inc <- vector()
   inc[1] <- round(1,digits=2)
   
   for (i in 1:(length(probs)-1))
    {
      inc[i+1] <- round(length(x$post.incl[x$post.incl[,i]>0.5,i]) / nrow(x$post.incl),digits=2)
    }
   rm(i)
   
   s <- cbind(E.coef,min.probs,probs,max.probs,inc)
   colnames(s) <- c("E(coeff.)","min pip", "mean pip", "max pip", "inc")
   rownames(s) <- colnames(x$models)
   if (! x$parameters[1,4] == "DMA") { s <- s[,c(1,5)] }
   
   cat("Model: ")
   cat(colnames(x$y))
   cat(" ~ ")
   cat("const ")
   cat(colnames(x$models))
   cat("\n")
   cat("\n")
     
   print(x$parameters)
   cat("\n")
   
   err <- rbind(x$MSE,x$MAE)
   colnames(err) <- c("model")
   err <- cbind(err,x$benchmarks)
   print(round(err,digits=4))

   cat("\n")
   cat("observations: ")
   cat(length(x$y.hat))
   cat("\n")
   cat("models: ")
   cat(nrow(x$models))
   cat("\n")
   cat("variables (incl. constant): ")
   cat(ncol(x$models))
   
   cat("\n")
   cat("\n")
   print(s)
   cat("\n")
   
   if (x$parameters[1,4] == "DMA") 
    { 
      cat("* pip = posteriori inclusion prob.")
      cat("\n")
      cat("* inc = relative frequency of a variable posteriori inclusion prob. > 1/2")
    }
   if (x$parameters[1,4] == "DMS" || x$parameters[1,4] == "MED")
    {
      cat("* inc = relative frequency of a variable inclusion")
    }
  }
