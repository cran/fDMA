
summary.altf2 <- function(object, ...)
  {
   x <- object
   
   nmods <- length(x$coeff.)
   
   cat("Mean coefficients: ")
   cat("\n")
   
   c <- matrix(,ncol=ncol(x$coeff.[[1]]),nrow=nmods)
   for (i in 1:nmods)
    {
      c[i,] <- colMeans(x$coeff.[[i]],na.rm=TRUE)
    }
   colnames(c) <- colnames(x$coeff.[[1]])
   rownames(c) <- names(x$coeff.)
   c <- round(c,digits=4)
   print(c)
   cat("\n")
   
   if (! all(names(x$coeff.)==c("av. TVP")))
    {
      cat("Frequency when p-values for t-test are less than: ")
    }
   cat("\n")
   
   for (k in 1:nmods)
    {
      if (! names(x$p.val.)[k]=="av. TVP")
       {
         vv <- rep.int(NA,3)
         for (i in 1:ncol(x$p.val.[[1]]))
          {
             v <- as.vector(na.exclude(x$p.val.[[k]][,i,drop=FALSE]))
             v1 <- v[v<0.01]
             v2 <- v[v<0.05]
             v3 <- v[v<0.10]
             v1 <- length(v1) / length(v)
             v2 <- length(v2) / length(v)
             v3 <- length(v3) / length(v)
             v <- c(v1,v2,v3)
             vv <- rbind(vv,v)
          }
         vv <- vv[-1,,drop=FALSE]
         vv <- t(vv)
         vv <- round(vv,digits=2)
         colnames(vv) <- colnames(x$coeff.[[1]])
         rownames(vv) <- c("0.01","0.05","0.10")
         temp <- names(x$p.val)[k]
         names(temp) <- c("")
         print(temp)
         cat("\n")
         print(vv)
         cat("\n")
       }
    }
   
   cat("\n")
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary)
  }
