
gNormalize <- function(data)
  {
    ### data - a matrix or xts object of, for example, Google Trends

    if (missing(data)) { stop("please, specify data") }
    if (! is.matrix(data)) 
      { 
        data <- as.matrix(data)
        warning("data should be a matrix, the function tried to convert data to a matrix") 
      }

    x <- data
    for (i in 1:nrow(data))
      {
        x[i,] <- data[i,] / sum(data[i,])
      }
    rownames(x) <- rownames(data)
    colnames(x) <- colnames(data)
    return(x)
  }
  