
upperBoundNull <- function(x, modeMethod="shorth", limits=c(-1,1), prob=0.99, ...){
  # function to compute an upper bound for a symmetric null distribution,
  #  given that the mode of the null distribution is close to zero
  #require("genefilter")
  stopifnot(is.numeric(x))
  x <- as.numeric(na.omit(x))
  modeMethod <- match.arg(modeMethod, c("shorth","half.range.mode","null"))
  xClose <- x[x >= limits[1] & x <= limits[2]]
  stopifnot(length(xClose)>0)
  xMode <- switch(modeMethod,
                  "shorth"=genefilter::shorth(xClose,...),
                  "half.range.mode"= genefilter::half.range.mode(xClose, ...),
                  "null"=0)
  nullSample <- x[x <= xMode]
  nullSample <- c(nullSample, xMode + (xMode - nullSample))
  nullUpperBound <- quantile(nullSample, probs=prob)
  names(nullUpperBound) <- NULL
  return(nullUpperBound)
}#upperBoundNull

