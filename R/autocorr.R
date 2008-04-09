# function to assess ChIP probe intensity autocorrelation

autocor <- function(x, probeAnno, chrom="1",samples=NULL, lag.max= 2000, lag.step=100, cor.method="pearson", channel=c("red","green","logratio"), verbose=TRUE)
{
  stopifnot(inherits(x,"ExpressionSet")|inherits(x,"RGList"),
            inherits(probeAnno, "probeAnno"), validObject(probeAnno))
  
  ## extract chromosomal information from probeAnno:
  probepos <- probeAnno[paste(chrom,"start",sep=".")]
  probeidx <- probeAnno[paste(chrom,"index",sep=".")]

  ## interpret arguments:
  lags=seq(0, lag.max, by=lag.step)
  match.winSize=round(lag.step/2)-1

  ## a. handle ExpressionSets (normalized data)
  if (inherits(x,"ExpressionSet")){
    if (is.null(samples))
      samples <- 1:ncol(exprs(x))
    dat <- rowMeans(exprs(x)[,samples,drop=FALSE])
  }

  ## b. hande RGLists (raw data)
  if (inherits(x,"RGList")){
    if (is.null(samples))
      samples <- 1:ncol(x$R)
    channel <- match.arg(channel, c("red","green","logratio"))
    dat <- switch(channel,
                  "red"=rowMeans(x$R[,samples,drop=FALSE]),
                  "green"=rowMeans(x$G[,samples,drop=FALSE]),
                  "logratio"=log2(rowMeans(x$R[,samples,drop=FALSE]))-log2(rowMeans(x$G[,samples,drop=FALSE])))
    names(dat) <- x$genes$'ID'
  }#if (inherits(x,"RGList"))
  
  if (verbose) cat("Lag: ")
  movedInt <- apply(matrix(lags),1, function(thismove){
    if (verbose) cat(thismove,"... ")
    movedpos   <- probepos + thismove
    matchedpos <- matchpt(movedpos, probepos)
    areInWindow <- matchedpos[,"distance"] <= match.winSize
    movedIntensities <- dat[probeidx[matchedpos[,"index"]]]
    movedIntensities[!areInWindow] <- NA
    return(movedIntensities)
  })#movedInt
  if (verbose) cat("\nComputing correlation...")
  ac <- cor(movedInt, use="pairwise.complete.obs", method=cor.method)[,1]
  names(ac) <- lags
  class(ac) <- "autocor.result"
  attr(ac, "chromosome") <- chrom
  return(ac)
}#autocor

#setOldClass("autocor.result")

plot.autocor.result <- function(x, plot.title="ChIP: Autocorrelation of Intensities", ...){
  stopifnot(inherits(x,"autocor.result"))
  ac <- x
  lags <- as.numeric(names(ac))
  chrom <- attr(ac, "chromosome")
  plot(x=lags,y=ac, type="h",lwd=5, xlab="Offset [bp]", ylab="Auto-Correlation", main=plot.title, sub=paste("Mean over Chromosome",chrom), col="darkblue",...)
  invisible(NULL)
}# plot.autocor.result

