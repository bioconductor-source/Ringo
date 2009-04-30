
### new class for quantiles over positions plots
setClass("qop", representation(data="list", genes="character", positions="numeric", samplenames="character", quantiles="numeric", mapping="character"), prototype=prototype(data=list(matrix(0, nrow=0, ncol=0))))

setValidity("qop", function(object){
  if (any(length(object@data)!=length(object@samplenames),
          !inherits(object@data[[1]], "matrix"),
          nrow(object@data[[1]])!=length(object@quantiles),
          ncol(object@data[[1]])!=length(object@positions)))
    {warning("Object does not fulfill requirements for class 'qop'.\n"); return(FALSE)}
  return(TRUE)
})# set validity of qop objects

setMethod("show",signature="qop", function(object){
  cat("Quantiles over relative positions for",
      length(object@data), "samples:\n")
  cat(object@samplenames)
  cat("\nPositions:\n", object@positions)
  cat("\nQuantiles:\n", object@quantiles)
  cat("\nbased on",length(object@genes),"genes.\n")
  invisible(NULL)
})

setMethod("plot",signature=c("qop","ANY"), function(x, y, xlab="Distance to feature start [bp]", ylab="Probe level [log2]", ylim, ...){
  if (missing(ylim))
    ylim <- range(do.call("rbind", x@data), na.rm=TRUE)
  plot(x=0,y=0, ylim=ylim, xlim=range(x@positions), xaxt="n", type="n", xlab=xlab, ylab=ylab,...)
  pos.kb <- sort(unique(round(x@positions/1000)))
  axis(side=1, at=pos.kb*1000,
       labels=paste(pos.kb,"kb",sep=""))
  if (missing(y))
    mycolors <- rainbow(length(x@data))
  else
    mycolors <- y
  for (i in 1:length(x@data)){
    for (j in 1:length(x@quantiles))
      lines(x=x@positions,y=x@data[[i]][j,], lwd=2, lty=j, col=mycolors[i])
  }
  legend(x="topleft", fill=mycolors, bty="n", 
         legend=x@samplenames)
  legend(x="topright", lty=seq(length(x@quantiles)), lwd=2, bty="n",
         legend=paste(round(x@quantiles*100),"% quantile", sep=""))
}) # plot method for qop objects


quantilesOverPositions <- function(xSet, selGenes, g2p,
                                   positions= seq(-5000, 10000, by=250),
                                   quantiles=c(0.1, 0.5, 0.9))
{
  stopifnot(inherits(xSet,"ExpressionSet"),
            is.list(g2p),
            all(selGenes %in% names(g2p)))
  ### get values per gene
  selGenesVal <- vector("list", length(selGenes))
  for (i in 1:length(selGenes)){
    if (i %% 1000 == 0) cat(i,"... ")
    ix <- g2p[[selGenes[i]]]
    if (length(ix)==0) next
    iy <- exprs(xSet)[names(ix),,drop=FALSE]
    if (sum(!is.na(iy[,1]))<2) next
    iy <- apply(iy, 2, function(a)
                approx(x=ix, y=a, xout=positions)$y)
    selGenesVal[[i]] <- iy
  }
  names(selGenesVal) <- selGenes

  ### get quantiles per sample
  valsBySample <- lapply(as.list(sampleNames(xSet)), function(thisSample){
    do.call("cbind", lapply(selGenesVal, function(theseValues){
      return(theseValues[,thisSample,drop=FALSE])}))})
  names(valsBySample) <- sampleNames(xSet)
  quantsBySample <- lapply(valsBySample, function(theseValues)
                           apply(theseValues, 1, quantile,
                                 probs=quantiles, na.rm=TRUE))
  
  ### convert to densities
  ## a. probe positions (for normalizing the smoothed signal)
  probePosDens <- density(unlist(g2p, use.names=FALSE))
  probePosY <- approx(x=probePosDens$x, y=probePosDens$y, 
                      xout=positions)$y
  probePosY <- probePosY/max(probePosY)
  ## b. apply normalization
  for (i in 1:length(quantsBySample)){
    quantsBySample[[i]] <- quantsBySample[[i]] * matrix(probePosY, ncol=length(probePosY), nrow=nrow(quantsBySample[[i]]), byrow=TRUE)
  }

  #  class(quantsBySample) <- c("qop", class(quantsBySample))
  # return(quantsBySample)
  res <- new("qop", data=quantsBySample, samplenames=sampleNames(xSet),
             positions=positions, genes=selGenes, quantiles=quantiles,
             mapping=deparse(substitute(g2p)))
  return(res)
} #quantilesOverPositions
