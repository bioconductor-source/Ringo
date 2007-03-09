
peakByThreshold <- function(positions, scores, threshold, distCutOff, minProbesInRow=3){
  # threshold: threshold for scores to be called peak
  # distCutOff: maximal distance between two positions before combining
  stopifnot(identical(length(positions), length(scores)), !is.unsorted(positions), is.numeric(minProbesInRow), minProbesInRow>0)
  highScores <- scores >= threshold
  highScores[is.na(highScores)] <- FALSE
  highScorePairs <- highScores[-length(highScores)] * highScores[-1]
  closePositions <- as.numeric(diff(positions)<=distCutOff)  
  peakPairs <- highScorePairs * closePositions
  peakClusters <- clusters(peakPairs, minLen=(minProbesInRow-1), doSelect=TRUE)
  nClusters <- nrow(peakClusters)
  if (nClusters==0) return(vector("list",0))
  res <- as.list(rep(NA, nClusters))
  peakMedScores <- rep(NA, nClusters)
  for (i in 1:nClusters){
    clustStart <- peakClusters[i,1]
    clustLen  <- peakClusters[i,2]
    clustScores <- scores[clustStart+(0:clustLen)]
    names(clustScores) <- positions[clustStart+(0:clustLen)]
    res[[i]]  <- clustScores
    ## peak score: "area" under the curve
    attr(res[[i]],"score") <- sum(clustScores - threshold)
    #peakMedScores[i] <- median(scores[clustStart+(0:clustLen)])
  }#  for (i in 1:nClusters)
  names(res) <- c(paste("Peak",1:nClusters,sep=""))#,"PeakMedianScores")
  return(res)
}#peakByThreshold


### main peak-finding funtion at the moment
findPeaksOnSmoothed <- function(smoothedX, probeAnno, thresholds, allChr=c(1:19,"X","Y"), distCutOff=600, minProbesInRow=3, cellType=NULL, verbose=TRUE)
{
  stopifnot(is.numeric(thresholds), length(thresholds)==ncol(smoothedX))
  resultPeaks <- vector("list",ncol(smoothedX))
  for (i in 1:ncol(smoothedX)){
    this.sample <- sampleNames(smoothedX)[i]
    if (verbose) cat("\n\nSample: ",this.sample,"...\n\nChr: ")
    thisModPeaks <- lapply(as.list(allChr), function(chr){
      if (verbose) cat(chr, "...")
      chrsta <- get(paste(chr,"start",sep="."), env=probeAnno)
      chrend <- get(paste(chr,"end",sep="."), env=probeAnno)
      stopifnot(all(chrend>chrsta))
      chrmid <- round((chrsta+chrend)/2)
      chridx <- get(paste(chr,"index",sep="."), env=probeAnno)
      chrrm <- exprs(smoothedX)[chridx,i]
      chr.peaks <- peakByThreshold(chrmid, chrrm, threshold=thresholds[i], distCutOff=distCutOff, minProbesInRow=minProbesInRow)
      ## new version: return objects of class peak instead
      names(chridx) <- chrmid
      chr.peaks <- lapply(as.list(1:length(chr.peaks)), function(z){
        x <- chr.peaks[[z]];
        peakID <- paste(cellType, this.sample, paste("chr",chr,sep=""), paste("peak",z,sep=""),sep=".")
        peakID <- gsub("(^\\.+)|(\\.+$)","",peakID)#remove any leading and trailing dots
        return(newPeak(peakID, chr=chr, start=as.integer(names(x)[1]), end=as.integer(names(x)[length(x)]), cellType=cellType, modification=this.sample, maxPeak=max(x), score=attr(x,"score"), probes=as.character(chridx[names(x)])))})
      return(chr.peaks)
    })
    #names(thisModPeaks) <- allChr
    resultPeaks[[i]] <- thisModPeaks
  }#for
  #names(resultPeaks) <- sampleNames(smoothedX)
  resultPeaks <- unlist(resultPeaks, recursive=FALSE, use.names=FALSE)
  resultPeaks <- unlist(resultPeaks, recursive=FALSE)
  class(resultPeaks) <- "peakList"
  return(resultPeaks)
}# findPeaksOnSmoothed


### Obsolote, remnant functions, maybe revived some day:

## find peaks at points where intensities most follow a Gaussian
### distribution around the genomic position
## obsolete function at the moment
getPeaks <- function(eSet, chrom, probeAnno, samples=1, est.sd=1000, condenseDist=2000, byStrand = FALSE, verbose=TRUE, ...)
{
  # 0. check arguments  
  stopifnot((inherits(eSet,"exprSet")|inherits(eSet,"ExpressionSet")))
  eSetProbeNames <- switch(class(eSet),
                           "ExpressionSet"=featureNames(eSet),
                           "exprSet"=geneNames(eSet))
  if (is.null(eSetProbeNames))
    stop("Could not determine probe identifiers from expression set.\nCheck 'featureNames' or 'geneNames' of expression set.\n")
  if (is.null(samples)) samples <- 1:ncol(exprs(eSet))
  thisCall <- match.call()

  # 1. get intensities in selected region
  if (verbose) cat("Getting probe intensities on chromosome...\n")
  sta = get(paste(chrom, "start", sep="."), envir=probeAnno)
  end = get(paste(chrom, "end",   sep="."), envir=probeAnno)
  mid <- round((sta+end)/2)
  ind = get(paste(chrom, "index", sep="."), envir=probeAnno)
  stopifnot(all(ind %in% eSetProbeNames))
  names(mid) <- ind
  nProbesOnChr <- length(mid)
  mid <- mid[order(mid)]
  usedProbes <- mid
  nProbes <- length(usedProbes)
  nSamples <- length(samples)
  usedProbesIdx <- match(names(usedProbes),eSetProbeNames)
  chromExprs <- exprs(eSet)[usedProbesIdx, samples, drop=FALSE]

  # 2. running compute residuals
  if (verbose) cat("Computing residuals...\n")
  probeResids <- sapply(1:nProbes, function(i) compResiduals(x1=usedProbes[i], y1=chromExprs[i], datx=usedProbes, daty=chromExprs, est.sd=est.sd))
  probeResids <- round(probeResids,digits=4)
  ordResids <- order(probeResids)

  ordPos <- usedProbes[ordResids]
  orderedResids <- probeResids[ordResids]
  # 3. condense results:
  if (verbose) cat("Condensing results...\n")
  ## a. compute distance between consecutive probes
  ordPosDiff <- abs(diff(ordPos))
  ordPosBinary <- c(as.numeric(ordPosDiff <= condenseDist),1)
  # at this stage isolated peak probes get dropped,
  #  since only at least two near probes can have a 1 in their distances, so build one cluster there
  ordPosClusters <- clusters(ordPosBinary)  
  nClusters <- nrow(ordPosClusters)
  # initialize result
  resStart <- integer(nClusters);resEnd <- integer(nClusters);resResid <- numeric(nClusters)
  for (j in 1:nClusters){
    clustStart <- ordPosClusters[j,1]
    clustLen <- ordPosClusters[j,2]
    clustPosRange <- range(ordPos[clustStart+(0:clustLen)])
    resStart[j] <- clustPosRange[1]
    resEnd[j] <- clustPosRange[2]
    resResid[j] <- mean(orderedResids[clustStart+(0:(clustLen-1))])
  }#for (j in 1:nClusters)

  # need to condense the final result again:
  ## drop clusters if approximately the same cluster is found earlier in the list
  resMid <- as.integer((resStart + resEnd) / 2)
  resMatchpt <- t(sapply(seq(nClusters), function(z) matchpt(resMid[z],resMid[1:(z-1)])))
  #resMatchpt <- matchpt(resMid)
  keepClusters <- (!(resMatchpt[,1] < seq(nClusters) & resMatchpt[,2]<=condenseDist))
  result <- data.frame(start=I(resStart[keepClusters]),end=I(resEnd[keepClusters]),residual=I(resResid[keepClusters]))
  return(result)
} # getPeaks

### Auxiliary functions:
compResiduals <- function(x1, y1, datx, daty, est.sd=4){
  # compute residuals from gaussian around point
  est.y   <- dnorm(datx, mean=x1, sd=est.sd)
  est.mod <- y1/dnorm(x1, mean=x1, sd=est.sd)
  est.y <- est.y * est.mod
  resid <- sum((est.y - daty)^2)
  return(resid)
}#compResiduals
