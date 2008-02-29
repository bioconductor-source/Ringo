
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
findPeaksOnSmoothed <- function(smoothedX, probeAnno, thresholds, allChr=c(1:19,"X","Y"), distCutOff=600, minProbesInRow=3, cellType=NULL, checkUnique=TRUE, uniqueCodes=c(0), verbose=TRUE)
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
      if (checkUnique){
        chruni <- get(paste(chr,"unique",sep="."), env=probeAnno)
        stopifnot(length(chruni)==length(chridx))
        chridx <- chridx[chruni %in% uniqueCodes]
      } #  if (checkUnique)
      chrrm <- exprs(smoothedX)[chridx,i]
      chr.peaks <- peakByThreshold(chrmid, chrrm, threshold=thresholds[i], distCutOff=distCutOff, minProbesInRow=minProbesInRow)
      ## new version: return objects of class peak instead
      names(chridx) <- chrmid
      chr.peaks <- lapply(as.list(1:length(chr.peaks)), function(z){
        x <- chr.peaks[[z]];
        peakID <- paste(cellType, this.sample, paste("chr",chr,sep=""), paste("peak",z,sep=""),sep=".")
        peakID <- gsub("(^\\.+)|(\\.+$)","",peakID)#remove any leading and trailing dots
        return(newPeak(peakID, chr=chr, start=as.integer(names(x)[1]), end=as.integer(names(x)[length(x)]), cellType=cellType, antibody=this.sample, maxPeak=max(x), score=attr(x,"score"), probes=as.character(chridx[names(x)])))})
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
