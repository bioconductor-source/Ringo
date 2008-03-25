

cherByThreshold <- function(positions, scores, threshold, distCutOff, minProbesInRow=3){
  # threshold: threshold for scores to be called cher
  # distCutOff: maximal distance between two positions before combining
  stopifnot(identical(length(positions), length(scores)), !is.unsorted(positions), is.numeric(minProbesInRow), minProbesInRow>0)
  highScores <- scores >= threshold
  highScores[is.na(highScores)] <- FALSE
  highScorePairs <- highScores[-length(highScores)] * highScores[-1]
  closePositions <- as.numeric(diff(positions)<=distCutOff)  
  cherPairs <- highScorePairs * closePositions
  cherClusters <- clusters(cherPairs, minLen=(minProbesInRow-1), doSelect=TRUE)
  nClusters <- nrow(cherClusters)
  if (nClusters==0) return(vector("list",0))
  res <- as.list(rep(NA, nClusters))
  cherMedScores <- rep(NA, nClusters)
  for (i in 1:nClusters){
    clustStart <- cherClusters[i,1]
    clustLen  <- cherClusters[i,2]
    clustScores <- scores[clustStart+(0:clustLen)]
    names(clustScores) <- positions[clustStart+(0:clustLen)]
    res[[i]]  <- clustScores
    ## cher score: "area" under the curve
    attr(res[[i]],"score") <- sum(clustScores - threshold)
    #cherMedScores[i] <- median(scores[clustStart+(0:clustLen)])
  }#  for (i in 1:nClusters)
  names(res) <- c(paste("cher",1:nClusters,sep=""))#,"cherMedianScores")
  return(res)
}#cherByThreshold


### main cher-finding funtion at the moment
findChersOnSmoothed <- function(smoothedX, probeAnno, thresholds, allChr=c(1:19,"X","Y"), distCutOff=600, minProbesInRow=3, cellType=NULL, checkUnique=TRUE, uniqueCodes=c(0), verbose=TRUE)
{
  stopifnot(is.numeric(thresholds), length(thresholds)==ncol(smoothedX))
  resultChers <- vector("list",ncol(smoothedX))
  for (i in 1:ncol(smoothedX)){
    this.sample <- sampleNames(smoothedX)[i]
    if (verbose) cat("\n\nSample: ",this.sample,"...\n\nChr: ")
    thisModChers <- lapply(as.list(allChr), function(chr){
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
        chrmid <- chrmid[chruni %in% uniqueCodes]
      } #  if (checkUnique)
      chrrm <- exprs(smoothedX)[chridx,i]
      chr.chers <- cherByThreshold(chrmid, chrrm, threshold=thresholds[i], distCutOff=distCutOff, minProbesInRow=minProbesInRow)
      if (length(chr.chers)==0) return(list())
      ## new version: return objects of class cher instead
      names(chridx) <- chrmid
      chr.chers <- lapply(as.list(1:length(chr.chers)), function(z){
        x <- chr.chers[[z]];
        cherID <- paste(cellType, this.sample, paste("chr",chr,sep=""), paste("cher",z,sep=""),sep=".")
        cherID <- gsub("(^\\.+)|(\\.+$)","",cherID)#remove any leading and trailing dots
        return(newCher(cherID, chr=chr, start=as.integer(names(x)[1]), end=as.integer(names(x)[length(x)]), cellType=cellType, antibody=this.sample, maxCher=max(x), score=attr(x,"score"), probes=as.character(chridx[names(x)])))})
      return(chr.chers)
    })
    #names(thisModChers) <- allChr
    resultChers[[i]] <- thisModChers
  }#for
  #names(resultChers) <- sampleNames(smoothedX)
  resultChers <- unlist(resultChers, recursive=FALSE, use.names=FALSE)
  resultChers <- unlist(resultChers, recursive=FALSE)
  class(resultChers) <- "cherList"
  return(resultChers)
}# findChersOnSmoothed
