

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
findChersOnSmoothed <- function(smoothedX, probeAnno, thresholds, allChr=NULL, distCutOff=600, minProbesInRow=3, cellType=NULL, antibodyColumn=NULL, checkUnique=TRUE, uniqueCodes=c(0), verbose=TRUE)
{
  stopifnot(is.numeric(thresholds), length(thresholds)==ncol(smoothedX),
            validObject(probeAnno), inherits(smoothedX,"ExpressionSet"))
  ## validate character vector of selected chromosome names and set default
  if (is.null(allChr))
      allChr <- chromosomeNames(probeAnno)
  else
      stopifnot(is.character(allChr),
                allChr %in% chromosomeNames(probeAnno))
  # look at the cellType definition
  if (!is.null(cellType)){
    stopifnot(is.character(cellType))
    if (length(cellType)==1){
      if (cellType %in% names(pData(smoothedX)))
        allCellTypes <- pData(smoothedX)[[cellType]]
      else
        allCellTypes <- rep(cellType, ncol(smoothedX))
    } else {
      if (length(cellType)!=ncol(smoothedX))
        stop("Argument 'cellType' must either be of length one or of length equal to the number of samples in the supplied ExpressionSet.\n")
      allCellTypes <- cellType
    }
  } else {
    allCellTypes <- vector("character", ncol(smoothedX))
  }
  if (!is.null(antibodyColumn)){
      antibodies <- as.character(pData(smoothedX)[[antibodyColumn]])
      if (is.null(antibodies))
          warning(paste("Column",antibodyColumn,"not defined not in the phenoData of object", deparse(substitute(xdf)),"."))
          antibodies <- sampleNames(smoothedX)
  } else {
      antibodies <- sampleNames(smoothedX)
  }
  resultChers <- vector("list",ncol(smoothedX))
  for (i in 1:ncol(smoothedX)){
    this.sample <- sampleNames(smoothedX)[i]
    thisCellType <- allCellTypes[i]
    if (verbose) cat("\n\nSample: ",this.sample,"...\n\nChr: ")
    thisModChers <- lapply(as.list(allChr), function(chr){
      if (verbose) cat(chr, "...")
      chrsta <- probeAnno[paste(chr,"start",sep=".")]
      chrend <- probeAnno[paste(chr,"end",sep=".")]
      stopifnot(all(chrend>chrsta))
      chrmid <- round((chrsta+chrend)/2)
      chridx <- probeAnno[paste(chr,"index",sep=".")]
      if (checkUnique){
        chruni <- probeAnno[paste(chr,"unique",sep=".")]
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
        cherID <- paste(thisCellType, this.sample, paste("chr",chr,sep=""), paste("cher",z,sep=""),sep=".")
        cherID <- gsub("(^\\.+)|(\\.+$)","",cherID)#remove any leading and trailing dots
        thisCher <- new("cher", name=cherID, chromosome=chr, start=as.integer(names(x)[1]), end=as.integer(names(x)[length(x)]), cellType=as.character(thisCellType), antibody=as.character(this.sample), maxLevel=max(x), score=attr(x,"score"), probes=as.character(chridx[names(x)]))
        return(thisCher)})
      return(chr.chers)
    })
    resultChers[[i]] <- thisModChers
  }#for
  resultChers <- unlist(resultChers, recursive=FALSE, use.names=FALSE)
  resultChers <- unlist(resultChers, recursive=FALSE)
  class(resultChers) <- c("cherList",class(resultChers))
  return(resultChers)
}# findChersOnSmoothed
