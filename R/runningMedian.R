
sliding.median <- function(positions, scores, half.width, return.counts=TRUE) {
  stopifnot(!is.unsorted(positions), length(positions) == length(scores), half.width >= 0)
  res <- .Call("sliding_median", as.integer(positions), as.numeric(scores), as.integer(half.width), PACKAGE="Ringo")
  if (return.counts){
    colnames(res) <- c("median","count")
    rownames(res) <- positions
  } else {
    res <- res[,1, drop=TRUE]
    names(res) <- positions
  }
  return(res)
}#sliding.median

sliding.quantile <- function(positions, scores, half.width, prob=0.5, return.counts=TRUE) {
  stopifnot(!is.unsorted(positions), length(positions) == length(scores), half.width >= 0, prob >= 0, prob <= 1)
  res <- .Call("sliding_quantile", as.integer(positions), as.numeric(scores), as.integer(half.width), as.numeric(prob), PACKAGE="Ringo")
  if (return.counts){
    colnames(res) <- c("quantile","count")
    rownames(res) <- positions
  } else {
    res <- res[,1, drop=TRUE]
    names(res) <- positions
  }
  return(res)
}#sliding.quantile

computeRunningMedians <- function(xSet, probeAnno, modColumn="Cy5", allChr=c(1:19,"X","Y"), winHalfSize=400, min.probes=5, quant=0.5, combineReplicates=FALSE, verbose=TRUE)
{
  stopifnot(inherits(xSet,"ExpressionSet"), all(is.character(allChr)),
            is.numeric(quant), (quant>=0)&(quant<=1), length(quant)==1)
  # initialize result matrix:
  if (combineReplicates)
    grouping <- factor(pData(xSet)[[modColumn]])
  else
    grouping <- factor(sampleNames(xSet))

  newExprs <- matrix(NA, nrow=nrow(exprs(xSet)), ncol=nlevels(grouping))
  rownames(newExprs) <- featureNames(xSet)
  for (chr in allChr){
    if (verbose) cat("\nChromosome",chr, "...\n")
    chrsta <- get(paste(chr,"start",sep="."), env=probeAnno)
    chrend <- get(paste(chr,"end",sep="."), env=probeAnno)
    chrmid <- round((chrsta+chrend)/2)
    chridx <- get(paste(chr,"index",sep="."), env=probeAnno)      
    
    for (i in 1:nlevels(grouping)){
      modSamples   <- which(grouping == levels(grouping)[i])
      if (verbose) cat(sampleNames(xSet)[modSamples],"... ")
      combined.dat <- as.vector(t(exprs(xSet)[chridx,modSamples,drop=FALSE]))
      # as.vector(t(X)) leads to columns (samples) being appended one value by one value into long vector
      combined.pos <- rep(chrmid, each=length(modSamples))
      ## old version: only median possible
      #slidingRes <- sliding.median(positions=combined.pos, scores=combined.dat, half.width=winHalfSize, return.counts=TRUE)
      slidingRes <- sliding.quantile(positions=combined.pos, scores=combined.dat, half.width=winHalfSize, prob=quant, return.counts=TRUE)
      slidingRes <- slidingRes[seq(1, nrow(slidingRes)+1-length(modSamples), by=length(modSamples)),,drop=FALSE]
      chrrm <- slidingRes[,"quantile"] #chrrm <- slidingRes[,"median"]
      slidingRes[,"count"] <- slidingRes[,"count"]/length(modSamples)
      areBelow <- slidingRes[,"count"] < min.probes
      if (any(areBelow)) chrrm[areBelow] <- NA
      stopifnot(length(chrrm) == length(chrmid))
      newExprs[chridx,i] <- chrrm
    }#for (i in 1:length(all.mods))
  } #for (chr in allChr)

  # cat construct ExpressionSet of results:
  sample.labels <- vector("character", nlevels(grouping))
  for (i in 1:nlevels(grouping)) sample.labels[i] <- as.character(levels(grouping)[i])
  newPD <- new("AnnotatedDataFrame", data=data.frame(label=sample.labels, row.names=sample.labels),
               varMetadata=data.frame("varLabel"=c("label"),row.names=c("label")))
  newEset <- new('ExpressionSet',exprs=newExprs,  phenoData = newPD)
      # experimentData = [MIAME], annotation = [character], ...
  featureNames(newEset) <- featureNames(xSet)
  sampleNames(newEset)  <- sample.labels
  return(newEset)
}#computeRunningMedians

