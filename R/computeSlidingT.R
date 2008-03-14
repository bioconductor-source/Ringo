computeSlidingT <- function(xSet, probeAnno, allChr=c(1:19,"X","Y"), test="one.sample", grouping=NULL, winHalfSize=400, min.probes=5, checkUnique=TRUE, uniqueCodes=c(0), verbose=TRUE)
{
  stopifnot(inherits(xSet,"ExpressionSet"), all(is.character(allChr)))
  test <- match.arg(test, c("one.sample","two.sample"))
  if (test == "two.sample"){
    grouping <- factor(grouping)
    if (any(is.null(grouping), length(grouping)!=ncol(xSet), nlevels(grouping)!=2))
      stop(paste("Argument 'factor' needs to be a factor of length",ncol(xSet),"with 2 levels.\n"))
  } else {
    grouping <- factor(1:ncol(xSet))
  }
  # initialize result matrix:
  probeMeans <- matrix(NA, nrow=nrow(exprs(xSet)), ncol=nlevels(grouping), dimnames=list(x=featureNames(xSet), y=as.character(grouping)))
  probeSds   <- matrix(NA, nrow=nrow(exprs(xSet)), ncol=nlevels(grouping), dimnames=list(x=featureNames(xSet), y=as.character(grouping)))
  probeCounts <- matrix(NA, nrow=nrow(exprs(xSet)), ncol=nlevels(grouping), dimnames=list(x=featureNames(xSet), y=as.character(grouping)))
  if (verbose) cat("\n computing probe-wise mean and standard deviation in sliding window.\n chr",chr, "...\n")
  for (chr in allChr){
    chrsta <- get(paste(chr,"start",sep="."), env=probeAnno)
    chrend <- get(paste(chr,"end",sep="."), env=probeAnno)
    chrmid <- round((chrsta+chrend)/2)
    chridx <- get(paste(chr,"index",sep="."), env=probeAnno)
    if (checkUnique){
      chruni <- get(paste(chr,"unique",sep="."), env=probeAnno)
      stopifnot(length(chruni)==length(chridx))
      chridx <- chridx[chruni %in% uniqueCodes]
      chrmid <- chrmid[chruni %in% uniqueCodes]
    } #  if (checkUnique)
    for (i in 1:nlevels(grouping)){
      modSamples   <- which(grouping == levels(grouping)[i])
      if (verbose) cat(sampleNames(xSet)[modSamples],"... ")
      combined.dat <- as.vector(t(exprs(xSet)[chridx,modSamples,drop=FALSE]))
      # as.vector(t(X)) leads to columns (samples) being appended one value by one value into long vector
      combined.pos <- rep(chrmid, each=length(modSamples))
      slidingRes <- sliding.meansd(positions=combined.pos, scores=combined.dat, half.width=winHalfSize)
      slidingRes <- slidingRes[seq(1, nrow(slidingRes)+1-length(modSamples), by=length(modSamples)),,drop=FALSE]
      haveSufficientProbes <- slidingRes[,"count"]>=min.probes
      probeMeans[chridx[haveSufficientProbes],i] <- slidingRes[haveSufficientProbes,"mean"]
      probeSds[chridx[haveSufficientProbes],i] <- slidingRes[haveSufficientProbes,"sd"]
      probeCounts[chridx[haveSufficientProbes],i] <- slidingRes[haveSufficientProbes,"count"]
    }#for (i in 1:length(all.mods))
  } #for (chr in allChr)

  ## do the probe-wise (regularized) t-testing:
  if (verbose) cat("\n computing t-statistics...\n")
  if (test == "one.sample"){
    ## regularization factor
    sampleS0s <- apply(probeSds,2,median, na.rm=TRUE)
    probeS0s  <- matrix(sampleS0s, nrow=nrow(xSet), ncol=nlevels(grouping),byrow=TRUE)
    ## one-sample t-statistic if mean is zero,using SE=SD/sqrt(n)
    probeTs <- probeMeans/(probeSds+probeS0s)*sqrt(probeCounts)
  }
  if (test == "two.sample"){
    ## compute combined SDs
    n1 <- which(grouping==levels(grouping)[1])
    n2 <- which(grouping==levels(grouping)[2])
    combSds <- as.vector(sqrt(probeSds^2 %*% matrix(c(1/n1, 1/n2))))
    stopifnot(length(combSds)==nrow(probeMeans))
    ## regularization factor
    combS0  <- median(combSds, na.rm=TRUE)
    ## two-sample regularized t-statistic (by Welch,i.e. unequal variances)
    probeTs <- matrix((probeMeans[,1]-probeMeans[,2])/(combSds+combS0))
  }
  # cat construct ExpressionSet of results:
  if (verbose) cat("preparing result...")
  newExprs <- probeTs
  rownames(newExprs) <- featureNames(xSet)
  if (test == "one.sample")
    sample.labels <- as.character(grouping)
  else if (test == "two.sample")
    sample.labels <- "t-stat"
  newPD <- new("AnnotatedDataFrame", data=data.frame(label=sample.labels, row.names=sample.labels), varMetadata=data.frame("varLabel"=c("label"),row.names=c("label")))
  newEset <- new('ExpressionSet',exprs=newExprs,  phenoData = newPD)
  featureNames(newEset) <- featureNames(xSet)
  sampleNames(newEset)  <- sample.labels
  if (verbose) cat("done.\n")
  return(newEset)
}#computeSlidingT
