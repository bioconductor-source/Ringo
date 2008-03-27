
newCher <- function(name, chr, start, end, cellType=NULL, antibody, maxLevel, score=NULL, probes=c(), ...) {
  ## element 'antibody' was named 'modification' before
  cher <- list(name=name, chr=chr, start=start, end=end, cellType=cellType, antibody=antibody, typeUpstream=list(), typeDownstream=list(), maxLevel=maxLevel, score=score, probes=probes)  
  ## add any further named arguments in "..."
  cher <- c(cher, lapply(as.list(match.call(expand.dots=FALSE)[["..."]]),eval))
  class(cher) = "cher"
  return(cher)
}#newCher

print.cher <- function(x, ...) {
  nUp = length(x$typeUpstream)
  nIn = length(x$typeInside)
  cat(x$name, "\nChr",x$chr,":",x$start, "-", x$end, "\n")
  with(x, if (!is.null(antibody)) cat("Antibody :",antibody,"\n"))
  with(x, if (!is.null(score)) cat("Score =",score,"\n"))
  with(x, if (!is.null(probes)) cat("Number of probes =",length(probes),"\n"))
  cat("nUp =", nUp, "nDown =", nIn, "\n")
  invisible(NULL)
}#print.cher

plot.cher <- function(x, dat, probeAnno, samples=NULL, extent=1000, gff=NULL,...){  
  stopifnot(inherits(x,"cher"), inherits(dat,"ExpressionSet"))
  if (!is.null(samples)&is.character(samples))
    samples <- match(samples, sampleNames(dat))
  if (!is.null(samples)&is.numeric(samples))
    stopifnot(all(samples) %in% 1:ncol(dat))
  cherVal <- chipAlongChrom(dat, chrom=x$chr, samples=samples, xlim=c(x$start-extent, x$end+extent), probeAnno=probeAnno, gff=gff, ...)
  rug(x=c(x$start,x$end), side=3, lwd=3, col="gold")
  if (!is.null(x$maxLevel))
    legend(x="topright", legend=paste("Max.Level:",round(x$maxLevel,digits=2)),fill="gold", bty="n")
  invisible(cherVal)
}#plot.cher

as.data.frame.cherList <- function(x, row.names=NULL, optional=FALSE,...) {
  stopifnot(is.list(x), inherits(x[[1]],"cher"))
  np <- length(x)
  cher.slots <- names(x[[1]])
  p.name <- vector("character",np)
  p.chr <- vector("character",np)
  p.start <- integer(np)
  p.end <- integer(np)
  p.cellType <- vector("character",np)
  p.antibody <- vector("character",np)
  p.features <- vector("character",np)
  p.maxLevel <- vector("numeric",np)
  p.score <- vector("numeric",np)
  for (i in 1:np){
     p <- x[[i]]
     p.name[i]  <- p$name
     p.chr[i]   <-  p$chr
     p.start[i] <-  p$start
     p.end[i] <-  p$end
     p.cellType[i] <-  ifelse(is.null(p$cellType), NA, p$cellType)
     p.antibody[i] <- ifelse(is.null(p$antibody), p$modification, p$antibody)
     p.features[i] <- paste(c(p$typeUpstream, p$typeInside), collapse=" ")
     p.maxLevel[i] <- ifelse(is.null(p$maxLevel), NA, p$maxLevel)
     p.score[i] <- ifelse(is.null(p$score), NA, p$score)
   }#for i
  df <- data.frame(name=p.name, chr=p.chr, start=p.start, end=p.end, cellType=p.cellType, antibody=p.antibody, features=p.features, maxLevel=p.maxLevel, score=p.score, stringsAsFactors=FALSE, row.names=row.names)
  return(df)
}#as.data.frame.cherList


generateCherList = function(chers,gff, g2t=NULL, allChr=c(1:19, "X", "Y"), tssCover=c(-25, 25), tssUpstream=5000, tssDownstream=20000, cellType="HL1", verbose=TRUE) {
  if (is.null(g2t)) {
    if (verbose) cat("generating GeneToTranscriptMap\n")
    g2t = getGeneToTranscriptMap(gff)
  }
  nChers = sum(sapply(chers, listLen))
  cherList = vector("list", nChers)
  index = 1
  for (antibody in names(chers)) {
    if (verbose) cat("\n| Processing antibody", antibody, "\n+-- chr")
    curMod = chers[[antibody]]
    for (chr in allChr) {
      if (verbose) cat("", chr)
      curChr = curMod[[chr]]
      if (length(curChr)==0) {next} # no chers on that chromosome
      curTss = gff[gff$chr == chr, ]
      generateCher = function(pName) {
        p = curChr[[pName]]
        n = names(p)
        this.cher.score = attr(p,"score")
        return(newCher(cherID=pName, chr, start=as.integer(n[1]), end=as.integer(n[length(n)]), cellType=cellType, antibody=antibody, maxLevel=max(p), score=this.cher.score))
      }# generateCher
      curChers = lapply(names(curChr), generateCher)
      curChers = addTypes(curTss, curChers, tssCover, tssUpstream, tssDownstream)
      # save the Chers
      cherList[index:(index+length(curChers)-1)] = curChers
      index = index + length(curChers)
    }
  }
  if (verbose) cat("\n")
  names(cherList) = sapply(cherList, function(p) p$name)
  cherList = typifyList(cherList, gff, g2t)
  class(cherList) = "cherList"
  return(cherList)
}#generateCherList

## older version of relateChers, see newer one below:
relateChers <- function(pl, gff, upstream=5000, verbose=TRUE){
  stopifnot(is.list(pl),inherits(pl[[1]],"cher"),
            inherits(gff,"data.frame"),
            all(c("strand","name","start","end","chr") %in% names(gff)),
            all(gff$start<gff$end))
  cherChr <- sapply(pl, function(x) x$chr)
  cherMid <- sapply(pl, function(x) round((x$start+x$end)/2))
  ## get borders of upstream region:
  gffUpStart <- ifelse(gff$strand==1, gff$start-upstream, gff$end+1)
  gffUpEnd   <- ifelse(gff$strand==1, gff$start-1, gff$end+upstream)
  realTSS <- ifelse(gff$strand==1, gff$start, gff$end)
  names(realTSS) <- gff$name
  if (verbose) cat("Relating",length(pl),"ChIP-enriched regions to GFF:\n")
  for (i in 1:length(pl)){
    if (verbose & i%%1000==0) cat(i," ")
    p <- pl[[i]]
    p$typeUpstream <- subset(gff, gff$chr==p$chr & cherMid[i]>=gffUpStart & cherMid[i] <= gffUpEnd, select="name", drop=TRUE)
    p$typeInside <- subset(gff, gff$chr==p$chr & cherMid[i]>=gff$start & cherMid[i] <= gff$end, select="name", drop=TRUE)
    p$distMid2TSS <- abs(realTSS[c(p$typeUpstream, p$typeInside, recursive=TRUE)]-cherMid[i])
    # the c(,recursive) is to prevent (illegal) indexing with 0-row data.frames
    if ("symbol" %in% names(gff)){
      p$upSymbol <-  subset(gff, gff$chr==p$chr & cherMid[i]>=gffUpStart & cherMid[i] <= gffUpEnd, select="symbol", drop=TRUE)
      p$inSymbol <- subset(gff, gff$chr==p$chr & cherMid[i]>=gff$start & cherMid[i] <= gff$end, select="symbol", drop=TRUE)
    }
    p$type <- paste(c("U","D")[c(length(p$typeUpstream)>0,length(p$typeInside)>0)],collapse="/")
    pl[[i]] <- p
  }#for
  return(pl)
}#relateChers


relateChers2 <- function(pl, gff, upstream=5000, verbose=TRUE){
  stopifnot(is.list(pl),inherits(pl[[1]],"cher"),
            inherits(gff,"data.frame"),
            all(c("strand","name","start","end","chr") %in% names(gff)),
            all(gff$start<gff$end))
  cherChr <- sapply(pl, function(x) x$chr)
  cherMid <- sapply(pl, function(x) round((x$start+x$end)/2))
  chersDf <- data.frame(chr=cherChr, start=cherMid, end=cherMid,
                        stringsAsFactors=FALSE)
  ## get borders of upstream region:
  realTSS <- ifelse(gff$strand==1, gff$start, gff$end)
  names(realTSS) <- gff$name

  gffUpDf <- data.frame(chr=gff$chr,
     start=ifelse(gff$strand==1, gff$start-upstream, gff$end+1),
     end=ifelse(gff$strand==1, gff$start-1, gff$end+upstream),
     stringsAsFactors=FALSE)
  ## overlap with upstream regions
  overU <- regionOverlap(chersDf, gffUpDf)
  ### inside features:
  overI <- regionOverlap(chersDf, gff)
  ## add information to Chers that have any overlap
  haveOverlap <- which(rowSums(overU+overI)>0)
  if (length(haveOverlap)==0) return(pl)
  if (verbose)
    cat("Relating",length(pl),"ChIP-enriched regions to GFF...\n")
  for (i in haveOverlap){
    #if (verbose & i%%1000==0) cat(i,"...")
    p <- pl[[i]]
    p$typeUpstream <- gff$name[overU[i,]>0]
    p$typeInside   <- gff$name[overI[i,]>0]
    p$distMid2TSS  <- abs(realTSS[c(p$typeUpstream, p$typeInside)]-cherMid[i])
    # the c(,recursive) is to prevent (illegal) indexing with 0-row data.frames
    if ("symbol" %in% names(gff)){
      p$upSymbol   <-  gff$symbol[match(p$typeUpstream, gff$name)]
      p$inSymbol <-  gff$symbol[match(p$typeInside, gff$name)]
    }
    p$type <- paste(c("U","D")[c(length(p$typeUpstream)>0,length(p$typeInside)>0)],collapse="/")
    pl[[i]] <- p
  }#for
  return(pl)
}#relateChers2
