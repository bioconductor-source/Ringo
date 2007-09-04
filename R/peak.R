#setOldClass("peak")
#setOldClass("peakList")

newPeak <- function(name, chr, start, end, cellType=NULL, antibody, maxPeak, score=NULL, probes=c(), ...) {
  ## element 'antibody' was named 'modification' before
  peak <- list(name=name, chr=chr, start=start, end=end, cellType=cellType, antibody=antibody, typeUpstream=list(), typeDownstream=list(), maxPeak=maxPeak, score=score, probes=probes)  
  ## add any further named arguments in "..."
  peak <- c(peak, lapply(as.list(match.call(expand.dots=FALSE)[["..."]]),eval))
  class(peak) = "peak"
  return(peak)
}#newPeak

print.peak <- function(x, ...) {
  nUp = length(x$typeUpstream)
  nDown = length(x$typeDownstream)
  cat(x$name, "\nChr",x$chr,":",x$start, "-", x$end, "\n")
  with(x, if (!is.null(antibody)) cat("Antibody :",antibody,"\n"))
  with(x, if (!is.null(score)) cat("Score =",score,"\n"))
  with(x, if (!is.null(probes)) cat("Number of probes =",length(probes),"\n"))
  cat("nUp =", nUp, "nDown =", nDown, "\n")
  invisible(NULL)
}#print.peak

plot.peak <- function(x, dat, probeAnno, samples=NULL, extent=1000, gff=NULL,...){  
  stopifnot(inherits(x,"peak"), inherits(dat,"ExpressionSet"))
  if (!is.null(samples)&is.character(samples))
    samples <- match(samples, sampleNames(dat))
  if (!is.null(samples)&is.numeric(samples))
    stopifnot(all(samples) %in% 1:ncol(dat))
  peakVal <- chipAlongChrom(dat, chrom=x$chr, samples=samples, xlim=c(x$start-extent, x$end+extent), probeAnno=probeAnno, gff=gff, ...)
  rug(x=c(x$start,x$end), side=3, lwd=3, col="gold")
  if (!is.null(x$maxPeak))
    legend(x="topright", legend=paste("Max.Level:",round(x$maxPeak,digits=2)),fill="gold", bty="n")
  invisible(peakVal)
}#plot.peak

as.data.frame.peakList <- function(x, row.names=NULL, optional=FALSE,...) {
  stopifnot(is.list(x), inherits(x[[1]],"peak"))
  np <- length(x)
  peak.slots <- names(x[[1]])
  p.name <- vector("character",np)
  p.chr <- vector("character",np)
  p.start <- integer(np)
  p.end <- integer(np)
  p.cellType <- vector("character",np)
  p.antibody <- vector("character",np)
  p.transcripts <- vector("character",np)
  p.maxPeak <- vector("numeric",np)
  p.score <- vector("numeric",np)
  for (i in 1:np){
     p <- x[[i]]
     p.name[i]  <- p$name
     p.chr[i]   <-  p$chr
     p.start[i] <-  p$start
     p.end[i] <-  p$end
     p.cellType[i] <-  ifelse(is.null(p$cellType), NA, p$cellType)
     p.antibody[i] <- ifelse(is.null(p$antibody), p$modification, p$antibody)
     p.transcripts[i] <- paste(c(p$typeUpstream, p$typeDownstream), collapse=" ")
     p.maxPeak[i] <- ifelse(is.null(p$maxPeak), NA, p$maxPeak)
     p.score[i] <- ifelse(is.null(p$score), NA, p$score)
   }#for i
   # df <- t(sapply(x, function(p) {c(p$name, p$chr, p$start, p$end, ifelse(is.null(p$cellType), NA, p$cellType), p$modification, length(p$typeUpstream), length(p$typeDownstream), paste(c(p$typeUpstream, p$typeDownstream), collapse=" "), p$maxPeak, p$score)}))
   #stopifnot(ncol(df)==11)
  #df <- data.frame(name=I(df[, 1]), chr=I(df[, 2]), start=as.integer(df[, 3]), end=as.integer(df[, 4]), cellType=I(df[, 5]), modification=I(df[, 6]), nUpstream=as.integer(df[, 7]), nDownstream=as.integer(df[, 8]),transcripts=I(df[, 9]), maxPeak=as.double(df[, 10]), score=I(round(as.double(df[, 11]),digits=4)), stringsAsFactors=FALSE, row.names=NULL)
  df <- data.frame(name=p.name, chr=p.chr, start=p.start, end=p.end, cellType=p.cellType, antibody=p.antibody,  transcripts=p.transcripts, maxPeak=p.maxPeak, score=p.score, stringsAsFactors=FALSE, row.names=row.names)
  return(df)
}#as.data.frame.peakList


generatePeakList = function(peaks,gff, g2t=NULL, allChr=c(1:19, "X", "Y"), tssCover=c(-25, 25), tssUpstream=5000, tssDownstream=20000, cellType="HL1", verbose=TRUE) {
  if (is.null(g2t)) {
    if (verbose) cat("generating GeneToTranscriptMap\n")
    g2t = getGeneToTranscriptMap(gff)
  }
  nPeaks = sum(sapply(peaks, listLen))
  peakList = vector("list", nPeaks)
  index = 1
  for (antibody in names(peaks)) {
    if (verbose) cat("\n| Processing antibody", antibody, "\n+-- chr")
    curMod = peaks[[antibody]]
    for (chr in allChr) {
      if (verbose) cat("", chr)
      curChr = curMod[[chr]]
      if (length(curChr)==0) {next} # no peaks on that chromosome
      curTss = gff[gff$chr == chr, ]
      generatePeak = function(pName) {
        p = curChr[[pName]]
        n = names(p)
        this.peak.score = attr(p,"score")
        return(newPeak(peakID=pName, chr, start=as.integer(n[1]), end=as.integer(n[length(n)]), cellType=cellType, antibody=antibody, maxPeak=max(p), score=this.peak.score))
      }# generatePeak
      curPeaks = lapply(names(curChr), generatePeak)
      curPeaks = addTypes(curTss, curPeaks, tssCover, tssUpstream, tssDownstream)
      # save the peaks
      peakList[index:(index+length(curPeaks)-1)] = curPeaks
      index = index + length(curPeaks)
    }
  }
  if (verbose) cat("\n")
  names(peakList) = sapply(peakList, function(p) p$name)
  peakList = typifyList(peakList, gff, g2t)
  class(peakList) = "peakList"
  return(peakList)
}#generatePeakList


relatePeaks <- function(pl, gff, upstream=5000, verbose=TRUE){
  stopifnot(is.list(pl),inherits(pl[[1]],"peak"),
            inherits(gff,"data.frame"),
            all(c("strand","name","start","end","chr") %in% names(gff)),
            all(gff$start<gff$end))
  peakChr <- sapply(pl, function(x) x$chr)
  peakMid <- sapply(pl, function(x) round((x$start+x$end)/2))
  ## get borders of upstream region:
  gffUpStart <- ifelse(gff$strand==1, gff$start-upstream, gff$end+1)
  gffUpEnd   <- ifelse(gff$strand==1, gff$start-1, gff$end+upstream)
  realTSS <- ifelse(gff$strand==1, gff$start, gff$end)
  names(realTSS) <- gff$name
  if (verbose) cat("Relating",length(pl),"peaks to GFF:\n...")
  for (i in 1:length(pl)){
    if (verbose & i%%1000==0) cat(i,"...")
    p <- pl[[i]]
    p$typeUpstream <- subset(gff, gff$chr==p$chr & peakMid[i]>=gffUpStart & peakMid[i] <= gffUpEnd, select="name", drop=TRUE)
    p$typeDownstream <- subset(gff, gff$chr==p$chr & peakMid[i]>=gff$start & peakMid[i] <= gff$end, select="name", drop=TRUE)
    p$distMid2TSS <- abs(realTSS[c(p$typeUpstream, p$typeDownstream, recursive=TRUE)]-peakMid[i])
    # the c(,recursive) is to prevent (illegal) indexing with 0-row data.frames
    if ("symbol" %in% names(gff)){
      p$upSymbol <-  subset(gff, gff$chr==p$chr & peakMid[i]>=gffUpStart & peakMid[i] <= gffUpEnd, select="symbol", drop=TRUE)
      p$downSymbol <- subset(gff, gff$chr==p$chr & peakMid[i]>=gff$start & peakMid[i] <= gff$end, select="symbol", drop=TRUE)
    }
    p$type <- paste(c("U","D")[c(length(p$typeUpstream)>0,length(p$typeDownstream)>0)],collapse="/")
    pl[[i]] <- p
  }#for
  return(pl)
}#relatePeaks
