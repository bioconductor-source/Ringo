
splitAndSimplify <- function(toSplit, bySplit, collapse.by="/", nameSep="[[:space:]]", doUnique=FALSE, doCollapse=TRUE){
  # function to
  # a: do split and make sure that no object has a combined name
  # b: return vector instead of list
  splitted <- split(toSplit,bySplit)
  if (doUnique) #should 2 or more MRs of the same class be treated as only one?
    collapsefun <- function(x) paste(unique(x),collapse=collapse.by)
  else
    collapsefun <- function(x) paste(sort(x),collapse=collapse.by)
  if (!doCollapse) collapsefun <- function(x) sort(x)
  # 'sort' since order in DF does not seem to be genomic order:
  splitted <- unlist(lapply(splitted,collapsefun),use.names=TRUE)
  singleNames <- strsplit(names(splitted), nameSep)
  splitted <- rep(splitted,listLen(singleNames))
  names(splitted) <- unlist(singleNames, use.names=FALSE)
  return(splitted)
}# splitAndSimplify

### Auxiliary function to convert NimbleGen's pair file format into
### .xys file. Only use this function, if actual xys-file is not available
pair2xys <- function(pair.file, path=getwd()){  
  stopifnot(length(pair.file)==1, is.character(pair.file), grep(".pair.*$",pair.file)==1, file.exists(pair.file))
  pair.header <- scan(pair.file,nlines=1,quiet=TRUE, what=character(0))
  xys.file <- gsub("_?pair.*$",".xys",pair.file)
  pair.data <- read.delim(pair.file, as.is=TRUE, comment.char="#")
  xys.data <- pair.data[,c("X", "Y", "PM", "PROBE_ID")]
  names(xys.data) <- c("X", "Y", "SIGNAL", "PROBE_ID")
  cat(pair.header,"\n", sep=" ", file=file.path(path,xys.file))
  write.table(xys.data, file=file.path(path, xys.file), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)
  invisible(NULL)
}#pair2xys

### Auxiliary function to convert NimbleGen's feature report file into
### .xys file. Only use this function, if actual xys-file is not available
ftr2xys <- function(ftr.file, path=getwd()){  
  stopifnot(length(ftr.file)==1, is.character(ftr.file), grep(".ftr$",ftr.file)==1, file.exists(ftr.file))
  ftr.header <- scan(ftr.file,nlines=1,quiet=TRUE, what=character(0))
  xys.file <- gsub("ftr$","xys",ftr.file)
  ftr.data <- read.delim(ftr.file, as.is=TRUE, comment.char="#")
  xys.data <- ftr.data[,c("X","Y","SIGNAL_MEAN","FGD_PIX")]
  names(xys.data) <- c("X","Y","SIGNAL","COUNT")
  cat(ftr.header,"\n", sep=" ", file=file.path(path,xys.file))
  write.table(xys.data, file=file.path(path, xys.file), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)
  invisible(NULL)
}#ftr2xys

## function to get runs of 1 in a binary vector, adopted from pairs project
clusters <- function(x, minLen=3, doSelect=FALSE) {
  x = c(0,x,0)
  diffx <- diff(x)
  start = which(diffx==1)
  end   = which(diffx==-1)
  len   = end-start
  if (doSelect) {
    sel   = len >= minLen
    start=start[sel]
    len=len[sel]
  } 
  cbind(start,len)
}#clusters

### function to take mean over sample groups in an ExpressionSet and
###  return an ExpressionSet holding those means
takeMeanOverGroups <- function(xSet, modColumn="Cy5")
{
  stopifnot(inherits(xSet,"ExpressionSet"), modColumn %in% names(pData(xSet)))

  grouping <- factor(pData(xSet)[[modColumn]])
  ngroups  <- nlevels(grouping)
  groupmat <- matrix(0, nrow=ncol(xSet), ncol=ngroups)
  for (i in 1:ngroups){
    thisLevel <- levels(grouping)[i]
    inGroup <- (grouping == thisLevel)
    nInGroup <- sum(inGroup, na.rm=TRUE)
    groupmat[inGroup,i] <- 1/nInGroup
  }
  # compute group-wise means for each probe
  datmat <- exprs(xSet)%*%groupmat
  colnames(datmat) <- as.character(levels(grouping))  
  ### TODO: combine pData as well meaningful
  new.df        <- data.frame(as.character(levels(grouping)))
  names(new.df) <- modColumn
  newPD <- new("AnnotatedDataFrame", data=new.df,
               varMetadata=data.frame("varLabel"=colnames(new.df),
                 row.names=colnames(new.df)))
  newEset <- new("ExpressionSet",  exprs=datmat, phenoData=newPD)
  featureNames(newEset) <- featureNames(xSet)
  return(newEset)
}#takeMeanOverGroups


compute.gc <- function(probe.sequences, digits=2){
  stopifnot(is.character(probe.sequences))
  splitted.seqs <- strsplit(toupper(probe.sequences),split="")
  round(sapply(splitted.seqs, function(x) length(grep("[GC]",x)))/
    listLen(splitted.seqs), digits=digits)
}#compute.gc

whichCsr <- function(X, arr.ind=TRUE){
  ## function to get a two-column matrix containing the indices of the
  ### non-zero elements in a "matrix.csr" class matrix
  stopifnot(inherits(X, "matrix.csr"))
  if (all(X@ra==0)) return(NULL)
  res <- cbind(rep(seq(dim(X)[1]),diff(X@ia)), # row indices
               X@ja )# column indices directly saved in matrix.csr format
  colnames(res) <- c("row","col")
  ## remove zero elements
  res <- res[X@ra != 0,,drop=FALSE]
  return(res)
}# whichCsr

getFeats <- function(cl){
  stopifnot(is.list(cl), inherits(cl[[1]],"cher"))
  return(unique(unlist(sapply(cl, function(cher) cher@extras[c("typeUpstream", "typeInside", "typeDownstream")]), use.names=FALSE)))
}# getFeats


exportCCData <- function(X, pA, method="GFF", outfile="myCCData.gff", samples, chrs, checkUnique=TRUE, uniqueCodes=c(0), verbose=TRUE)
{
  stopifnot(inherits(X, "ExpressionSet"), inherits(pA, "probeAnno"),
            validObject(pA))
  if (missing(chrs))  chrs <- chromosomeNames(pA)
  if (missing(samples)) samples <- 1:ncol(X)
  ## take mean over selected samples
  exprs(X)[,1] <- rowMeans(exprs(X)[,samples, drop=FALSE])
      
  method <- match.arg(method, c("GFF"))

  if (method=="GFF"){
    if (verbose) cat("Preparing GFF...\n")
    # init gff file
    cat("##gff-version 2",
        paste("## source-version","Ringo",package.version("Ringo")),
        paste("## date",Sys.Date()),
        "## format:",
        paste("## seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", sep="\t"),
        sep="\n", #collapse="\n",
        file=outfile, append=FALSE)
     
    for (chr in chrs){
      if (verbose) cat("Chromosome",chr, "...\n")
      chrsta <- pA[paste(chr,"start",sep=".")]
      chrend <- pA[paste(chr,"end",sep=".")]
      chridx <- pA[paste(chr,"index",sep=".")]
      if (checkUnique){
        chruni <- pA[paste(chr,"unique",sep=".")]
        stopifnot(length(chruni)==length(chridx))
        chridx <- chridx[chruni %in% uniqueCodes]
        chrsta <- chrsta[chruni %in% uniqueCodes]
        chrend <- chrend[chruni %in% uniqueCodes]
        if (length(chridx)==0){
          warning(paste("No reporters with unique hits",
                        "on chromosome",chr,".\n")); next}
      } #  if (checkUnique)
      stopifnot(all(chrend >= chrsta))
      #GFF format: <seqname>\t<source>\t<feature>\t<start>\t<end>\t<score>\t<strand>\t<frame>\t<attribute>
      n <- length(chridx)
      dat <- round(exprs(X)[chridx, 1], digits=3)
      chrdat <- data.frame(
            "seqname"=rep(paste("chr",gsub("^chr","",chr),sep=""), n),
            "source"=rep("ChIP", n), #as.character(chridx),
            "feature"=rep("reporter level", n),
            "start"=as.integer(chrsta), "end"=as.integer(chrend),
            "score"=dat,
            "strand"=rep(".", n), "frame"=rep(".", n),
            "attribute"=paste(chridx,";","chr",chr,":",
              chrsta,"-",chrend,sep=""),stringsAsFactors=FALSE)
      write.table(chrdat, sep="\t", append=TRUE, row.names=FALSE,
                  col.names=FALSE, file=outfile, quote=FALSE)
    }# for chr in chrs
  }# if method=="GFF"
  if (verbose)
    cat("Written to file '", outfile, "'.\n", sep="")
}# exportCCData
