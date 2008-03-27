addTypes = function(curTss, curRegions, tssCover=c(-25, 25), tssUpstream=5000, tssDownstream=20000, useDist2NextUpGene = FALSE) {
  # check, which transcripts are covered
  theTSS = ifelse(curTss$strand==1, curTss$start, curTss$end)
  start = theTSS + tssCover[1]
  end = theTSS + tssCover[2]
  curRegions = lapply(curRegions, function(p) { p$typeTSSCover = c(p$typeTSSCover, curTss$name[(start >= p$start) & (end <= p$end)]); return(p) })
  # cers in upstream area
  if (useDist2NextUpGene) {
    delta = pmin(curTss$dist2NextUpGene, tssUpstream)
  }
  else {
    delta = tssUpstream
  }
  stopifnot(all(delta > 0))
  #delta = tssUpstream
  start = ifelse(curTss$strand==1, curTss$start - delta, curTss$end)
  end = ifelse(curTss$strand==1, curTss$start, curTss$end + delta)
  curRegions = lapply(curRegions, function(p) { p$typeUpstream = c(p$typeUpstream, curTss$name[(start <= p$start) & (end >= p$end)]); return(p) })
  # cers in downstream area
  strandDelta = ifelse(curTss$strand==1,
    curTss$intron1start - curTss$start + 1,
    curTss$end - curTss$intron1end + 1)
  strandDelta[is.na(strandDelta)] = (curTss$end - curTss$start + 1)[is.na(strandDelta)]
  stopifnot(all(strandDelta > 0))
  delta = pmin(strandDelta, tssDownstream)
  start = ifelse(curTss$strand==1, curTss$start, curTss$end - delta)
  end = ifelse(curTss$strand==1, curTss$start + delta, curTss$end)
  curRegions = lapply(curRegions, function(p) { p$typeDownstream = c(p$typeDownstream, curTss$name[(start <= p$start) & (end >= p$end)]); return(p) })
  # cers in intron area
  start = curTss$intron1start
  end = curTss$intron1end
  start[is.na(start)] = 0
  end[is.na(end)] = 0
  stopifnot(length(start) == length(end))
  curRegions = lapply(curRegions, function(p) { p$typeIntron = c(p$typeIntron, curTss$name[(start <= p$start) & (end >= p$end)]); return(p) })
  return(curRegions)
}

typifyList = function(theList, geneList, g2t) {
  t2g = getTranscriptToGeneMap(geneList)
  allTypes = c("typeUpstream", "typeDownstream", "typeIntron", "typeTSSCover")
  isUnique = function(p, type) {
    cur = p[[type]]
    if (length(cur) == 0) {
      return(FALSE)
    }
    else {
      if (length(cur) > 1 & !all(cur %in% g2t[[t2g[[cur[[1]]]]]])) {
        return(FALSE)
      }
      for (t in allTypes) {
        if (t == type) {
          next
        }
        if (length(p[[t]]) > 0) {
          return(FALSE)
        }
      }
      return(TRUE)
    }
    stop
  }
  isLostCover = function(p) {
    return(length(p$typeUpstream) == 0 &
           length(p$typeDownstream) == 0 &
           length(p$typeIntron) == 0 &
           length(p$typeTSSCover) == 0)
  }
  typify = function(p) {
    if (isLostCover(p)) {
      p$type = paste(p$type, "L/C", sep="")
    }
    else if (isUnique(p, "typeUpstream")) {
      p$type = paste(p$type, "U", sep="")
    }
    else if (isUnique(p, "typeDownstream")) {
      p$type = paste(p$type, "D", sep="")
    }
    else if (isUnique(p, "typeIntron")) {
      p$type = paste(p$type, "I", sep="")
    }
    else if (isUnique(p, "typeTSSCover")) {
      p$type = paste(p$type, "C", sep="")
    }
    else {
      p$type = paste(p$type, "Non unique", sep="")
    }
    return(p)
  }
  theList = lapply(theList, typify)
  return(theList)
}


getGeneToTranscriptMap = function(geneList) {
  allGenes = unique(geneList$gene)
  g2t = sapply(allGenes, function(g) { return(geneList$name[geneList$gene == g]) })
  names(g2t) = allGenes
  return(g2t)
}

getTranscriptToGeneMap = function(geneList) {
  t2g = as.list(geneList$gene)
  names(t2g) = geneList$name
  return(t2g)
}

combine.tables = function(tables) {
  stopifnot(all(sapply(tables, is.table)))
  allCols = unique(c(lapply(tables, names), recursive=TRUE))
  mat = matrix(0, length(tables), length(allCols))
  colnames(mat) = sort(allCols)
  rownames(mat) = names(tables)
  for (n in names(tables)) {
    tab = tables[[n]]
    for (m in names(tab)) {
      mat[n, m] = tab[m]
    }
  }
  return(mat)
}

calcOverlapDfVsDf = function(df.check, df.target, overlap = 0.75, returnCode=TRUE, verbose=TRUE, nonSymmetric=FALSE) {
  vals = list()
  for (chr in c(as.character(1:19), "X", "Y")) {
    if (verbose) cat(chr, "")
    toCheck = which(df.check$chr == chr)
    if (length(toCheck) == 0) {
      next
    }
    targetInd = which((df.target$chr == chr))
    if (length(targetInd) == 0) {
      next
    }
    target = df.target[targetInd, ]
    calcOverlap = function(ind) {
      val = pmin(target$end, df.check$end[ind]) - pmax(target$start, df.check$start[ind]) + 1
      if (nonSymmetric) {
        # this was the non-symmetric version
        theOverlap = val / (df.check$end[ind] - df.check$start[ind] + 1)
      }
      else {
        theOverlap = val / (pmin(target$end - target$start, df.check$end[ind] - df.check$start[ind]) + 1)
      }
      overlapIndex = (theOverlap >= overlap)
      nOverlap = sum(overlapIndex)
      if (returnCode) {
        ret = cbind(rep(ind, nOverlap), targetInd[overlapIndex],
          target$code[overlapIndex], theOverlap[overlapIndex])
      }
      else {
        ret = cbind(rep(ind, nOverlap), targetInd[overlapIndex],
          theOverlap[overlapIndex])
      }
      return(ret)
    }
    if (length(toCheck) > 0) {
      newVals = lapply(toCheck, calcOverlap)
      newVals = newVals[listLen(newVals) > 0]
      if (verbose) cat(sprintf("(%d)", length(newVals)))
      vals = c(vals, newVals)
    }
  }
  if (verbose) cat("\n")
  # sapply doesn't do the job... so we have to loop explicitly... very ugly
  nColumns = ifelse(returnCode, 4, 3)
  ret = matrix("", sum(listLen(vals) / nColumns), nColumns)
  if (length(vals) != 0) {
    ind = 1
    for (i in 1:length(vals)) {
      cur = vals[[i]]
      ret[ind:(ind + nrow(cur) - 1), ] = cur
      ind = ind + nrow(cur)
    }
  }
  if (returnCode) {
    ret = data.frame(checkIndex=as.integer(ret[, 1]), targetIndex=as.integer(ret[, 2]), targetCode=I(ret[, 3]), overlap=as.double(ret[, 4]))
  }
  else {
    ret = data.frame(checkIndex=as.integer(ret[, 1]), targetIndex=as.integer(ret[, 2]), overlap=as.double(ret[, 3]))
  }
  return(ret)
}


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

## Auxiliary function to convert NimbleGen's feature report file into
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
  return(unique(unlist(sapply(cl, function(cher) cher[c("typeUpstream", "typeInside", "typeDownstream")]), use.names=FALSE)))
}# getFeats

          
