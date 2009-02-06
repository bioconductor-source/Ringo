##---------------------------------------------------------------------------
## Use either
## a.) NimbleGen's POS file or the output from MUMmer or BLAT or Exonerate
##    (querying all probes on the chip vs genome) OR
## b.) a data frame that provides the mapping of probe identifiers to
##     genomic match positions
## to construct probe annotaton datastructures
##
## IMPORTANT: After running MUMmer or BLAT condense their various output files
## into one table using the scripts 'condenseExonerateOutput.pl' or
## 'condenseBlatOutput.pl' respectively. This condensed table is simlilar to
## a NimbleGen POS file and usable by this script!
##
## This script has the following parts
## 1. Read the table containing the POS file or MUMmer or BLAT hits
## 2. Construct along-chromosome vectors for probeAnno environment:
##    start, end: (integer)
##    index: positions of the PM on the chip (integer)
##    unique: whether the probe hits multiple places (logical)
##-------------------------------------------------------------------------

## Note: This script has been adopted from the 'davidTiling' package.
## However, since ChIP-chip is not strand-specific,
##  the strand details were omitted.

## a. path to the result table after running the condense*Output.pl script

## FORMAT of that position file should be in NimbleGen's pos format:
# here's an example of their typical head (comment char and one whitespace at start of each line have been added only just here)

# SEQ_ID   CHROMOSOME     PROBE_ID        POSITION        LENGTH  COUNT
# chr1:4483138-4484138    CHR01   CHR01P004483140 4483140 50      1
# chr1:4483138-4484138    CHR01   CHR01P004483245 4483245 50      1

posToProbeAnno <- function(pos, chrNameColumn="CHROMOSOME",probeColumn="PROBE_ID", chrPositionColumn="POSITION", lengthColumn="LENGTH", sep="\t", genome="unknown", microarrayPlatform="unknown", verbose=TRUE, ...){
    if (!is.data.frame(pos) && !is.matrix(pos)){
        stopifnot(is.character(pos), file.exists(pos))
        hits <- read.delim(pos, header=TRUE, as.is=TRUE, ...)
        con <- file(pos, open="r")
        headerLine <- readLines(con, 1)
        headers <- unlist(strsplit(headerLine, split=sep))
        columnsAreIn <- c(chrNameColumn, probeColumn, chrPositionColumn, lengthColumn) %in% headers
        if (!all(columnsAreIn))
            stop(paste("\nColumn(s)",paste(c(chrNameColumn, probeColumn, chrPositionColumn, lengthColumn)[!columnsAreIn],sep=", ", collapse=", "),"are not found in file or data.frame",deparse(substitute(pos)),".\n"))
        signat <- as.list(rep("", length(headers)))
        names(signat) <- headers
        signat[[chrPositionColumn]] <- 0L
        signat[[lengthColumn]] <- 0L
        dat <- scan(con, what=signat, sep=sep, quiet=TRUE, ...)
        close(con)
        hits <- as.data.frame(dat, stringsAsFactors=FALSE)
    } else {
        hits <- as.data.frame(pos)
        columnsAreIn <- c(chrNameColumn, probeColumn, chrPositionColumn, lengthColumn) %in% names(hits)
        if (!all(columnsAreIn))
            stop(paste("\nColumn(s)",paste(c(chrNameColumn, probeColumn, chrPositionColumn, lengthColumn)[!columnsAreIn],sep=", ", collapse=", "),"are not found in file or data.frame",deparse(substitute(pos)),".\n"))
    }
    hits$END    <- hits[[chrPositionColumn]] + hits[[lengthColumn]] - 1
    hits$MIDDLE <- round((hits[[chrPositionColumn]] + hits$END)/2)
    hitOrder <- order(hits[[chrNameColumn]], hits$MIDDLE)
    hits <- hits[hitOrder,]
    if (length(grep("chr",hits[[chrNameColumn]], ignore.case=TRUE)))
    hits[[chrNameColumn]] <- gsub("chr(0)?","", hits[[chrNameColumn]], ignore.case=TRUE)
    probeindex <- 1:nrow(hits)
    probeMultiplicity <- table(hits[[probeColumn]])
    allProbeNames <- names(probeMultiplicity)
    if (verbose) cat("Creating probeAnno mapping for chromosome ")
    sp <- split(probeindex, as.factor(hits[[chrNameColumn]]))
    thisMapping = new.env(hash=TRUE, size=length(sp)*4)
    # for every chromosome, assign the corresponding vectors
    for(i in seq(along=sp)) {
        chromprobes <- sp[[i]]
        nm <- names(sp)[i]
        if (verbose) cat(nm, "")
        assign(paste(nm, "start", sep="."),  as.integer(hits[chromprobes, chrPositionColumn]), envir=thisMapping)
        assign(paste(nm, "end", sep="."),  as.integer(hits[chromprobes,"END"]) , envir=thisMapping)
        assign(paste(nm, "index", sep="."),  hits[chromprobes, probeColumn],  envir=thisMapping)
        idx <- match(hits[chromprobes, probeColumn], allProbeNames)
        assign(paste(nm, "unique", sep="."), as.integer(as.numeric(probeMultiplicity[idx]>1)*3), envir=thisMapping)
        ## uniqueness codes:  0 = probe has one unique hit;   3= probe has multiple hits
    }# for(i in seq(along=sp))
    if (verbose) cat("Done.\n")
    myProbeAnno <- new("probeAnno", map=thisMapping)
    genome(myProbeAnno) <- genome
    arrayName(myProbeAnno) <- microarrayPlatform
    return(myProbeAnno)
}#posToProbeAnno

validProbeAnno <- function(probeAnno){
  stopifnot(inherits(probeAnno, "probeAnno")||is.environment(probeAnno))
  if (!is.environment(probeAnno))
    validObject(probeAnno)
  else {
    envElements <- ls(probeAnno)
    chromElements <- grep("[[:digit:]]+\\..*",envElements, value=TRUE)
    if (length(chromElements)==0){
      warning("Enviroment does not contain any chromosome/strand probe mappings, expecting elements named for example 1.start\n"); return(FALSE)}
    isStrandSpecific <- length(grep("\\.\\+\\.",chromElements))>0 || length(grep("\\.-\\.",chromElements))>0
    sp <- strsplit(chromElements,split="\\.")
    if (isStrandSpecific){
      uniChromNames <- sapply(sp, function(x) paste(x[1:2],collapse="."))
    } else {
      uniChromNames <- sapply(sp, function(x) x[1])
    }
    for (thisName in uniChromNames){
      theseElemNames <- paste(thisName,c("start","end","index","unique"),sep=".")
      theseAreIn <- theseElemNames %in% envElements
      if (!all(theseAreIn)){
        warning(paste("Environment ",deparse(substitute(probeAnno))," seems to hold information for chromosome/strand '",thisName,"', but does not contain elements '",paste(theseElemNames[!theseAreIn],collapse="', '"),"'.\n",sep=""))
        return(FALSE)}
      theseElements <- mget(theseElemNames, env=probeAnno)
      ## check if the objects related to the same chromosome/strand match up
      if (length(unique(sapply(theseElements,length)))!=1){
        print(sapply(theseElements, length))
        warning("Lengths of per chromosome/strand vectors differ!\n")
        return(FALSE)
      }
      ## check if match positions are correctly coded, end positions >= start position
      if (!all(theseElements[[2]] >= theseElements[[1]])){
        warning(paste("Some match positions end before they actually start.\nPlease check elements",paste(thisName,"start",sep="."),"and",paste(thisName,"end",sep="."),".\n"))
        return(FALSE)
      }# if (!all(theseElements[[2]] >= theseElements[[1]]))
      mids <- round((theseElements[[1]]+theseElements[[2]])/2)
      if (is.unsorted(mids)){
        warning(paste("Probe matches are not sorted (in increasing order) by their middle position on chromosome",thisName,"and possible others.\n"))
        return(FALSE)
      }
    }# for (thisName in uniChromNames)
    return(TRUE)
  }
}# validProbeAnno


features2Probes <- function(gff, probeAnno, upstream=5000, checkUnique=TRUE, uniqueCodes=c(0), mem.limit=1e8, verbose=TRUE){
  stopifnot(inherits(gff,"data.frame"), validObject(probeAnno),
            all(c("strand","name","start","end","chr") %in% names(gff)),
            all(gff$strand %in% c(-1,1)), all(gff$start<gff$end))
  ## get borders of upstream region:
  gff$start2 <- ifelse(gff$strand==1, pmax(gff$start-upstream,1), gff$start)
  gff$end2   <- ifelse(gff$strand==1, gff$end, gff$end+upstream)
  realTSS    <- ifelse(gff$strand==1, gff$start, gff$end)
  names(realTSS) <- gff$name
  ## prepare result
  f2p <- vector("list", nrow(gff))
  names(f2p) <- gff$name
  allChr <- intersect(unique(gff$chr), chromosomeNames(probeAnno))
  if (verbose) cat("Chromosome ")
  for (chr in allChr){
    if (verbose) cat(chr,"... ")
    chrsta <- probeAnno[paste(chr,"start",sep=".")]
    chrend <- probeAnno[paste(chr,"end",sep=".")]
    chrmid <- round((chrsta+chrend)/2)
    chridx <- probeAnno[paste(chr,"index",sep=".")]
    if (checkUnique){
      chruni <- probeAnno[paste(chr,"unique",sep=".")]
      stopifnot(length(chruni)==length(chridx))
      chridx <- chridx[chruni %in% uniqueCodes]
      chrmid <- chrmid[chruni %in% uniqueCodes]
    }
    if (length(chrmid)==0) next
    chrProbeDf <- data.frame(chr=rep(chr, length(chrmid)), start2=chrmid, end2=chrmid, stringsAsFactors=FALSE)
    chrGff <- subset(gff, chr==chr)
    chrOverlap <- regionOverlap(chrGff, chrProbeDf, startColumn="start2", endColumn="end2", mem.limit=mem.limit)
    idxOverlap <- nonzero(chrOverlap)
    if (nrow(idxOverlap)==0) next
    ## for each overlapping Feature:
    for (i in unique(idxOverlap[,1])){
      fprobes <- (chrmid[idxOverlap[idxOverlap[,1]==i,2]]-realTSS[chrGff$name[i]])*chrGff$strand[i]
      names(fprobes) <- chridx[idxOverlap[idxOverlap[,1]==i,2]]
      f2p[[chrGff$name[i]]] <- fprobes
    }
  }# for (chr in allChr)
  return(f2p)
}# features2Probes

## function to extract the probe match positions from the 'genes' slot of an
##  RGList as is the case with some ChIP-chip microarray platforms, e.g.
##  with Agilent after reading in the data with read.maimages(...,"agilent")
extractProbeAnno <- function(object, format="agilent", ...){
  stopifnot(inherits(object, "RGList"), "genes" %in% names(object))
  format <- match.arg(format, c("agilent"))
  G <- object$genes
  if (format=="agilent"){
    idColumn <- "ProbeName"
    positionColumn <- "SystematicName"
    stopifnot(all(c(idColumn, positionColumn) %in% names(G)))
    idxWithPos <- grep("^chr.+\\:[[:digit:]]+\\-[[:digit:]]+$",G[[positionColumn]])
    if (!all(seq(nrow(G)) %in% idxWithPos))
      warning(paste("Some reporters had no or an unrecognized genome position in ",deparse(substitute(object)),"$genes$", positionColumn,".\n", sep=""))
    G <- G[idxWithPos, , drop=FALSE]
    splitted <- do.call("rbind", strsplit(G[[positionColumn]], split="[:-]"))
    ## TO DO:make G[[idColumn]] unique
    probePos <- data.frame("PROBE_ID" = G[[idColumn]],
                           CHROMOSOME = as.character(splitted[,1]),
                           COORD1 = as.integer(splitted[,2]),
                           COORD2 = as.integer(splitted[,3]),
                           stringsAsFactors=FALSE)
    probePos$POSITION <- with(probePos, pmin(COORD1, COORD2))
    probePos$LENGTH   <- with(probePos, abs(COORD2 - COORD1)+1)
    thisProbeAnno <- posToProbeAnno(probePos, ...)
  }
  return(thisProbeAnno)
}#extractProbeAnno

