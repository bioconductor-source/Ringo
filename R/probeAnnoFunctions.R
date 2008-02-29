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

posToProbeAnno <- function(pos, chrNameColumn="CHROMOSOME",probeColumn="PROBE_ID", chrPositionColumn="POSITION", lengthColumn="LENGTH", verbose=TRUE, ...){
  if (!exists(deparse(substitute(pos)))){
    stopifnot(is.character(pos), file.exists(pos))
    hits <- read.delim(pos, header=TRUE, as.is=TRUE, ...)
  } else {
    hits <- as.data.frame(pos)
  }
  columnsAreIn <- c(chrNameColumn, probeColumn, chrPositionColumn, lengthColumn) %in% names(hits)
  if (!all(columnsAreIn))
    stop(paste("\nColumn(s)",paste(c(chrNameColumn, probeColumn, chrPositionColumn, lengthColumn)[!columnsAreIn],sep=", ", collapse=", "),"are not found in file or data.frame",deparse(substitute(pos)),".\n"))
  hits$END    <- hits[[chrPositionColumn]] + hits[[lengthColumn]] - 1
  hits$MIDDLE <- round((hits[[chrPositionColumn]] + hits$END)/2)
  hitOrder <- order(hits[[chrNameColumn]], hits$MIDDLE)
  hits <- hits[hitOrder,]
  if (length(grep("chr",hits[[chrNameColumn]], ignore.case=TRUE)))
    hits[[chrNameColumn]] <- gsub("chr(0)?","", hits[[chrNameColumn]], ignore.case=TRUE)
  probeindex <- 1:nrow(hits)
  probeMultiplicity <- table(hits[[probeColumn]])
  if (verbose) cat("Creating probeAnno environment for chromosome")
  probeAnno = new.env()
  sp <- split(probeindex, as.factor(hits[[chrNameColumn]]))
  # for every chromosome, assign the hiting probes to the corresponding vector:
  for(i in seq(along=sp)) {
    chromprobes = sp[[i]]
    nm  = names(sp)[i]
    cat(nm, "")
    assign(paste(nm, "start", sep="."),  as.integer(hits[chromprobes, chrPositionColumn]), envir=probeAnno)
    assign(paste(nm, "end", sep="."),  as.integer(hits[chromprobes,"END"]) , envir=probeAnno)
    assign(paste(nm, "index", sep="."),  hits[chromprobes, probeColumn],  envir=probeAnno)
    assign(paste(nm, "unique", sep="."), as.integer(probeMultiplicity[hits[[probeColumn]][chromprobes]]>1)*3, envir=probeAnno)
    ## uniqueness codes:  0 = probe has one unique hit;   3= probe has multiple hits
  }# for(i in seq(along=sp))
  if (verbose) cat("Done.\n")
  return(probeAnno)
}#posToProbeAnno

validProbeAnno <- function(probeAnno){
  stopifnot(is.environment(probeAnno))
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
}# validProbeAnno
