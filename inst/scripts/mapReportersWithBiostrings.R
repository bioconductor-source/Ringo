## June 2008                            Joern Toedling

## use Biostrings (version >= 2.0) for mapping reporters to the genome
##  shown on the reporters of the example data.

### these functions were adopted from the GenomeSearching vignette of
###  Biostrings and slightly modified

library("Ringo")
library("Biostrings")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg18")

#### FUNCTIONS:

## 1. for mapping with mismatches
doMapping <- function(stringSet, chromosomes, max.mismatch=1)
{
    seqnames <- seqnames(Hsapiens)
    if (!missing(chromosomes))
      seqnames <- seqnames[match( chromosomes, seqnames)]
    stopifnot(length(seqnames)!=0)
    res <- vector("list", length(seqnames))
    names(res) <- seqnames
    for (seqname in seqnames) {
      subject <- Hsapiens[[seqname]]
      cat("\n>>> Finding all hits in chromosome", seqname, "...\n")
      chrRes <- vector("list", length(stringSet))
      for (i in seq_len(length(stringSet))) {
        if (i %% 100 == 0) cat(i, " ")
        patternID <- names(stringSet)[i]
        pattern <- stringSet[[i]]
        plus.matches <- matchPattern(pattern, subject, max.mismatch=max.mismatch)
        rcpattern <- reverseComplement(pattern)
        minus.matches <- matchPattern(rcpattern, subject, max.mismatch=max.mismatch)
        nMatches <- length(plus.matches) + length(minus.matches)
        patRes   <- data.frame("CHROMOSOME"=rep.int(seqname, nMatches), 
                               "PROBE_ID"=rep.int(patternID, nMatches),
                               "POSITION"=c(start(plus.matches), start(minus.matches)),
                               "END"=c(end(plus.matches), end(minus.matches)),
                               "STRAND"=rep(c("+","-"), c(length(plus.matches), length(minus.matches))),
                               stringsAsFactors=FALSE)
        chrRes[[i]] <- patRes
      }#for
      chrRes <- do.call("rbind", chrRes)
      res[[seqname]] <- chrRes
      unload(Hsapiens, seqname)
    }# for (seqname in seqnames)
    res <- do.call("rbind",res)
    stopifnot(all(res$END >= res$POSITION))
    res$LENGTH <- res$END - res$POSITION + 1
    res$END <- NULL
    return(res)
}# doMapping

### 2. for only exact reporter to genome mappings. This function requires
##  that all reporter sequences are of the same length and it does not
##  allow for mismatches. Is is MUCH faster than the previous one.
doExactMapping <- function(stringSet, chromosomes)
{
    seqnames <- seqnames(Hsapiens)
    if (!missing(chromosomes))
      seqnames <- seqnames[match( chromosomes, seqnames)]
    stopifnot(length(seqnames)!=0)
    res <- vector("list", length(seqnames))
    names(res) <- seqnames
    
    pdict <- PDict(stringSet)
    mdict <- PDict(reverseComplement(stringSet))
    for (seqname in seqnames) {
      subject <- Hsapiens[[seqname]]
      cat("\n>>> Finding all hits in chromosome", seqname, "...\n")
      plus.mindex <- matchPDict(pdict, subject)
      plus.matches <- extractAllMatches(subject=subject, plus.mindex)
      minus.mindex <- matchPDict(mdict, subject)
      minus.matches <- extractAllMatches(subject=subject, minus.mindex)
      nMatches <- length(plus.matches) + length(minus.matches)
      chrRes   <- data.frame("CHROMOSOME"=rep.int(seqname, nMatches), 
                             "PROBE_ID"=c(names(plus.matches), names(minus.matches)),
                             "POSITION"=c(start(plus.matches), start(minus.matches)),
                             "END"=c(end(plus.matches), end(minus.matches)),
                             "STRAND"=rep(c("+","-"), c(length(plus.matches), length(minus.matches))),
                             stringsAsFactors=FALSE)
    }#for
    res[[seqname]] <- chrRes
    unload(Hsapiens, seqname)
    res <- do.call("rbind",res)
    stopifnot(all(res$END >= res$POSITION))
    res$LENGTH <- res$END - res$POSITION + 1
    res$END <- NULL
    return(res)
}# doExactMapping

exDir <- system.file("exData", package="Ringo")

ndf <- read.delim(file.path(exDir,"MOD_2003-12-05_SUZ12_1in2.ndf"), header=TRUE, as.is=TRUE)

## get the reporter sequences and set their unique identifiers:
repseqs <- DNAStringSet(ndf$"PROBE_SEQUENCE")
stopifnot(!any(duplicated(ndf$"PROBE_ID")))
names(repseqs) <- ndf$"PROBE_ID"

### have a look at seqnames in BSgenome package for H. sapiens
seqnames(Hsapiens)

## this takes some hours:
# hits <- doMapping(repseqs, "chr9")

### alternative for only exact matches and length of all reporters equal:
hits <- doExactMapping(repseqs, "chr9")

### compare these hits to the ones supplied in the POS file
pos <- read.delim(file.path(ringoExampleDir,"MOD_2003-12-05_SUZ12_1in2.pos"), header=TRUE, as.is=TRUE)

chits <- with(hits, paste(PROBE_ID, CHROMOSOME, POSITION, LENGTH, sep="|"))
cpos <- with(pos, paste(PROBE_ID, CHROMOSOME, POSITION, LENGTH, sep="|"))

length(intersect(chits, cpos))
# [1] 660  ## not that great

## build a probeAnno object out of this data.frame
exProbeAnno <- posToProbeAnno(hits, genome="H.sapiens (hg18)", microarrayPlatform="NimbleGen MOD_2003-12-05_SUZ12_1in2")
