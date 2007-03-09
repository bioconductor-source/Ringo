##--------------------------------------------------------------------------------
## Use the output from MUMmer or BLAT or Exonerate
## (querying all probes on the chip vs genome)
## to construct probe annotaton datastructures
##
## IMPORTANT: After running MUMmer or BLAT condense their various output files
## into one table using the scripts 'condenseMUMmerOutput.pl' or
## 'condenseBlatOutput.pl' respectively. This condensed table is then used by this
## script!
##
## This script has the following parts
## 1. Read the table containing the MUMmer or BLAT hits
## 2. Construct along-chromosome vectors for probeAnno environment:
##    start, end: (integer) 
##    index: index (1...6553600) of the PM on the chip (integer) 
##    unique: whether the probe hits multiple places (logical)
## 3. Construct GFF dataframe
## 4. Construct probeAnnoReverse and probeAnnoDirect lists
## 5. Write 'probeAnno.rda'
##--------------------------------------------------------------------------------

## Note: This script has been adopted from the 'davidTiling' package. However, since ChIP-chip
##       is not strand-specific, the strand details were omitted.

### PREPARATION #########################################
### check these settings before running the rest of the script!

## a. path to the result table after running the condense*Output.pl script

## FORMAT of that position file should be in NimbleGen's pos format:
# here's an example of their typical head (comment char and one whitespace at start of each line have been added only just here)

# SEQ_ID   CHROMOSOME     PROBE_ID        POSITION        LENGTH  COUNT
# chr1:4483138-4484138    CHR01   CHR01P004483140 4483140 50      1
# chr1:4483138-4484138    CHR01   CHR01P004483245 4483245 50      1
# ...

# what to use:
hitResultFile   <- "~/projects/tfchip/DesignFiles/2006-11-15_Sperling_mm8_array1.pos"
hitResultFile2   <- "~/projects/tfchip/DesignFiles/2006-11-15_Sperling_mm8_array2.pos"

transcriptsFile <- "~/groupsvn/projects/Sperling/chipDesign/doc/mm8EnsemblTranscriptAnno.txt"
genesFile       <- "~/groupsvn/projects/Sperling/chipDesign/doc/ensembl_all_genes_mm8.tsv"

# what to build:
gffFile       <- "~/groupsvn/projects/Sperling/Jenny/chip/data/mm8gff.rda"
probeAnnoFile <- "ngchip2probeAnno_array2.rda"
# probeAnnoFile <- "~/projects/tfchip/data/ngchip2probeAnno.rda"
chrLocObjName <- "ngchip2chrloc"

# how far upstream to consider probes to be located in a genes upstream region
maxUpstream <- 10000   # 10 kb  , only necessary for probeAnnoExtension

# how are the chromosomes referred to?
allChr <- c(1:19,"X","Y")  # no probes match to MT

# what can this script do:
possibleProducts <- c("probeAnno","probeAnnoExtension","chromLocationObj")

# IMPORTANT  what do we want to do now:
what <- possibleProducts[1]

### READ IN ###################################################
library(Ringo)

hits <- read.delim(hitResultFile, header=TRUE, as.is=TRUE)
# do some sorting again to be on the safe side:

if (exists("hitResultFile2")){
  hits2 <- read.delim(hitResultFile2, header=TRUE, as.is=TRUE)
  hits <- rbind(hits, hits2)
}

if (!all(c("SEQ_ID","CHROMOSOME","PROBE_ID","POSITION","LENGTH")%in%names(hits)))
  stop(paste("File",hitResultFile,"does not seem to be a POS file, since at least one of the
required column headers (SEQ_ID,CHROMOSOME,PROBE_ID,POSITION,LENGTH) is missing.\n"))


hitOrder <- order(hits$CHROMOSOME, hits$POSITION)
hits <- hits[hitOrder,]
if (length(grep("chr",hits$CHROMOSOME, ignore.case=TRUE)))
  hits$CHROMOSOME <- gsub("chr(0)?","", hits$CHROMOSOME, ignore.case=TRUE)
probeindex <- 1:nrow(hits)
probeMultiplicity <- table(hits$PROBE_ID)

##--------------------------------------------------
## The along-chromosome data in 'probeAnno'
##--------------------------------------------------
## group PMs by chromosome
if("probeAnno" %in% what) { 
  cat("Making probeanno: ")
  probeAnno = new.env()
  sp <- split(probeindex, as.factor(hits$CHROMOSOME))
  # for every chromosome, assign the hiting probes to the corresponding vector:
  for(i in seq(along=sp)) {
    chromprobes = sp[[i]]
    nm  = names(sp)[i]
    cat(nm, "")
    assign(paste(nm, "start", sep="."),  as.integer(hits[chromprobes, "POSITION"]), envir=probeAnno)
    assign(paste(nm, "end", sep="."),  as.integer(hits[chromprobes, "POSITION"])+as.integer(hits[chromprobes, "LENGTH"])-1 , envir=probeAnno)
    assign(paste(nm, "index", sep="."),  hits[chromprobes, "PROBE_ID"],  envir=probeAnno)
    assign(paste(nm, "unique", sep="."), as.integer(probeMultiplicity[hits$PROBE_ID[chromprobes]]>1)*3,
           envir=probeAnno)
    ## uniqueness codes:  0 = probe has one unique hit;   3= probe has multiple hits
  }
  cat("\n")
}# if("probeAnno" %in% what)



##--------------------------------------------------
## gff
##--------------------------------------------------
if("gff" %in% what) {
  # augment list of transcript details with further info and save as gff
  tss <- read.delim(file=transcriptsFile, header=TRUE, as.is=TRUE, comment.char="", quote="")
  names(tss)[c(5,6)] <- c("start","end")
  tss  <- tss[-grep("^NT_",tss$chr),]
  ntss <- nrow(tss)
  tss$seqname <- paste("chr",tss$chr,sep="")
  tss$feature <- rep("transcript",ntss)

  # for each TSS compute its distance to next upstream gene:
  # load gene information
  geneInfo <- read.delim(genesFile, header=TRUE, as.is=TRUE, comment.char="", quote="")
  geneInfo <- geneInfo[order(geneInfo$chr, geneInfo$start),]
  ## first get strand-dependend real TSS of each transcript:
  realTSS <- ifelse(tss$strand==1, tss$start, tss$end)
  geneInfoRealStart <-  ifelse(geneInfo$strand==1, geneInfo$start, geneInfo$end)
  names(geneInfoRealStart) <- geneInfo$gene
  geneInfoInter <- ifelse(geneInfo$strand==1, c(Inf, geneInfo$start[-1]-geneInfo$end[-nrow(geneInfo)]),
                          c(geneInfo$start[-1]-geneInfo$end[-nrow(geneInfo)],Inf))
  # handle first and last transcript on each strand:
  geneInfochromDiff <- diff(as.numeric(as.factor(geneInfo$chr)))
  geneInfoInter[c(FALSE, geneInfochromDiff!=0) & (geneInfo$strand=="+")] <- Inf
  geneInfoInter[c(geneInfochromDiff!=0,FALSE) & (geneInfo$strand=="-")] <- Inf
  
  # what to do with overlapping transcripts: YET TO DECIDE!
  geneInfoInter[geneInfoInter<0] <- Inf

  # distance between transcript start and gene start:
  transcript2Gene <- abs(realTSS-geneInfoRealStart[tss$gene])
  allChr <- unique(tss$chr)
  tssUpDist <- rep(NA, nrow(tss))
  for (i in 1:nrow(tss)){
    if (i %% 1000 == 0) cat(i, "... ")
    thisGene <- tss$gene[i]
    geneInfoPos <- which(geneInfo$gene==thisGene)
    tssUpDist[i] <- geneInfoInter[geneInfoPos]+transcript2Gene[i]
  }# tssUpDist
  
  tss <- cbind(tss, tssUpDist)
  names(tss)[ncol(tss)] <- "dist2NextUpGene"

  exStarts <- strsplit(tss$exonstarts,",[[:space:]]?")
  exEnds   <- strsplit(tss$exonends,",[[:space:]]?")

  intron1starts <- integer(ntss)
  intron1ends   <- integer(ntss)

  for (i in 1:ntss){
    if (i %% 1000 == 0) cat(i,"... ")
    theseStarts <- as.numeric(exStarts[[i]])
    theseEnds   <- as.numeric(exEnds[[i]])
    if (length(theseStarts)<2){
      intron1starts[i] <- NA; intron1ends[i] <- NA; next}
    stopifnot(!is.unsorted(theseStarts), !is.unsorted(theseEnds))
    if (tss$strand[i] > 0){
      intron1starts[i] <- theseEnds[1] + 1
      intron1ends[i]   <- theseStarts[2] - 1
    } else { # minus strand
      intron1starts[i] <- rev(theseEnds)[2] + 1
      intron1ends[i]   <- rev(theseStarts)[1]-1
    }
  }#for

  gff <- tss
  
  gff$intron1start <- intron1starts
  gff$intron1end   <- intron1ends

  geneSyms <- vector("character", ntss)
  geneDescs <- vector("character",ntss)
  ## get additional gene annotation using biomaRt:
  mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
  martRes <- getGene(tss$gene,type="ensembl",mart=mart)

  geneDescs <- martRes$description[!duplicated(martRes$ID)]
  names(geneDescs) <- martRes$ID[!duplicated(martRes$ID)]
  gff$description <- geneDescs[gff$gene]

  martResSyms <- split(martRes$symbol, martRes$ID)
  geneSyms    <- sapply(martResSyms, function(thisSym)
                        return(thisSym[which.min(nchar(thisSym)[nchar(thisSym)>0])]))
  gff$symbol  <- sapply(geneSyms[gff$gene], function(x) if (length(x)==0) return(NA) else return(x[1]))

  names(gff)[names(gff)=="transcript"] <- "name"
  
  # GO annotation
  library(biomaRt)
  
  mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
  getGO("ENSMUSG00000038193",type="ensembl",mart=mart)
  goRes <- getGO(unique(tss$gene),type="ensembl",mart=mart)

  save(goRes, file="~/temp/goRes.RData")
  goRes2 <- split(goRes$go_id, goRes$ensembl_gene_id)
  goRes3 <- lapply(goRes2, function(x) paste(x, collapse=","))
  sum(! gff$gene %in% names(goRes2))
  
  gene2go <- goRes3[gff$gene]
  gff$go <- gene2go

  ## different setup: use library 'annQueR'
  if (require(annQueR)){
    mm8Genes <- getAllGenes(host="ensembldb.ensembl.org",dbname = "mus_musculus_core_39_36")
    str(mm8Genes)
    mm8Transcripts <- getTranscripts(mm8Genes, host="ensembldb.ensembl.org",dbname = "mus_musculus_core_39_36")
    str(mm8Transcripts)
    mm8Ext <- getFurtherAnnotation(mm8Genes)  
    mm8Transcripts$entrezgene <- as.character(mm8Ext[mm8Transcripts$gene,"entrezgene"])
    mm8Transcripts$refseq <- as.character(mm8Ext[mm8Transcripts$gene,"refseq_dna"])
    mm8Transcripts$symbol <- as.character(mm8Ext[mm8Transcripts$gene,"mgi_symbol"])
    mm8GO <- getGOids(mm8Genes)
    mm8Transcripts$go <- unlist(mm8GO[mm8Transcripts$gene])
    gff <- mm8Transcripts
    gff$feature <- rep("transcript",nrow(gff))
    gff$seqname <- paste("chr",gff$chr,sep="")
    names(gff) <- gsub("mgi_symbol","symbol",names(gff))
    names(gff) <- gsub("transcript","name",names(gff))
  }# if (require(annQueR))

  #gffFile = "~/temp/mm8gff.rda"
  
  save(gff, file=gffFile, compress=TRUE)

  
  # make a condensed form of the gff with one row for each gene only:
  uniGene <- unique(gff$gene)
  nGenes <- length(uniGene)
  cat(nGenes, "genes to process:\n")
  gff2   <- list(chr=rep("NA",nGenes), start=integer(nGenes), end=integer(nGenes),
                 strand=rep("-",nGenes), gene=uniGene)
  for (i in 1:nGenes){
    if (i %% 1000 == 0) cat(i,"... ")
    g1 <- gff[gff$gene == uniGene[i],]
    gff2$chr[i] <- g1$chr[1]
    gff2$start[i] <- min(g1$start)
    gff2$end[i] <- max(g1$end)
    gff2$strand[i] <- g1$strand[1]
  }#for all unique genes
  cat("\n")
}

##--------------------------------------------------
## Part 5: the per-probe data in 'probeAnno'
##--------------------------------------------------
### takes quite long unfortunately
if("probeAnnoExtension" %in% what) {
  uniProbe <- unique(hits$PROBE_ID); nProbes <- length(uniProbe)
  probeInGene <- vector("list",nProbes)
  probeUpstreamGene <- vector("list",nProbes)
  names(probeInGene) = names(probeUpstreamGene) <- uniProbe
  
  #featNames = c("transcript", "CDS", "ncRNA", "nc_primary_transcript","rRNA", "snRNA", "snoRNA", "tRNA")
  #selGff  = which(gff$feature %in% featNames)  
  sgff    = gff2
  
  stopifnot(nrow(sgff)>1000)
  stopifnot(all(sgff$strand %in% c("+", "-")))  

  for(chr in allChr) {
    ifeat  = which(sgff$chr == chr)
    nfeat <- length(ifeat)
    cat(chr, ": ", nfeat, " features\n", sep="")
    chrfeatstart   = sgff$start[ifeat]
    chrfeatend   = sgff$end[ifeat]
    stopifnot(all(chrfeatend>chrfeatstart))
    chrfeatstrand = sgff$strand[ifeat]
    Name = sgff$gene[ifeat]
    stopifnot(all(chrfeatstart<=chrfeatend))
    # left(1) and right(2) border of transcripts' upstream region
    chrfeatupstream1 <- ifelse(chrfeatstrand=="+",chrfeatstart-maxUpstream, chrfeatend+1)
    chrfeatupstream2 <- ifelse(chrfeatstrand=="+",chrfeatstart-1, chrfeatend+maxUpstream)
    
    sta = get(paste(chr, "start", sep="."), envir=probeAnno)
    end = get(paste(chr, "end",   sep="."), envir=probeAnno)
    mid <- round((sta+end)/2)
    ind = get(paste(chr, "index", sep="."), envir=probeAnno)
    nProbesOnChr <- length(mid)

    combRank <- rank(c(mid,chrfeatstart,chrfeatend,chrfeatupstream1,chrfeatupstream2),tie="first")
    midRank    <- combRank[1:nProbesOnChr];
    names(ind) <- as.character(midRank)
    featstartRank <- combRank[(nProbesOnChr)+(1:nfeat)];
    featendRank <- combRank[(nProbesOnChr+nfeat)+(1:nfeat)];
    featup1Rank <- combRank[(nProbesOnChr+nfeat*2)+(1:nfeat)];
    featup2Rank <- combRank[(nProbesOnChr+nfeat*3)+(1:nfeat)];

    selfeatIn <- which(featstartRank + 1 != featendRank)
    selfeatUp <- which(featup1Rank + 1 != featup2Rank)

    ## first find those probes with  ranks inside each feature
    ranksInsideFeatsList <- sapply(selfeatIn, function(i){
      inFeatRanks <- (featstartRank[i]+1):(featendRank[i]-1)
      inFeatRanks <- inFeatRanks[inFeatRanks %in% midRank]
      return(as.character(inFeatRanks))})
    names(ranksInsideFeatsList) <- Name[selfeatIn]
    ranksInsideFeatsList <- ranksInsideFeatsList[listLen(ranksInsideFeatsList)>0]
    # reverse to get a ranks-to-features relation
    if (length(ranksInsideFeatsList)>0){
      featsIncludingRanks <- reverseSplit(ranksInsideFeatsList)
      probeInGene[ind[names(featsIncludingRanks)]] <-
        mapply(c, probeInGene[ind[names(featsIncludingRanks)]], featsIncludingRanks, SIMPLIFY=FALSE)
    }
    ## second find those probes with  ranks upstream each feature
    ranksUpstreamFeatsList <- sapply(selfeatUp, function(i){
      upFeatRanks <- (featup1Rank[i]+1):(featup2Rank[i]-1)
      upFeatRanks <- upFeatRanks[upFeatRanks %in% midRank]
      return(as.character(upFeatRanks))})
    names(ranksUpstreamFeatsList) <- Name[selfeatUp]
    ranksUpstreamFeatsList <- ranksUpstreamFeatsList[listLen(ranksUpstreamFeatsList)>0]
    # reverse to get a ranks-to-features relation
    if (length(ranksUpstreamFeatsList)>0){
      featsUpstreamRanks <- reverseSplit(ranksUpstreamFeatsList)
      probeInGene[ind[names(featsUpstreamRanks)]] <-
        mapply(c, probeUpstreamGene[ind[names(featsUpstreamRanks)]], featsUpstreamRanks, SIMPLIFY=FALSE)
    }
  } ## for chr
  
  cat("\n")
  assign("probeInGene", probeInGene, envir=probeAnno)
  assign("probeUpstreamGene", probeUpstreamGene, envir=probeAnno)
  
}# if "probeAnnoExtended"

# ---------------------------------
#  Object of class 'chromLocation'
# ---------------------------------
if ("chromLocationObj" %in% what){
  probemiddle <- as.integer(paste(hits$Strand, rowMeans(hits[,c("Start","End")]),sep=""))
  names(probemiddle) <- hits$PROBE_ID
  chrLocs <- split(probemiddle, as.factor(hits$CHROMOSOME))
  names(chrLocs) <- gsub("chr","", names(chrLocs))
  probeList <- split(hits$CHROMOSOME,as.factor(hits$PROBE_ID))
  probesToChrom <- new.env(hash=TRUE)
  multiassign(probeList,  env=probesToChrom)
  # get chromosome lengths:
  library("MmusculusGenome.mm7")
  source(system.file("R","MmusculusGenome.mm7",package="MmusculusGenome.mm7"),echo=T)
  uniChr <- unique(hits$CHROMOSOME)
  if (!("chr" %in% uniChr)) uniChr <- paste("chr",uniChr, sep="")
  chrLengths <- sapply(uniChr, function(x) get(x,mm7@data_env)@last)
  rm(mm7); detach(package:MmusculusGenome.mm7);  gc()
  allID <- names(probeList)
  geneSym <- new.env(hash=TRUE)
  multiassign(allID, value=allID, env=geneSym)
  chrlocobj <- new("chromLocation",
                   organism="Mus musculus", dataSource="ngchipchip1",
                   chromLocs=chrLocs, probesToChrom=probesToChrom,
                   chromInfo=chrLengths, geneSymbols=geneSym)
  assign(chrLocObjName, chrlocobj)
}# if ("chromLocationObj" %in% what)

##
##  save
##
toSave <- (c(possibleProducts,chrLocObjName))
save(list=toSave[toSave %in% ls()], file=probeAnnoFile, compress=TRUE)

