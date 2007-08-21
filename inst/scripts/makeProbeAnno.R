##--------------------------------------------------------------------------------
## Use NimbleGen's POS file or the output from MUMmer or BLAT or Exonerate
## (querying all probes on the chip vs genome)
## to construct probe annotaton datastructures
##
## IMPORTANT: After running MUMmer or BLAT condense their various output files
## into one table using the scripts 'condenseMUMmerOutput.pl' or
## 'condenseBlatOutput.pl' respectively. This condensed table is simlilar to
## a NimbleGen POS file and usable by this script!
##
## This script has the following parts
## 1. Read the table containing the POS file or MUMmer or BLAT hits
## 2. Construct along-chromosome vectors for probeAnno environment:
##    start, end: (integer) 
##    index: index (1...6553600) of the PM on the chip (integer) 
##    unique: whether the probe hits multiple places (logical)
## 3. Construct GFF dataframe
## 4. optional: compute probes' GC content from NDF file supplied by NimbleGen
## 5. Write 'probeAnno.rda'
##--------------------------------------------------------------------------------

## Note: This script has been adopted from the 'davidTiling' package.
## However, since ChIP-chip is not strand-specific,
##  the strand details were omitted.

## a. path to the result table after running the condense*Output.pl script

## FORMAT of that position file should be in NimbleGen's pos format:
# here's an example of their typical head (comment char and one whitespace at start of each line have been added only just here)

# SEQ_ID   CHROMOSOME     PROBE_ID        POSITION        LENGTH  COUNT
# chr1:4483138-4484138    CHR01   CHR01P004483140 4483140 50      1
# chr1:4483138-4484138    CHR01   CHR01P004483245 4483245 50      1
# ...

# what to use:
hitResultFile   <- "MOD_2003-12-05_SUZ12_1in2.pos"
#hitResultFile2   <- "MOD_2003-12-05_SUZ12_2in2.pos"

# what to build:
gffFile       <- "ensembl_gff.rda"
probeAnnoFile <- "exampleProbeAnno.rda"

# how are the chromosomes referred to?
# allChr <- c(1:22,"X","Y")  # homo sapiens
allChr <- "9" # the example data only contains information concerning reporters mapped to chromosome 9

# what can this script do:
possibleProducts <- c("probeAnno","gff","probegc")

# IMPORTANT  what do we want to do now:
what <- possibleProducts[c(1,2,3)]

### READ IN ###################################################

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

  if (probegc %in% what) {
    ### add the GC content of each probe on the used arrays to the probeAnno environment
    list.files(pattern="\\.ndf")
    ndf.file.name <- "MOD_2003-12-05_SUZ12_1in2.ndf"
    ndf <- read.delim(ndf.file.name, header=TRUE, as.is=TRUE)
    stopifnot(all(c("PROBE_ID","PROBE_SEQUENCE") %in% names(ndf)))
    are.unique <- !duplicated(ndf[["PROBE_ID"]])
    probes.gc <- compute.gc(ndf[["PROBE_SEQUENCE"]][are.unique])
    # my take some minutes
    names(probes.gc) <- ndf[["PROBE_ID"]][are.unique]
    assign("probes.gc", value=probes.gc, env=probeAnno)
    rm(ndf, ndf2, ndf.file.name1, ndf.file.name2);gc()
  }# if (probegc %in% what)
  
  save(probeAnno, file=probeAnnoFile)
}# if("probeAnno" %in% what)

##--------------------------------------------------
## gff
##--------------------------------------------------
if("gff" %in% what) {
  # augment list of transcript details with further info and save as gff

  library("biomaRt")
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

  gene.ids <- unique(unlist(lapply(as.list(allChr), function(this.chr) getFeature(type="ensembl_gene_id", chr=this.chr, mart=ensembl)[,2]), use.names=FALSE))

  # what information to get for each transcript:
  sel.attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "strand", "start_position","end_position", "transcript_start", "transcript_end", "description")
  # if 'hgnc_symbol' is not a defined attribute for gene symbols in your data species, try to find an equivalent, using commands such as this one:
  grep("symbol",listAttributes(ensembl)[,1], value=TRUE)

  # retreive information:
  gff <- getBM(attributes=sel.attributes, filters="ensembl_gene_id", value=gene.ids, mart=ensembl)

  martDisconnect(ensembl)

  ## replace attribute names by standardized names
  gff$gene <- gff$"ensembl_gene_id"
  gff$name <- gff$ensembl_transcript_id
  gff$chr <- gff$chromosome_name
  gff$symbol <- gff$"hgnc_symbol"
  gff$start <- gff$transcript_start
  gff$end <- gff$transcript_end
  gff$"gene_start" <- gff$"start_position"
  gff$"gene_end" <- gff$"end_position"
  gff$feature <- rep("transcript",nrow(gff))
  gff$chromosome_name <- gff$markersymbol <- NULL
  gff$transcript_start <- gff$transcript_end <- NULL
  gff$start_position <- gff$end_position <- NULL
  gff$ensembl_transcript_id <- NULL
  gff$ensembl_gene_id <- NULL

  
  ### some transcripts names my occurr in multiples, usually this is because an Ensembl transcript can have more than one Symbol defined for it:
  if (any(duplicated(gff$name))){
    dupl.trans <- unique(gff$name[duplicated(gff$name)])
    G <- lapply(as.list(dupl.trans), function(this.trans){
      this.gff <- subset(gff,name == this.trans)
      if (nrow(unique(this.gff[,c("gene","name","chr","start","end","description")]))>1) return(this.gff[1,,drop=FALSE])
      non.zero.gff <- subset(this.gff, nchar(symbol)>0)
      this.other.sym <- NULL
      if (nrow(non.zero.gff)> 0){
        this.new.sym <- non.zero.gff$symbol[which.min(nchar(non.zero.gff$symbol))]
        if (nrow(non.zero.gff)>1)
          this.other.sym <- paste("Synonyms",paste(non.zero.gff$symbol[-which.min(nchar(non.zero.gff$symbol))],collapse=","),sep=":")
      } else { this.new.sym <- "" }
      this.gff$symbol[1] <- this.new.sym
      if (!is.null(this.other.sym))
        this.gff$description[1] <- paste(this.gff$description[1],this.other.sym,sep=";")
      return(this.gff[1,,drop=FALSE])
    })
    gff <- rbind(gff[-which(gff$name %in% dupl.trans),],do.call("rbind",G))
    rm(G, dupl.trans);gc()
  }#if (any(duplicated(gff$name)))

  
  # reorder by chromosome and start position
  gff <- gff[order(gff$chr, gff$gene_start, gff$start),c("name","gene","chr","strand","gene_start","gene_end","start","end","symbol", "description","feature")]
  rownames(gff) <- NULL
  
  save(gff, file=gffFile, compress=TRUE)

} # make gff
