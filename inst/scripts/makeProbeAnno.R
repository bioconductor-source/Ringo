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
## 4. Construct probeAnnoReverse and probeAnnoDirect lists
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
hitResultFile   <- "~/projects/tfchip/DesignFiles/2006-11-15_Sperling_mm8_array1.pos"
hitResultFile2   <- "~/projects/tfchip/DesignFiles/2006-11-15_Sperling_mm8_array2.pos"

# what to build:
gffFile       <- "~/groupsvn/projects/Sperling/Jenny/chip/data/mm8gff.rda"
probeAnnoFile <- "ngchip2probeAnno_array2.rda"
# probeAnnoFile <- "~/projects/tfchip/data/ngchip2probeAnno.rda"

# how are the chromosomes referred to?
allChr <- c(1:19,"X","Y")  # mouse musculus! our arrays have no probes matching mtDNA

# what can this script do:
possibleProducts <- c("probeAnno","gff")

# IMPORTANT  what do we want to do now:
what <- possibleProducts[1]

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

  save(probeAnno, file=probeAnnoFile, save=TRUE)
  
}# if("probeAnno" %in% what)



##--------------------------------------------------
## gff
##--------------------------------------------------
if("gff" %in% what) {
  # augment list of transcript details with further info and save as gff

  library("biomaRt")
  ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")  
  trans.ids <- getFeature(type="ensembl_transcript_id", chr=c(1:19,"X","Y"), mart=ensembl)[,2]

  # what information to get for each transcript:
  sel.attributes=c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "strand", "transcript_start", "transcript_end","markersymbol", "description")

  # retreive information:
  gff <- getBM(attributes=sel.attributes, filters="ensembl_transcript_id", value=trans.ids, mart=ensembl)

  martDisconnect(ensembl)
  
  gff$name <- gff$ensembl_transcript_id
  gff$chr <- gff$chromosome_name
  gff$symbol <- gff$marker_symbol
  gff$feature <- rep("transcript",nrow(gff))

  gff$chromosome_name <- gff$marker_symbol <- NULL

  save(gff, file=gffFile, compress=TRUE)

} # make gff
