### July 2007                  Joern Toedling

## this is a small convenience script that should provide users with ideas how
##  one can easily construct the tab-separated "targets" file that is needed
##  by readNimblegen from the file SampleKey.txt that is often provided by
##  NimbleGen.
##  You are welcome and hereby encouraged to modify this code
##  according to your own needs.


## how the targets file could look like:
exDir <- system.file("exData",package="Ringo")
read.delim(file.path(exDir,"example_files.txt"), header=TRUE)
### Important are the columns "FileNameCy3", "FileNameCy5", "Cy3" and "Cy5"

## convert SampleKey.txt into targets.txt
sample.key = read.delim("SampleKey.txt",header=TRUE, as.is=TRUE)

uni.slides <- unique(sample.key$"CHIP_ID")
cy3part <- subset(sample.key, DYE=="Cy3")
cy5part <- subset(sample.key, DYE=="Cy5")

cy3order <- match(uni.slides, cy3part$"CHIP_ID")
cy5order <- match(uni.slides, cy5part$"CHIP_ID")
stopifnot(all.equal(cy3part$"CHIP_ID"[cy3order], cy5part$"CHIP_ID"[cy5order]))

slides <- data.frame(SlideNumber=uni.slides, FileNameCy3=paste(uni.slides,"_532.pair.txt", sep=""), FileNameCy5=paste(uni.slides,"_635.pair.txt", sep=""), Species=make.names(cy3part$"SAMPLE_SPECIES"[cy3order]), Cy3=make.names(cy3part$"SAMPLE_DESCRIPTION"[cy3order]),  Cy5=make.names(cy5part$"SAMPLE_DESCRIPTION"[cy5order]))

### write slides
write.table(slides, file="slide_files.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=TRUE)
###  check if readNimblegen could read it:
targets <- readTargets("slide_files.txt")
