
exportCherList <- function(object, filename="chers.gff", format="gff3",...)
{
  require("rtracklayer")
  stopifnot(inherits(object, "cherList"))
  chersDf <- as.data.frame.cherList(object)
  chersDf$strand <- NA
  chersDf$chrom <- paste("chr", sub("^chr", "", chersDf$chr), sep = "")
  chersDf$start <- as.integer(as.character(chersDf$start))
  chersDf$end <- as.integer(as.character(chersDf$end))
  chersAdf <- as(chersDf, "AnnotatedDataFrame")
  descs <- c(chrom = "Chromosome ID",
             start = "Start position on chromosome",
             end = "End position on chromosome",
             strand = "DNA strand, sense (+) or antisense (-)")
  md <- data.frame(labelDescription = descs, stringsAsFactors = FALSE)
  varMetadata(chersAdf)[rownames(md),] <- md
  chersTS <- trackSet(chersAdf, ...)
  export(chersTS, filename, format=format)
  invisible(NULL)
}
