## reimplement class "cher", or more precisely "ChIP-enriched region" from package Ringo in S4:
setClass("cher", representation(name="character", chromosome="character", start="integer", end="integer", cellType="character", antibody="character", maxLevel="numeric", score="numeric", probes="vector", extras="list"))

setMethod("initialize", "cher", function(.Object, name, chromosome, start, end, cellType=NA, antibody, maxLevel=NA, score=NA, probes=vector("character",0), extras=list(),...){
  .Object@name <- as.character(name)
  .Object@chromosome <- as.character(chromosome)
  .Object@start <- as.integer(start)
  .Object@end <- as.integer(end)
  .Object@cellType <- as.character(cellType)
  .Object@antibody <- as.character(antibody)
  .Object@maxLevel <- as.numeric(maxLevel)
  .Object@score <- as.numeric(score)
  .Object@probes <- probes
  further.args <- lapply(as.list(match.call(expand.dots=FALSE)[["..."]]),eval)
  if (length(further.args)>0)
    extras <- c(extras, further.args)
  .Object@extras <- extras
  validObject(.Object)
  .Object
})# initialize cher

setValidity("cher", function(object){
  if (any(length(object@name)!=1, length(object@chromosome)!=1, length(object@start)!=1, length(object@end)!=1, length(object@antibody)!=1, length(object@maxLevel)!=1, length(object@score)!=1)){
    warning("Slots 'name','chromosome','start','end','antibody','maxLevel',and 'score' may each only contain one element.\n"); return(FALSE)}
  if (object@start > object@end){
    warning("Region 'end' is smaller than 'start'.\n");return(FALSE)}
  return(TRUE)
})# set validity of cher objects

setMethod("show",signature="cher", function(object){
  cat(object@name, "\nChr",object@chromosome,":",object@start, "-", object@end, "\n")
  if (!is.null(object@antibody)) cat("Antibody :",object@antibody,"\n")
  if (!is.na(object@maxLevel)) cat("Maximum level =",object@maxLevel,"\n")
  if (!is.na(object@score)) cat("Score =",object@score,"\n")
  if (length(object@probes)>0) cat("Spans",length(object@probes),"probes.\n")
  if (length(object@extras)>0) cat("Defined extras:", paste(names(object@extras), collapse=", "),"\n")
  invisible(NULL)
})

setMethod("plot",signature=c("cher","missing"), function(x){
  stop("For plotting objects of class 'cher', you need to provide at least\n 1) the 'cher' object\n 2) an ExpressionSet that holds the data the region was found on\n 3) a probeAnno object.\n")
})

setMethod("plot",signature=c("cher","ExpressionSet"), function(x, y, probeAnno, samples=NULL, extent=1000, gff=NULL, ...){
  if (!is.null(samples)&is.character(samples))
    samples <- match(samples, sampleNames(y))
  if (!is.null(samples)&is.numeric(samples))
    stopifnot(all(samples) %in% 1:ncol(y))
  cherVal <- chipAlongChrom(y, chrom=x@chromosome, samples=samples, xlim=c(x@start-extent, x@end+extent), probeAnno=probeAnno, gff=gff, ...)
  rug(x=c(x@start,x@end), side=3, lwd=3, col="gold")
  if (!is.null(x@maxLevel))
    legend(x="topright", legend=paste("Max.Level:",round(x@maxLevel,digits=2)),fill="gold", bty="n")
  invisible(cherVal)
})

setMethod("update", signature=c("cher"), function(object, ...){
  further.args <- lapply(as.list(match.call(expand.dots=FALSE)[["..."]]),eval, envir=parent.frame(1))
  for (this.name in intersect(names(further.args), slotNames(object))){
    slot(object, this.name) <- further.args[[this.name]]
  }
  if (is.element("extras",slotNames(object)))
    extras <- object@extras
  else
    extras <- list()
  for (this.name in setdiff(names(further.args), slotNames(object))){
    extras[[this.name]] <- further.args[[this.name]]
  }
  object@extras <- extras
  return(object)
})
