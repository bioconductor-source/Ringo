
setClass("probeAnno",representation(map="environment", arrayName="character", genome="character"))

setMethod("initialize", "probeAnno", function(.Object, map=new.env(), arrayName="", genome=""){
  if (missing(map)) map <- new.env()
  stopifnot(is.environment(map))
  .Object@map <- map
  stopifnot(validObject(.Object))
  .Object
})# initialize probeAnno

setValidity("probeAnno", function(object){
  map <- object@map
  envElements <- ls(map)
  if (length(envElements)==0){
    return(TRUE)}
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
      warning(paste("Environment ",deparse(substitute(map))," seems to hold information for chromosome/strand '",thisName,"', but does not contain elements '",paste(theseElemNames[!theseAreIn],collapse="', '"),"'.\n",sep=""))
      return(FALSE)}
    theseElements <- mget(theseElemNames, env=map)
    ## check if the maps related to the same chromosome/strand match up
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
})# validity of probeAnno


setMethod("[",signature(x="probeAnno"),
   function(x, i, j="missing", ..., drop="missing"){
     stopifnot(is.character(i), length(i)==1)
     if (!exists(i, envir=x@map))
       stop(paste("No mapping '",i,"' in this 'probeAnno' object.\n", sep=""))
     get(i, env=x@map)
})

setReplaceMethod("[",signature(x="probeAnno"),
   function(x, i, j="missing", value){
     stopifnot(is.character(i), length(i)==1)
     assign(x=i, value=value, envir=x@map)
     x
})

setMethod("get", signature(x="character", pos="missing", envir="probeAnno", mode="missing", inherits="missing"),  function( x, pos="missing", envir, mode="missing", inherits="missing"){
  stopifnot(is.character(x), length(x)==1)
  if (!exists(x, envir=envir@map))
    stop(paste("No mapping '",x,"' in this 'probeAnno' object.\n", sep=""))
  get(x, envir=envir@map)
})

setMethod("ls", signature(name="probeAnno", pos="missing", envir="missing", all.names="missing", pattern="missing"), function(name, pos="missing", envir="missing", all.names="missing", pattern="missing") ls(name@map))

setGeneric("chromosomeNames", function(x) standardGeneric("chromosomeNames"))

setMethod("chromosomeNames", signature(x="probeAnno"), function(x){
  envElements <- ls(x@map)
  chromElements <- grep("[[:digit:]]+\\..*",envElements, value=TRUE)
  if (length(chromElements)==0){
    warning("Enviroment does not contain any chromosome/strand probe mappings,\n  expecting elements named, e.g., '1.start'\n"); return(vector("character",0))}
  sp <- strsplit(chromElements,split="\\.")
  uniChromNames <- unique(sapply(sp, function(x) x[1]))
  return(uniChromNames)
})

setGeneric("arrayName", function(x) standardGeneric("arrayName"))

setMethod("arrayName", signature(x="probeAnno"), function(x){
  x@arrayName})

setGeneric("arrayName<-", function(x, value) standardGeneric("arrayName<-"))

setReplaceMethod("arrayName", signature(x="probeAnno", value="character"),
  function(x, value){x@arrayName <- value; return(x)})

setGeneric("genome", function(x) standardGeneric("genome"))

setMethod("genome", signature(x="probeAnno"), function(x){
  x@genome})

setGeneric("genome<-", function(x, value) standardGeneric("genome<-"))

setReplaceMethod("genome", signature(x="probeAnno", value="character"),
  function(x, value){x@genome <- value; return(x)})


setMethod("show", signature(object="probeAnno"), function(object){
  stopifnot(validObject(object))
  cat("A 'probeAnno' object holding the mapping between\nreporters and genomic positions.\nChromosomes: ")
  cat(chromosomeNames(object)); cat("\n")
  cat("Microarray platform:", object@arrayName,"\nGenome:",object@genome)
  cat("\n")
  invisible(NULL)
})

setAs(from="environment", to="probeAnno", def=function(from){
  new("probeAnno", map=from)
})
