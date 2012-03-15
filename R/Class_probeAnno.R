
setClass("probeAnno",representation(map="environment", arrayName="character", genome="character"))

setMethod("initialize", "probeAnno", function(.Object, map, arrayName="", genome=""){
  if (missing(map)) map <- new.env(hash=TRUE)
  stopifnot(is.environment(map))
  .Object@map <- map
  .Object@arrayName <- arrayName
  .Object@genome <- genome
  stopifnot(validObject(.Object))
  .Object
})# initialize probeAnno

setValidity("probeAnno", function(object){
  fails <- character(0)
  map <- object@map
  envElements <- ls(map)
  if (length(envElements)==0){
    return(TRUE)}
  chromElements <- grep("^.+\\.start$",envElements, value=TRUE)
  if (length(chromElements)==0){
    fails <- c(fails, "Enviroment does not contain any chromosome/strand probe mappings, expecting elements named for example 1.start\n")}
  uniChromNames <- gsub("\\.start$","", chromElements)
  #isStrandSpecific <- length(grep("\\.\\+\\.",chromElements))>0 || length(grep("\\.-\\.",chromElements))>0
  #sp <- strsplit(chromElements,split="\\.")
  #if (isStrandSpecific){
  #  uniChromNames <- sapply(sp, function(x) paste(x[1:2],collapse="."))
  #} else {
  #  uniChromNames <- sapply(sp, function(x) x[1])
  #}
  for (thisName in uniChromNames){
    theseElemNames <- paste(thisName,c("start","end","index","unique"),sep=".")
    theseAreIn <- theseElemNames %in% envElements
    if (!all(theseAreIn)){
      fails <- c(fails, paste("Environment ",deparse(substitute(map))," seems to hold information for chromosome/strand '",thisName,"', but does not contain elements '",paste(theseElemNames[!theseAreIn],collapse="', '"),"'.\n",sep=""))
    }
    theseElements <- mget(theseElemNames, envir=map)
    ## check if the maps related to the same chromosome/strand match up
    if (length(unique(sapply(theseElements,length)))!=1){
      print(sapply(theseElements, length))
      fails <- c(fails, "Lengths of per chromosome/strand vectors differ!")
    }
    ## check if probe-name vectors are factors which can be a problem
    if (any(sapply(theseElements, is.factor))){
      fails <- c(fails, 
         paste("Factor vectors in probeAnno object pose problems.\n",
               "Replace factors in the object or the original probe",
               "positions data.frame by character or integer",
               "vectors, respectively.\n") )
      return(fails)
    }
    ## check if match positions are correctly coded, end positions >= start position
    if (!all(theseElements[[2]] >= theseElements[[1]])){
      fails <- c(fails, paste("Some match positions end before they actually start.\nPlease check elements",paste(thisName,"start",sep="."),"and",paste(thisName,"end",sep="."),".\n"))
    }# if (!all(theseElements[[2]] >= theseElements[[1]]))
    mids <- round((theseElements[[1]]+theseElements[[2]])/2)
    if (is.unsorted(mids)){
      fails <- c(fails, paste("Probe matches are not sorted (in increasing order) by their middle position on chromosome",thisName,"and possibly others.\n"))
    }
  }# for (thisName in uniChromNames)
  if (length(fails) > 0) return(fails)
  return(TRUE)
})# validity of probeAnno


setMethod("[",signature(x="probeAnno"),
   function(x, i, j="missing", ..., drop="missing"){
     stopifnot(is.character(i), length(i)==1)
     if (!exists(i, envir=x@map))
       stop(paste("No mapping '",i,"' in this 'probeAnno' object.\n", sep=""))
     get(i, envir=x@map)
})

setReplaceMethod("[",signature(x="probeAnno"),
   function(x, i, j="missing", value){
     stopifnot(is.character(i), length(i)==1)
     assign(x=i, value=value, envir=x@map)
     x
})

setMethod("get", signature(x="character", pos="missing", envir="probeAnno"),
          function(x, pos="missing", envir, mode="missing", inherits="missing"){
  stopifnot(is.character(x), length(x)==1)
  if (!exists(x, envir=envir@map))
    stop(paste("No mapping '",x,"' in this 'probeAnno' object.\n", sep=""))
  get(x, envir=envir@map)
})

setMethod("ls", signature(name="probeAnno", pos="missing", envir="missing", all.names="missing", pattern="missing"), function(name, pos="missing", envir="missing", all.names="missing", pattern="missing") ls(name@map))

setMethod("chromosomeNames", signature(x="probeAnno"), function(x){
  envElements <- ls(x@map)
  chromElements <- grep("^.+\\.start$",envElements, value=TRUE)
  if (length(chromElements)==0){
    warning("Enviroment does not contain any chromosome/strand probe mappings,\n  expecting elements named, e.g., '1.start'\n"); return(vector("character",0))}
  uniChromNames <- gsub("\\.start$","", chromElements)
  return(uniChromNames)
})

setMethod("arrayName", signature(x="probeAnno"), function(x){
  x@arrayName})

setReplaceMethod("arrayName", signature(x="probeAnno", value="character"),
  function(x, value){x@arrayName <- value; return(x)})

setMethod("genome", signature(x="probeAnno"), function(x){
  x@genome})

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
