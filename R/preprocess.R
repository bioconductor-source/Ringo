preprocess <- function(myRG, method="vsn", returnMAList=FALSE,
                       verbose=TRUE, ...){
  stopifnot(inherits(myRG, "RGList"))
  method <- match.arg(method, choices=c("vsn","loess","median","Gquantile", "nimblegen","none"))
  if (method %in% c("loess","median")){
    if (verbose) cat("Background correction...\n")
    myRG <- backgroundCorrect(myRG, method="normexp", offset=50)
  }
  if (verbose) cat("Normalizing...\n")
  myMA <- switch(method,
                 "loess"=normalizeWithinArrays(myRG, method="loess",...),
                 "median"=normalizeWithinArrays(myRG, method="median",...),
                 "vsn"=normalizeBetweenArrays(myRG, method="vsn",...),
                 "Gquantile"=normalizeBetweenArrays(myRG, method="Gquantile",...),
                 "nimblegen"=nimblegenNorm(myRG,...),
                 "none"=normalizeWithinArrays(myRG, method="none",...)
                 )
  if (returnMAList){return(myMA)}
  else {return(asExprSet(myMA))}
}# preprocess

asExprSet <- function(myMA){
  stopifnot(inherits(myMA,"MAList"), !is.null(myMA$targets), nrow(myMA$targets)==ncol(myMA$M))
  if (is.null(rownames(myMA$targets)))
    rownames(myMA$targets) <- as.character(myMA$targets[[1]])
  ## include some matching between colnames(myMA$M) and
  ##  rownames of 'targets' here
  colnames(myMA$M) <- rownames(myMA$targets)
  ## old version using "exprSet":
  #myPD <- new("phenoData", pData=myMA$targets, varLabels=as.list(colnames(myMA$targets)))
  #myEset <- new("exprSet", exprs=myMA$M, phenoData=myPD)
  #myEset <- as(myEset, "ExpressionSet")  
  myPD <- new("AnnotatedDataFrame", data=myMA$targets,
              varMetadata=data.frame("varLabel"=colnames(myMA$targets),
                row.names=colnames(myMA$targets)))
  myEset <- new("ExpressionSet",  exprs=myMA$M, phenoData=myPD)
  if (!is.null(myMA$genes$PROBE_ID))
    featureNames(myEset) <- myMA$genes$PROBE_ID
  return(myEset)
}#asExprSet

nimblegenNorm <- function(myRG, ...){
  # function to compute Nimblegen's scaled log ratios
  require(affy)
  stopifnot(inherits(myRG, "RGList"), all(c("genes","R","G","targets") %in% names(myRG)))
  srat <- log2(myRG$R) - log2(myRG$G)
  srat.tbw <- apply(srat, 2, tukey.biweight, ...)
  stopifnot(length(srat.tbw)==ncol(myRG$R))
  srat <- srat - matrix(srat.tbw, nrow=nrow(srat), ncol=ncol(srat), byrow=TRUE)
  resList <- list(M=srat, A=(log2(myRG$R)+log2(myRG$G))/2, genes=myRG$genes, targets=myRG$targets)
  resMA <- new("MAList", resList)
  return(resMA)
}#nimblegenNorm
