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
                 #"vsn"=normalizeBetweenArrays(myRG, method="vsn",...),
                 "vsn"=normalizeBetweenArraysVSN(myRG,...),
                 "Gquantile"=normalizeBetweenArrays(myRG, method="Gquantile",...),
                 "nimblegen"=nimblegenScale(myRG,...),
                 "none"=normalizeWithinArrays(myRG, method="none",...)
                 )
  if (returnMAList){return(myMA)}
  else {return(asExprSet(myMA))}
}# preprocess

asExprSet <- function(myMA){
  stopifnot(inherits(myMA,"MAList"), !is.null(myMA$targets))
  myMA$M <- as.matrix(myMA$M) # if only one sample, M is a vector
  stopifnot(nrow(myMA$targets)==ncol(as.matrix(myMA$M)))
  if (is.null(rownames(myMA$targets)))
    rownames(myMA$targets) <- as.character(myMA$targets[[1]])
  ## include some matching between colnames(myMA$M) and
  ##  rownames of 'targets' here
  colnames(myMA$M) <- rownames(myMA$targets)
  myPD <- new("AnnotatedDataFrame", data=myMA$targets,
              varMetadata=data.frame("varLabel"=colnames(myMA$targets),
                row.names=colnames(myMA$targets)))
  myEset <- new("ExpressionSet",  exprs=myMA$M, phenoData=myPD)
  if (!is.null(myMA$genes$PROBE_ID))
    featureNames(myEset) <- myMA$genes$PROBE_ID
  return(myEset)
}#asExprSet

nimblegenScale <- function(myRG, ...){
  # function to compute Nimblegen's scaled log ratios
  stopifnot(inherits(myRG, "RGList"), all(c("genes","R","G","targets") %in% names(myRG)))

  ## copied from 'affy', no need to include the whole package because of this
  tukey.biweight <- function (x, c = 5, epsilon = 1e-04)
    {
      m <- median(x)
      s <- median(abs(x - m))
      u <- (x - m)/(c * s + epsilon)
      w <- rep(0, length(x))
      i <- abs(u) <= 1
      w[i] <- ((1 - u^2)^2)[i]
      t.bi <- sum(w * x)/sum(w)
      return(t.bi)
    }

  srat <- log2(myRG$R) - log2(myRG$G)
  srat.tbw <- apply(srat, 2, tukey.biweight, ...)
  stopifnot(length(srat.tbw)==ncol(myRG$R))
  srat <- srat - matrix(srat.tbw, nrow=nrow(srat), ncol=ncol(srat), byrow=TRUE)
  resList <- list(M=srat, A=(log2(myRG$R)+log2(myRG$G))/2, genes=myRG$genes, targets=myRG$targets)
  resMA <- new("MAList", resList)
  return(resMA)
}#nimblegenScale


### modified version of normalizeBetweenArrays to handle new clases in VSN
normalizeBetweenArraysVSN <- function(object, targets=NULL, ...) {
  require("vsn")
  stopifnot(is(object,"RGList")|is(object,"MAList"))
  y <- NULL
  if(!is.null(object$G) && !is.null(object$R)) {
    y <- cbind(object$G,object$R)
    object$G <- object$R <- NULL
  }
  if(!is.null(object$M) && !is.null(object$A)) y <- 2^cbind(object$A-object$M/2,object$A+object$M/2)
  if(is.null(y)) stop("object doesn't appear to be RGList or MAList")
  y <- vsnMatrix(y,...)
  n2 <- ncol(exprs(y))/2
  G <- exprs(y)[,1:n2]/log(2)
  R <- exprs(y)[,n2+(1:n2)]/log(2)
  object$M <- R-G
  object$A <- (R+G)/2
  if(!is(object,"MAList")) object <- new("MAList",unclass(object))
  return(object)
}#normalizeBetweenArraysVSN
