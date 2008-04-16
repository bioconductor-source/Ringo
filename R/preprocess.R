preprocess <- function(myRG, method="vsn", returnMAList=FALSE,
                       idColumn="PROBE_ID", verbose=TRUE, ...){
  stopifnot(inherits(myRG, "RGList"))
  if (!returnMAList)
    stopifnot(is.data.frame(myRG$genes), idColumn %in% names(myRG$genes))
  method <- match.arg(method, choices=c("vsn","loess","median","Gquantile", "Rquantile", "nimblegen", "Gvsn", "Rvsn", "none"))
  if (method %in% c("loess","median")){
    if (verbose) cat("Background correction...\n")
    myRG <- backgroundCorrect(myRG, method="normexp", offset=50)
  }
  if (verbose) cat("Normalizing...\n")
  myMA <- switch(method,
                 "loess"=normalizeWithinArrays(myRG, method="loess",...),
                 "median"=normalizeWithinArrays(myRG, method="median",...),
                 "vsn"=normalizeBetweenArraysVSN(myRG,...),
                 "Gquantile"=normalizeBetweenArrays(myRG, method="Gquantile",...),
                 "Rquantile"=normalizeBetweenArrays(myRG, method="Rquantile",...),
                 "nimblegen"=nimblegenScale(myRG,...),
                 "Gvsn"=oneChannelVSN(myRG, "G", ...),
                 "Rvsn"=oneChannelVSN(myRG, "R", ...),
                 "none"=normalizeWithinArrays(myRG, method="none",...)
                 )
  if (returnMAList){return(myMA)}
  else {return(asExprSet(myMA, idColumn=idColumn))}
}# preprocess

asExprSet <- function(from, idColumn="PROBE_ID"){
  stopifnot(inherits(from,"MAList"), !is.null(from$targets),
            is.data.frame(from$genes), idColumn %in% names(from$genes))
  from$M <- as.matrix(from$M) # if only one sample, M is a vector
  stopifnot(nrow(from$targets)==ncol(as.matrix(from$M)))
  if (is.null(rownames(from$targets)))
    rownames(from$targets) <- as.character(from$targets[[1]])
  colnames(from$M) <- rownames(from$targets)
  myPD <- new("AnnotatedDataFrame", data=from$targets,
              varMetadata=data.frame("varLabel"=colnames(from$targets),
                row.names=colnames(from$targets)))
  myEset <- new("ExpressionSet",  exprs=from$M, phenoData=myPD)
  featureNames(myEset) <- make.names(from$genes[[idColumn]], unique=TRUE)
  featureData(myEset) <- as(from$genes,"AnnotatedDataFrame")
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
  G <- exprs(y)[,1:n2]
  R <- exprs(y)[,n2+(1:n2)]
  object$M <- R-G
  object$A <- (R+G)/2
  if(!is(object,"MAList")) object <- new("MAList",unclass(object))
  return(object)
}#normalizeBetweenArraysVSN

### do VSN normalization based on one channel only:
oneChannelVSN <- function(object, channel="G", ...) {
  stopifnot(is(object,"RGList"))
  channel <- match.arg(channel, c("G","R"))
  y <- object[[channel]]
  yfit <- vsnMatrix(y,...)
  n2 <- ncol(exprs(y))/2
  G <- exprs(predict(yfit, newdata=object[["G"]]))
  R <- exprs(predict(yfit, newdata=object[["R"]]))
  object$M <- R-G
  object$A <- (R+G)/2
  if(!is(object,"MAList")) object <- new("MAList",unclass(object))
  return(object)
}#oneChannelVSN
