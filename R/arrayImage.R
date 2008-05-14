"image.RGList" <-
function(x, arrayno, channel=c("red","green","logratio"),
         mycols=NULL, mybreaks=NULL, dim1="X", dim2="Y", ...){
  stopifnot(inherits(x, "RGList"), "genes" %in% names(x),
            is.character(dim1), length(dim1)==1,
            is.character(dim2), length(dim2)==1,
            all(c(dim1, dim2) %in% names(x$genes)))
  myRG <- x
  maxX <- max(myRG$genes[[dim1]]); maxY <- max(myRG$genes[[dim2]])
  channel <- match.arg(channel, c("red","green","logratio"))
  arrayMat <- matrix(0, nrow=maxX, ncol=maxY)
  myDat <- switch(channel,
                  "red"=myRG$R[,arrayno],
                  "green"=myRG$G[,arrayno],
                  "logratio"=log2(myRG$R[,arrayno])-log2(myRG$G[,arrayno]))
  stopifnot(length(myDat)==nrow(myRG$genes))
  for (i in 1:nrow(myRG$genes))
    arrayMat[myRG$genes[[dim1]][i], myRG$genes[[dim2]][i]] <- myDat[i]
  imageTitle <- switch(channel,
                  "red"=colnames(myRG$R)[arrayno],
                  "green"=colnames(myRG$G)[arrayno],
                  "logratio"=paste("log-ratio",colnames(myRG$R)[arrayno]))
  if (is.null(mycols)){
    paletteName <- switch(channel, "red"="Reds", "green"="Greens",
                          "logratio"="Blues")
    if (is.null(mybreaks))
      mycols <- colorRampPalette(c("#FFFFFF", brewer.pal(9,paletteName)),space="rgb")(25)
    else
      mycols <- colorRampPalette(c("#FFFFFF", brewer.pal(9,paletteName)),space="rgb")(length(mybreaks)-1)
  }#if (is.null(mycols))
  
  if (is.null(mybreaks)) mybreaks=c(-1,quantile(myDat, probs=seq(0,1,length.out=length(mycols)), na.rm=TRUE))
  
  stopifnot(length(mybreaks)- 1 == length(mycols))
  image(arrayMat, breaks=mybreaks, col=mycols, main=imageTitle, xaxt="n", yaxt="n",...)
  invisible(arrayMat)
}#image.RGList

"arrayImage" <-
  function(x,...) { image.RGList(x,...) }
