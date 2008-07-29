"image.RGList" <- function (x, arrayno, channel = c("red", "green", "logratio"),
    mycols = NULL, mybreaks = NULL, dim1 = "X", dim2 = "Y", ppch=20, pcex=0.3, verbose=TRUE, ...)
{
    stopifnot(inherits(x, "RGList"), "genes" %in% names(x), is.character(dim1),
        length(dim1) == 1, is.character(dim2), length(dim2) ==
            1, all(c(dim1, dim2) %in% names(x$genes)))
    myRG <- x
    maxX <- max(myRG$genes[[dim1]])
    maxY <- max(myRG$genes[[dim2]])
    if (verbose) cat("Dimensions of array:", maxX,"x", maxY,".\n")
    channel <- match.arg(channel, c("red", "green", "logratio"))
    myDat <- switch(channel, red = myRG$R[, arrayno], green = myRG$G[,
        arrayno], logratio = log2(myRG$R[, arrayno]) - log2(myRG$G[,
        arrayno]))
    stopifnot(length(myDat) == nrow(myRG$genes))
    if (is.null(mycols)) {
        paletteName <- switch(channel, red = "Reds", green = "Greens",
            logratio = "Blues")
        if (is.null(mybreaks))
            mycols <- colorRampPalette(c("#FFFFFF", brewer.pal(9,
                paletteName)), space = "rgb")(25)
        else mycols <- colorRampPalette(c("#FFFFFF", brewer.pal(9,
            paletteName)), space = "rgb")(length(mybreaks) -
            1)
    }
    if (is.null(mybreaks))
        mybreaks = c(-1, quantile(myDat, probs = seq(0, 1, length.out = length(mycols)), na.rm = TRUE))
    stopifnot(length(mybreaks) - 1 == length(mycols))
    myCut <- cut(myDat, mybreaks)
    plot(x=1, y=1, xlim=c(1, maxX), ylim=c(1,maxY), type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA, pch=ppch, cex=pcex, ...)
    points(x=myRG$genes[[dim1]], y=myRG$genes[[dim2]], col=mycols[as.numeric(myCut)], pch=ppch, cex=pcex, ...)
    invisible(NULL)  
}#image.RGList

"arrayImage" <-
  function(x,...) { image.RGList(x,...) }
