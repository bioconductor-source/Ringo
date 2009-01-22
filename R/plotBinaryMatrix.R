plotBM <- function(x, boxCol="darkblue", reorder=FALSE, frame=TRUE, ...){
    stopifnot(is.matrix(x))
    if (reorder) {
        ## treat them 
        whichGroup <- x %*% 2^((ncol(x)-1):0)
        numTimes <- table(whichGroup)
        ## to avoid breaks in case >=2 categories occur equally often:
        #numTimes <- numTimes + cumsum(rep(0.1, length(numTimes)))
        ord <- order(numTimes[as.character(whichGroup)], whichGroup, decreasing=FALSE)
    } else {
        ord <- nrow(x):1
    }
    x <- x[ord,]
    blockBorders <- apply(x,2,function(x) diff(c(FALSE,x,FALSE)))
    plot(y=0, xlim=c(0,ncol(x)),  x=0, ylim=c(0,nrow(x)), type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA, frame.plot=FALSE, ...)
    for (j in 1:ncol(blockBorders)){
        theseBlocks <- rbind(c(0,0), cbind(which(blockBorders[,j]==1)-1,which(blockBorders[,j]==-1)-1))
        for (i in 2:nrow(theseBlocks))
            polygon(x=j+c(-1,-1,0, 0), y=theseBlocks[i,c(1,2,2,1)], col=boxCol, border=NA)
    }
    ## draw frame around plotted matrix:
    if (frame){
        arrows(y0=0, y1=0, x0=0, x1=ncol(x), length=0)
        arrows(y0=nrow(x), y1=nrow(x), x0=0, x1=ncol(x), length=0)
        arrows(y0=0, y1=nrow(x), x0=ncol(x), x1=ncol(x), length=0)
        arrows(y0=0, y1=nrow(x), x0=0, x1=0, length=0)
    }# if (frame)
    ## annotate axes
    if (!is.null(rownames(x)))
        axis(side=2, at=seq(nrow(x))-0.5, labels=rownames(x), tick=FALSE, line=0, las=1)
    if (!is.null(colnames(x)))
        axis(side=1, at=seq(ncol(x))-0.5, labels=colnames(x), tick=FALSE, line=0)
    invisible(x)
}#plotBM
