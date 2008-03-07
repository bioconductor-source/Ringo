# Function to compute overlaps between (enriched) genomic regions in one list and those in another list
peakOverlap <- function(xdf, ydf, chrColumn="chr",startColumn="start",endColumn="end", mem.limit=1e8) {
  stopifnot(is.data.frame(xdf), is.data.frame(ydf),
            all(c(chrColumn, startColumn, endColumn) %in% names(xdf)),
            all(c(chrColumn, startColumn, endColumn) %in% names(ydf)))
  if (!all(xdf[[startColumn]]<=xdf[[endColumn]]))
    stop("Some regions in",deparse(substitute(xdf)),"have end postions that are smaller than their start positions.\n")
  if (!all(ydf[[startColumn]]<=ydf[[endColumn]]))
    stop(paste("Some regions in",deparse(substitute(ydf)),"have end postions that are smaller than their start positions.\n"))
  ### avoid creation of too big matrices here:
  if (nrow(xdf)*nrow(ydf)<mem.limit){
    res <- .Call(ringoPeakOverlap, as.character(xdf[[chrColumn]]), as.integer(as.character(xdf[[startColumn]])), as.integer(as.character(xdf[[endColumn]])), as.character(ydf[[chrColumn]]), as.integer(as.character(ydf[[startColumn]])), as.integer(as.character(ydf[[endColumn]])))
    ## convert the result to matrix.csr
    res <- as.matrix.csr(res)
  } else {
    ## iterative procedure:
    if (nrow(xdf)> mem.limit)
      stop("Matrix ", deparse(substitute(xdf)), " has too many rows. Work with subsets of this matrix or increase parameter 'mem.limit'.\n") ## maybe implement a workaround that as well later
    maxYRows <- floor(mem.limit/nrow(xdf))
    yRowSeq  <- unique(as.integer(c(seq(maxYRows, nrow(ydf), by=maxYRows), nrow(ydf))))
    partRess <- vector("list", length(yRowSeq))
    firstYCol <- as.integer(1)
    for (j in seq(length(yRowSeq))){
      lastYCol <-  yRowSeq[j]
      partYdf <- ydf[firstYCol:lastYCol,,drop=FALSE]
      thisRes <- .Call(ringoPeakOverlap, as.character(xdf[[chrColumn]]), as.integer(as.character(xdf[[startColumn]])), as.integer(as.character(xdf[[endColumn]])), as.character(partYdf[[chrColumn]]), as.integer(as.character(partYdf[[startColumn]])), as.integer(as.character(partYdf[[endColumn]])))
      thisRes <- as.matrix.csr(thisRes)
      partRess[[j]] <- thisRes
      firstYCol <-  lastYCol+1
    }# for (j in seq(length(yRowSeq)))
    res <- do.call("cbind", partRess)
  }
  return(res)
}#peakOverlap
