# Function to compute overlaps between (enriched) genomic regions in one list and those in another list
peakOverlap <- function(xdf, ydf, chrColumn="chr",startColumn="start",endColumn="end") {
  stopifnot(is.data.frame(xdf), is.data.frame(ydf),
            all(c(chrColumn, startColumn, endColumn) %in% names(xdf)),
            all(c(chrColumn, startColumn, endColumn) %in% names(ydf)))
  if (!all(xdf[[startColumn]]<=xdf[[endColumn]]))
    stop("Some regions in",deparse(substitute(xdf)),"have end postions that are smaller than their start positions.\n")
  if (!all(ydf[[startColumn]]<=ydf[[endColumn]]))
    stop(paste("Some regions in",deparse(substitute(ydf)),"have end postions that are smaller than their start positions.\n"))
  res <- .Call(ringoPeakOverlap, as.character(xdf[[chrColumn]]), as.integer(as.character(xdf[[startColumn]])), as.integer(as.character(xdf[[endColumn]])), as.character(ydf[[chrColumn]]), as.integer(as.character(ydf[[startColumn]])), as.integer(as.character(ydf[[endColumn]])))
  return(res)
}#peakOverlap
