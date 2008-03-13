
sliding.meansd <- function(positions, scores, half.width) {
  stopifnot(!is.unsorted(positions), length(positions) == length(scores), half.width >= 0)
  res <- .Call(ringoMovingMeanSd, as.integer(positions), as.numeric(scores), as.integer(half.width))
  colnames(res) <- c("mean","sd","count")
  rownames(res) <- positions
  return(res)
}#sliding.meansd
