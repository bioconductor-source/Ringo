\name{plotBM}
\Rdversion{1.1}
\alias{plotBM}
\alias{plotBinaryMatrix}
\title{
  Visualization of a binary matrix
}
\description{
  This function produces simple, heatmap-like visualizations of binary
  matrices.
}
\usage{
plotBM(x, boxCol = "darkblue", reorder = FALSE, frame = TRUE, ...)
}
\arguments{
  \item{x}{Binary matrix to visualize}
  \item{boxCol}{Color to use for boxes of '1's}
  \item{reorder}{logical; states whether the rows
    shall be reordered according to the size of the category}
  \item{frame}{logical; states whether a frame should be drawn around
    the visualization. In contrast to the frame drawn in
    \code{plot.default}, there is no gap between the visualization and
    this frame.}
  \item{\dots}{
    further arguments passed on to \code{plot.default}}
}
\details{
  For reordering, each row is interpreted as a binary matrix, for
  example a row z=(1,0,0,1) would be interpreted as the binary number
  1001 = 9 in the decimal system. Rows are then reordered by the
  frequency of each binary number with the rows that correspond to the
  most frequent binary number shown at the top in the visualization.
}
\note{
  An alternative way to display such matrices are given by
  \code{heatmap} or, the simpler version thereof, \code{image}.
  However, image files produced with this functions tend to be very
  large. This function uses \code{plot.default} and \code{polygon} which
  results in much smaller file sizes and is sufficient for binary
  matrices. 
}
\value{
  The function invisibly returns the (reordered)
  matrix \code{x}, but its mainly called for its side effect of
  producing the visualization.
}
\author{
  Joern Toedling
}
\seealso{\code{\link{polygon}},\code{\link{colors}}}
\examples{
  A <- matrix(round(runif(80)), ncol=4, byrow=TRUE)
  dimnames(A)=list(letters[seq(nrow(A))],
                   as.character(as.roman(seq(ncol(A)))))
  show(A)
  plotBM(A, reorder=FALSE)
  plotBM(A, reorder=TRUE)
}
\keyword{hplot}
