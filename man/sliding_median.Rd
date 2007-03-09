\name{sliding.median}
\alias{sliding.median}
\alias{slidingmedian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute median score over a running window}
\description{
  Compute median score over a sliding window.
}
\usage{
sliding.median(positions, scores, half.width, return.counts = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{positions}{numeric; sorted vector of (genomic) positions of scores}
  \item{scores}{numeric; scores to be smoothed associated to the
    \code{positions}}
  \item{half.width}{numeric, half the window size of the sliding window}
  \item{return.counts}{logical; should the number of points, e.g. probes,
  that were used for computing the median in each sliding window also
  returned?}
}
\details{
  For each position of argument \code{positions}, the median score over
  a window of size 2*\code{half.width} is computed as a smoothed score.
  For the smoothed score at position i, scores at adjacent positions j are
  included in the median computation as long as
  they are not further than \code{half.width} away from \code{positions[i]}.
}
\value{
  Matrix with two columns:
  \item{median}{medians over running window centered at the
    positions that were specified in argument \code{positions}.}
  \item{count}{number of points that were considerd for computing
    the median at each position}
  These positions are given as \code{row.names} of the resulting vector.

  If argument \code{return.counts} is \code{FALSE}, only a vector
  of the medians is returned, with the positions as \code{names}.
}
\author{Oleg Sklyer and Joern Toedling}
\seealso{\code{\link[stats]{median}}}
\examples{
  sampleSize <- 1000
  ap <- cumsum(1+round(runif(sampleSize)*10))
  as <- c(rnorm(floor(sampleSize/3)),
          rnorm(ceiling(sampleSize/3),mean=1.5),
          rnorm(floor(sampleSize/3)))
  arm <- sliding.median(ap, as, 20)
  plot(ap, as, pch=20, xlab="position",ylab="level")
  points(ap, arm[,1], type="l", col="red", lwd=2)
}
\keyword{manip}% __ONLY ONE__ keyword per line
