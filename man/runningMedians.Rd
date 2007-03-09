\name{runningMedians}
\alias{runningMedians}
\title{Compute median score over a running window}
\description{
  Compute median score over a running window
}
\usage{
runningMedians(positions, scores, cutoff)
}
\arguments{
  \item{positions}{numeric; sorted vector of (genomic) positions of scores}
  \item{scores}{numeric; scores to be smoothed associated to the
    \code{positions}}
  \item{cutoff}{numeric, half the window size of the sliding window}
}
\details{
  For each position of argument \code{positions}, the median score over
  a window of size 2*\code{cutoff} is computed. That means that for the
  smoothed score at position i, scores at adjacent positions j are
  included in the median computation as long as
  they are not further than \code{cutoff} away from \code{positions[i]}.
}
\value{
  numeric vector of medians over running window centered at the
  positions that were specified in argument \code{positions}.
  These positions are given as \code{names} of the resulting vector.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link[stats]{median}}}
\examples{
 sampleSize <- 1000
 ap <- cumsum(1+round(runif(sampleSize)*10))
 as <- rnorm(sampleSize)
 arm <- runningMedians(ap, as, 20)
 plot(ap, as)
 points(ap, arm, type="l", col="red", lwd=2)
}
\keyword{manip}% at least one, from doc/KEYWORDS
