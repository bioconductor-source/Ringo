\name{sliding.quantile}
\alias{sliding.quantile}
\alias{slidingquantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute quantile of scores in a sliding window}
\description{
  This functions is used to slide a window of specified size over
  scores at given positions.
  Computed is the quantile over the scores in the window.
}
\usage{
sliding.quantile(positions, scores, half.width, prob = 0.5,
                 return.counts = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{positions}{numeric; sorted vector of (genomic) positions of scores}
  \item{scores}{numeric; scores to be smoothed associated to the
    \code{positions}}
  \item{half.width}{numeric, half the window size of the sliding window}
  \item{prob}{numeric specifying which quantile is to be computed over
    the scores in the window; default 0.5 means compute the median over
    the scores.}
  \item{return.counts}{logical; should the number of points, e.g. probes,
  that were used for computing the median in each sliding window also
  returned?}
}
\details{
  This function is meant as a more general alternative to the
  \code{sliding.median} function. In most, but not in all, situations,
  we found it to outperform \code{sliding.median}.
}
\value{
  Matrix with two columns:
  \item{quantile}{medians over running window centered at the
    positions that were specified in argument \code{positions}.}
  \item{count}{number of points that were considerd for computing
    the median at each position}
  These positions are given as \code{row.names} of the resulting vector.
  If argument \code{return.counts} is \code{FALSE}, only a vector
  of the medians is returned, with the positions as \code{names}.
}
\author{Oleg Sklyar and Joern Toedling}
\seealso{\code{\link[stats]{quantile}},\code{\link{sliding.median}}}
\examples{
  sampleSize <- 1000
  ap <- cumsum(1+round(runif(sampleSize)*10))
  as <- c(rnorm(floor(sampleSize/3)),
          rnorm(ceiling(sampleSize/3),mean=1.5),
          rnorm(floor(sampleSize/3)))
  arm <- sliding.quantile(ap, as, 20)
  arq <- sliding.quantile(ap, as, 20, prob=0.25)
  plot(ap, as, pch=20, xlab="position",ylab="level")
  points(ap, arm[,1], type="l", col="red", lwd=2)
  points(ap, arq[,1], type="l", col="green", lwd=2)
  legend(x="topleft", legend=c("median","1st quartile"),
         col=c("red","green"), lty=1, lwd=2)
}
\keyword{manip}% __ONLY ONE__ keyword per line
