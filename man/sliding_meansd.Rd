\name{sliding.meansd}
\alias{sliding.meansd}
\title{Compute mean and standard deviation of scores in a sliding window}
\description{
  This functions is used to slide a window of specified size over
  scores at given positions.
  Computed is the mean and standard deviation over the scores
  in the window.
}
\usage{
sliding.meansd(positions, scores, half.width)
}
\arguments{
  \item{positions}{numeric; sorted vector of (genomic) positions of scores}
  \item{scores}{numeric; scores to be smoothed associated to the
    \code{positions}}
  \item{half.width}{numeric, half the window size of the sliding window}
}
\value{
  Matrix with three columns:
  \item{mean}{means over scores in running window centered at the
    positions that were specified in argument \code{positions}.}
  \item{sd}{standard deviations over scores in running window centered
    at the positions that were specified in argument \code{positions}.}
  \item{count}{number of points that were considered for computing
    the mean and standard deviation at each position}
}
\author{Joern Toedling and Oleg Sklyar}
\seealso{\code{\link{sliding.quantile}}}
\examples{
   set.seed(123)
   sampleSize <- 10
   ap <- cumsum(1+round(runif(sampleSize)*10))
   as <- c(rnorm(floor(sampleSize/3)),
           rnorm(ceiling(sampleSize/3),mean=1.5),
           rnorm(floor(sampleSize/3)))
   sliding.meansd(ap, as, 20)
   ap
   mean(as[1:3])
   sd(as[1:3])
}
\keyword{manip}
