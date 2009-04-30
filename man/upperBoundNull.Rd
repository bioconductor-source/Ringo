\name{upperBoundNull}
\alias{upperBoundNull}
\alias{nullUpperBound}
\title{function to estimate upper limit of null distribution}
\description{
  The data is assumed to arise from a mixture of two distributions,
  a symmetric null distribution with its mode close to zero,
  and an alternative distribution that is stochastically larger than the
  null. This function tries to pinpoint the minimum of data values that
  are more likely to arise from the alternative distribution, i.e. an
  upper bound for values following the null distribution.
}
\usage{
upperBoundNull(x, modeMethod = "shorth", limits = c(-1, 1), prob = 0.99, ...)
}
\arguments{
  \item{x}{numeric vector of data values}
  \item{modeMethod}{character string; which method to use for estimating
  the mode of the null distribution; see details}
  \item{limits}{numeric of length 2; data values within this range are
    used for estimating the mode of the null distribution}
  \item{prob}{quantile of the null distribution that will be used as an
    upper bound}
  \item{\dots}{additional arguments that are passed on to the function
    for mode estimation}
}
\details{
  For estimating the mode of the null distribution, current
  options are
  \describe{
    \item{"shorth"}{the function \code{\link[genefilter]{shorth}}}
    \item{"half.range.mode"}{the function
      \code{\link[genefilter]{half.range.mode}}}
    \item{"null"}{does not estimate the mode from the data, but sets it
      to 0}
  }
}
\value{
  a single numeric value which is the estimated upper bound for the null
  distribution.
}
\author{Joern Toedling, based on suggestions
  by Richard Bourgon}
\note{
  This way of estimating the null distribution is mentioned in the PhD
  thesis of Richard Bourgon.
}
\seealso{\code{\link[genefilter]{shorth}}, \code{\link[genefilter]{half.range.mode}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
       modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  apply(exprs(smoothX), 2, upperBoundNull)
}
\keyword{manip}
