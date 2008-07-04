\name{nimblegenScale}
\alias{nimblegenScale}
\title{Function to compute scaled log-ratios}
\description{
  This function compute the scaled log-ratios from raw probe
  intensities, as done by Nimblegen for scaling of ChIP-chip data.
}
\usage{
nimblegenScale(myRG, ...)
}
\arguments{
  \item{myRG}{Object of class \code{RGList}}
  \item{\dots}{further arguments passed on to \code{tukey.biweight}}
}
\details{
  Nimblegen provides scaled  log-ratios as normalized values for the
  probes.
  log.ratio = log2(R)-log2(G)
  scaled.ratio = log.ratio - tukey.biweight(log.ratio)
}
\value{
  Return an \code{MAList}, with the M slot of the list holding the
  scaled log ratios.
}
\references{for more details on the Tukey biweight estimator:
  Statistical Algorithms Description Document, 2002,
  Affymetrix. }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
%\note{ ~~further notes~~ }
%\seealso{}
%\examples{}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{internal}% __ONLY ONE__ keyword per line
