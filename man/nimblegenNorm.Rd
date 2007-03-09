\name{nimblegenNorm}
\alias{nimblegenNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute scaled log-ratios}
\description{
  This function compute the scaled log-ratios from raw probe
  intensities, as done by Nimblegen for normalization.
}
\usage{
nimblegenNorm(myRG, ...)
}
%- maybe also 'usage' for other objects documented here.
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
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
%\note{ ~~further notes~~ }
\seealso{\code{\link[affy]{tukey.biweight}}}
%\examples{}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{internal}% __ONLY ONE__ keyword per line
