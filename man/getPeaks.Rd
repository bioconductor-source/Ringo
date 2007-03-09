\name{getPeaks}
\alias{getPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Peaks in ChIP data}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPeaks(eSet, chrom, probeAnno, samples = 1, est.sd = 1000, condenseDist = 2000, byStrand = FALSE, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eSet}{ ~~Describe \code{eSet} here~~ }
  \item{chrom}{ ~~Describe \code{chrom} here~~ }
  \item{probeAnno}{ ~~Describe \code{probeAnno} here~~ }
  \item{samples}{ ~~Describe \code{samples} here~~ }
  \item{est.sd}{ ~~Describe \code{est.sd} here~~ }
  \item{condenseDist}{ ~~Describe \code{condenseDist} here~~ }
  \item{byStrand}{ ~~Describe \code{byStrand} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{  ~~ If necessary, more details than the description above ~~}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
%\note{ ~~further notes~~ }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
    # r19 <- getPeaks(hl1V, chrom="19", probeAnno=probeAnno)
}
\keyword{manip}% at least one, from doc/KEYWORDS
