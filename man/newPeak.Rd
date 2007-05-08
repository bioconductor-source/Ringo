\name{newPeak}
\alias{newPeak}
\title{Create a list object of class peak}
\description{
  This is just a simple convenience function to create a \code{list}
  of class \code{peak}.
}
\usage{
newPeak(name, chr, start, end, cellType = NULL, antibody, maxPeak, score = NULL, probes = c(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{character; (if possible) unique identifier for the peak}
  \item{chr}{character; chromosome the peak is located on}
  \item{start}{integer; genomic start coordinate of the peak}
  \item{end}{integer; genomic end coordinate of the peak}
  \item{cellType}{optional character describing the cell type in which
    the peak was identified.}
  \item{antibody}{required character vector describing the
    antibody or other characteristic for which fragments were supposedly
    enriched in immuno-precipitation step}
  \item{maxPeak}{numeric; maximal probe level in the peak}
  \item{score}{optional numeric score of the peak}
  \item{probes}{optional character vector of probe identifiers of probes
  that make up the peak}
  \item{\dots}{further arguments that will be additional elements of the
  \code{peak} list object}
}
\value{
  a list object of class \code{peak}, see \code{\link{peak-class}}
}
\author{Tammo Krueger, Joern Toedling}
\note{
  this function is called by other Ringo functions, such as
  \code{\link{findPeaksOnSmoothed}} and normally need not be called by
  the user.
}
\seealso{\code{\link{peak-class}},  \code{\link{findPeaksOnSmoothed}}}
\examples{
  my.peak <- newPeak("H3ac.HAND2.upstream.peak",
                     chr="8", start=60210400, end=60213000,
                     antibody="H3ac",cellType="cellLine.HL1",
                     maxPeak=3.5, score=100)
  print(my.peak)
}
\keyword{internal}
