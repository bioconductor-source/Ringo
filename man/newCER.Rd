\name{newCher}
\alias{newCher}
\title{Create a list object of class cher}
\description{
  This is just a simple convenience function to create a \code{list}
  of class \code{cher} (ChIP-enriched region).
}
\usage{
newCher(name, chr, start, end, cellType = NULL, antibody, maxCher, score = NULL, probes = c(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{character; (if possible) unique identifier for the cher}
  \item{chr}{character; chromosome the cher is located on}
  \item{start}{integer; genomic start coordinate of the Cher}
  \item{end}{integer; genomic end coordinate of the Cher}
  \item{cellType}{optional character describing the cell type in which
    the cher was identified.}
  \item{antibody}{required character vector describing the
    antibody or other characteristic for which fragments were supposedly
    enriched in immuno-precipitation step}
  \item{maxCher}{numeric; maximal probe level in the cher}
  \item{score}{optional numeric score of the cher}
  \item{probes}{optional character vector of probe identifiers of probes
  that make up the cher}
  \item{\dots}{further arguments that will be additional elements of the
  \code{cher} list object}
}
\value{
  a list object of class \code{cher}, see \code{\link{cher-class}}
}
\author{Tammo Krueger, Joern Toedling}
\note{
  this function is called by other Ringo functions, such as
  \code{\link{findChersOnSmoothed}} and normally need not be called by
  the user.
}
\seealso{\code{\link{cher-class}},  \code{\link{findChersOnSmoothed}}}
\examples{
  cher1 <- newCher("Suz12.Nudt2.upstream.cher",
                     chr="9", start=34318900, end=34320100,
                     antibody="Suz12",cellType="human",
                     maxCher=2, score=99)
  print(cher1)
}
\keyword{internal}
