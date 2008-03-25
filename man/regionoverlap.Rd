\name{regionOverlap}
\alias{regionOverlap}
\alias{region.overlap}
\title{Function to compute overlap of genomic regions}
\description{
  Given two data frames of genomic regions, this function computes the
  base-pair overlap, if any, between every pair of regions from the two
  lists.
}
\usage{
regionOverlap(xdf, ydf, chrColumn = "chr", startColumn = "start",
endColumn = "end", mem.limit=1e8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdf}{\code{data.frame} that holds the first set of genomic regions}
  \item{ydf}{\code{data.frame} that holds the first set of genomic regions}
  \item{chrColumn}{character; what is the name of the column that holds
    the chromosome name of the regions in \code{xdf} and \code{ydf}}
  \item{startColumn}{character; what is the name of the column that holds
    the start position of the regions in \code{xdf} and \code{ydf}}
  \item{endColumn}{character; what is the name of the column that holds
    the start position of the regions in \code{xdf} and \code{ydf}}
  \item{mem.limit}{integer value; what is the maximal allowed size of
    matrices during the computation}
}
\value{
  Originally, a matrix with \code{nrow(xdf)} rows and
  \code{nrow(ydf)} columns, in which entry \code{X[i,j]} specifies the
  length of the overlap between region \code{i} of the first list
  (\code{xdf}) and region \code{j} of the second list (\code{ydf}).
  Since this matrix is very sparse, we use the \code{matrix.csr}
  representation from the \code{SparseM} package for it.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\note{
  The function only return the absolute length of overlapping regions in
  base-pairs. It does not return the position of the overlap or the
  fraction of region 1 and/or region 2 that overlaps the other regions.

  The argument \code{mem.limit} is not really a limit to used RAM, but
  rather the maximal size of matrices that should be allowed during the
  computation. If larger matrices would arise, the second regions list
  is split into parts and the overlap with the first list is computed
  for each part. During computation, matrices of size
  \code{nrow(xdf)} times \code{nrow(ydf)} are created.
}
\seealso{\code{\link[SparseM]{matrix.csr-class}}}
\examples{
  ## toy example:
  regionsH3ac <- data.frame(chr=c("chr1","chr7","chr8","chr1","chrX","chr8"), start=c(100,100,100,510,100,60), end=c(200, 200, 200,520,200,80))
  regionsH4ac <- data.frame(chr=c("chr1","chr2","chr7","chr8","chr9"),
start=c(500,100,50,80,100), end=c(700, 200, 250, 120,200))

  ## compare the regions first by eye
  ##  which ones do overlap and by what amount?
  regionsH3ac
  regionsH4ac

  ## compare it to the result:
  as.matrix(regionOverlap(regionsH3ac, regionsH4ac))
  whichCsr(regionOverlap(regionsH3ac, regionsH4ac))
}
\keyword{manip}
