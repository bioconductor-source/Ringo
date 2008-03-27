\name{preprocess}
\alias{preprocess}
\title{Preprocess Nimblegen Raw Intensities}
\description{
  Calls one of various \code{limma} functions to transform raw probe
  intensities into (background-corrected) normalized log ratios
  (M-values). 
}
\usage{
preprocess(myRG, method = "vsn", returnMAList=FALSE, verbose=TRUE, ...)
}
\arguments{
  \item{myRG}{object of class \code{RGList}}
  \item{method}{string; denoting which normalization method to choose,
    see below for details}
  \item{returnMAList}{logical; should an MAList object be returned?
    Default is to return an ExpressionSet object.}
  \item{verbose}{logical; progress output to STDOUT?}
  \item{\dots}{further arguments to be passed on
    \code{normalizeWithinArrays} and \code{normalizeBetweenArrays}}
}
\details{
 The procedure and called \code{limma} functions depend on the choice of
 method.
 \describe{
   \item{loess}{First performs background correction of intensity values
     using \code{limma}'s function \code{backgroundCorrect} with
     \code{method="normexp"} and \code{offset=50}. Then calls
     \code{normalizeWithinArrays}.}
   \item{vsn}{Does not do any background correction, but rather directly
     calls \code{normalizeBetweenArrays} with \code{method="vsn"}.}
  \item{Gquantile}{Does not do any background correction, but rather directly
     calls \code{normalizeBetweenArrays} with \code{method="Gquantile"}.}
  \item{Rquantile}{Does not do any background correction, but rather directly
     calls \code{normalizeBetweenArrays} with \code{method="Rquantile"}.}
   \item{median}{First performs background correction of intensity values
     using \code{limma}'s function \code{backgroundCorrect} with
     \code{method="normexp"} and \code{offset=50}. Then calls
     \code{normalizeWithinArrays} with \code{method="median".}}
    \item{nimblegen}{Scaling procedure used by Nimblegen. Yields
      scaled log-ratios by a two step procedure:
      srat = log2(R) - log2(G)
      srat = srat - tukey.biweight(srat)}
    \item{none}{No normalization of probe intensities,
      takes raw \code{log2(R)-log2(G)} as component \code{M}
      and \code{(log2(R)+log2(G))/2} as component \code{A};
      uses \code{normalizeWithinArrays} with \code{method="none"}.}
  }
}
\value{
  Returns normalized, transformed values as an object of class
  \code{ExpressionList} or \code{MAList}.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link[limma]{backgroundCorrect}},
  \code{\link[limma]{normalizeWithinArrays}},
  \code{\link[limma]{normalizeBetweenArrays}},
  \code{\link[limma]{malist}},\code{\link[Biobase]{ExpressionSet}},
  \code{\link[affy]{tukey.biweight}}}
\examples{
   exDir <- system.file("exData",package="Ringo")
   exRG <- readNimblegen("example_targets.txt","spottypes.txt",path=exDir)
   exampleX <- preprocess(exRG)
   sampleNames(exampleX) <- 
     make.names(paste(exRG$targets$Cy5,"vs",exRG$targets$Cy3,sep="_"))
   print(exampleX)

   ### compare VSN to NimbleGen's tukey-biweight scaling
   exampleX.NG <- preprocess(exRG, method="nimblegen")
   sampleNames(exampleX.NG) <- sampleNames(exampleX)
   if (interactive())
     corPlot(cbind(exprs(exampleX),exprs(exampleX.NG)),
       grouping=c("VSN normalized","Tukey-biweight scaled"))
}
\keyword{manip}% at least one, from doc/KEYWORDS
