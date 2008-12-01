\name{preprocess}
\alias{preprocess}
\title{Preprocess Raw ChIP-chip Intensities}
\description{
  Calls one of various (\code{limma}) functions to transform raw probe
  intensities into (background-corrected) normalized log ratios
  (M-values).
}
\usage{
preprocess(myRG, method="vsn", ChIPChannel="R", inputChannel="G",
           returnMAList=FALSE, idColumn="PROBE_ID", verbose=TRUE, ...)
}
\arguments{
  \item{myRG}{object of class \code{RGList}}
  \item{method}{string; denoting which normalization method to choose,
    see below for details}
  \item{ChIPChannel}{string; which element of the \code{RGList} holds
    the  ChIP result, see details}
  \item{inputChannel}{string; which element of the \code{RGList} holds
    the  untreated \emph{input} sample; see details}
  \item{returnMAList}{logical; should an MAList object be returned?
    Default is to return an ExpressionSet object.}
  \item{idColumn}{string; indicating which column of the \code{genes}
    data.frame of the RGList holds the identifier for reporters on the
    microarray. This column, after calling
    \code{\link[base]{make.names}} on it, will make up the unique
    \code{featureNames} of the resulting \code{ExpressionSet}.
    If argument \code{returnMAList} is \code{TRUE}, this argument is
    ignored.}
  \item{verbose}{logical; progress output to STDOUT?}
  \item{\dots}{further arguments to be passed on
    \code{normalizeWithinArrays}
    and \code{normalizeBetweenArrays}}
}
\details{
 The procedure and called \code{limma} functions depend on the choice of
 method.
 \describe{
   \item{loess}{Calls \code{normalizeWithinArrays} with
     \code{method="loess"}.}
   \item{vsn}{Calls \code{normalizeBetweenArrays} with
     \code{method="vsn"}.}
   \item{Gquantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="Gquantile"}.}
   \item{Rquantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="Rquantile"}.}
   \item{median}{Calls \code{normalizeWithinArrays} with
     \code{method="median".}}
   \item{nimblegen}{Scaling procedure used by Nimblegen. Yields
     scaled log-ratios by a two step procedure:
     srat = log2(R) - log2(G)
     srat = srat - tukey.biweight(srat)}
   \item{Gvsn}{Learns \code{vsn} model on green channel intensities
     only and applies that transformation to both channels before
     computing fold changes.}
   \item{Rvsn}{Learns \code{vsn} model on red channel intensities
     only and applies that transformation to both channels before
     computing fold changes.}
   \item{none}{No normalization of probe intensities,
     takes raw \code{log2(R)-log2(G)} as component \code{M}
     and \code{(log2(R)+log2(G))/2} as component \code{A};
     uses \code{normalizeWithinArrays} with \code{method="none"}.}
 }

  Mostly with two-color ChIP-chip, the ChIP sample is marked with the
  red Cy5 dye and for the untreated \emph{input} sample the green Cy3
  dye is used. In that case the RGList\code{myRG}'s element \code{R}
  holds the ChIP data, and element \code{G} holds the input data.
  If this is not the case with your data, use the  arguments
  \code{ChIPChannel} and \code{inputChannel} to specify the respective
  elements of \code{myRG}.
}
\value{
  Returns normalized, transformed values as an object of class
  \code{ExpressionList} or \code{MAList}.
}
\note{
  Since Ringo version 1.5.6, this function does not call limma's
  function \code{\link[limma]{backgroundCorrect}} directly any
  longer. If wanted by the user, background correction should be
  indicated as additional arguments passed on to
  \code{\link[limma]{normalizeWithinArrays}} or
  \code{\link[limma]{normalizeBetweenArrays}}, or alternatively call
  \code{\link[limma]{backgroundCorrect}} on the \code{RGList} before
  \code{preprocess}ing.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link[limma]{normalizeWithinArrays}},
  \code{\link[limma]{normalizeBetweenArrays}},
  \code{\link[limma]{malist}},\code{\link[Biobase]{ExpressionSet}},
  \code{\link[vsn]{vsnMatrix}}}
\examples{
   exDir <- system.file("exData",package="Ringo")
   exRG <- readNimblegen("example_targets.txt","spottypes.txt",
                         path=exDir)
   exampleX <- preprocess(exRG)
   sampleNames(exampleX) <- make.names(paste(exRG$targets$Cy5,"vs",
                                        exRG$targets$Cy3,sep="_"))
   print(exampleX)
   ### compare VSN to NimbleGen's tukey-biweight scaling
   exampleX.NG <- preprocess(exRG, method="nimblegen")
   sampleNames(exampleX.NG) <- sampleNames(exampleX)
   if (interactive())
     corPlot(cbind(exprs(exampleX),exprs(exampleX.NG)),
       grouping=c("VSN normalized","Tukey-biweight scaled"))
}
\keyword{manip}% at least one, from doc/KEYWORDS
