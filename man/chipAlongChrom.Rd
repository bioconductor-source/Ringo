\name{chipAlongChrom}
\alias{chipAlongChrom}
\alias{plot,ExpressionSet,probeAnno-method}
\title{Visualize ChIP intensities along the chromosome}
\description{
  This function can visualize the array intensities from a ChIP chip
  experiment for a chromosomal region or the whole chromosome. It's
  based on the \code{plotAlongChrom} function from the package
  \code{tilingArray}, but provides a different visualization.
}
\usage{
\S4method{plot}{ExpressionSet,probeAnno}(x, y, ...)

chipAlongChrom(eSet,  probeAnno, chrom, xlim, ylim, 
 samples = NULL, paletteName = "Set2", colPal = NULL,
 ylab = "Fold change [log]", ipch = 16, ilwd = 3, ilty = 1,
 icex = 3, gff = NULL,
 featureExclude=c("chromosome", "nucleotide_match","insertion"),
 zeroLine = TRUE, sampleLegend = TRUE, sampleLegendPos = "topleft",
 featureLegend = FALSE, maxInterDistance = 200,  coord = NULL,
 highlight, main, ...)
}
\arguments{
  \item{eSet}{An expression set containing the (normalized) ChIP
    intensities, e.g. the result objects from functions
    \code{preprocess} and \code{computeRunningMedians}.}
  \item{x}{Corresponds to argument \code{eSet} when calling the S4
    method} 
  \item{probeAnno}{An object of class \code{probeAnno} holding genomic
    position, index and gene association of probes on array.}
  \item{y}{Corresponds to argument \code{probeAnno} when calling the S4
    method}
  \item{chrom}{character; the chromosome to visualize}
  \item{xlim}{start and end genomic coordinates on the chromosome to
    visualize}
  \item{ylim}{minimum and maximum probe intensities for the plot, if
    missing (default) set to \code{range(exprs(eSet))}}
  \item{samples}{numeric; which samples from the \code{eSet} are to be
    shown. Default is to show all samples in the \code{eSet},}
  \item{paletteName}{character; Name of the RColorBrewer palette to use
    for sample colors. If the number of samples is greater than the palette
    size, random colors are taken.}
  \item{colPal}{vector of colors to use for the sample intensities.
    This is alternative to the argument \code{paletteName}
    for specifying which colors to use.}
   \item{ylab}{character; label for the y-axis, passed on to the
    plotting function as \code{ylab}}
  \item{ipch}{plot character to use}
  \item{icex}{character expansion to use for plotting symbol}
  \item{ilwd}{line width of plotted data lines}
  \item{ilty}{line type of plotted data lines; passed on
    to \code{par(lty)}.}
  \item{gff}{data frame containing annotation for genomic feature to be
    used to further annotate the plot.}
  \item{featureExclude}{character vector specifying the feature types
    in the data.frame \code{gff} that should not be shown in the plot}
  \item{zeroLine}{logical; should a dashed horizontal line at y=0 be
    put into the plot?}
  \item{sampleLegend}{logical; should a sample legend be put into the
    plot?} 
  \item{sampleLegendPos}{character; where to put the sample legend; one
    of \sQuote{topleft} (default), \sQuote{bottomleft},
    \sQuote{topright}, or \sQuote{bottomrigth}}
  \item{featureLegend}{logical; should a feature legend be put beneath
    the plot?}
  \item{maxInterDistance}{numeric; only used when \code{itype} is either
    "r" or "u"; specifies the maximal distance up to which
    adjacent probe positions should be connected by a line.}
  \item{coord}{optional integer of length 2; can be used instead of
    \code{xlim} to specify the start and end coordinates of the genomic
    region to plot}
 \item{highlight}{optional list specifying a genomic region to be
   highlighted in the shown plot}
 \item{main}{optional main title for the plot; if not specified: the
   default is \sQuote{Chromosome coordinate [bp]}}
  \item{\dots}{further parameters passed on to \code{grid.polyline} and
    \code{grid.points}}
}
\value{
  \code{\link{invisible}} list of probe positions (element \code{x}) and
  probe levels (element \code{y}) in the selected genomic region.
}
%\details{}
\author{Joern Toedling \email{joern.toedling@curie.fr}}
\note{
  The S4 method is provided as a mere convenience wrapper.

  When plotting a new \sQuote{grid} plot in an active x11 window that
  already contains a plot, remember to call \code{grid.newpage()}
  before.
}
\seealso{
  \code{\link[Biobase]{ExpressionSet-class}},
  \code{\link{probeAnno-class}},
  \code{\link[grid]{grid.points}},
  \code{plotAlongChrom} in package \code{tilingArray}}
\examples{
  ### load data
  ringoExampleDir <- system.file("exData",package="Ringo")
  load(file.path(ringoExampleDir,"exampleProbeAnno.rda"))
  load(file.path(ringoExampleDir,"exampleX.rda"))

  ### show a gene that is well represented on this microarray
  plot(exampleX, exProbeAnno, chrom="9",
       xlim=c(34318000,34321000), ylim=c(-2,4), gff=exGFF)

  ### this should give you the same result as:
  chipAlongChrom(exampleX, chrom="9", xlim=c(34318000,34321000),
                 ylim=c(-2,4), probeAnno=exProbeAnno, gff=exGFF)

}
\keyword{hplot}% at least one, from doc/KEYWORDS
