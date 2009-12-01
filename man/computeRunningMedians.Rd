\name{computeRunningMedians}
\alias{computeRunningMedians}
\title{Function to compute running medians on a tiling expression set}
\description{
  Function to compute running medians (or other quantiles) on a tiling
  expression set.
}
\usage{
computeRunningMedians(xSet, probeAnno, modColumn = "Cy5",
 allChr, winHalfSize = 400, min.probes = 5, quant = 0.5,
 combineReplicates = FALSE, nameSuffix = ".sm", checkUnique=TRUE,
 uniqueCodes=c(0), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xSet}{Object of class \code{ExpressionSet} holding the
    normalized probe intensity data}
  \item{probeAnno}{Environment holding the genomic positions of probes
    in the ExpressionSet}
  \item{modColumn}{Column of the ExpressionSet's \code{phenoData}
    holding the samples' difference of interest}
  \item{allChr}{Character vector of all chromosomes in genome; if not
    specified (defaul) all chromosomes annotated in the supplied
    \code{probeAnno} are used.}
  \item{winHalfSize}{Half the size of the window centered at a probe
  position, in which all other probes contribute to the calculation of
  the median.}
  \item{min.probes}{integer; if less probes are in the sliding window,
    NA instead of the median is returned. This meant to avoid to
    computing non-meaningful medians. If unwanted, set this to 1 or
    less} 
  \item{quant}{numeric; which quantile to use for the smoothing.
    The default \code{0.5} means compute the median over the values in
    the sliding window.}
  \item{combineReplicates}{logical; should the median not be computed
    over individual samples in the ExpressionSet, but should samples be
    combined according to the column \code{modColumn} of the
    \code{phenoData}. The median is then computed across all probe
    levels and samples of the same type in the window. The resulting
    ExpressionSet has so many columns as are there different entries in
    the column \code{modColumn}}
  \item{nameSuffix}{character; suffix attached to the sample labels of
    the supplied ExpressionSet \code{xSet} for the sample names of the
    resulting ExpressionSet.}
  \item{checkUnique}{logical; indicates whether the uniqueness
    indicator of probe matches from the probeAnno environment should be
    used.}
  \item{uniqueCodes}{numeric; which numeric codes in the chromosome-wise
    match-uniqueness elements of the probeAnno environment indicate
    uniqueness?}
  \item{verbose}{logical; detailed progress output to STDOUT?}
}
\value{
  An object of class \code{ExpressionSet}, holding smoothed intensity
  values for the probes of the supplied ExpressionSet.
  The number of results samples is either the number of levels in the
  supplied \code{modColumn} of the supplied ExpressionSet's phenoData if
  \code{combineReplicates} is set to TRUE or equal to the number of
  samples in the supplied ExpressionSet \code{xSet} otherwise.
}
\author{Joern Toedling \email{joern.toedling@curie.fr}}
\seealso{\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link{sliding.quantile}}, \code{\link{probeAnno-class}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
                                   winHalfSize = 400)
  combX <- combine(exampleX, smoothX)
  if (interactive()){
   grid.newpage()
   plot(combX, exProbeAnno, chrom="9", xlim=c(34318000,34321000),
        ylim=c(-2,4), gff=exGFF)
  }
}
\keyword{manip}% at least one, from doc/KEYWORDS
