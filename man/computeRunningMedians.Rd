\name{computeRunningMedians}
\alias{computeRunningMedians}
\title{Function to compute running medians on a tiling expression set}
\description{
  Function to compute running medians on a tiling expression set 
}
\usage{
computeRunningMedians(xSet, probeAnno, modColumn = "Cy5",
allChr = c(1:19, "X", "Y"), winHalfSize = 400, min.probes = 5,
combineReplicates = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xSet}{Object of class \code{ExpressionSet} holding the
    normalized probe intensity data}
  \item{probeAnno}{Environment holding the genomic positions of probes
    in the ExpressionSet}
  \item{modColumn}{Column of the ExpressionSet's \code{phenoData}
    holding the samples' difference of interest}
  \item{allChr}{Character vector of all chromosomes in genome}
  \item{winHalfSize}{Half the size of the window centered at a probe
  position, in which all other probes contribute to the calculation of
  the median.}
  \item{min.probes}{integer; if less probes are in the sliding window,
    NA instead of the median is returned. This meant to avoid to
    computing non-meaningful medians. If unwanted, set this to 1 or
    less}
  \item{combineReplicates}{logical; should the median not be computed
    over individual samples in the ExpressionSet, but should samples be
    combined according to the column \code{modColumn} of the
    \code{phenoData}. The median is then computed across all probe
    levels and samples of the same type in the window. The resulting
    ExpressionSet has so many columns as are there different entries in
    the column \code{modColumn}}
  \item{verbose}{logical; detailed progress output to STDOUT?}
}
\value{
  An object of class \code{ExpressionSet}, holding smoothed intentesity
  values for the probes of the supplied ExpressionSet. The number of
  results samples is the number of levels in the supplied
  \code{modColumn} of the supplied ExpressionSet's phenoData.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link[Biobase]{ExpressionSet}},
  \code{\link{sliding.median}}}
\examples{
  current.dir <- getwd()
  setwd(system.file("exData",package="Ringo"))
  load("exampleProbeAnno.rda")
  load("exampleX.rda")
  setwd(current.dir)
  smoothX = computeRunningMedians(exampleX, probeAnno=exProbeAnno,
modColumn = "Cy5", allChr = c("8"), winHalfSize = 400, verbose = TRUE)
  par(mfrow=c(1,1))
  chipAlongChrom(exampleX, chrom="8", xlim=c(60208000,60216000),
                 ylim=c(-1,6), probeAnno=exProbeAnno, gff=exGFF)  
  chipAlongChrom(smoothX, chrom="8", xlim=c(60208000,60216000),
                 ylim=c(-1,6), probeAnno=exProbeAnno, gff=exGFF,
                 itype="l", ilwd=3, add=TRUE)

  ## combine over replicates:
  smoothX2 = computeRunningMedians(exampleX, probeAnno=exProbeAnno,
       modColumn = "Cy5", allChr = c("8"), winHalfSize = 400,
       combineReplicates=TRUE, verbose = TRUE)
  chipAlongChrom(exampleX, chrom="8", xlim=c(60208000,60216000),
                 ylim=c(-1,6), probeAnno=exProbeAnno, gff=exGFF)  
  chipAlongChrom(smoothX2, chrom="8", xlim=c(60208000,60216000),
                 ylim=c(-1,6), probeAnno=exProbeAnno, gff=exGFF,
                 itype="l", ilwd=3, paletteName="Spectral", add=TRUE)
}
\keyword{manip}% at least one, from doc/KEYWORDS
