\name{findPeaksOnSmoothed}
\alias{findPeaksOnSmoothed}
\title{Find peaks on smoothed ExpressionSet}
\description{
  Given an ExpressionSet of smoothed probe intensities, an environment
  with the mapping of probes to chromosomes, and a vector of thresholds
  for calling genomic sites enriched, this function finds the peaks
  consisting of enriched genomic positions, with probes mapped to them.
  'Adjacent' enriched positions are condensed into a single peak.
}
\usage{
findPeaksOnSmoothed(smoothedX, probeAnno, thresholds, allChr = c(1:19,
"X", "Y"), distCutOff = 600, minProbesInRow = 3, cellType=NULL,
verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smoothedX}{Object of class \code{ExpressionSet} holding the
    smoothed probe intensities, e.g. the result of function
    \code{computeRunningMedians}.}
  \item{probeAnno}{environment containing the probe to genome mapping}
  \item{thresholds}{numeric vector of threshold above which smoothed
    probe intensities are considered to correspond to enriched
    probes. The vector has to be of length equal the number of samples
    in \code{smoothedX}, with a single threshold for each sample.}
  \item{allChr}{character vector of all chromosomes, on which peaks are
    to be found. Every chromosome here has to have probes mapped to it
    in the \code{probeAnno} environment.}
  \item{distCutOff}{integer; maximum amount of base pairs at which
    enriched probes are condensed into one peak.}
  \item{minProbesInRow}{integer; minimum number of enriched
    probes required for a peak; see \code{details} for further
    explanation.}
  \item{cellType}{character; name of cell type the data stems from,
    maybe could be taken from \code{pData(smoothedX)}}
  \item{verbose}{logical; extended output to STDOUT?}
}
\value{
  A list of class \code{peakList}, holding objects of class \code{peak}
  that were found on the supplied data.
}
\details{
   Specifying a minimum number of probes for a peak (argument
   \code{minProbesInRow})  guarantees that a peak is supported by a
   reasonable number of measurements in probe-sparse regions.
   For example, if there's only one enriched probe within a certain
   genomic 1kb region and no other probes can been mapped to that
   region, this single probe does arguably not provide enough evidence
   for calling this genomic region enriched.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{peakByThreshold}},\code{\link{computeRunningMedians}},
  \code{\link{relatePeaks}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
     modColumn = "Cy5", allChr = c("8"), winHalfSize = 400, combineReplicates=TRUE)
  peaksX <- findPeaksOnSmoothed(smoothX, probeAnno=exProbeAnno,
thresholds=1.5, allChr = c("8"),distCutOff=600, cellType="HL1")
  if (interactive())
    plot(peaksX[[2]], smoothX, probeAnno=exProbeAnno, gff=exGFF)
  peaksX <- relatePeaks(peaksX, exGFF)
  as.data.frame.peakList(peaksX)
}
\keyword{manip}% at least one, from doc/KEYWORDS
