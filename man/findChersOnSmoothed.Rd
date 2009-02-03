\name{findChersOnSmoothed}
\alias{findChersOnSmoothed}
\title{Find ChIP-enriched regions on smoothed ExpressionSet}
\description{
  Given an ExpressionSet of smoothed probe intensities, an environment
  with the mapping of probes to chromosomes, and a vector of thresholds
  for calling genomic sites enriched, this function finds the
  'chers' (ChIP-enriched regions)
  consisting of enriched genomic positions, with probes mapped to them.
  'Adjacent' enriched positions are condensed into a single Cher.
}
\usage{
findChersOnSmoothed(smoothedX, probeAnno, thresholds, allChr = NULL,
   distCutOff = 600, minProbesInRow = 3, cellType = NULL,
   antibodyColumn=NULL, checkUnique = TRUE, uniqueCodes = c(0),
   verbose = TRUE)
}
\arguments{
  \item{smoothedX}{Object of class \code{ExpressionSet} holding the
    smoothed probe intensities, e.g. the result of function
    \code{computeRunningMedians}.}
  \item{probeAnno}{environment containing the probe to genome mapping}
  \item{thresholds}{numeric vector of threshold above which smoothed
    probe intensities are considered to correspond to enriched
    probes. The vector has to be of length equal the number of samples
    in \code{smoothedX}, with a single threshold for each sample.}
  \item{allChr}{character vector of all chromosomes on which enriched
    regions are sought. Every chromosome here has to have probes mapped
    to it in the \code{probeAnno} environment. By default (\code{NULL})
    the \code{chromosomeNames} of the probeAnno object are used.}
  \item{distCutOff}{integer; maximum amount of base pairs at which
    enriched probes are condensed into one Cher.}
  \item{minProbesInRow}{integer; minimum number of enriched
    probes required for a Cher; see \code{details} for further
    explanation.}
  \item{cellType}{character; name of cell type the data comes from, is
    either a. of length one indicating the column of
    \code{pData(smoothedX)} that holds the cell type OR b. of length one
    indicating the common cell type for all samples in the
    \code{ExpressionSet} OR c. of length equal to \code{ncol(smoothedX)}
    specifying the cell type of each sample individually.}
  \item{antibodyColumn}{the name or number of the column of the
    \code{pData(smoothedX)} that holds the description of the antibody
    used for each sample. This information is used to annotate found
    ChIP-enriched regions accordingly. If \code{NULL} (default), the
    \code{sampleNames} of \code{smoothedX} are used.}
  \item{checkUnique}{logical; indicates whether the uniqueness
    indicator of probe matches from the probeAnno environment should be
    used.}
  \item{uniqueCodes}{numeric; which numeric codes in the chromosome-wise
    match-uniqueness elements of the probeAnno environment indicate
    uniqueness?}
  \item{verbose}{logical; extended output to STDOUT?}
}
\value{
  A list of class \code{cherList}, holding objects of class \code{cher}
  that were found on the supplied data.
}
\details{
   Specifying a minimum number of probes for a Cher (argument
   \code{minProbesInRow})  guarantees that a Cher is supported by a
   reasonable number of measurements in probe-sparse regions.
   For example, if there's only one enriched probe within a certain
   genomic 1kb region and no other probes can been mapped to that
   region, this single probe does arguably not provide enough evidence
   for calling this genomic region enriched.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{cherByThreshold}},\code{\link{computeRunningMedians}},
  \code{\link{relateChers}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
       modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  chersX <- findChersOnSmoothed(smoothX, probeAnno=exProbeAnno,
       thresholds=0.45, allChr="9", distCutOff=600, cellType="human")
  if (interactive())
    plot(chersX[[1]], smoothX, probeAnno=exProbeAnno, gff=exGFF)
  chersX <- relateChers(chersX, exGFF)
  as.data.frame.cherList(chersX)
}
\keyword{manip}% at least one, from doc/KEYWORDS
