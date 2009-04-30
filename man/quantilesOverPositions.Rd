\name{quantilesOverPositions}
\Rdversion{1.1}
\alias{quantilesOverPositions}
\title{
  show ChIP-chip data aligned over genome features, e.g. TSSs
}
\description{
   Function to show the ChIP-chip data aligned over certain genome
   features, for example transcription start sites (TSSs).
}
\usage{
quantilesOverPositions(xSet, selGenes, g2p,
                       positions = seq(-5000, 10000, by = 250),
                       quantiles = c(0.1, 0.5, 0.9))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xSet}{an \code{ExpressionSet} holding the ChIP-chip data}
  \item{selGenes}{character; vector of genome features,
    e.g. transcripts, to use for the plot}
  \item{g2p}{
    A list object containing the mapping between genome positions and
    probes on the microarray. Created with the function
    \code{features2Probes}.}
  \item{positions}{
    Numeric vector of positions related to the coordinates of the genome
    features, such as in which distances of the TSS the values should be
    computed over the aligned data
  }
  \item{quantiles}{
    numeric; which quantiles to compute over the aligned data}
}
%\details{}
\value{
  An object of class \code{qop}, which can be visualized by its plot
  method.
}
\author{Joern Toedling}
%\note{}
\seealso{
  \code{\link{features2Probes}}, \code{\link{qop-class}}
}
\examples{
  ringoExampleDir <- system.file("exData",package="Ringo")
  load(file.path(ringoExampleDir,"exampleProbeAnno.rda"))
  trans2Probe <- features2Probes(exGFF, exProbeAnno)
  load(file.path(ringoExampleDir,"exampleX.rda"))
  exampleSX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
     modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  exampleC <- findChersOnSmoothed(exampleSX, probeAnno=exProbeAnno,
     thresholds=0.2, allChr="9", distCutOff=600, cellType="human")
  exampleC <- relateChers(exampleC, exGFF)
  exampleQop <- quantilesOverPositions(exampleSX,
     selGenes=getFeats(exampleC), quantiles=c(0.5, 0.9),
     g2p=trans2Probe, positions=seq(-4000, 1000, by=250))
  show(exampleQop)
  plot(exampleQop, ylim=c(-0.5, 2.1))
}
\keyword{manip}
