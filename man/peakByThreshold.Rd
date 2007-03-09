\name{peakByThreshold}
\alias{peakByThreshold}
\title{Function to identify peaks based on thresholds}
\description{
  Given a vector of probe positions on the chromosome, a vector of
  smoothed intensities on these positions, and a threshold for
  intensities to indicated enrichment, this function identifies enrichment
  peaks on this chromosome.
  
  This function is called by the function \code{findPeaksOnSmoothed}.
}
\usage{
peakByThreshold(positions, scores, threshold, distCutOff,
  minProbesInRow = 3)
}
\arguments{
  \item{positions}{numeric vector of genomic positions of probes}
  \item{scores}{scores (intensities) of probes on those positions}
  \item{threshold}{threshold for scores to be called a peak}
  \item{distCutOff}{maximal positional distance between two probes to be
    part of the same peak}
 \item{minProbesInRow}{integer; minimum number of enriched
    probes required for a peak; see \code{details} for further
    explanation.}
}
\value{
  A LIST with n components, where the first n components are the peak
  clusters, each one holding the scores and, as their names, the genomic
  positions of probes in that cluster.
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
\seealso{\code{\link{findPeaksOnSmoothed}}}
\examples{
 ## example with random generated data:
 rpos <- cumsum(round(runif(200)*5))
 rsco <- rnorm(200)+0.2
 plot(rpos, rsco, type="l", col="seagreen3", lwd=2)
 rug(rpos, side=1, lwd=2); abline(h=0, lty=2)              
 rpeaks <- peakByThreshold(rpos, rsco, threshold=0, distCutOff=2)
 sapply(rpeaks[-length(rpeaks)], function(thisClust){
  points(x=as.numeric(names(thisClust)), y=thisClust, type="h", lwd=2,
 col="gold")})
}
\keyword{manip}% at least one, from doc/KEYWORDS
