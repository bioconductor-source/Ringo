\name{relatePeaks}
\alias{relatePeaks}
\title{Relate found peaks to genomic features}
\description{
  This function relates found peaks to annotated genomic features,
  such as transcripts.
}
\usage{
relatePeaks(pl, gff, upstream = 5000, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pl}{Object of class \code{peakList}}
  \item{gff}{\code{data.frame} holding genomic feature annotation}
  \item{upstream}{up to how many bases upstream of annotated genomic
    features should peaks be counted as related to that feature (see
    details)}
  \item{verbose}{logical; extended output to STDOUT?}
}
\details{
  Peaks will be counted as related to genomic features, if
  \itemize{
    \item their middle position is located between start and end
    position of the feature
    \item their middle position is located not more than argument
    \code{upstream} bases upstream of the feature start
  }.
  
  One can visualize such peak-feature relations as a graph using
  the Bioconductor package \code{Rgraphviz}.
  See the script 'graphPeaks2Transcripts.R' in Ringo's scripts
  directory for an example. 
}
\value{
  An object of class \code{peakList} with for each peak the elements
  \code{typeUpstream} and \code{typeInside} filled in with the names
  of the features that have been related to.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\examples{
  # see findPeaksOnSmoothed for an example
}
\keyword{manip}
