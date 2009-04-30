\name{relateChers}
\alias{relateChers}
\title{Relate found Chers to genomic features}
\description{
  This function relates found 'cher's (ChIP-enriched regions) to
  annotated genomic features, such as transcripts.
}
\usage{
relateChers(pl, gff, upstream = 5000, verbose = TRUE)
}
\arguments{
  \item{pl}{Object of class \code{cherList}}
  \item{gff}{\code{data.frame} holding genomic feature annotation}
  \item{upstream}{up to how many bases upstream of annotated genomic
    features should chers be counted as related to that feature (see
    details)}
  \item{verbose}{logical; extended output to STDOUT?}
}
\details{
  chers will be counted as related to genomic features, if
  \itemize{
    \item their middle position is located between start and end
    position of the feature
    \item their middle position is located not more than argument
    \code{upstream} bases upstream of the feature start
  }.
  
  One can visualize such cher-feature relations as a graph using
  the Bioconductor package \code{Rgraphviz}.
  See the script 'graphChers2Transcripts.R' in Ringo's scripts
  directory for an example. 
}
\value{
  An object of class \code{cherList} with for each cher the elements
  \code{typeUpstream} and \code{typeInside} filled in with the names
  of the features that have been related to.
}
\author{Joern Toedling}
\examples{
  # see findChersOnSmoothed for an example
}
\keyword{manip}
