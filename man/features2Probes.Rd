\name{features2Probes}
\alias{features2Probes}
\title{Function for mapping genomic features to probes}
\description{
  This function creates a mapping between annotated genomic features and
  probes on the array whose matching genomic positions are stored in a
  probeAnno environment.
}
\usage{
features2Probes(gff, probeAnno, upstream = 5000, checkUnique = TRUE, uniqueCodes = c(0), mem.limit=1e8, verbose = TRUE)
}
\arguments{
  \item{gff}{\code{data.frame} holding genomic feature annotation}
  \item{probeAnno}{Object of class \code{environment} holding the
    genomic positions of probes in the ExpressionSet}
  \item{upstream}{up to how many bases upstream of annotated genomic
    features should probes be counted as related to that feature (see
    details)}
  \item{checkUnique}{logical; indicates whether the uniqueness
    indicator of probe matches from the probeAnno environment should be
    used.}
  \item{uniqueCodes}{numeric; which numeric codes in the chromosome-wise
    match-uniqueness elements of the probeAnno environment indicate
    uniqueness?}
  \item{mem.limit}{integer value; what is the maximal allowed size of
    matrices during the computation; see \code{regionOverlap}}
  \item{verbose}{logical; detailed progress output to STDOUT?}
}
\value{
  The results is a list of length equal to the number of rows in the
  provided \code{gff}, the data.frame of genomic features.
  The \code{names} of the list are the names specified in the
  \code{gff}. Each element of the list is specified by the probes
  mapping into the genomic region from \code{upstream} bases upstream of
  the feature's start site to the feature's end site. The entries itself
  are either \code{NULL}, if no probe was mapped into this region, or a
  named numeric vector, with its values being the distances of the
  probes' middle positions to the feature's start site (which depends on
  the strand the feature is on) and its names being the identifiers of
  these probes.
}
\author{Joern Toedling}
\note{
  This resulting mapping is not used excessively by other Ringo
  functions, so creating this mapping is optional at this time,
  but it may simplify subsequent gene/transcript-based analyses.
  
  Here, the term \emph{feature} describes a genomic entity such as a
  gene, transcript, non-coding RNA or a similar feature annotated to a
  genome. It does NOT refer to oligo-nucleotide or cDNA probes on the
  microarray.
}
\seealso{\code{\link{regionOverlap}}}
\examples{
   ringoExampleDir <- system.file("exData",package="Ringo")
   load(file.path(ringoExampleDir,"exampleProbeAnno.rda"))
   trans2Probe <- features2Probes(exGFF, exProbeAnno)
   trans2Probe[exGFF$name[match("NUDT2", exGFF$symbol)]]
   exGFF[match(names(trans2Probe)[listLen(trans2Probe)>0],exGFF$name),]
   trans2Probe[listLen(trans2Probe)==1]
}
\keyword{manip}
