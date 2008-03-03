\name{peak}
\alias{peak}
\alias{peak-class}
\alias{print.peak}
\alias{peakList-class}
\alias{peakList}
\docType{class}
\title{Class peak - A list holding ChIP-chip peak characteristics}
\description{
  An object of class \code{peak} is a list holding characteristics of an
  enrichment peak from ChIP chip data. While this is only a S3 rather
  than a formal \code{S4} class, the presence of certain list
  elements enables a standardized approach to working with enriched
  regions identified from ChIP-chip data.
  Note that the term "peak", while commonly used in ChIP-chip context,
  is slightly misleading and the term "ChIP-enriched region" would be
  more appropriate. Within such regions the actual signal could show two
  or more actual signal peaks or none at all (long plateau). For
  succinctness sake, we use the term "peak" for a ChIP-enriched region,
  though.
}
\section{Slots/List Components}{
  \describe{
    \item{name}{required character vector of length 1
    unequivocally describing the peak, e.g. "Suz12.Nudt2.upstream.peak"}
    \item{chr}{required character vector of length one, naming the
    chromosome of the peak, e.g. "9" }
    \item{start}{required integer of peak start position on the
    chromosome, e.g. 34318900 }
    \item{end}{required integer of peak start position on the
    chromosome, e.g. 34320100 }
    \item{cellType}{optional character vector describing the cell
    type the ChIP chip experiment has been done on, e.g. "human"}
    \item{antibody}{required character vector describing the
    antibody or characteristic for which fragments were supposedly
    enriched in immuno-precipitation step,
    e.g. "Suz12" for the protein Suz12}
    \item{typeUpstream}{optional character vector of features that this
    peak is located upstream of,
    e.g. the transcriptional start site of "ENST00000379158".
    See \code{\link{relatePeaks}} for details.}
    \item{typeInside}{ optional character vector of features
    that this peak is located inside of }
    \item{maxPeak}{optional numeric of maximal (smoothed) probe level
    in the peak, e.g. 2.00}
    \item{score}{optional numeric of a peak score, currently we use the
    sum of smoothed probe levels (log fold changes), e.g. 69.16}
    \item{probes}{optional character vector of probe identifiers of
    all probes with match positions in the peak }
    \item{distMid2TSS}{optional named numeric vector of distances of
    the peak middle position to features,
    e.g. TSSs of features upstream and inside; names are
    the features to which the distances are given; only
    meaningful in combination with \code{typeUpstream} and
    \code{typeInside}; e.g. 55 with name "ENST00000379158"}
    \item{upSymbol}{optional character vector of gene symbols of
    features the peak is located upstream of; supplements
    \code{typeUpstream}; e.g. "Nudt2" }
    \item{inSymbol}{optional character vector of gene symbols of
    features the peak is located upstream of; supplements
    \code{typeInside}.}
    \item{\ldots}{further list elements can be added like normal
      list elements.}
 }
}
\section{Creation}{
  Objects of class \code{peak} are created by other functions,
  such as \code{\link{findPeaksOnSmoothed}}; or by creating a list with
  the respective components and assigning it the class \code{peak},
  which is essentially, what the function \code{\link{newPeak}} does.
}
\section{Methods}{
  Methods defined for lists of class \code{peak} are
  \code{print.peak} and \code{\link{plot.peak}}.
}
\section{peakList}{
  A list in which each element is of class \code{peak}, is called
  a \code{peakList}. This class, however, is not much in use yet.
}
\seealso{
  \code{\link{plot.peak}}, \code{\link{findPeaksOnSmoothed}},
  \code{\link{relatePeaks}}
}
\author{Joern Toedling, Tammo Krueger, PT Bote}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
   modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  peaksX <- findPeaksOnSmoothed(smoothX, probeAnno=exProbeAnno,
   thresholds=0.45, allChr="9", distCutOff=600)
  peaksX <- relatePeaks(peaksX, exGFF)
  str(peaksX[[1]])
}
\keyword{classes}
\keyword{data}
