\name{computeSlidingT}
\alias{computeSlidingT}
\title{Function to compute sliding T statistics on a tiling expression set}
\description{
  Function to compute sliding (regularized) one- or two-sample T
  statistics on a tiling expression set.
}
\usage{
computeSlidingT(xSet, probeAnno, allChr = c(1:19, "X", "Y"), test = "one.sample", grouping = NULL, winHalfSize = 400, min.probes = 5, checkUnique = TRUE, uniqueCodes = c(0), verbose = TRUE)
}
\arguments{
  \item{xSet}{Object of class \code{ExpressionSet} holding the
    normalized probe intensity data}
  \item{probeAnno}{Environment holding the genomic positions of probes
    in the ExpressionSet}
  \item{allChr}{Character vector of all chromosomes in genome}
  \item{test}{character; one of \code{one.sample} or \code{two.sample}}
  \item{grouping}{factor vector of length equal to number of samples,
    currently not required}
  \item{winHalfSize}{Half the size of the window centered at a probe
    position, in which all other probes contribute to the calculation of
    the mean and standard deviation.}
  \item{min.probes}{integer; if less probes are in the sliding window,
    NA instead of the mean and sd is returned. This is meant to avoid to 
    computing non-meaningful means and standard deviations. If unwanted,
    set this to 1 or less}
  \item{checkUnique}{logical; indicates whether the uniqueness
    indicator of probe matches from the probeAnno environment should be
    used.}
  \item{uniqueCodes}{numeric; which numeric codes in the chromosome-wise
    match-uniqueness elements of the probeAnno environment indicate
    uniqueness?}
  \item{verbose}{logical; detailed progress output to STDOUT?}
}
%\details{}
\value{
  An object of class \code{ExpressionSet}, holding the T statistics
  values for the probes of the supplied ExpressionSet. The number of
  results samples is the number of levels in the supplied factor
  \code{grouping}.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\note{the option \code{two.sample} is not implemented yet} 
\seealso{\code{\link{sliding.meansd}}}
\examples{
 exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  tX <- computeSlidingT(exampleX, probeAnno=exProbeAnno,
             allChr=c("9"), winHalfSize=400)
  sampleNames(tX) <- "t-Stat_Suz12vsTotal"
  if (interactive()){
    grid.newpage()
    plot(cbind2(exampleX, tX), exProbeAnno, chrom="9",
         xlim=c(34318000,34321000), ylim=c(-2,8.5), gff=exGFF,
         paletteName="Paired")
  }
}
\keyword{manip}
