\name{posToProbeAnno}
\alias{posToProbeAnno}
\alias{posToProbeAnnoEnvironment}
\alias{createProbeAnno}
\title{Function for creating a probeAnno environment}
\description{
  This function allows the user to create a probeAnno environment that
  holds the mapping between probes on the array and their genomic match
  position(s). As input, the function takes either
  a.) one of NimbleGen's POS file or a similar file that holds the
  mapping of probes to the genome. OR
  b.) a \code{data.frame} holding this information
}
\usage{
posToProbeAnno(pos, chrNameColumn = "CHROMOSOME",
   probeColumn = "PROBE_ID", chrPositionColumn = "POSITION",
   lengthColumn = "LENGTH", sep="\t", genome="unknown",
   microarrayPlatform="unknown", verbose = TRUE, ...)
}
\arguments{
  \item{pos}{either a file-name that specifies the path to the POS or
    other mapping file OR a data.frame holding the mapping}
  \item{chrNameColumn}{name of the column in the file or
    \code{data.frame} that holds the chromosome name of the match}
  \item{probeColumn}{name of the column that holds the matching probe's
    unique identifier}
  \item{chrPositionColumn}{name of the column that holds the match
    genomic position/coordinate on the chromosome}
  \item{lengthColumn}{name of the column that holds the length of the 
    match position, in case of perfect match should correspond to the
    sequence length of the probe}
  \item{sep}{string; denotes the seperator between elements in the
    supplied mappings file \code{pos}; passed on to function
    \code{scan}; ignored if \code{pos} is not a filename.}
  \item{genome}{string; denotes genome (and assembly) the reporters have
    been mapped to for this probeAnno object, e.g. "M. musculus (mm9)"}
  \item{microarrayPlatform}{string; denotes the commercial or custom
    microarray platform/design that holds the reporters whose mapping is
    stored in this probeAnno object, e.g. "NimbleGen MOD SUZ12"}
  \item{verbose}{logical; should progress be written to STDOUT?}
  \item{\dots}{further arguments passed on to function \code{scan},
    which is used for reading in the file \code{pos}.}
}
\details{
  The default column names correspond to the column names in a NimbleGen
  POS file.

  For custom mappings, using the tools Exonerate, BLAT or MUMmer, the
  scripts directory of this package holds Perl scripts to generate such
  a POS file from the respective output files.
}
\value{
  The results is an object of class \code{probeAnno}.
}
\author{Joern Toedling}
\seealso{\code{\link{probeAnno-class}}, \code{\link{scan}}}
\examples{
 exPos <- read.delim(file.path(system.file("exData",package="Ringo"),
                     "MOD_2003-12-05_SUZ12_1in2.pos"),
                     header=TRUE,as.is=TRUE)
 str(exPos)
 exProbeAnno <- posToProbeAnno(exPos,
     genome="M. musculus (assembly mm8)",
     microarrayPlatform="NimbleGen 2005-06-17_Ren_MM5Tiling_Set1")
 ## is equivalent to 
 exProbeAnno2 <- posToProbeAnno(file.path(
    system.file("exData",package="Ringo"),"MOD_2003-12-05_SUZ12_1in2.pos"),
    genome="M. musculus (assembly mm8)",
    microarrayPlatform="NimbleGen 2005-06-17_Ren_MM5Tiling_Set1")
 ls(exProbeAnno)
 chromosomeNames(exProbeAnno2)
}
\keyword{manip}
