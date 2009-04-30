\name{image.RGList}
\alias{arrayImage}
\alias{image.RGList}
\title{Function to visualize spatial distribution of raw intensities}
\description{
  Function to visualize spatial distribution of raw intensities on
  NimbleGen Oligoarrays. Requires \code{RGList} with component
  \code{genes} complete with \code{genes$X} and \code{genes$X}
  coordinates of probes on array.
  \code{arrayImage} is a synonym of \code{image.RGList}.
}
\usage{
\method{image}{RGList}(x,arrayno,channel=c("red","green","logratio"),
         mycols=NULL, mybreaks=NULL, dim1="X", dim2="Y",
         ppch=20, pcex=0.3, verbose=TRUE, ...)
}
\arguments{
  \item{x}{object of class \code{RGList} containing red and green
    channel raw intensities; possibly result of \code{readNimblegen}.}
  \item{arrayno}{integer; which array to plot; one of \code{1:ncol(x$R)}}
  \item{channel}{character; which channel to plot, either \code{red},
    \code{green} or the \code{logratio} \code{log2(red)-log2(green)}}
  \item{mycols}{vector of colors to use for image; if \code{NULL}
    defaults to \code{colorRampPalette(c("White", "Yellow", "Red"))(10)}}
  \item{mybreaks}{optional numeric vector of breaks to use as argument
    \code{breaks} in \code{image.default}; default \code{NULL} means take
    \code{length(mycols)+1} quantiles of the data as breaks.}
  \item{dim1}{string; which column of the 'genes' element of the
    supplied \code{RGList} indicates the first dimension of the reporter
    position on the microarray surface; for example this column is
    called 'X' with some NimbleGen arrays and 'Row' with some Agilent
    arrays.}
  \item{dim2}{string; which column of the 'genes' element of the
    supplied \code{RGList} indicates the second dimension of the
    reporter position on the microarray surface; for example this column
    is called 'Y' with some NimbleGen arrays and 'Col' with some Agilent
    arrays.}
  \item{ppch}{which symbol to use for intensities; passed on as
    \code{pch} to \code{points.default}}
  \item{pcex}{enlargement factor for intensity symbols; passed on as
    \code{cex} to \code{points.default}}
  \item{verbose}{logical; extended output to STDOUT?}
  \item{\dots}{further arguments passed on to \code{plot.default} and
    \code{points.default}}
}
\value{
  invisibly returns NULL; function is called for its side effect, this
  is producing the plot
}
\author{Joern Toedling}
\seealso{\code{\link{readNimblegen}},\code{\link[graphics]{plot.default}},
  \code{\link[graphics]{points.default}}}
\examples{
   exDir <- system.file("exData",package="Ringo")
   exRG <- readNimblegen("example_targets.txt","spottypes.txt",path=exDir)
   image(exRG, 1, channel="red", mycols=c("black","darkred","red"))
   ## this example looks strange because the example data files only
   ## includes the probe intensities of probes mapped to the forward
   ## strand of chromosome 9.
   ## you can see these probes are distributed all over the array
}
\keyword{hplot}
