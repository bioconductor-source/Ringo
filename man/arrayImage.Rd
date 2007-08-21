\name{image.RGList}
\alias{arrayImage}
\alias{image.RGList}
%%- Also NEED an '\alias' for EACH other topic documented here.
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
         mycols=NULL, mybreaks=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{\dots}{further arguments passed on to \code{image.default}}
}
\value{
  invisibly returns matrix of raw intensities, in which rows and columns
  correspond to X and Y coordinates on arrays. Can be directly plotted
  with function \code{image.default}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{readNimblegen}},\code{\link[graphics]{image}}}
\examples{
   exDir <- system.file("exData",package="Ringo")
   exRG <- readNimblegen("example_targets.txt","spottypes.txt",path=exDir)
   image(exRG,arrayno=1,channel="red")
   ## this example looks strange because the example data files only
   ## includes the probe intensities of probes mapped to the forward
   ## strand of chromosome 9.
   ## you can see these probes are distributed all over the array
}
\keyword{hplot}% at least one, from doc/KEYWORDS
