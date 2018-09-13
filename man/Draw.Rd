\name{Draw}
\alias{Draw}
\title{
Draw signals or their PCA scores/loadings.
}
\description{
Draws FIDs, spectra or their PCA scores/loadings.
}
\usage{
Draw(Signal_data, type.draw = c("signal","pca"),
     output = c("default","window","png","pdf"), 
     dirpath = ".",filename = "\%003d", height = 480, 
     width = 640, pdf.onefile = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Signal_data}{Matrix containing the FIDs or spectra, one line per FID/spectrum.}
  \item{type.draw}{Either "signal" or "pca", which calls respectively \code{\link{DrawSignal}} or \code{\link{DrawPCA}} to do the drawing.}
  \item{output}{Specifies how to display the drawings:
  \describe{
  \item{\code{default}}{The output is the default one.}
  \item{\code{window}}{Create a new window for each page.}
  \item{\code{png}}{Create and save a new png image for each image.}
  \item{\code{pdf}}{Create and save a new pdf image for each image.}
  }}
  \item{dirpath}{The path to the directory where the png or pdf are outputted.}
  \item{filename}{The filenames of the png and pdf, see argument \code{filename} in \code{grDevices::png} for more details.}
  \item{height}{Height of the png and pdf in pixels.}
  \item{width}{Width of the png and pdf in pixels.}
  \item{pdf.onefile}{Wen \code{output} is set to \code{"pdf"} and there are mutliples pages,
  if \code{pdf.onefile} is \code{TRUE}, all the pages are in the same file and
  if it is \code{FALSE} all the pages are in a different pdf file.}
  \item{\dots}{The remaining arguments are passed either to \code{\link{DrawSignal}} or \code{\link{DrawPCA}}.}
}
\details{
Depending on the \code{type.draw} value, it can draw each row of \code{Signal_data} in a way described by \code{subtype} or the PCA scores or loadings (depending on the \code{type.pca} value) of all the FIDs/spectra in \code{Signal_data}. 
}
\author{
Benoît Legat & Manon Martin
}

\seealso{
See Also \code{\link{DrawSignal}} and \code{\link{DrawPCA}}.
}
\examples{
# Draw each signal Real part and Mod in separate png with name end001 end002, ...
# Draw the spectra
require(PepsNMRData)
Draw(FinalSpectra_HS, type.draw = "signal", 
         output="window",subtype="together")

# Draw a PCA
Draw(FinalSpectra_HS, type.draw="pca",output="window")
}
\keyword{ hplot }
