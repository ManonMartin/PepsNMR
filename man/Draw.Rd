\name{Draw}
\alias{Draw}
\title{
Draw a Signal or its PCA
}
\description{
Draws signals, spectra or their PCA.
}
\usage{
Draw(Signal_data, type.draw = c("signal","pca"),
     output = c("default","window","png","pdf"), dirpath = ".",
     filename = "\%003d", height = 480, width = 640, pdf.onefile = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Signal_data}{Matrix containing the FIDs, RawSpectra or Spectra, as outputted by \code{\link{ReadFids}}.}
  \item{type.draw}{Either "signal" or "pca", which calls respectively \code{\link{DrawSignal}} or \code{\link{DrawPCA}} to do the drawing.}
  \item{output}{Specifies how to display the drawings.
  \describe{
  \item{default}{The output is the default.}
  \item{window}{Create a new window for each page.}
  \item{png}{Create and save a new png image for each image.}
  \item{pdf}{Create and save a new pdf image for each image.}
  }}
  \item{dirpath}{The path to the directory where the png or pdf are outputted.}
  \item{filename}{The filenames of the png and pdf, see the argument \code{filename} in \code{\link{png}} doc.}
  \item{height}{height of the png and pdf in pixels.}
  \item{width}{width of the png and pdf in pixels.}
  \item{pdf.onefile}{Wen \code{output} is set to \code{"pdf"} and there are mutliples pages,
  if \code{pdf.onefile} is \code{TRUE}, all the pages are on the same file and
  if it is \code{FALSE} all the pages are on a different pdf file.}
  \item{\dots}{The remaining arguments are passed to \code{\link{DrawSignal}} or \code{\link{DrawPCA}}.}
}
\details{
Depending on the \code{type.draw}, it can draw each row of \code{Signal_data} in a way described by \code{subtype} or the PCA of all the signals described by the rows. It keeps the default output but it can also be a window, a png or a pdf.
}
\author{
Benoît Legat
}

\seealso{
See Also \code{\link{DrawSignal}} and \code{\link{DrawPCA}}.
}
\examples{
# Draw each signal Real part and Mod in separate png with name end001 end002, ...
# Draw the spectra
Draw(FinalSpectra_HS, type.draw = "signal", output="window", subtype="stacked")

# Draw a PCA
Draw(FinalSpectra_HS, type.draw="pca",output="window")
}
\keyword{ hplot }
