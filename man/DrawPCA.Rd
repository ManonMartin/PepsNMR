\name{DrawPCA}
\alias{DrawPCA}
\title{
Draw the PCA of the signals
}
\description{
Draw the PCA of the signals given as row of the matrix \code{Signal_data}.

Do not call this function directly but rather call \code{\link{Draw}} to specify how the plot will be outputted.
}
\usage{
DrawPCA(Signal_data, drawNames = TRUE, createWindow=F, class = NULL, axes =c(1,2))
}
\arguments{
  \item{Signal_data}{Matrix containing the the FIDs, RawSpectra or Spectra, as outputted by \code{\link{ReadFids}}.}
  \item{drawNames}{Specifies whether the names of the spectra have to be shown alonside the points on the graph.}
  \item{createWindow}{If \code{TRUE}, create a new window for each page.}
  \item{class}{Additional numeric vector indicating the class of each spectra.}
  \item{axes}{Two axes to draw the scores plot.}
}
\author{
Benoît Legat
}
\seealso{
See Also \code{\link{Draw}} and \code{\link{DrawSignal}}.
}
\examples{
DrawPCA(FinalSpectra_HS, drawNames = TRUE)
}
\keyword{ hplot }