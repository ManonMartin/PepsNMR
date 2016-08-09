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
DrawPCA(Signal_data, drawNames=TRUE, createWindow=F, main = "PCA score plot", class = NULL, axes =c(1,2),
                     type =c("scores", "loadings"), loadingstype=c("l", "p"), num.stacked=4, xlab="rowname")

}
\arguments{
  \item{Signal_data}{Matrix containing the the FIDs, RawSpectra or Spectra, as outputted by \code{\link{ReadFids}}.}
  \item{drawNames}{Specifies whether the names of the spectra have to be shown alonside the points on the scores plot.}
  \item{main}{Specify the plot title.}
  \item{class}{Additional vector (numeric or character) indicating the class of each spectra. Used for scores plot only}
  \item{axes}{Two or more axes to draw the scores or loadings. If for drawing scores, only the first two axes are used.}
  \item{type}{The type of plot, either \code{"scores"} or \code{"loadings"}}
  \item{loadingstype}{The type of Loadings plot, either a line plot (\code{"l"}) or points with histogram-like vertical lines (\code{"p"}).}
  \item{num.stacked}{Number of stacked plots.}
  \item{xlab}{Label of the x-axis of loadings plot.}
  
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