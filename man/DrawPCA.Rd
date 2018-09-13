\name{DrawPCA}
\alias{DrawPCA}
\title{
Draw the PCA scores or loadings of the signals
}
\description{
The function draws the PCA scores or loadings of the FIDs/spectra given in the matrix \code{Signal_data}.

Do not call this function directly but rather call \code{\link{Draw}} to specify how the plot will be returned.
}
\usage{
DrawPCA(Signal_data, drawNames = TRUE,  main = "PCA score plot", Class = NULL, 
              axes = c(1,2), type.pca = c("scores", "loadings"), 
              loadingstype=c("l", "p"), num.stacked = 4, xlab = "rowname",
              createWindow)

}
\arguments{
  \item{Signal_data}{Matrix containing the FIDs or spectra, one line per FID/spectrum.}
  \item{drawNames}{If \code{TRUE}, the names of the spectra have to be shown alonside the points on the scores plot.}
  \item{main}{Plot title.}
  \item{Class}{Vector (numeric or character) indicating the class of each spectra. Used for scores plot only.}
  \item{axes}{Vector of score or loading numbers to be plotted. If it represents the score's numbers, only the first two elements are used.}
  \item{type.pca}{The type of plot, either \code{"scores"} or \code{"loadings"}}
  \item{loadingstype}{The type of loadings plot, either a line plot (\code{"l"}) or points with histogram-like vertical lines (\code{"p"}).}
  \item{num.stacked}{Number of stacked plots for the loadings plots.}
  \item{xlab}{Label of the x-axis of loadings plots.}
  \item{createWindow}{If \code{TRUE}, will open a new window to display the graphs.}
}
\author{
Benoît Legat & Manon Martin
}
\seealso{
See also \code{\link{Draw}} and \code{\link{DrawSignal}}.
}
\examples{
require(PepsNMRData)
# Draw loadings
DrawPCA(FinalSpectra_HS, main = "PCA loadings plot", 
      Class = NULL, axes =c(1,3, 5), type ="loadings", loadingstype="l", 
      num.stacked=4, xlab="ppm", createWindow = TRUE)

# Draw scores
class = substr(rownames(FinalSpectra_HS),5,5)
DrawPCA(FinalSpectra_HS, drawNames = TRUE, main = "PCA scores plot", 
        Class = class,  axes = c(1,2), type = "scores", createWindow = TRUE)

}
\keyword{ hplot }