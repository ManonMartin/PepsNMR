\name{DrawSignal}
\alias{DrawSignal}
\title{
Draw Signals
}
\description{
Depending on the \code{subtype}, will draw the different parts of the complex FIDs/spectra. 
}
\usage{
DrawSignal(Signal_data, subtype = c("stacked", "together", "separate", 
           "diffmean", "diffmedian", "diffwith"), ReImModArg = c(TRUE, 
           FALSE, FALSE, FALSE), vertical = T , xlab = "index", 
           RowNames = NULL, row = 1, num.stacked = 4, main = NULL, 
           createWindow)
} 

\arguments{
  \item{Signal_data}{Matrix containing the FIDs or spectra, one line per FID/spectrum.}
  \item{subtype}{ Specifies the drawing array:
  \describe{
  \item{\code{together}}{Plots all the signals in the same plot.}
  \item{\code{separate}}{Plots each signal on a different page.}
  \item{\code{stacked}}{Plots \code{num.stacked} signals on stacked plots with the same x-axis.}
  \item{\code{diffmean}}{Plots all the signals in the same plot but subtracted by their mean at each point.}
  \item{\code{diffmedian}}{Plots all the signals in the same plot but subtracted by their median at each point.}
  \item{\code{diffwith}}{Plots all the signals in the same plot but subtracted by the \code{row}\eqn{^{th}}{^th} signal at each point.}
  }
  }
  \item{ReImModArg}{Specifies which of the real, imaginary, modulus, or argument part of the complex signal has to be plotted. Those plots are on the same page.}
  \item{vertical}{Specifies whether the parts of the complex signal have to be put vertically or horizontally on the page if there are only 2 parts.
  If more, there will be 2 horizontally and 2 vertically anyway.}
  \item{xlab}{Label of the x-axis.}
  \item{RowNames}{Strings to use instead of the rownames as labels for the plots if \code{subtype = "separate"}. It should be a vector of the same length than the number of FIDs.}
  \item{row}{\code{row} to be compared to if the \code{subtype} is \code{"diffwith"}.}
  \item{num.stacked}{Number of stacked plots if \code{subtype} is \code{"stacked"}.}
  \item{main}{If not \code{NULL}, the main title when \code{subtype} is different from  \code{"separate"}.}
  \item{createWindow}{If \code{TRUE}, will open a new window.}
}
\details{

Don't call this function directly but rather call \code{\link{Draw}} to specify how the plot will be outputted.
}
\author{
Benoît Legat & Manon Martin
}
\seealso{
See also \code{\link{Draw}} and \code{\link{DrawPCA}}.
}
\examples{

plots <- DrawSignal(FinalSpectra_HS[1:4,], subtype = "together",
            ReImModArg = c(TRUE, TRUE, FALSE, FALSE), createWindow = TRUE)

grid::grid.draw(plots)

}
\keyword{ hplot }