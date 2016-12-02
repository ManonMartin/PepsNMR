\name{DrawSignal}
\alias{DrawSignal}
\title{
Draw Signals
}
\description{
Depending on the \code{subtype}, will draw the different parts of complex signals. 
}
\usage{
DrawSignal(Signal_data, subtype = c("stacked", "together", "separate", 
           "diffmean", "diffmedian", "diffwith"), ReImModArg = c(TRUE, 
           FALSE, FALSE, FALSE), vertical = T , xlab = "rowname", RowNames = NULL, 
           nticks = 42,row = 1, num.stacked = 4, main.title = NULL, createWindow)
} 

\arguments{
  \item{Signal_data}{Matrix containing the signals or spectra, one line per signal/spectrum.}
  \item{subtype}{ Specifies the drawing array:
  \describe{
  \item{together}{Plots all the signals in the same plot.}
  \item{separate}{Plots each signal on a different page.}
  \item{stacked}{Plots \code{num.stacked} signals on stacked plots with the same x-axis.}
  \item{diffmean}{Plots all the signals in the same plot but substracted by their mean at each point.}
  \item{diffmedian}{Plots all the signals in the same plot but substracted by their median at each point.}
  \item{diffwith}{Plots all the signals in the same plot but substracted by the \code{row}\eqn{^{th}}{^th} signal at each point.}
  }
  }
  \item{ReImModArg}{Specifies which of the real, imaginary, modulus, or argument part of the complex signal has to be outputted. Those plots are on the same page.}
  \item{vertical}{Specifies whether the parts of the complex signal have to be put vertically or horizontally on the page if there are only 2 parts.
  If more, there will be 2 horizontally and 2 vertically anyway.}
  \item{xlab}{Label of the x-axis.}
  \item{RowNames}{Strings to use instead of the rownames as labels for the plots if\code{subtype = "separate"}. It should be a vector of the same length than the number of FIDs.}
  \item{nticks}{Number of coordinates to display on the x-axis.}
  \item{row}{\code{row} to be compared to if the \code{subtype} is \code{"diffwith"}.}
  \item{num.stacked}{Number of stacked plots if \code{subtype} is \code{"stacked"}.}
  \item{main.title}{If not \code{NULL}, the main title when \code{subtype} is different from  \code{"separate"}.}
  \item{createWindow}{If \code{TRUE}, will open a new window.}
}
\details{

Each page contains one or several plot if several part of the complex signal are to be outputted.

Don't call this function directly but rather call \code{\link{Draw}} to specify how the plot will be outputted.
}
\author{
Benoît Legat
}
\seealso{
See Also \code{\link{Draw}} and \code{\link{DrawPCA}}.
}
\examples{

DrawSignal(FinalSpectra_HS[1:4,], subtype = "together",
            ReImModArg = c(TRUE, TRUE, FALSE, FALSE), createWindow = TRUE)

}
\keyword{ hplot }