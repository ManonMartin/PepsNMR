\name{WindowSelection}
\alias{WindowSelection}
\title{
Spectral window selection
}
\description{
Selects an interval in the ppm scale and only returns the value of the spectrum in that interval.


}
\usage{
WindowSelection(Spectrum_data, from = 0.2, to = 10, reverse.axis = TRUE)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{from}{The lower ppm value of the interval. A typical value is 0.2. If NULL, default value is the lower index without NA}
  \item{to}{The higher ppm value of the interval, it should be greater or equal to from. A typical value is 10. If NULL, default value is the higher index without NA}
  \item{reverse.axis}{If \code{TRUE}, the ppm interval and the spectrum are reversed.}
}
\details{
If \code{from} and/or \code{to} are not specified we calculate it so that we have the largest window without NAs.
Those NAs are typicaly produced by PPMconversion with misaligned TMSP peaks.
}
\value{
\item{Spectrum_data}{The matrix of the value of the spectra in the specified interval.}
}
\author{
Benoît Legat
}
\examples{
# The interval is chosen so that we have the largest interval without NAs
Ws.spec <- WindowSelection(Data_HS$RawSpectData_HS_9)
# from is chosen so that we have the largest interval possible without NAs
Ws.spec <- WindowSelection(Data_HS$RawSpectData_HS_9, to=10)
Ws.spec <- WindowSelection(Data_HS$RawSpectData_HS_9, from=0.2, to=10)
}
\keyword{ manip }

