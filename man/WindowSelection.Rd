\name{WindowSelection}
\alias{WindowSelection}
\title{
Spectral window selection
}
\description{
Selects an interval in the ppm scale and only returns the value of the spectrum in that interval.


}
\usage{
WindowSelection(Spectrum_data, from.ws = 0.2, to.ws = 10)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{from.ws}{The lower ppm value of the interval. A typical value is 0.2. If NULL, default value is the lower index without NA.}
  \item{to.ws}{The higher ppm value of the interval, it should be greater or equal to from. A typical value is 10. If NULL, default value is the higher index without NA.}
}
\details{
If \code{from.ws} and/or \code{to.ws} are not specified we calculate it so that we have the largest window without NA.
Those NAs are typicaly produced by PPMconversion with misaligned TMSP peaks.
}
\value{
\item{Spectrum_data}{The matrix of the value of the spectra in the specified interval.}
}
\author{
Benoît Legat
}
\examples{
# The interval is chosen so that we have the largest interval without NA
Ws.spec <- WindowSelection(Data_HS_sp$Spectrum_data_HS_9)

# or
Ws.spec <- WindowSelection(Data_HS_sp$Spectrum_data_HS_9, from.ws=0.2, to.ws=10)
}
\keyword{ manip }

