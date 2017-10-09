\name{WindowSelection}
\alias{WindowSelection}
\title{
Spectral window selection
}
\description{
Selects an interval in the ppm scale and returns the value of the spectra in that interval.


}
\usage{
WindowSelection(Spectrum_data, from.ws = 10, to.ws = 0.2)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{from.ws}{The left ppm value of the interval. A typical value is 10. If NULL, default value is the first index without NA.}
  \item{to.ws}{The right ppm value of the interval. A typical value is 0.2. If NULL, default value is the last index without NA.}
}
\details{
If \code{from.ws} and/or \code{to.ws} are not specified we calculate it so that we have the largest window without NA.
Those NAs are typically produced by the \code{InternalReferencing} function.
}
\value{
\item{Spectrum_data}{The matrix of the value of the spectra in the specified interval.}
}
\author{
Benoît Legat & Manon Martin
}
\examples{
# The interval is chosen so that we have the largest interval without NA
Ws.spec <- WindowSelection(Data_HS_sp$Spectrum_data_HS_9)

# or
Ws.spec <- WindowSelection(Data_HS_sp$Spectrum_data_HS_9, from.ws=10, to.ws=0.2)
}
\keyword{ manip }

