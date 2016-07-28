\name{ZoneAggregation}
\alias{ZoneAggregation}
\title{
Aggregates the values in a given ppm interval.
}
\description{
Replace the values given in specified intervals by triangular shaped peaks with the same area than the original peaks.

}
\usage{
ZoneAggregation(Spectrum_data, fromto.za = list(Citrate =c(2.5, 2.7)))
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{fromto.za}{List containing the borders in ppm of the intervals to aggregate.}

}
\details{
The interval is specified in the unit of the column names (which should be ppm).
This aggregation is usually performed with urine samples that contain citrate.
}
\value{
\item{Spectrum_data}{The matrix of spectra with their zone aggregated.}
}
\author{
Benoît Legat
}
\examples{
\dontrun{Spectrum_data <- ZoneAggregation(Data_HU_sp, fromto.za = list(Citrate =c(2.5, 2.7)))}

}
\keyword{ manip }

