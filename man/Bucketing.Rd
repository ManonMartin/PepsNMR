\name{Bucketing}
\alias{Bucketing}
\title{
Spectral data reduction
}
\description{
Reduces the number of data points by aggregating frequencies into \eqn{m}{m} buckets.

The input spectrum can have its ppm axis in increasing or decreasing order and it does not
have to be equispaced.
}
\usage{
Bucketing(Spectrum_data, m = 500, width = FALSE, ppmleft = NULL, 
            ppmright = NULL, intmeth = c("r", "t"), tolbuck = 10^-4)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{m}{The number of buckets or the buckets width. If \code{m} represents the number of buckets, should be an integer smaller or equal to the number of frequencies in \code{Spectrum_data}.}
  \item{width}{If \code{width} is not \code{FALSE}, then \code{m} represents the buckets width, otherwise, it represents the number of buckets.}
  \item{ppmleft}{Left boundary for ppm integration.}
  \item{ppmright}{Right boundary for ppm integration.}
  \item{intmeth}{Type of bucketing: rectangular (\code{"r"}) or trapezoidal (\code{"t"}). See details below.}
  \item{tolbuck}{Tolerance threshold to check if the buckets of the original spectra are of constant length.}
}
\details{
Bucketing will decrease data precision but has two main interests:
\itemize{
  \item{Facilitate the statistical analysis}
  \item{Decrease the impact of peaks misalignment between different spectra that should be aligned; assuming we are in the ideal case where they fall in the same bucket.
  Of course, the better the prior warping, the larger can be \eqn{m}{m} and the more informative the spectra will be.}
}

  The ppm interval of \code{Spectrum_data} - let's say is \eqn{[a,b]}{[a,b]} where \eqn{a > b}{a > b} - is divided into \eqn{m}{m} buckets of size \eqn{(a-b)/m}{(a-b)/m}.
  The new ppm scale contains the \eqn{m}{m} centers of these intervals.
  The value of the spectrum at these centers is the integral of the initial spectral intensity on
  this bucket using either trapezoidal (trapezium rule) or rectangular integration.


}
\value{
\item{Spectrum_data}{The matrix of spectra with their new ppm axis.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat, Bernadette Govaerts.
}
\note{
}

\seealso{
}
\examples{
Bucket.spec <- Bucketing(Data_HS_sp$RawSpectData_HS_10, m = 500)
}
\keyword{ manip }
\keyword{ array }
