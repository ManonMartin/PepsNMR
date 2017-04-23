\name{FidInfo_HU}
\alias{FidInfo_HU}
\docType{data}
\title{Information about the FIDs for the Human Urine database.}
\description{Matrix containing acquisition parameters of the Human Urine FIDs.}
\usage{data("FidInfo_HU")}
\format{
  A matrix with 24 observations and 10 variables:
    \describe{
   \item{TD}{Time domain size}
   \item{BYTORDA}{Determine the endiness of stored data. If 0 -> Little Endian; if 1 -> Big Endian}
   \item{DIGMOD}{Digitization mode}
   \item{DECIM}{Decimation rate of digital filter}
   \item{DSPFVS}{DSP firmware version}
   \item{SW_h}{Sweep width in Hz}
   \item{SW}{Sweep width in ppm}
   \item{O1}{Spectrometer frequency offset}
   \item{GPRDLY}{Group Delay}
   \item{DT}{Dwell time in microseconds}
   }
}
\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
}

\seealso{
See also \code{\link{FidData_HU}} for the FIDs and \code{\link{Group_HU}} for the class membership of FID samples.}

\examples{
data(FidInfo_HU)
}
\keyword{datasets}
