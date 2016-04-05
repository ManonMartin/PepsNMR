\name{FidInfo_HU}
\alias{FidInfo_HU}
\docType{data}
\title{Information about the FIDs.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing acquisition parameters of Human Urine FIDs.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FidInfo_HU")}
\format{
  A matrix with 24 observations and 6 variables.
    \describe{
   \item{TD}{Time domain size}
   \item{BYTORDA}{Determine the endiness of stored data. If 0 -> Little Endian; if 1 -> Big Endian}
   \item{DIGMOD}{Digitization mode}
   \item{DECIM}{Decimation rate of digital filter}
   \item{DSPFVS}{DSP firmware version}
   \item{SW_h}{Sweep width in Hz}
   \item{SW}{Sweep width in ppm}
   \item{GPRDLY}{Group Delay}
   \item{DT}{Dwell time in microseconds}
   }
}
\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
}
\references{
}
\examples{
data(FidInfo_HU)
}
\keyword{datasets}
