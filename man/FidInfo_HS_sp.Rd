\name{FidInfo_HS_sp}
\alias{FidInfo_HS_sp}
\docType{data}
\title{Information about the Human Serum FIDs.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing acquisition parameters of the Human Serum FIDs.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FidInfo_HS_sp")}
\format{
A matrix with 4 observations and 8 variables. }

\details{
Variables are:
\describe{
   \item{TD}{Time domain size}
   \item{BYTORDA}{The endiness of stored data. If 0 -> Little Endian; if 1 -> Big Endian}
   \item{DIGMOD}{Digitization mode}
   \item{DECIM}{Decimation rate of digital filter}
   \item{DSPFVS}{DSP firmware version}
   \item{SW_h}{Sweep width in Hz}
   \item{SW}{Sweep width in ppm}
   \item{DT}{Dwell time in microseconds}
   }
}
\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
}

\examples{
data(FidInfo_HS)
}
\keyword{datasets}
