\name{FidData_HU}
\alias{FidData_HU}
\docType{data}
\title{Raw FIDs for the Human Urine database.}
\description{Matrix containing the raw Free Indiction Decays of the Human Urine database. The experimental design is as follow: urine was collected from 3 donors (their class membership is available in \code{Group_HU}).}
\usage{data("FidData_HU")}
\format{
  A complex matrix of 24 spectra with 29411 time points.
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.}

\seealso{
See also \code{\link{FidInfo_HU}} for acquisition parameters and \code{\link{Group_HU}} for the class of FID samples.
}


\examples{
data(FidData_HU)
plot(Re(FidData_HU[1,]), type = "l") 
}
\keyword{datasets}
