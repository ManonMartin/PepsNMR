\name{FidData_HS}
\alias{FidData_HS}
\docType{data}
\title{Raw FIDs for the Human Serum database.}
\description{Matrix containing the raw Free Induction Decays of the Human Serum database.
The experimental design is as follow: serum was collected from 4 blood donors (their class membership is available in \code{Group_HS}). Measurements were performed on 8 different days with replicates. The order of the measurements were permuted according to a latin hypercube sampling method.}
\usage{data("FidData_HS")}
\format{
  A complex matrix of 32 spectra with 32768 time points.
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.}

\seealso{
See also \code{\link{FidInfo_HS}} for acquisition parameters and \code{\link{Group_HS}} for the class of FID samples.
}

\examples{
data(FidData_HS)
plot(Re(FidData_HS[1,]), type = "l") 
}
\keyword{datasets}
