\name{FidData_HU}
\alias{FidData_HU}
\docType{data}
\title{Raw FIDs for the Human Urine database.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing the raw Free Indiction Decays of the Human Urine database.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FidData_HU")}
\format{
  A complex matrix of 24 spectra with 29411 time points.
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(FidData_HU)
plot(Re(FidData_HU[1,]), type="l") 
}
\keyword{datasets}
