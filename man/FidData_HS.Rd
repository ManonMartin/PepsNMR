\name{FidData_HS}
\alias{FidData_HS}
\docType{data}
\title{Raw FIDs for the Human Serum database.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing the raw Free Indiction Decays of the Human Serum database.
The experimental design is as follow: serum was collected from 4 blood donors (Class membership available in \code{Group_HS}). Measurements were performed on 8 different days with replicates. The order of measurements were permuted accordingly to a latin hypercube sampling method.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FidData_HS")}
\format{
  A complex matrix of 32 spectra with 32768 time points.
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(FidData_HS)
plot(Re(FidData_HS[1,]), type="l") 
}
\keyword{datasets}
