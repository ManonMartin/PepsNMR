\name{FidData_HS}
\alias{FidData_HS}
\docType{data}
\title{Raw FIDs for the Human Serum database.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing the raw Human Serum Free Indiction Decays.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FidData_HS")}
\format{
  A complex matrix with 32 observations and 32768 time points.

}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(FidData_HS)
plot(Re(FidData_HS[1,]), type="l") 
}
\keyword{datasets}
