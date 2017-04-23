\name{FinalSpectra_HS}
\alias{FinalSpectra_HS}
\docType{data}
\title{Spectra for the Human Serum database after the advised preprocessing workflow.
}
\description{Matrix containing the Human Serum spectra after the full pre-treatment process and with the default parameters. At this stage the spectra are fully pre-processed. 
For more details on this database, see \code{help(FidData_HS)}.
}
\usage{data("FinalSpectra_HS")}

\format{A complex matrix with 32 observations and 500 frequency data points in a ppm scale.
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
}

\seealso{
See also \code{\link{FidData_HS}} for the raw FIDs and \code{\link{Group_HS}} for the samples class.
}


\examples{
data(FinalSpectra_HS)
plot(Re(FinalSpectra_HS[1,]), type="l") 
}
\keyword{datasets}
