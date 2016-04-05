\name{FinalSpectra_HS}
\alias{FinalSpectra_HS}
\docType{data}
\title{Spectra for the Human Serum database after Scaling.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing the Human Serum spectra after the full pre-treatment process and with the default parameters. At this stage the spectra are fully pre-processed.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FinalSpectra_HS")}

\format{A complex matrix with 32 observations and 500 frequency data points (ppm scale).
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(FinalSpectra_HS)
plot(Re(FinalSpectra_HS[1,]), type="l") 
}
\keyword{datasets}
