\name{FinalSpectra_HU}
\alias{FinalSpectra_HU}
\docType{data}
\title{Spectra for the Human Urine database after Scaling.
%%   ~~ data name/kind ... ~~
}
\description{Matrix containing the Human Urine spectra after the full pre-treatment process and with the default parameters. At this stage the spectra are fully pre-processed.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("FinalSpectra_HU")}

\format{A complex matrix with 24 observations and 500 frequency data points (ppm scale).
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
data(FinalSpectra_HU)
plot(Re(FinalSpectra_HU[1,]), type="l") 
}
\keyword{datasets}
