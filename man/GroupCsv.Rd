\name{GroupCsv}
\alias{GroupCsv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the numerical vector of spectra classes 
}
\description{
Returns the classes of spectra from a csv file.
}
\usage{
GroupCsv(Spectrum_data, gr.file, csv2 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra, one row per spectrum.}
  \item{gr.file}{If non \code{NULL}, character string with the complete path to the .csv file with the class information.}
  \item{csv2}{If \code{TRUE}, with read with the comma as decimal and the semicolon as field separator.}
}

\value{
A numerical vector with classes of spectra.
}

\author{
Manon Martin
}

\examples{
path <-  system.file("extdata/Group_HS.csv", package = "SOAP")
GroupCsv(FinalSpectra_HS, gr.file = path)
}

\keyword{ manip}


