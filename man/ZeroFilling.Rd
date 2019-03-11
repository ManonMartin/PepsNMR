\name{ZeroFilling}
\alias{ZeroFilling}
\title{Zero Filling}
\description{The function applies zero filling to the FIDs.}
\usage{
ZeroFilling(Fid_data, fn = ncol(Fid_data), verbose = FALSE)
}
\arguments{
  \item{Fid_data}{Matrix containing the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{fn}{Number of 0 to be added.}
  \item{verbose}{If \code{"TRUE"}, will print processing information.}
}
\details{

Zero filling does not improve the spectral resolution but lead to better visually defined lines in the spectra.
During zero filling, \code{fn} zeros are appended at the end of the FIDs. This number is rounded to the nearest 2^x value to ease the upcoming Fourier Transform of the FIDs.
}
\value{
  \item{Fid_data}{The zero-filled FIDs.}
}

\author{
Manon Martin
}

\examples{
require(PepsNMRData)
ZF_fid <- ZeroFilling(Data_HS_sp$FidData_HS_3, fn = ncol(Data_HS_sp$FidData_HS_3))
}

\keyword{ manip }
