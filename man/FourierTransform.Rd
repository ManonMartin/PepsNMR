\name{FourierTransform}
\alias{FourierTransform}
\title{
Provides the Fourier transform of the FID.
}
\description{
The function takes as input the values of the FIDs in the time domain and computes these values in the frequency domain.

}
\usage{
FourierTransform(Fid_data, Fid_info = NULL, SW_h = NULL, reverse.axis = TRUE)
}
\arguments{
  \item{Fid_data}{Matrix containing the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{Fid_info}{Matrix containing the info about the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{SW_h}{Sweep Width in Hertz, \emph{i.e.} the period of the Fourier Transform. It is equal to \eqn{1/DT} where \eqn{DT} is the Dwell Time, the time between two points in the time domain. If given, the value in \code{Fid_info} is ignored.}
  \item{reverse.axis}{If \code{TRUE}, the frequency scale is reversed.}
  
}
\details{
The number of points \eqn{m} doesn't change and the frequency interval is from \eqn{-SW/2} to \eqn{SW/2 - SW/m} (the \eqn{-SW/m} is due to the fact that we only have \eqn{m} points, not \eqn{m+1}, the Fourier Transform is periodic with period \eqn{SW} so it is the same at \eqn{-SW/2} and \eqn{SW/2} anyway).

\eqn{SW} is taken in \code{Fid_info} in the column with name \code{SW_h} of the first row.
It is assumed to be the same for every FID since their column names are shared.
}
\value{
\item{RawSpect_data}{The matrix of spectra in Hertz, that is, their coordinates in Hertz are written as column names of the matrix.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat
}

\examples{
FT.spec <- FourierTransform(Data_HS_sp$FidData_HS_3, FidInfo_HS, SW_h=12019.23)
}
\keyword{ manip }

