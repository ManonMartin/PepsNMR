\name{FourierTransform}
\alias{FourierTransform}
\title{
Applies the fourier transforationm to the FIDs.
}
\description{
The function takes the FIDs in the time domain and translate it into the frequency domain. It also converts the frequency scale from hertz to part per million (ppm).

}
\usage{
FourierTransform(Fid_data, Fid_info = NULL, SW_h = NULL, SW = NULL, 
                 O1 = NULL, reverse.axis = TRUE, verbose = FALSE)
}
\arguments{
  \item{Fid_data}{Matrix containing the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{Fid_info}{Matrix containing the info about the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{SW_h}{Sweep Width in hertz. If given, the value in \code{Fid_info} is ignored.}
  \item{SW}{Sweep width in ppm. If given, the value in \code{Fid_info} is ignored.}
  \item{O1}{Spectrometer frequency offset. If given, the value in \code{Fid_info} is ignored.}
  \item{reverse.axis}{If \code{TRUE}, the frequency scale is reversed.}
  \item{verbose}{If\code{"TRUE"}, will print processing information.}
  
}
\details{
The number of points \eqn{m} doesn't change and the frequency interval is from \eqn{-SW/2} to \eqn{SW/2 - SW/m} (the \eqn{-SW/m} is due to the fact that we only have \eqn{m} points, not \eqn{m+1} and the fourier transform is periodic with period \eqn{SW} so it is the same at \eqn{-SW/2} and \eqn{SW/2} anyway).

\eqn{SW}, \eqn{SW_h} and \eqn{O1} are usually taken from the \code{Fid_info} matrix.
\eqn{SW} and \eqn{SW_h} are assumed to be the same for every FID since their column names are shared.


The frequency scale is dependent on the kind of spectrometer used, more precisely on its external magnetic field. We therefore translate it to a ppm (part per million) scale which is independent of this external magnetic field thanks to the recovered transmitter frequency offset value (\eqn{O1}).

}
\value{
\item{RawSpect_data}{The matrix of spectra in ppm.}
}
\references{
Martin, M., Legat, B., Leenders, J., Vanwinsberghe, J., Rousseau, R., Boulanger, B., & Govaerts, B. (2018). PepsNMR for 1H NMR metabolomic data pre-processing. \emph{Analytica chimica acta}, 1019, 1-13.

  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat & Manon Martin
}

\examples{
require(PepsNMRData)
FT.spec <- FourierTransform(Data_HS_sp$FidData_HS_3,FidInfo_HS_sp, SW_h = 12019.23)
}
\keyword{ manip }

