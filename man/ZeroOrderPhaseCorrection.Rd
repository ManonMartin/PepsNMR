\name{ZeroOrderPhaseCorrection}
\alias{ZeroOrderPhaseCorrection}
\title{
Zero Order Phase Correction 
}
\description{
Corrects the spectra in order to have the more absorptive mode of their real part.


}
\usage{
ZeroOrderPhaseCorrection(RawSpect_data, plot_rms=NULL, returnAngle = FALSE, createWindow = TRUE, Angle = NULL, p.zo=0.8, plot_spectra = FALSE)
}
\arguments{
  \item{RawSpect_data}{Matrix containing the spectra in Hertz, one row per spectrum.}
  \item{plot_rms}{Contains a vector of row names for which a debug plot should be made showing the value of the function we try to minimize as a function of the phase.}
  \item{returnAngle}{If \code{TRUE}, will return the rotation angle used for phase correction.}
  \item{createWindow}{If \code{TRUE}, will open a new window to draw the rms plots, if \code{FALSE}, plots are drawn in the current device.}
 \item{Angle}{If not \code{NULL}, is a numeric vector with angles specified in radian to maually rotate the spectra. One angle is given per spectrum.}
  \item{p.zo}{Positive and negative intensities \emph{resp.} above and bellow the \emph{p}th quantile are averaged to detect an incorrect rotation due to a high water signal. p is the probability for the sample quantile used as threshold.}
 \item{plot_spectra}{If \code{TRUE}, will draw real and imaginary parts of the rotated spectra.}
 \item{rotation}{If \code{TRUE}, an automatic correction on rotation is applied to spectra.}
}
\details{
After the First Order Phase Correction, every part of the signal, \emph{i.e.} every frequency, theorically starts synchronized at 0.
All the values of the Fourier Transform should therefore theorically have the same argument (or angle).
The purpose is to find this argument in a robust way and to rotate the signal of this angle in the complex plane so that the spectrum becomes real positive, \emph{i.e.} has an absoptive mode.

With an ideal FID and after the rotation, we will have a zero imaginary part and a positive real part.
However, in practical FIDs, all the values of the Fourier Transform do not have the same argument.

In the ideal case, we would discard the imaginary part and keep the positive real part
which can be interpreted as the density of some molecules since it is positive.

In the practical case, since we discard the imaginary part, we will focus our optimization on the positiveness of the real part.
You will observe in the output that the imaginary part is not so small but the real part should almost have only positive values.
\code{\link{BaselineCorrection}} and \code{\link{NegativeValuesZeroing}} will take care of the last negative values of the real part. See the reference for more details.

If a very high water signal is present in the spectrum, the positiveness criterion can erroneously rotate the spectrum of 180 degrees. An automatic procedure to detect this issue will average positive and negative intensities \emph{resp.} above and bellow the \emph{p}th quantile. If the mean value is below 0, an automatic 180 degrees rotation is applied to the spectrum.

}
\value{
\item{RawSpect_data}{The matrix of rotated spectra.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat
}

\examples{
Zopc.spec <- ZeroOrderPhaseCorrection(Data_HS_sp$RawSpectData_HS_4)
}

\keyword{ manip }
