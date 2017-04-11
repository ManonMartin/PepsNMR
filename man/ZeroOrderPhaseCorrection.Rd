\name{ZeroOrderPhaseCorrection}
\alias{ZeroOrderPhaseCorrection}
\title{
Zero Order Phase Correction 
}
\description{
The function corrects the spectra in order to have the more absorptive mode of their real part.


}
\usage{
ZeroOrderPhaseCorrection(Spectrum_data, method = c("rms", "manual", "max"), 
                        plot_rms = NULL, returnAngle = FALSE, 
                        createWindow = TRUE, 
                        angle = NULL, plot_spectra = FALSE,  
                        quant = 0.95, freq = TRUE, fromto.0OPC = NULL)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in hertz, one row per spectrum.}
  \item{method}{Method used to select the angles to rotate the spectra. See details.}
  \item{plot_rms}{Contains a vector of row names for which a debug plot should be made showing the value of the function we try to minimize as a function of the phase.}
  \item{returnAngle}{If \code{TRUE}, will return the rotation angle used for phase correction.}
  \item{createWindow}{If \code{TRUE}, will open a new window to draw the rms or spectra plots, if \code{FALSE}, plots are drawn in the current device.}
 \item{angle}{If not \code{NULL}, a numeric vector with angles specified in radian to maually rotate the spectra. One angle is given per spectrum.}
 \item{plot_spectra}{If \code{TRUE}, will draw real and imaginary parts of the rotated spectra.}
 \item{quant}{Quantile probability for the positiveness criterion. See details.}
 \item{freq}{If \code{TRUE}, the values in \code{fromto.0OPC} represent frequencies (column names of spectra), if \code{FALSE} these values are indices.}
 \item{fromto.0OPC}{If not \code{NULL}, a list containing the extremities of the intervals used to compute the positiveness criterion, either in ppm OR column index.}
}
\details{

When \code{method} is \code{"rms"}, a positiveness criterion is applied on the spectrum. \code{"manual"} is used when a vector of angles are specified in \code{Angle} and \code{"max"} will optimize a maximum spectral intensity. Beware that if \code{fromto.0OPC} is not \code{NULL}, the optimization will only consider the selected spectral window(s).


After the first order phase correction, every part of the signal, \emph{i.e.} every frequency, theorically starts synchronized at 0.
All the values of the Fourier transform should therefore theorically have the same argument (or angle).
The purpose is to find this argument in a robust way and to rotate the signal of this angle in the complex plane so that the spectrum becomes real positive, \emph{i.e.} has an absoptive mode.

With an ideal FID and after the rotation, we will have a zero imaginary part and a positive real part.
However, in practical FIDs, all the values of the Fourier Transform do not have the same argument.

In the ideal case, we would discard the imaginary part and keep the positive real part
which can be interpreted as the density of some molecules since it is positive.

In the practical case, since we discard the imaginary part, we will focus our optimization on the positiveness of the real part.
You will observe in the output that the imaginary part is not so small but the real part should almost have only positive values.
\code{\link{BaselineCorrection}} and \code{\link{NegativeValuesZeroing}} will take care of the last negative values of the real part. See the reference for more details.

%If a very high water signal is present in the spectrum, the positiveness criterion can erroneously rotate the spectrum of 180 degrees. An automatic procedure to detect this issue will average positive and negative intensities \emph{resp.} above and bellow the \emph{p}th quantile. If the mean value is below 0, an automatic 180 degrees rotation is applied to the spectrum.



}
\value{
\item{Spectrum_data}{The matrix of rotated spectra.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat, Manon Martin
}

\examples{
Zopc.res <- ZeroOrderPhaseCorrection(Data_HS_sp$Spectrum_data_HS_4, 
            freq = FALSE, fromto.0OPC = list(c(5000,15000)))
}

\keyword{ manip }
