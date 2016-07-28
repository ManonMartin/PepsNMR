\name{Apodization}
\alias{Apodization}
\title{Apodization of the FID}
\description{Multiplies the FID with a decreasing function in order to increase the Signal-to-Noise Ratio (SNR).
This is based on the fact that the signal intensity is decreasing over time unlike the noise that keeps a constant amplitude, leaving a noisy tail at the end of the FID signal.}
\usage{
Apodization(Fid_data, Fid_info=NULL, DT=NULL,type.apod = c("exp", "cos2", 
  "blockexp", "blockcos2", "gauss", "hanning", "hamming"),phase=0,
  rectRatio=1/2, gaussLB=1, expLB=1, plotWindow=F, returnFactor=F)
}
\arguments{
  \item{Fid_data}{Matrix containing the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{Fid_info}{Matrix containing the info about the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{DT}{If given, used instead of \code{Fid_info} to give the Dwell Time, the time between 2 points of the FID.}
  \item{type.apod}{Type of apodization, see details.}
  \item{phase}{Pphase at which the apodization window is maximum for \code{cos2}, \code{hanning} and \code{hamming} types. For example, if phase is 0.2, the maximum is at 20\% of the signal.}
  \item{rectRatio}{If there is a rectangular window, ratio between the width of the window and the width of the signal.}
  \item{gaussLB}{Line Broadening for the gaussian window, see details.}
  \item{expLB}{Line Broadening for the exponential window, see details.}
  \item{plotWindow}{Plot the window which is applied to the signal.}
    \item{returnFactor}{Return a list with the final FIDs and the apodization function, see the examples.}
}
\details{
The different types of apodization are:
\describe{
  \item{exp}{The signal is multiplied by a decreasing exponential \eqn{\exp(-t/\textsf{expLB)}}{exp(-t/expLB)}.}
  \item{cos2}{The signal is multiplied by the value of a \eqn{\cos^2}{cosinus squared} from 0 (where its value is 1) until \eqn{\pi/2}{pi/2} (where its value is 0).}
  \item{blockexp}{The first part of the signal (described by \code{rectRatio}) is left unchanged and the second is multiplied by \eqn{\exp(-t/\textsf{expLB)}}{exp(-t/expLB)} starting at value 1.}
  \item{blockcos2}{the first part is left unchanged as with \code{blockexp} and the second part is multiplied by a \eqn{\cos^2}{cosinus squared} where its value starts at 1 at the end of the block and ends at 0 at the end of the signal.}
  \item{gauss}{The signal is multiplied by a gaussian window centered at the beginning of the FID and with \eqn{\sigma = 1/\textsf{gaussLB}}{sigma=1/gaussLB}.}
  \item{hanning}{The signal is multiplied by a hanning window : \eqn{0.5 + 0.5\cos}{0.5 + 0.5cos}.}
  \item{hamming}{The signal is multiplied by a hamming window : \eqn{0.54 + 0.46\cos}{0.54 + 0.46cos}.}
}

A FID of the form \eqn{s_0\exp(i2\pi\nu t)\exp(-t/\tau)}{s0 exp(i*2*pi*nu*t) exp(-t/tau)} has a peak
in its spectrum at the frequency \eqn{\nu}{nu} of width that is inversely proportional to \eqn{\tau}{tau}.
This peak is called a \emph{spectral line} and its width a \emph{spectral width}.

With the exponential multiplication, the decaying exponential becomes
\deqn{\exp(-t(1/\tau + LB))}{exp(-t(1/tau + LB))}
The new decay \eqn{\tau_2}{tau2} which satisfied \eqn{1/\tau_2 = 1/\tau + LB}{1/tau2 = 1/tau + LB} is therefore smaller so the spectral line is \emph{broader}.
That is why we call this parameter the Line Broadening.
}
\value{
If \code{returnFactor = TRUE}, will return a list with the following elements: \code{Fid_data} and \code{Factor}. Otherwise, the function will just return \code{Fid_data}. 
  \item{Fid_data}{The apodized FIDs.}
  \item{Factor}{The apodization function.}
}
\references{
  Inspired from the matNMR library.

  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat
}

\examples{
Apod.fid <- Apodization(Data_HS_sp$FidData_HS_2, FidInfo_HS, plotWindow=FALSE)

#or
Apod.res <- Apodization(Data_HS_sp$FidData_HS_2, FidInfo_HS, plotWindow=FALSE, returnFactor=TRUE)
Apod.fid = Apod.res[["Fid_data"]]
plot(Apod.res[["Factor"]], type="l")
}
\keyword{ manip }
