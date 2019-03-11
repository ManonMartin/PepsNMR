\name{ZeroOrderPhaseCorrection}
\alias{ZeroOrderPhaseCorrection}
\title{
Zero Order Phase Correction 
}
\description{
The function corrects the spectra in order to have their real part in an absorptive mode.


}
\usage{
ZeroOrderPhaseCorrection(Spectrum_data, type.zopc = c("rms", "manual", "max"), 
                        plot_rms = NULL, returnAngle = FALSE, 
                        createWindow = TRUE, angle = NULL, plot_spectra = FALSE,  
                        ppm.zopc = TRUE, exclude.zopc = list(c(5.1,4.5)), verbose = FALSE)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{type.zopc}{Method used to select the angles to rotate the spectra. See details.}
  \item{plot_rms}{Contains a vector of row names for which a debug plot should be made showing the value of the function we try to minimize as a function of the phase.}
  \item{returnAngle}{If \code{TRUE}, will return the rotation angle used for phase correction.}
  \item{createWindow}{If \code{TRUE}, will open a new window to draw the rms or spectra plots, if \code{FALSE}, plots are drawn in the current device.}
 \item{angle}{If not \code{NULL}, a numeric vector with angles specified in radian to manually rotate the spectra, one angle per spectrum. By convention, the spectra are rotated with the correction angle \code{- angle}.}
 \item{plot_spectra}{If \code{TRUE}, will draw real and imaginary parts of the rotated spectra.}
 \item{ppm.zopc}{If \code{TRUE}, the values in \code{exclude.zopc} represent frequencies in ppm value (column names of spectra), if \code{FALSE} these values are column indices.}
 \item{exclude.zopc}{If not \code{NULL}, a list containing the extremities of the intervals excluded for the computation of the positiveness criterion, either expressed in ppm (decreasing values) OR in column indices (increasing values), e.g. \code{exclude.zopc = list(c(0,10000))} if \code{ppm.zopc == FALSE} or \code{exclude.zopc = list(c(1,-1))} if \code{ppm.zopc == TRUE}.}
\item{verbose}{If \code{"TRUE"}, will print processing information.}
}
\details{

We focus our optimization on the positiveness of the real part which should be in an absoptive mode.

When \code{type.zopc} is \code{"rms"}, a positiveness criterion is measured for each spectrum. \code{"manual"} is used when a vector of angles are specified in \code{angle} and \code{"max"} will optimize the maximum spectral intensity in the non-excluded window(s). Beware that if \code{exclude.zopc} is not \code{NULL}, the optimization will only consider the non-excluded spectral window(s). 

By default the water region (5.1 - 4.5) is ignored.

\code{\link{BaselineCorrection}} and \code{\link{NegativeValuesZeroing}} will take care of the last negative values of the real part of the spectra. See the reference for more details.
}
\value{
\item{Spectrum_data}{The matrix of rotated spectra.}
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
Zopc.res <- ZeroOrderPhaseCorrection(Data_HS_sp$Spectrum_data_HS_4, 
            ppm.zopc = FALSE, exclude.zopc = list(c(5000,15000)))
}

\keyword{ manip }
