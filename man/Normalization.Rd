\name{Normalization}
\alias{Normalization}
\title{
Normalizes the spectrum
}
\description{

Spectra normalization to correct for the dilution factor common to all biofuid samples.
}
\usage{
Normalization(Spectrum_data, type.norm = c("mean", "pqn", "median", 
              "firstquartile", "peak"), from.norm = 3.05, to.norm = 4.05, ref)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{type.norm}{Different types of normalization are available. See the details section.}
  \item{from.norm}{Used if \code{type.norm} is \code{peak}. See details.}
  \item{to.norm}{Used if \code{type.norm} is \code{peak}. See details.}
  \item{ref}{Row index for the reference spectrum if \code{type.norm} is \code{pqn}. See details.}
}
\details{

Normalization of spectra before their warping or their statistical analysis is necessary in order to be able to efficiently compare their relative peak intensities. 

It is therefore appropriate to call this filter at the end of the preprocessing.

Types can be:
\describe{
\item{\code{mean}}{Each spectrum is divided by its mean so that its mean becomes 1.}
\item{\code{median}}{Each spectrum is divided by its median so that its median becomes 1.}
\item{\code{firstquartile}}{Each spectrum is divided by its first quartile so that its first quartile becomes 1.}
\item{\code{peak}}{Each spectrum is divided by the value of the peak of the spectrum contained between \code{from.norm} and \code{to.norm} inclusive (\emph{i.e.} the maximum value of the spectrum in that interval).}
\item{\code{pqn}}{Probabilistic Quotient Normalization from Dieterle et al. (2006). If \code{ref} is missing, will use the median spectrum as the reference spectrum, otherwise will use the spectrum specified by the row index \code{ref}.}
}
}
\value{
\item{Spectrum_data}{The matrix of normalized spectra.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
  
  
  Dieterle, F., Ross, A. , Schlotterbeck, G.,and Senn, H (2006). Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. \emph{Analytical Chemistry} 78 (13), 4281-4290 
}
\author{
Benoît Legat
}
\examples{
Norm.spec <- Normalization(Data_HS$RawSpectData_HS_12)
}
\keyword{ manip }

