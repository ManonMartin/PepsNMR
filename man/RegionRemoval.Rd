\name{RegionRemoval}
\alias{RegionRemoval}
\title{
Removal of non-informative regions
}
\description{
Removes the non-informative regions by setting the values of the spectra in these intervals to zero.

}
\usage{
RegionRemoval(Spectrum_data, typeofspectra = c("manual", "serum", "urine"), 
                          type.rr = c( "zero", "NA"), 
                          fromto.rr = list(Water = c(4.5, 5.1)), verbose = FALSE)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{typeofspectra}{Type of spectra, if not \code{"manual"}, will automatically remove unwanted regions depending on the nature of spectra.}
  \item{type.rr}{Type of region removal method. If \code{type.rr = "zero"}, intensities are set to 0; if \code{type.rr = "NA"}, intensities are set to NA.}
  \item{fromto.rr}{List containing the extremities of the intervals to be removed.}
  \item{verbose}{If \code{"TRUE"}, will print processing information.}
}
\details{
The presence of non-informative regions can strongly bias the subsequent statistical analysis.

The inclusive ppm interval \code{fromto.rr} is set to zero or completed with NAs for every spectrum.
The ppm scale can be increasing or decreasing (i.e. \code{from < to} or \code{from > to}).

The type of spectra can be NULL to manually specify the area to be removed otherwise it is specified as \code{typeofspectra = "serum"} or \code{typeofspectra = "urine"} and the removed area are for \code{typeofspectra = "serum"}: water (4.5 - 5.1 ppm) and for \code{typeofspectra = "urine"}: water, uree and maleic acid (4.5 - 6.1 ppm).



}
\value{
\item{Spectrum_data}{The matrix of spectra with the removed regions.}
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

# Remove the lactate and water regions for serum spectra
require(PepsNMRData)
fromto <- list(Water =c(4.5, 5.1), Lactate=c(1.32, 1.36))
Rr.spec <- RegionRemoval(Data_HS_sp$Spectrum_data_HS_11,fromto.rr = fromto) 



}
\keyword{ manip }
