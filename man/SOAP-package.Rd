\name{SOAP-package}
\alias{SOAP-package}
\alias{SOAP}
\docType{package}
\title{Metabolomic data pre-treatment strategy for 1H-NMR spectroscopic data
}
\description{
This package provides R functions for common pre-treatment steps that are applied
on 1H-NMR data.
It also provides the function \code{\link{ReadFids}} to read the FID directly from the Bruker output.
Those pre-treatments are cited below in the advised order of their use:
\describe{
  \item{\code{\link{FirstOrderPhaseCorrection}}}{Correct for the first order phase correction.}
  \item{\code{\link{SolventSuppression}}}{Remove solvent signal from the FID.}
  \item{\code{\link{Apodization}}}{Increase the Signal-to-Noise ratio of the FID.}
  \item{\code{\link{FourierTransform}}}{Transform the FID into a spectrum.}
  \item{\code{\link{ZeroOrderPhaseCorrection}}}{Correct for the zero order phase correction.}
  \item{\code{\link{BaselineCorrection}}}{Remove the spectral baseline.}
  \item{\code{\link{NegativeValuesZeroing}}}{Set negatives values to 0.}
  \item{\code{\link{PPMConversion}}}{Convert the frequency scale (Hertz -> ppm).}
  \item{\code{\link{Warping}}}{Warp the samples according to a reference spectrum.}
  \item{\code{\link{WindowSelection}}}{Select the informative part of the spectrum.}
  \item{\code{\link{Bucketing}}}{Proceed to data reduction.}
  \item{\code{\link{RegionRemoval}}}{Set a desired region to 0.}
  \item{\code{\link{ZoneAggregation}}}{Aggregate a region to a single peak.}
  \item{\code{\link{Normalization}}}{Normalize the spectra.}
}
}
\details{
\tabular{ll}{
Package: \tab SOAP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-24\cr
License: \tab GPLv2\cr
}
The FIDs are read using \code{\link{ReadFids}} which also gives a matrix with meta information about each FID.
The other functions apply different pre-treatment steps on these FIDs, some needing the info matrix.
During this pre-treatment, the FID is transformed through its Fourier Transform and the frequency scale is expressed in PPM.
For more detailed and illustrated explanation about these pre-treatment steps, see the documentation of each function and/or the chapter 1 of the reference below.
}
\author{
Benoît Legat

Maintainer: Benoît Legat <benoit.legat@gmail.com>
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\keyword{ package }
\seealso{
}
\examples{
path <-  system.file("extdata", package = "SOAP")
dir(path)

fidList <- ReadFids(file.path(path, "HumanSerum"))
Fid_data <- fidList[["Fid_data"]]
Fid_info <- fidList[["Fid_info"]]
Fid_data <- FirstOrderPhaseCorrection(Fid_data, Fid_info)
Fid_data <- SolventSuppression(Fid_data)
Fid_data <- Apodization(Fid_data, Fid_info)
RawSpect_data <- FourierTransform(Fid_data, Fid_info)
RawSpect_data <- ZeroOrderPhaseCorrection(RawSpect_data)
RawSpect_data <- BaselineCorrection(RawSpect_data)
RawSpect_data <- NegativeValuesZeroing(RawSpect_data)
Spectrum_data <- PPMConversion(RawSpect_data, Fid_info)
Spectrum_data <- Warping(Spectrum_data)
Spectrum_data <- WindowSelection(Spectrum_data)
Spectrum_data <- Bucketing(Spectrum_data)
Spectrum_data <- RegionRemoval(Spectrum_data)
Spectrum_data <- ZoneAggregation(Spectrum_data)
Spectrum_data <- Normalization(Spectrum_data)
}
