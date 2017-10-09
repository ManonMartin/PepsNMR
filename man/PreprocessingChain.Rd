\name{PreprocessingChain}
\alias{PreprocessingChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Proprocessing workflow for 1H-NMR data
}
\description{
The function is a wrapper for all the preprocessing steps available in PepsNMR.
}
\usage{
PreprocessingChain(Fid_data = NULL, Fid_info = NULL, data.path = NULL, readFids = TRUE, 
    groupDelayCorr = TRUE, solventSuppression = TRUE, apodization = TRUE, 
    fourierTransform = TRUE, zeroOrderPhaseCorr = TRUE, internalReferencing = TRUE,
    baselineCorrection = TRUE, negativeValues0 = TRUE, warping = TRUE, 
    windowSelection = TRUE, bucketing = TRUE, regionRemoval = TRUE, 
    zoneAggregation = TRUE,normalization = TRUE, ...,  export = FALSE, 
    format = c("Rdata", "csv", "txt"), out.path = ".",  filename = "filename", 
    writeArg = c("none", "return", "txt"))
}

\arguments{
  \item{Fid_data}{
If non \code{NULL}, matrix containing the complex FIDs, one row per FID.
}
  \item{Fid_info}{
If non \code{NULL}, matrix containing the information for each spectrum, one row per spectrum, as returned by \code{\link{ReadFids}}.
}
  \item{data.path}{
A character string specifying the directory where the FIDs are searched.
}
  \item{readFids}{
If \code{TRUE}, applies the \code{\link{ReadFids}} function to the data.
}
  \item{groupDelayCorr}{
If \code{TRUE}, applies the \code{\link{GroupDelayCorrection}} function to the data.
}
  \item{solventSuppression}{
If \code{TRUE}, applies the \code{\link{SolventSuppression}} function to the data.
}
  \item{apodization}{
If \code{TRUE}, applies the \code{\link{Apodization}} function to the data.
}
  \item{fourierTransform}{
If \code{TRUE}, applies the \code{\link{FourierTransform}} function to the data.
}
  \item{zeroOrderPhaseCorr}{
If \code{TRUE}, applies the \code{\link{ZeroOrderPhaseCorrection}} function to the data.
}
  \item{internalReferencing}{
If \code{TRUE}, applies the \code{\link{InternalReferencing}} function to the data.
}
  \item{baselineCorrection}{
If \code{TRUE}, applies the \code{\link{BaselineCorrection}} function to the data.
}
  \item{negativeValues0}{
If \code{TRUE}, applies the \code{\link{NegativeValuesZeroing}} function to the data.
}
  \item{warping}{
If \code{TRUE}, applies the \code{\link{Warping}} function to the data.
}
  \item{windowSelection}{
If \code{TRUE}, applies the \code{\link{WindowSelection}} function to the data.
}
  \item{bucketing}{
If \code{TRUE}, applies the \code{\link{Bucketing}} function to the data.
}
  \item{regionRemoval}{
If \code{TRUE}, applies the \code{\link{RegionRemoval}} function to the data.
}
  \item{zoneAggregation}{
If \code{TRUE}, applies the \code{\link{ZoneAggregation}} function to the data.
}
  \item{normalization}{
If \code{TRUE}, applies the \code{\link{Normalization}} function to the data.
}
  \item{\dots}{
Other optionnal arguments of the above pre-processing functions.
}
  \item{export}{
If \code{TRUE}, will export the spectral intensities and the aquisition parameters matrices.
}
  \item{format}{
Format chosen to export the spectral intensities and the aquisition parameters matrices. 
}
  \item{out.path}{
Path used to export the spectral intensities and the aquisition parameters matrices if \code{export == TRUE} and the function argument if \code{writeArg == "txt"}.
}
  \item{filename}{
Name given to exported files.
}
  \item{writeArg}{
If not \code{"none"}, will export the function arguments, either in the return of the function (\code{"return"}) or as a text file (\code{"txt"}).
}

}

\value{

The function will return a list with the spectral intensities and the aquisition parameters matrices. If \code{writeArg == "return"}, an additionnal list element is returned (\code{arguments}).

  \item{Spectrum_data}{The pre-processed spectra.}
  \item{Fid_info}{The acquisition parameters.}
  \item{arguments}{The function arguments.}
  
}
\references{
 Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Manon Martin
}


\examples{

path <-  system.file("extdata", package = "PepsNMR")
data.path <-  file.path(path, "HumanSerum")
res <-  PreprocessingChain(Fid_data = NULL, Fid_info = NULL, data.path = data.path, ReadFids = TRUE,
                            export = FALSE, writeArg = "return")

}

\keyword{ manip }

