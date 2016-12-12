\name{NegativeValuesZeroing}
\alias{NegativeValuesZeroing}
\title{
Zeroing of negative values.
}
\description{
The function sets negative intensities to zero.
}
\usage{
NegativeValuesZeroing(Spectrum_data)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in Hertz, one row per spectrum..}
}
\value{
\item{Spectrum_data}{The matrix of spectrums with the negative values set to zero.}
}
\details{
As explained in \code{\link{BaselineCorrection}},
negative values does not make sense and can have bad impacts on our statistical analyses.
\code{\link{BaselineCorrection}} do its best to remove negative values but there might be some remaining.

This filter simply sets them to zero.
After the \code{\link{BaselineCorrection}} they should be close to zero anyway because of the high penalty given to negative values of the signal after the correction.
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat
}
\examples{
Nvz.spec <- NegativeValuesZeroing(Data_HS_sp$RawSpectData_HS_6)
}
\keyword{ manip }

