\name{PPMConversion}
\alias{PPMConversion}
\title{
PPM scale conversion.
}
\description{
Converts the frequency scale from Hertz to part per million (PPM).
}
\usage{
PPMConversion(RawSpect_data, RawSpect_info,
              shiftHandling = c( "cut", "zerofilling",
              "NAfilling", "circular"),c = 2)
}
\arguments{
\item{RawSpect_data}{Matrix containing the spectra in Hertz, one row per spectrum.}
\item{RawSpect_info}{Matrix containing the information for each spectrum, one row per spectrum, as returned by \code{\link{ReadFids}}.}
\item{shiftHandling}{See the details section.}
\item{c}{Parameter used to fix the threshold for the TMSP peak.}
}
\details{

The frequency scale is dependent on the kind of spectrometer used, more precisely on its
external magnetic field. We therefore translate it to a PPM (Part Per Million) scale which is independent of this factor.

TMSP peak is defined as the first peak in the spectrum higher than a predefined threshold. This threshold is computed as: \code{c*(cumulated_mean/cumulated_sd)}.

Since the spectra can be shifted differently, we need to handle misalignment of the left and right of the spectrum.

This can be illustrated here:

\preformatted{
before
| : TMSP peak
 1  2  3  |  5  6  7  8  9
 1  2  3  4  5  |  7  8  9
 1  2  3  4  |  6  7  8  9

shifted
-5 -4 -3 -2 -1  0  1  2  3  4  5 : ppm scale
       1  2  3  |  5  6  7  8  9
 1  2  3  4  5  |  7  8  9
    1  2  3  4  |  6  7  8  9
}

The different shift handlings are the following:
\describe{
\item{\code{NAfilling}}{The extremities at which a spectrum is not defined are replaced by \code{NA}. It is detected by \code{\link{WindowSelection}} which produces a warning if there are NAs in the selected window.
\preformatted{
-5 -4 -3 -2 -1  0  1  2  3  4  5  ppm scale
NA NA  1  2  3  |  5  6  7  8  9
 1  2  3  4  5  |  7  8  9 NA NA
NA  1  2  3  4  |  6  7  8  9 NA
}}
\item{\code{zerofilling}}{The extremities at which a spectrum is not defined are replaced by \code{0}. It makes sense since in practice the spectrum is close to zero at the extremities.
\preformatted{
-5 -4 -3 -2 -1  0  1  2  3  4  5  ppm scale
 0  0  1  2  3  |  5  6  7  8  9
 1  2  3  4  5  |  7  8  9  0  0
 0  1  2  3  4  |  6  7  8  9  0
}}
\item{\code{circular}}{The spectra are shifted circularly which means that the end of a spectrum is reproduced at the beginning.
It makes sense since the spectrum is periodic since it is the result of FFT.
\preformatted{
-5 -4 -3 -2 -1  0  1  2  3        ppm scale
 8  9  1  2  3  |  5  6  7
 1  2  3  4  5  |  7  8  9
 9  1  2  3  4  |  6  7  8
}}
\item{\code{cut}}{The ppm values for which some spectra are not defined are removed.
\preformatted{
      -3 -2 -1  0  1  2  3        ppm scale
       1  2  3  |  5  6  7
       3  4  5  |  7  8  9
       2  3  4  |  6  7  8
}}
The difference between these shift handlings should not be critical in practice since the extremities of the spectra are not used most of the time and are removed in \code{\link{WindowSelection}.
}}}
\value{
\item{Spectrum_data}{The matrix of the spectral value in the ppm scale.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat
}
\examples{
PpmConv.spec <- PPMConversion(Data_HS_sp$RawSpectData_HS_7, FidInfo_HS, shiftHandling ="zerofilling")
}
\keyword{ manip }
