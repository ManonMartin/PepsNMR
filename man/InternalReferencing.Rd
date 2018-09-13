\name{InternalReferencing}
\alias{InternalReferencing}
\title{
Chemical shift referencing.
}
\description{
Chemical shifts are referenced against a Reference Compound (RC, \emph{e.g.} TMSP).
}
\usage{
InternalReferencing(Spectrum_data, Fid_info, method = c("max", 
               "thres"), range = c("nearvalue", "all", "window"), 
               ppm.value = 0, direction = "left", 
               shiftHandling = c("zerofilling",  "cut", "NAfilling", 
               "circular"), c = 2, pc = 0.02, fromto.RC = NULL, 
               ppm.ir = TRUE, rowindex_graph = NULL)
}
\arguments{
\item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
 
\item{Fid_info}{Matrix containing the information for each spectrum, one row per spectrum, as returned by \code{\link{ReadFids}}.}

\item{method}{Method used to find the RC peak in the spectra, See the details section.}

\item{range}{How the search zone is defined. Either accross the whole ppm axis (\code{"all"}), near the 0 ppm location (\code{nearvalue}) with parameter \code{pc}, or in a manually specified area of the ppm axis (\code{"window"}) with the non-null parameter \code{fromto.RC}.}

\item{ppm.value}{By default, the ppm value of the reference compound is set to 0, but any arbitrary value in the ppm interval of spectra can be used instead.}


\item{direction}{If \code{method = "thres"}, the direction towards which to search for the RC peak.}

\item{shiftHandling}{See the details section.}

\item{c}{If \code{method = "thres"}, parameter used to fix the threshold for the RC peak.}

\item{pc}{If \code{range = "nearvalue"}, percentage of the ppm axis around the \code{ppm.value} ppm value to look for the RC peak (\emph{e.g.} for \code{pc = 0.02}, intensities whose index values are 0.01\% below and above 0 ppm are investiguated).}

\item{fromto.RC}{If \code{range = "window"}, a list containing numerical vectors indicating the extremities of the intervals within which to search for the RC peak. These extremities are either frequencies in ppm (decreasing values) OR in column indices (increasing values) depending on the \code{ppm.ir} value (e.g. \code{fromto.RC = list(c(0,10000))} if \code{ppm.ir == FALSE} or \code{fromto.RC = list(c(1,-1))} if \code{ppm.ir == TRUE}).}
 
\item{ppm.ir}{If \code{TRUE}, the values in \code{fromto.RC} represent frequencies in ppm (column names of spectra), if \code{FALSE} these values are column indices.}

\item{rowindex_graph}{If not \code{NULL}, a numeric vector with the row numbers of spectra that need to be plotted for inspection.}

}
\details{

Once the search zone is defined with \code{range}, the RC is found depending on the \code{method}. If \code{method = "thres"}, RC is the first peak in the spectrum higher than a predefined threshold which is computed as: \code{c*(cumulated_mean/cumulated_sd)}. If \code{method = "max"}, the maximum intensity in the search zone is defined as the RC.


Since the spectra can be shifted differently, we need to handle misalignment of the left and right of the spectrum.

This can be illustrated here:

\preformatted{
| : TMSP peak

before
 1  2  3  |  5  6  7  8  9
 1  2  3  4  5  |  7  8  9
 1  2  3  4  |  6  7  8  9

shifted
-5 -4 -3 -2 -1  0  1  2  3  4  5 : ppm scale
       1  2  3  |  5  6  7  8  9
 1  2  3  4  5  |  7  8  9
    1  2  3  4  |  6  7  8  9
}


The different \bold{shift handlings} (\code{shiftHandling}) are the following:
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
if \code{rowindex_graph} is \code{NULL}:
\item{Spectrum_data}{The matrix of the spectral value in the ppm scale.}


if \code{rowindex_graph} is not \code{NULL}:
\item{Spectrum_data}{The matrix of the spectral value in the ppm scale.}
\item{plots}{The spectra that need to be plotted for inspection.}


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
PpmConv.spec <- InternalReferencing(Data_HS_sp$Spectrum_data_HS_5, 
                             FidInfo_HS, shiftHandling = "zerofilling")
}
\keyword{ manip }
