\name{InternalReferencing}
\alias{InternalReferencing}
\title{
Chemical shift referencing.
}
\description{
Chemical shift are referenced against a Reference Compound (RC, \emph{e.g.} TMSP).
}
\usage{
InternalReferencing(RawSpect_data, RawSpect_info, method = c("max", 
               "thres"), range = c("near0", "all", "window"), 
               ppm.ref = 0, shiftHandling = c("zerofilling", 
               "cut", "NAfilling", "circular"), c = 2, pc = 0.02, 
               fromto.TMSP = NULL, ppmxaxis = TRUE, 
               rowindex_graph = NULL)
}
\arguments{
\item{RawSpect_data}{Matrix containing the spectra in Hertz, one row per spectrum, as returned by \code{\link{ReadFids}}.}
 
\item{RawSpect_info}{Matrix containing the information for each spectrum, one row per spectrum, as returned by \code{\link{ReadFids}}.}

\item{method}{Method used to find the RC peak in the spectra, \emph{cf.} Details section.}

\item{range}{How the search zone is defined. Either accross the whole ppm axis (\code{"all"}), near the 0 ppm location (\code{near0}) with parameter \code{pc}, or in a manually specified area of the ppm axis (\code{"window"}) with the non-null parameter \code{fromto.TMSP}.}

\item{ppm.ref}{By default, the ppm value of the reference compound is set to 0, but any arbitrary value in the ppm interval can be used instead.}

\item{shiftHandling}{See the details section.}

\item{c}{If \code{method = "thres"}, parameter used to fix the threshold for the RC peak.}

\item{pc}{If \code{range == "near0"}, percentage of the ppm axis around 0 ppm to look for the RC peak (\emph{e.g.} for \code{pc = 0.02}, intensities whose index values are 0.01\% below and above 0 ppm are investiguated).}

\item{fromto.TMSP}{If \code{range == "window"}, a list containing numerical vectors indicating the extremities of the intervals within which to search for the TMSP peak (e.g. \code{fromto.TMSP = list(c(0,10000))}). These extremities are either frequencies in ppm value OR column indices depending on the \code{ppmxaxis} value.}
 
\item{ppmxaxis}{If \code{TRUE}, the values in \code{fromto.TMSP} represent frequencies (column names of spectra), if \code{FALSE} these values are column indices.}

\item{rowindex_graph}{If not \code{NULL}, a numeric vector with the row numbers of spectra that need to be plotted for inspection.}

}
\details{

Once the search zone is defined with \code{range}, the Referencing Compound (RC) is found depending on the \code{method}. If \code{method = "thres"}, RC is the first peak in the spectrum higher than a predefined threshold which is computed as: \code{c*(cumulated_mean/cumulated_sd)}. If \code{method = "max"}, the maximum intensity in the search zone is defined as the RC.


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
PpmConv.spec <- InternalReferencing(Data_HS_sp$Spectrum_data_HS_7, FidInfo_HS, 
                              shiftHandling = "zerofilling")
}
\keyword{ manip }
