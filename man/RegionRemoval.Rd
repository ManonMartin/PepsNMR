\name{RegionRemoval}
\alias{RegionRemoval}
\title{
Removal of non-informative regions
}
\description{
Removes the non-informative regions by setting the values of the spectrum in these intervals to zero.

}
\usage{
RegionRemoval(Spectrum_data, typeofspectra = c("Serum", "Urine", NULL), 
                          type = c( "zero", "NA"), fromto=NULL)
}
\arguments{
  \item{Spectrum_data}{The matrix containing the spectra in ppm, one row per spectrum.}
  \item{typeofspectra}{Type of spectra.}
  \item{type}{The type of region removal method.}
  \item{fromto}{List containing the extremities of the intervals to be removed.}

}
\details{
The presence of non-informative regions can strongly bias the subsequent statistical analysis.

The interval from \code{from} ppm to \code{to} ppm inclusive is set to zero or completed with NAs for every spectrum.
The ppm scale can be increasing or decreasing and \code{from} and \code{to} can be given in any order (i.e. \code{from < to} or \code{from > to}).

The type of spectra can be NULL to manually specify the area to be removed otherwise it is specified as "Serum" or "Urine" and the default removed area are: 

Serum spectra: Water (4.5 - 5.04 ppm) and lactate (1.36 - 2.8 ppm) peaks
Urine spectra: Water (4.5 - 5.04 ppm) and uree (4.5 - 6 ppm).


}
\value{
\item{Spectrum_data}{The matrix of spectra with the region removed.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat
}
\examples{

# Remove the lactate and water regions for serum spectra
fromto=list(Water =c(4.5, 5), Lactate=c(1.28, 1.36))
Rr.spec <- RegionRemoval(Data_HS$RawSpectData_HS_11,fromto = fromto) 



}
\keyword{ manip }
