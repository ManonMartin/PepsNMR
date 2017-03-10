\name{Warping}
\alias{Warping}
\title{Warping of the spectra}
\description{
Warps the frequency \emph{x}-axis to minimize the pairwise distance between a sample spectrum and a reference spectrum.
}
\usage{
Warping(Spectrum_data, normalization.type=c("median","mean",
        "firstquartile","peak","none"), fromto.normW = c(3.05, 4.05), 
        reference.choosing=c("fixed", "before", 
        "after"), reference=1, optim.crit=c("RMS", "WCC"), 
        ptw.wp=F, K=3, L=40, lambda.smooth=0, deg=3, 
        lambda.bspline=0.01, kappa=0.0001, max_it_Bspline=10, 
        returnReference = FALSE, returnWarpFunc = FALSE)
}
\arguments{
\item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
\item{normalization.type}{Type of normalization applied to the spectrum prior to warping.
See \code{\link{Normalization}} for details about the different types.
\code{none} means that no normalization is applied.
It is advised to use \code{median} instead of the default \code{mean} normalization.}
\item{fromto.normW}{Used by \code{\link{Normalization}} when \code{normalization.type} is \code{peak}.}
\item{reference.choosing}{Specifies how the reference will be chosen:
\describe{
\item{\code{"fixed"}}{The reference is specified by the rowname given in \code{reference}.}
\item{\code{"before"}}{The reference is taken as the spectrum with the minimum sum of square difference with the other spectra.}
\item{\code{"after"}}{Each spectrum is taken as the reference and the sum of square difference with the other spectra is calculated after the warping. See details below.}
}}
\item{reference}{The row number or name of the reference spectrum when \code{reference.choosing} is \code{fixed}.}
\item{optim.crit}{If \code{ptw.wp} is set to \code{TRUE}, \code{WCC} can also be considered as a criterion for optimization, see \code{ptw::ptw} for details.}
\item{ptw.wp}{If set to \code{TRUE}, it uses \code{ptw::ptw} to do the warping.
In this case the arguments \code{K}, \code{L}, \code{deg}, \code{lambda.bspline}, \code{kappa} and \code{max_it_Bspline} are ignored.
The warping does not use Bspline and \code{K} is fixed to 2, \emph{i.e.} a second order polynomial.}
\item{K}{It is the degree of the polynomial used for the warping (see details).}
\item{L}{This is the number of B-Splines that are used for the warping.
It should be either 0 or greater than \code{deg}.}
\item{lambda.smooth}{Nonnegative coefficient for the smoothing \code{lambda.smooth = 0} means no smoothing. See \code{ptw::difsm} for more details.}
\item{deg}{Degree of the B-Splines.}
\item{lambda.bspline}{Nonnegative roughness penalty coefficient for the B-Spline warping. See the reference for more details.}
\item{kappa}{Nonnegative ridge penalty coefficient for the B-Spline warping.
See the reference for more details.}
\item{max_it_Bspline}{Maximum number of iterations for the B-Spline warping.}
\item{returnReference}{If \code{TRUE}, will return the name of the reference spectrum.}
\item{returnWarpFunc}{If \code{TRUE}, will return the warping functions.}

}
\details{

When \code{reference.choosing} is \code{"after"}, the reference with the minimum sum is taken as the reference and the warped spectra according to this reference (that have already been calculated at this stage) are returned. This is \eqn{n}{n} times slower than the 2 others where \eqn{n}{n} is the number of spectra.


\bold{Principle}: 

We try to find a warping function between a reference spectrum and a sample.
This function is a sum of polynomial of degree \code{K} and \code{L} B-Splines of degree \code{deg}.
The unknowns are the polynomial and B-Splines coefficients.

No warping is equivalent to warping with a, the polynomial identity and all the coefficients of the B-Splines with value 0.
See the reference for details.

First, the polynomial is estimated on the reference and the sample both smoothed with parameter \code{lambda.smooth}.
The B-Splines are estimated on the non-smoothed reference and sample using the polynomial just found.
}
\value{
If \code{returnReference = TRUE}, the function will return the name of the reference spectrum and if \code{returnWarpingfunc = TRUE}, it will also return the warping function.

  \item{Spectrum_data}{The warped spectra.}
  \item{Reference}{The name of the reference spectrum.}
  \item{Warpingfunc}{The warping function.}

}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).

  
}
\author{
Benoît Legat
}

\examples{
Warp.spec <- Warping(Data_HS_sp$Spectrum_data_HS_8, reference.choosing="fixed", 
                             reference = row.names(Data_HS_sp$Spectrum_data_HS_8)[1],
                             returnReference = FALSE)
#or
Warp.res <- Warping(Data_HS_sp$Spectrum_data_HS_8, reference.choosing="fixed", 
                             reference = row.names(Data_HS_sp$Spectrum_data_HS_8)[1],
                             returnReference = TRUE)
Warp.spec <- Warp.res[["Spectrum_data"]]
Warp.res[["Reference"]]


}
\keyword{ manip }
