\name{PreprocessingChain}
\alias{PreprocessingChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pre-treatment chain.
%%  ~~function to do ... ~~
}
\description{Pre-treatment facilities organised in the advised sequence.
}
\usage{
PreprocessingChain(dataname = "Dataset", data.path = getwd(), out.path = getwd(), 
                        nspectr = 1, save = FALSE, saveall = FALSE, ImpG= FALSE, RetArgs = TRUE, 
                        export = c(NULL, "csv", "rdata"), 
                        Fopc = TRUE, Ss = TRUE, A = TRUE, Zopc = TRUE, Bc = TRUE, 
                        Zsnv = TRUE, W = TRUE, B = TRUE, Zs = TRUE, Za=FALSE, N = TRUE,
                        l=1, subdirs = FALSE, 
                        group_delay=NULL, 
                        lambda.ss=1e6, ptw.ss=TRUE,  
                        DT=NULL,type.apod = "exp",phase=0, rectRatio=1/2, gaussLB=1, expLB=1, 
                        SW_h=NULL, 
                        Angle = NULL,  p.zo=0.8, 
                        ptw.bc=TRUE, maxIter = 42,lambda.bc=1e7, p.bc=0.05, eps=1e-8, 
                        shiftHandling="cut", c = 2, 
                        normalization.type="median", from.normW=3.05, 
                        to.normW=4.05,reference.choosing="fixed", 
                        reference=1,optim.crit="RMS", ptw.wp=F, K=3, L=40,
                        lambda.smooth=0, deg=3, lambda.bspline=0.01, kappa=0.0001,max_it_Bspline=10, 
                        from.ws = 0.2, to.ws = 10, reverse.axis = TRUE, 
                        m = 500, width = FALSE, ppmleft=NULL, ppmright = NULL, intmeth = c("r", "t"), tolbuck = 10^-4, 
                        typeofspectra = "serum",type.rr =  "zero", fromto.rr=list(Water =c(4.5, 5.1)), 
                        fromto.za = list(Citrate =c(2.5, 2.7)), 
                        type.norm="mean", from.norm=3.05, to.norm=4.05, ref.norm=1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{dataname}{Includes the dataset name and additional infos.}
  \item{data.path}{Path to the main directory containing the data in Bruker format.}
  \item{out.path}{Path for optional outputs.}
  \item{nspectr}{Row number of the spectrum selected for each step illustrations.}
  \item{save}{If \code{TRUE}, the spectra and acquisition parameters are saved.}
  \item{saveall}{If \code{TRUE}, datasets from each step are saved into a list.}
  \item{ImpG}{If \code{TRUE}, provides graphical illustration for each step.}
  \item{RetArgs}{If \code{TRUE}, the function arguments are saved.}
  \item{export}{If not \code{NULL}, the spectra and acquisition parameters are saved either as .csv or .rdata.}  
  \item{Fopc}{If \code{TRUE}, performs the FirstOrderPhaseCorrection routine}
  \item{Ss}{If \code{TRUE}, performs the SolventSuppression routine.}
  \item{A}{If \code{TRUE}, performs the Apodization routine.}
  \item{Zopc}{If \code{TRUE}, performs the ZeroOrderPhaseCorrection routine.}
  \item{Bc}{If \code{TRUE}, performs the BaselineCorrection routine.}
  \item{Zsnv}{If \code{TRUE}, performs the NegativeValuesZeroing routine.}
  \item{W}{If \code{TRUE}, performs the Warping routine.}
  \item{B}{If \code{TRUE}, performs the Bucketing routine.}
  \item{Zs}{If \code{TRUE}, performs the RegionRemoval routine.}
  \item{Za}{If \code{TRUE}, performs the ZoneAggregation routine.}
  \item{N}{If \code{TRUE}, performs the Normalization routine.}

  
Other subfunction parameters are:

  \item{l}{A positive number indicating which line of the title file to use as spectra names.}
  \item{subdirs}{If \code{TRUE}, will search inside subdirectories for FIDs and will merge them to have one FID matrix and one info matrix.}
  \item{group_delay}{If given, it is used instead of \code{Fid_info} to give how much the FIDs must be shifted to the left. It can be non-integer, in that case the values are interpolated. However it has to be non-negative since in our practical case, it would make no sense to add a part of the end of the FID at the beginning.}
  \item{lambda.ss}{Penalty on roughness used to calculate the smoothed version of the FID.
  The higher lambda is, the smoother the detected solvent signal will be.}
  \item{ptw.ss}{If \code{TRUE}, calculates the solvent signal in C using the \code{ptw} package which is a lot faster.}
  \item{DT}{If given, used instead of \code{Fid_info} to give the Dwell Time, the time between 2 points of the FID.}
  \item{type.apod}{Type of apodization, see details.}
  \item{phase}{Pphase at which the apodization window is maximum for \code{cos2}, \code{hanning} and \code{hamming} types. For example, if phase is 0.2, the maximum is at 20\% of the signal.}
  \item{rectRatio}{If there is a rectangular window, ratio between the width of the window and the width of the signal.}
  \item{gaussLB}{Line Broadening for the gaussian window, see details.}
  \item{expLB}{Line Broadening for the exponential window, see details.}
  \item{SW_h}{Sweep Width in Hertz, \emph{i.e.} the period of the Fourier Transform. It is equal to \eqn{1/DT} where \eqn{DT} is the Dwell Time, the time between two points in the time domain. If it is given, the value in \code{Fid_info} is ignored.}
  \item{Angle}{If not \code{NULL}, is a numeric vector with angles specified in radian to maually rotate the spectra. One angle is given per spectrum.}
  \item{p.zo}{Positive and negative intensities \emph{resp.} above and bellow the \emph{p}th quantile are averaged to detect an incorrect rotation due to a high water signal. p is the probability for the sample quantile used as threshold.}
  \item{ptw.bc}{If \code{TRUE}, calculates the baseline in C using the \code{ptw} library which is a lot faster.
  The R version is only kept because it is easier to understand than C and in case of problems with the installation of \code{ptw}.}
  \item{maxIter}{Maximum of iterations if \code{ptw.bc} is set to \code{FALSE}.}
  \item{lambda.bc}{Smoothing parameter (generally 1e5 -- 1e8).}
  \item{p.bc}{Asymmetry parameter.}
  \item{eps}{Numerical precision for convergence when estimating the baseline.}
  \item{shiftHandling}{See the details section.}
  \item{c}{Parameter used to fix the threshold for the TMSP peak.}
  \item{normalization.type}{Type of normalization applied to the spectrum prior to warping.
See \code{\link{Normalization}} for details about the different types.
\code{none} means that no normalization needs to be applied (\emph{e.g.} you have already done it prior to calling \code{Warping}).
It is advised to use \code{median} instead of \code{mean}, see the reference below.}
  \item{from.normW}{Used by \code{\link{Normalization}} when \code{normalization.type} is \code{peak}.}
  \item{to.normW}{Used by \code{\link{Normalization}} when \code{normalization.type} is \code{peak}.}
  \item{reference.choosing}{Specifies how the reference will be chosen.
   \describe{
     \item{\code{fixed}}{The reference is specified by the rowname given in \code{reference}.}
     \item{\code{before}}{The reference is taken as the spectrum with the minimum sum of square difference with the other spectra.}
     \item{\code{after}}{Each spectrum is taken as the reference and the sum of square difference with the other spectra is calculated after the warping.
The reference that has the minimum sum is taken as the reference and the warped spectra according to this reference (that have already been calculated at this stage) are returned. This is \eqn{n}{n} times slower than the 2 others where \eqn{n}{n} is the number of spectra.}
}}
  \item{reference}{The row number of the reference when \code{reference.choosing} is \code{fixed}.}
  \item{optim.crit}{If \code{ptw.wp} is set to \code{TRUE}, \code{WCC} can also be considered as a criterion for optimization, see \code{ptw::ptw} for details.}
  \item{ptw.wp}{If set to \code{TRUE}, it uses \code{ptw::ptw} to do the warping.
In this case the arguments \code{K}, \code{L}, \code{deg}, \code{lambda.bspline}, \code{kappa} and \code{max_it_Bspline} are ignored.
The warping does not use Bspline and \code{K} is fixed to 2, \emph{i.e.} a second order polynomial.}
  \item{K}{It is the degree of the polynomial used for the warping (see details).}
  \item{L}{This is the number of B-Splines that are used for the warping.
It should be either 0 or greater than \code{deg}.}
  \item{lambda.smooth}{Nonnegative coefficient for the smoothing \code{lambda_smooth}. \code{0} means no smoothing.}
  \item{deg}{Degree of the B-Splines.}
  \item{lambda.bspline}{Nonnegative roughness penalty coefficient for the B-Spline warping. See the reference for more details.}
  \item{kappa}{Nonnegative ridge penalty coefficient for the B-Spline warping.
See the reference for more details.}
  \item{max_it_Bspline}{Maximum number of iteration for the B-Spline warping.}
  \item{from.ws}{The lower ppm value of the interval. A typical value is 0.2. If NULL, default value is the lower index without NA.}
  \item{to.ws}{The higher ppm value of the interval, it should be greater or equal to from. A typical value is 10. If NULL, default value is the higher index without NA}
  \item{reverse.axis}{If \code{TRUE}, the ppm interval and the spectrum are reversed.}
  \item{m}{The number of buckets, should be an positive integer smaller or equal to the number of frequencies in \code{Spectrum_data}.}
  \item{width}{If \code{width} is not \code{FALSE}, then \code{m} represents the buckets width, otherwise, it represents the number of buckets.}
  \item{ppmleft}{Left boundary for ppm integration.}
  \item{ppmright}{Right boundary for ppm integration.}
  \item{intmeth}{Type of bucketing: rectangular (\code{"r"}) or trapezoidal (\code{"t"}). See details below.}
  \item{tolbuck}{Tolerance threshold to check if the buckets of the original spectra are of constant length.}
  
  \item{typeofspectra}{Type of spectra, if not \code{NULL}, will automatically remove unwanted regions wrt the nature of the spectra.}
  \item{type.rr}{The type of region removal method.}
  \item{fromto.rr}{List containing the extremities of the intervals to be removed.}  
  \item{fromto.za}{List containing the borders in ppm of the intervals to aggregate.}
  \item{type.norm}{Different types of normalization are available. See the details section.}
  \item{from.norm}{Used if \code{type.norm} is \code{peak}. See details.}
  \item{to.norm}{Used if \code{type.norm} is \code{peak}. See details.}
  \item{ref.norm}{Row index of the reference spectrum for normalization if \code{type.norm} is \code{pqn}. See details.}   
}
\details{
Note that few steps cannot be skipped : \code{\link{FourierTransform}}, \code{\link{PPMConversion}} \& \code{\link{WindowSelection}}.
Here is a reminder of all the steps involved:
\describe{
  \item{\code{\link{ReadFids}}}{Read FIDs in Bruker format from a directory.}
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
\value{
If \code{RetArgs = TRUE}, will return a list with the following elements: \code{Spectra}, \code{Arguments} and \code{Fid_info}. Otherwise, the function will just return \code{Spectra}.
\item{Spectra}{The matrix of pre-treated spectra.}
\item{Arguments}{The function parameters.}
\item{Fid_info}{The matrix containing the info about the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}

}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}(Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{Manon Martin
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data.path <-  system.file("extdata/HumanSerum", package = "SOAP")
dir(data.path)

Res <- PreprocessingChain(dataname="HumanSerum",data.path = data.path)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

