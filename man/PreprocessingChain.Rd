\name{PreprocessingChain}
\alias{PreprocessingChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pre-treatment chain.
%%  ~~function to do ... ~~
}
\description{Pre-treatment facilities organised in the advised sequence.
}
\usage{
PreprocessingChain(title = "Run\%003d", dataname="Dataset", data.path = getwd(),
            out.path = getwd(), nspectr = 1, save = FALSE, saveall = FALSE, 
            ImpG= FALSE, RetArgs = TRUE, Fopc = TRUE, Ss = TRUE, A = TRUE, 
            Zopc = TRUE, Bc = TRUE, Zsnv = TRUE, W = TRUE, B = TRUE, Zs = TRUE, 
            Za = FALSE, N = TRUE, ...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{title}{Run title.}
  \item{dataname}{Name of the data folder contained in the \code{data.path} directory.}
  \item{data.path}{Path to the data in Bruker format from the specified directory.}
  \item{out.path}{Path for optional outputs.}
  \item{nspectr}{Row number of the spectrum selected for each step illustrations.}
  \item{save}{If \code{TRUE}, Fourier transformed and final data are saved.}
  \item{saveall}{If \code{TRUE}, datasets from each step are saved into a list.}
  \item{ImpG}{If \code{TRUE}, provides graphical illustration for each step.}
  \item{RetArgs}{If \code{TRUE}, the function arguments are saved.}
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
  \item{...}{Other arguments to be passed to inner routines.}
  
}
\details{
Note that few steps cannot be skipped which are: FourierTransform, PPMConversion \& WindowSelection
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
\item{Spectra}{The matrix of pre-treated spectra.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}(Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{Manon Martin
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data.path <-  system.file("extdata", package = "SOAP")
dir(data.path)

Res <- PreprocessingChain(title = "Run\%003d", dataname="HumanSerum", 
                                          data.path = data.path)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

