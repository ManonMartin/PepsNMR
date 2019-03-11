\name{Normalization}
\alias{Normalization}
\title{
Normalizes the spectra
}
\description{
Spectra normalization to correct for the dilution factor common to all biofuid samples.
}
\usage{
Normalization(Spectrum_data, type.norm, fromto.norm = c(3.05, 4.05), ref.norm = "median",
              returnFactor = FALSE, verbose = FALSE)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{type.norm}{Different types of normalization are available: \code{"mean"}, \code{"pqn"}, \code{"median"}, \code{"firstquartile"} or \code{"peak"}. No default value is provided. See the details section for more info.}
  \item{fromto.norm}{Used if \code{type.norm} is \code{"peak"}. See details.}
  \item{ref.norm}{The reference spectrum if \code{type.norm} is \code{"pqn"}. See details.}
  \item{returnFactor}{If \code{TRUE}, returns a vector with the normalization factors.}
  \item{verbose}{If \code{"TRUE"}, will print processing information.}
}
\details{

Normalization of spectra before their warping or their statistical analysis is necessary in order to be able to efficiently compare their relative peak intensities. 

It is therefore appropriate to call this filter at the end of the preprocessing workflow.

Normalization types can be:
\describe{
\item{\code{mean}}{Each spectrum is divided by its mean so that its mean becomes 1.}
\item{\code{median}}{Each spectrum is divided by its median so that its median becomes 1.}
\item{\code{firstquartile}}{Each spectrum is divided by its first quartile so that its first quartile becomes 1.}
\item{\code{peak}}{Each spectrum is divided by the value of the peak of the spectrum contained between \code{"fromto.norm"} inclusive (\emph{i.e.} the maximum value of spectral intensities in that interval).}
\item{\code{pqn}}{Probabilistic Quotient Normalization from Dieterle et al. (2006). If \code{ref.norm} is \code{"median"} or \code{"mean"}, will use the median or the mean spectrum as the reference spectrum ; if it is a single number, will use the spectrum located at that row in the spectral matrix; if \code{ref.norm} is a numeric vertor of length equal to the number of spectral variables, it defines manually the reference spectrum.}
}

The choice of a proper normalisation method is a crucial although not straightforward step in a metabolomic analysis.


Applying CSN is accurate in the following situations: 
\itemize{
\item	when working on human/animal sera in the case of not serious pathology, given the homeostasis principle and since no dilution effect is present. 
\item		When working on biopsies, the “metabolome quantity” is set constant across the samples by adding a varying volume of a buffer and the same applies when working with cell media, where the quantity of cells is made constant.
}

To counteract all the dilution effects and the excretion differences between urine samples, the PQN approach is often recommended in the literature (Dieterle et al., 2006).

For any other situation (large difference between the groups, other kind of sample, etc.), the choice of the normalisation method is not straightforward. A solution is to refer to endogenous stable metabolites that are present in a constant quantity across samples and use them as standards to normalize all spectral profiles. For the urine samples, the creatinine has been considered as such standard (this option is also implemented in PepsNMR), even though it has been shown that the creatinine concentration could fluctuate given specific parameters (Tang et al., 2015). A review on normalization techniques for mass spectroscopy metabolomics from Wu \& Li (2015) provides some guidance in the choice on the normalization approach regarding the type of sample analysed and can be transposed to the NMR spectra normalisation.

}
\value{
\item{Spectrum_data}{The matrix of normalized spectra.}
}
\references{

Martin, M., Legat, B., Leenders, J., Vanwinsberghe, J., Rousseau, R., Boulanger, B., & Govaerts, B. (2018). PepsNMR for 1H NMR metabolomic data pre-processing. \emph{Analytica chimica acta}, 1019, 1-13.
  
  
  Yiman Wu, Liang Li. (2016). \emph{Sample normalization methods in quantitative metabolomics}, Journal of Chromatography A, Volume 1430, Pages 80-95, ISSN 0021-9673
 
 
 Tang KWA, Toh QC, Teo BW. (2015). \emph{Normalisation of urinary biomarkers to creatinine for clinical practice and research – when and why.} Singapore Medical Journal. 56(1):7-10. 
 
   Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
 
  
  Dieterle, F., Ross, A. , Schlotterbeck, G.,and Senn, H (2006). Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. \emph{Analytical Chemistry} 78 (13), 4281-4290 
  
}
\author{
Benoît Legat & Manon Martin
}
\examples{
require(PepsNMRData)
Norm.spec <- Normalization(Data_HS_sp$Spectrum_data_HS_12, 
                                type.norm = "mean")
}
\keyword{ manip }

