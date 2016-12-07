\name{Data_HU_sp}
\alias{Data_HU_sp}
\docType{data}
\title{FIDs or raw spectra from the Human Urine database.
%%   ~~ data name/kind ... ~~
}
\description{Four first FIDs or raw spectra of the Human Urine database given after each preprocessing step.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("Data_HU_sp")}
\format{
  A list with 15 objects that are complex matrices with 4 observations/object
  \describe{
    \item{\code{FidData_HU_0}}{Raw FIDs of the Human Urine database.}
    \item{\code{FidData_HU_1}}{FIDs of the Human Urine database after First Order Phase Correction.}
    \item{\code{FidData_HU_2}}{FIDs of the Human Urine database after Solvent residuals suppression.}
    \item{\code{FidData_HU_3}}{FIDs of the Human Urine database after Apodization.}     \item{\code{RawSpectData_HU_4}}{Spectra of the Human Urine database after Fourier transformation.}   
    \item{\code{RawSpectData_HU_5}}{Spectra of the Human Urine database after Zero Order Phase Correction.}   
    \item{\code{RawSpectData_HU_6}}{Spectra of the Human Urine database after Baseline correction.}   
    \item{\code{RawSpectData_HU_7}}{Spectra of the Human Urine database after Negative Values Zeroing.}   
    \item{\code{RawSpectData_HU_8}}{Spectra of the Human Urine database after ppm conversion.}   
    \item{\code{RawSpectData_HU_9}}{Spectra of the Human Urine database after Warping.}   
    \item{\code{RawSpectData_HU_10}}{Spectra of the Human Urine database after window selection.}   
    \item{\code{RawSpectData_HU_11}}{Spectra of the Human Urine database after Bucketing}   
    \item{\code{RawSpectData_HU_12}}{Spectra of the Human Urine database after Region removal.}  
    \item{\code{RawSpectData_HU_13}}{Spectra of the Human Urine database after Zone aggregation.}      
    \item{\code{RawSpectData_HU_14}}{Spectra of the Human Urine database after Normalization.}   
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.

}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

plot(Re(Data_HU_sp$FidData_HU_0[1,]), type="l") 
comment(Data_HU_sp$FidData_HU_0)

}
\keyword{datasets}
