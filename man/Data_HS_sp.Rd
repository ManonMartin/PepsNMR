\name{Data_HS_sp}
\alias{Data_HS_sp}
\docType{data}
\title{FIDs and spectra from the Human Serum database.}
\description{4 first FIDs and spectra of the Human Serum database after each preprocessing step.
For more details on this database, use \code{help(FidData_HS)}.}
\usage{data("Data_HS_sp")}
\format{
  A list with 14 objects that are complex matrices. 4 observations/object:
  \describe{
    \item{\code{FidData_HS_0}}{Raw FIDs.}
    \item{\code{FidData_HS_1}}{FIDs after first order phase correction.}
    \item{\code{FidData_HS_2}}{FIDs after solvent residuals suppression.}
    \item{\code{FidData_HS_3}}{FIDs after apodization.}     
    \item{\code{Spectrum_data_HS_4}}{Spectra after fourier transformation.} 
    \item{\code{Spectrum_data_HS_5}}{Spectra after zero order phase correction.}  
    \item{\code{Spectrum_data_HS_6}}{Spectra after internal calibration.}
    \item{\code{Spectrum_data_HS_7}}{Spectra after baseline correction.}   
    \item{\code{Spectrum_data_HS_8}}{Spectra after negative values zeroing.}   
    \item{\code{Spectrum_data_HS_9}}{Spectra after warping.}   
    \item{\code{Spectrum_data_HS_10}}{Spectra after window selection.}   
    \item{\code{Spectrum_data_HS_11}}{Spectra after bucketing.}   
    \item{\code{Spectrum_data_HS_12}}{Spectra after region removal.}   
    \item{\code{Spectrum_data_HS_13}}{Spectra after normalization.}   
  }
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
}

\examples{

plot(Re(Data_HS_sp$FidData_HS_0[1,]), type="l") 
comment(Data_HS_sp$FidData_HS_0)

}
\keyword{datasets}
