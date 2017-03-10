\name{Data_HS_sp}
\alias{Data_HS_sp}
\docType{data}
\title{FIDs or raw spectra from the Human Serum database.
%%   ~~ data name/kind ... ~~
}
\description{4 first FIDs or spectra of the Human Serum database given after each preprocessing step.
For more details on this database, see \code{help(FidData_HS)}.
}
\usage{data("Data_HS_sp")}
\format{
  A list with 14 objects that are complex matrices. 4 observations/object.
  \describe{
    \item{\code{FidData_HS_0}}{Raw FIDs of the Human Serum database.}
    \item{\code{FidData_HS_1}}{FIDs of the Human Serum database after First Order Phase Correction.}
    \item{\code{FidData_HS_2}}{FIDs of the Human Serum database after Solvent residuals suppression.}
    \item{\code{FidData_HS_3}}{FIDs of the Human Serum database after Apodization.}     \item{\code{Spectrum_data_HS_4}}{Spectra of the Human Serum database after Fourier transformation.}   
    \item{\code{Spectrum_data_HS_5}}{Spectra of the Human Serum database after Zero Order Phase Correction.}   
    \item{\code{Spectrum_data_HS_6}}{Spectra of the Human Serum database after Baseline correction.}   
    \item{\code{Spectrum_data_HS_7}}{Spectra of the Human Serum database after Negative Values Zeroing.}   
    \item{\code{Spectrum_data_HS_8}}{Spectra of the Human Serum database after ppm conversion.}   
    \item{\code{Spectrum_data_HS_9}}{Spectra of the Human Serum database after Warping.}   
    \item{\code{Spectrum_data_HS_10}}{Spectra of the Human Serum database after window selection.}   
    \item{\code{Spectrum_data_HS_11}}{Spectra of the Human Serum database after Bucketing}   
    \item{\code{Spectrum_data_HS_12}}{Spectra of the Human Serum database after Region removal.}   
    \item{\code{Spectrum_data_HS_13}}{Spectra of the Human Serum database after Normalization.}   
  }
}

\source{Data provided by Dr. Pascal de Tullio and coworkers of the Pharmaceutical Chemistry Laboratory in the Pharmacy Department of the University of Liege (ULg), Belgium.
}

\examples{

plot(Re(Data_HS_sp$FidData_HS_0[1,]), type="l") 
comment(Data_HS_sp$FidData_HS_0)

}
\keyword{datasets}
