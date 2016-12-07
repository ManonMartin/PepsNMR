\name{ReadFids}
\alias{ReadFids}
\title{Read FIDs in Bruker format from a directory}
\description{Finds all directories of \code{path} which contain a valid FID (\emph{i.e.} contain the files \code{fid}, \code{acqu} and \code{acqus}) and loads
them in a matrix.
}
\usage{
ReadFids(path, l = 1, subdirs = FALSE)
}
\arguments{
  \item{path}{A character string specifying the directory where the FIDs are searched.}
  \item{l}{A positive number indicating which line of the title file to use as spectra names.}
  \item{subdirs}{If \code{TRUE}, will search inside subdirectories for FIDs and will merge them to have unique FID and info matrices.}
}
\details{
The row names are the first line of the file "pdata/1/title" in the directory or the directory name if the title file doesn't exists or the line "l" is blank. The column names are the time coordinates of the FID.
All the FIDs therefore need to have the same length and time interval between points.

\code{SW_h} is used in \code{Fid_info} as the Sweep Width in hertz.
This is half of the width (in hertz, not in angular frequency) of the Fourier Transform \emph{i.e.} \eqn{1/2DT}{1/2DT}, where DT is the Dwell Time, the time between two data points in the FID.
}
\value{
Returns a list with the FIDs and their related information.
\item{Fid_data}{The matrix containing the FIDs.}
\item{Fid_info}{A matrix containing the information about the FIDs.
The naming of the row is the same than for \code{Fid_data}.

The columns are:
\describe{
   \item{TD}{Time domain size}
   \item{BYTORDA}{Determine the endiness of stored data. If 0 -> Little Endian; if 1 -> Big Endian}
   \item{DIGMOD}{Digitization mode}
   \item{DECIM}{Decimation rate of digital filter}
   \item{DSPFVS}{DSP firmware version}
   \item{SW_h}{Sweep width in Hz}
   \item{SW}{Sweep width in ppm}
   \item{DT}{Dwell time in microseconds}
   \item{O1}{Irradiation frequency offset}
}
}
}

\author{
Benoît Legat
}


\examples{
path <-  system.file("extdata", package = "SOAP")
dir(path)

fidList_HS <- ReadFids(file.path(path, "HumanSerum"))
FidData_HS_0 <- fidList_HS[["Fid_data"]]
FidInfo_HS <- fidList_HS[["Fid_info"]]
}
\keyword{ IO }

