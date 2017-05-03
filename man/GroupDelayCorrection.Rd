\name{GroupDelayCorrection}
\alias{GroupDelayCorrection}
\title{
Perform a first order phase correction.
}
\description{
The function removes the group delay at the beginning of the FIDs. 

}
\usage{
GroupDelayCorrection(Fid_data, Fid_info = NULL, group_delay = NULL)
}
\arguments{
  \item{Fid_data}{Matrix containing the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{Fid_info}{Matrix containing the info about the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{group_delay}{If given, it is used instead of \code{Fid_info} to decide how much the FIDs must be shifted to the left. It can be non-integer, in that case the values are interpolated. However it has to be non-negative since in our practical case, it would make no sense to add a part of the end of the FID at the beginning.}
}
\details{
First Order Phase Correction step could also called "removal of Bruker digital filter".

Due to Bruker's digital filter and to other technical reasons a first order phase shift caused by a group delay is present in the FID and needs to be removed.
Luckily, information about this delay is available when loading the FID with \code{\link{ReadFids}} and is written in \code{Fid_info}.

This function shifts circularly each FID in order to cancel this delay.
By circularly, we mean that the starting portion of the FID becomes its ending portion when applied.

Each FID is shifted by the same amount since it can be non-integer and the columns names which are the time coordinates are shared between all the FIDs.
}
\value{
\item{Fid_data}{The matrix of FIDs corrected for the first order phase shift.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat & Manon Martin
}

\examples{
Fopc.fid <- GroupDelayCorrection(FidData_HS_sp$FidData_HS_0, FidInfo_HS)
}
\keyword{manip }
