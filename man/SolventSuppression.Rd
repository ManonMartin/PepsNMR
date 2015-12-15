\name{SolventSuppression}
\alias{SolventSuppression}
\title{
Suppress the Solvent signal present in each FID.
}
\description{
Signal smooting for water residuals resonance removal.

}
\usage{
SolventSuppression(Fid_data, lambda.ss = 1e6, ptw.ss=TRUE, 
                              plotSolvent=F, returnSolvent=F)
}
\arguments{
  \item{Fid_data}{Matrix containing the FIDs, one row per signal, as outputted by \code{\link{ReadFids}}.}
  \item{lambda.ss}{Penalty on roughness used to calculate the smoothed version of the FID.
  The higher lambda is, the smoother the detected solvent signal will be.}
  \item{ptw.ss}{If \code{TRUE}, calculates the solvent signal in C using the \code{ptw} package which is a lot faster.
  The R version is only kept because it is easier to understand than C and in case of problems with the installation of \code{ptw}.}
  \item{plotSolvent}{If \code{TRUE}, plots the solvent signal on the current device. If you want a specific device or device parameter, set it before calling \code{SolventSuppression}.}
    \item{returnSolvent}{If \code{TRUE}, returns a list with the resulting FID, the real part of the solvent signal and the imaginary part of the solvent signal, see the Examples.}
}
\details{
FIDs usually present a wavy shape. Under the assumption that water is the main component of the analyzed solution, its signal can be modelled by the smoothing of the composite FID.
We then substract this wave, \emph{i.e.} the solvent residuals resonance, from the original FID.


The smoothing is done with a Whittaker smoother which is obtained by the minimization of
\deqn{V + \lambda R}{V + \lambda R}
where
\itemize{
  \item \eqn{V} is the sum of the squared differences between the original and the smoothed signal.
  \item \eqn{R} measures the roughness of the smoothed signal. The smoother it is the lower \eqn{R} will be.
}
}
\value{
\item{Fid_data}{The matrix of FIDs with the solvent residuals signal removed.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).

  Eilers, PHC. (2003). \emph{A perfect smoother}. Analytical Chemistry,  75(14), 3631-3636.

}
\author{
Benoît Legat
}
\note{
}

\seealso{
See also \code{\link{BaselineCorrection}} which also uses the Whittaker smoother.
}
\examples{
Ss.fid <- SolventSuppression(Data_HS$FidData_HS_1, plotSolvent=FALSE, returnSolvent=FALSE)

#or
Ss.res <- SolventSuppression(Data_HS$FidData_HS_1, plotSolvent=FALSE, returnSolvent=TRUE)
Ss.fid = Ss.res[["Fid_data"]]
SolventRe = Ss.res[["SolventRe"]]
plot(SolventRe[1,], type="l")

}
\keyword{ manip }
