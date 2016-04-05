\name{Scaling}
\alias{Scaling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scaling of the spectra.
}
\description{
Unit Variance or Pareto scaling for spectra.
}
\usage{
Scaling(Spectrum_data, type.scaling=c("pareto", "uv"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{type.scaling}{Two types of scaling are available: Unit Variance \code{"uv"} or Pareto scaling \code{"pareto"}.}
}
\details{
If the scaling factor is 0, then it is set to 1 and no scaling is performed on these variables.
}
\value{
\item{Spectrum_data}{The matrix of scaled spectra.}
}
\references{
For more information about scaling methods:  Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & Werf, M. J. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. \emph{BMC genomics}, 7(1), 1.
}

\author{
Manon Martin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

Scaled.spec <- Scaling(Data_HS$RawSpectData_HS_12)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~manip }