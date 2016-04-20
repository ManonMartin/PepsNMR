\name{RepValid}
\alias{RepValid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA and inertia measures for spectral repeatability
}
\description{
This function enables to check with PCA results and inertia measures the spectral repeatability given a spectral dataset where the group corresponds to repeated measures of the same biological samples. 
}
\usage{
RepValid(Spectrum_data, group=NULL, read = c("csv2", "csv", NULL), gr.file=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra, one row per spectrum.}
  \item{group}{If provided, numeric vector containing the class information.}
    \item{gr.file}{If non \code{NULL}, character string with the complete path to the .csv file with the class information.
}
  \item{read}{If non \code{NULL}, specifies }

}
\details{
If group is \code{NULL}, the function will read a .csv file from \code{gr.file} tyhe functions \code{read.csv} or \code{read.csv2} from the \code{utils} package. The csv file needs to be formatted as a two-colunm matrix with column 1: the spectra IDs matching with the rownames of \code{Spectrum_data} and column 2: the classes defined as 1, 2, 3 etc.

  \preformatted{
  
  "ID1"  1
  "ID2"  1
  "ID3"  3
  "ID4"  2
   ...  ...
  
  }

Note that if provided, classes indicated in \code{group} need to have the same order that the row order of the spectra.

}
\value{
\item{plot.pca}{Scores plots of the 4 first PCs.}
\item{res.pca}{List with PCA results.}
\item{res.inertia}{List with inertia results.}
}

\author{
Manon Martin
}


\examples{
res.valid = RepValid(FinalSpectra_HS, group = Group_HS)

do.call(grid.arrange, res.valid[["plot.pca"]])
print(summary(res.valid[["res.pca"]]))
print(res.valid[["res.inertia"]])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
