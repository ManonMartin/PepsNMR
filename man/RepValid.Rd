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


The csv file should be like this:
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
path <-  system.file("extdata/Group_HS.csv", package = "SOAP")
res.valid = RepValid(FinalSpectra_HS, read="csv2", gr.file = path)
#or res.valid = RepValid(FinalSpectra_HS, read="csv2", gr.file=path)

do.call(gridExtra::grid.arrange, res.valid[["plot.pca"]])
print(summary(res.valid[["res.pca"]]))
print(res.valid[["res.inertia"]])

}

\keyword{ manip }