\name{Bucketing}
\alias{Bucketing}
\title{
Spectral data reduction
}
\description{
Reduces the number of data points by aggregating frequencies into \eqn{m}{m} buckets.

The input spectrum can have its ppm axis in increasing or decreasing order and it does not
have to be equispaced.
}
\usage{
Bucketing(Spectrum_data, m = 500)
}
\arguments{
  \item{Spectrum_data}{Matrix containing the spectra in ppm, one row per spectrum.}
  \item{m}{The number of buckets, should be an positive integer smaller or equal to the number of frequencies in \code{Spectrum_data}.}
}
\details{
Bucketing will decrease data precision but has two main interests:
\itemize{
  \item{Facilitate the statistical analysis}
  \item{Decrease the impact of peaks misalignment between different spectra that should be aligned; assuming we are in the ideal case where they fall in the same bucket.
  Of course, the better the prior warping, the larger can be \eqn{m}{m} and the more informative the spectra will be.}
}

  The ppm interval of \code{Spectrum_data} - let's say is \eqn{[a,b]}{[a,b]} where \eqn{a > b}{a > b} - is divided into \eqn{m}{m} buckets of size \eqn{(a-b)/m}{(a-b)/m}.
  The new ppm scale contains the \eqn{m}{m} centers of these intervals.
  The value of the spectrum at these centers is the integral of the initial spectral intensity on
  this bucket using the trapezium rule divided by the length of the bucket.

  This way, the integral of the initial signal is the same
  as the integral of the bucketted signal using the trapezium rule (integral in \code{#} in the figure below), plus if the integral (area in \code{$} in the figure below) between \eqn{a}{a} and the first center and the integral between
  the last center and \eqn{b}{b} where we consider the signal constant and with the same
  value as the value of the signal at the first center (resp. last center).

  It can be explained by the fact that the trapezium rule for equi-spaced points (space \eqn{h}{h}) is the sum of the middle points times \eqn{h}{h} plus the sum of the 2 extremities times \eqn{h/2}{h/2}.
  Adding the integral between \eqn{a}{a} and  the first point and between the last point and \eqn{b}{b} we therefore have the sum of all our points times \eqn{h}{h} which is by construction the trapezium rule applied to the initial signal.

  \preformatted{
        __
   a_  /##\
    $\/####\
    $#######\  _b
    $########\/$
  }
}
\value{
\item{Spectrum_data}{The matrix of spectra with their new ppm axis.}
}
\references{
  Rousseau, R. (2011). \emph{Statistical contribution to the analysis of metabonomics data in 1H NMR spectroscopy}
  (Doctoral dissertation, PhD thesis. Institut de statistique, biostatistique et sciences actuarielles, Université catholique de Louvain, Belgium).
}
\author{
Benoît Legat, inspired from a code written by Bernadette Govaerts.
}
\note{
}

\seealso{
}
\examples{
Bucket.spec <- Bucketing(Data_HS$RawSpectData_HS_10, m = 500)
}
\keyword{ manip }
\keyword{ array }
