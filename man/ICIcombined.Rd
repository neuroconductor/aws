\name{ICIcombined}
\alias{ICIcombined}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptive smoothing by Intersection of Confidence Intervals (ICI) using multiple windows
}
\description{
The function performs adaptive smoothing by Intersection of Confidence Intervals (ICI)
using multiple windows as described in Katkovnik et al (2006)
}
\usage{
ICIcombined(y, hmax, hinc = 1.45, thresh = NULL, kern = "Gaussian", m = 0,
            sigma = NULL, nsector = 1, symmetric = FALSE, presmooth = FALSE,
            combine = "weighted", unit = c("SD","FWHM"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Object of class \code{"array"} containing the original (response) data
  on a grid
}
  \item{hmax}{
maximum bandwidth
}
  \item{hinc}{
factor used to increase the bandwidth from scale to scale
}
  \item{thresh}{
threshold used in tests to determine the best scale
}
  \item{kern}{Determines the kernel function.
Object of class \code{"character"} kernel, can be any of
     \code{c("Gaussian","Uniform","Triangle","Epanechnicov","Biweight","Triweight")}.
Defaults to \code{kern="Gaussian"}.
}
  \item{m}{
Object of class \code{"integer"} vector of length \code{length(dy)}
       determining the order of derivatives specified for the coordinate directios.}
  \item{sigma}{
  error standard deviation
}
  \item{nsector}{
     number of sectors to use.
}
  \item{symmetric}{
Object of class \code{"logical"} determines if sectors are symmetric with respect to
the origin.
}
  \item{presmooth}{
    Object of class \code{"logical"} determines if bandwidths are  smoothed
    for more stable results.
}
  \item{combine}{
    Either \code{"weighted"} or \code{"minvar"}. Determines how whether to combine
    sectorial results a weighted (with inverse variance) mean or to chose the sectorial
    estimate with minimal variance.
}
\item{unit}{
 How should the bandwidth be interpreted in case of a Gaussian kernel.
 For \code{"SD"} the bandwidth refers to the standard deviation of the
 kernel while \code{"FWHM"} interprets the banwidth in terms of Full Width Half
  Maximum of the kernel.
}
}
\details{
This mainly follows Chapter 6.2 in Katkovnik et al (2006).
}
\value{
An object of class \code{ICIsmooth}
}
\references{
             J. Polzehl, K. Papafitsoros, K. Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
             Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06.

V. Katkovnik, K. Egiazarian and J. Astola,
\emph{Local Approximation Techniques in Signal And Image Processing},
 SPIE Society of Photo-Optical Instrumentation Engin., 2006, PM157
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

\seealso{
\code{\link{ICIsmooth}}, \code{\link{ICIsmooth-class}}, \code{\link{kernsm}}
}
\keyword{ smooth }
\keyword{ nonparametric }
