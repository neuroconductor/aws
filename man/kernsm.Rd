\name{kernsm}
\alias{kernsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel smoothing on a 1D, 2D or 3D grid
}
\description{
Performs Kernel smoothing on a 1D, 2D or 3D grid by fft
}
\usage{
kernsm(y, h = 1, kern = "Gaussian", m = 0, nsector = 1, sector = 1,
       symmetric = FALSE, unit = c("SD","FWHM"))
}
\arguments{
  \item{y}{
Object of class \code{"array"} containing the original (response) data
  on a grid
}
  \item{h}{
bandwidth
}
  \item{kern}{
Determines the kernel function.
Object of class \code{"character"} kernel, can be any of
     \code{c("Gaussian","Uniform","Triangle","Epanechnicov","Biweight","Triweight")}.
Defaults to \code{kern="Gaussian"}
}
  \item{m}{
Object of class \code{"integer"} vector of length \code{length(dy)}
       determining the order of derivatives specified for the coordinate directios.}
  \item{nsector}{
    number of sectors to use. Positive weights are restricted to the sector
     selected by \code{sector}
}
  \item{sector}{
Object of class \code{"integer"} between \code{1} and \code{nsector}. sector used.
}
  \item{symmetric}{
Object of class \code{"logical"} determines if sectors are symmetric with respect to the origin.
}
\item{unit}{
 How should the bandwidth be interpreted in case of a Gaussian kernel.
 For \code{"SD"} the bandwidth refers to the standard deviation of the
 kernel while \code{"FWHM"} interprets the banwidth in terms of Full Width Half
  Maximum of the kernel.
}
}
\details{
In case of \code{any(m>0)} derivative kernels are generated and applied for the
corresponding coordinate directions. If \code{nsector>1} the support of the kernel
is restricted to a circular sector determined by \code{sector}.
}
\value{
An object of class \code{kernsm}
}
\references{
             J. Polzehl, K. Papafitsoros, K. Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
             Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06 (URL: http://doi.org/10.18637/jss.v095.i06).

V. Katkovnik, K. Egiazarian and J. Astola,
\emph{Local Approximation Techniques in Signal And Image Processing},
 SPIE Society of Photo-Optical Instrumentation Engin., 2006, PM157
}
\author{
%%  ~~who you are~~
Joerg Polzehl \email{polzehl@wias-berlin.de}
}


\seealso{
 \code{\link{kernsm-class}}, \code{\link{ICIsmooth}},\code{\link{ICIcombined}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ nonparametric }
