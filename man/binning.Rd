\name{binning}
\alias{binning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Binning in 1D, 2D or 3D}
\description{
  The function performs a binning in 1D, 2D or 3D.
}
\usage{
binning(x, y, nbins, xrange = NULL)
}
\arguments{
  \item{x}{ design matrix, dimension n x d, \code{d \%in\% 1:3}. }
  \item{y}{ either a response vector of length n or NULL }
  \item{nbins}{ vector of length d containing number of bins for each dimension,
    may be set to NULL }
  \item{xrange}{ range for endpoints of bins for each dimension, either matrix
   of dimension 2 x d or NULL. \code{xrange} is increased if the cube defined  does not contain all design points. }
}
\value{ A list with components
\item{x}{matrix of coordinates of non-empty bin centers }
\item{x.freq}{number of observations in nonempty bins}
\item{midpoints.x1}{Bin centers in dimension 1}
\item{midpoints.x2}{if \code{d>1} Bin centers in dimension 2}
\item{midpoints.x3}{if \code{d>2} Bin centers in dimension 3}
\item{breaks.x1}{Break points dimension 1}
\item{breaks.x2}{if \code{d>1} Break points dimension 2}
\item{breaks.x3}{if \code{d>2} Break points dimension 3}
\item{table.freq}{number of observations per bin}
\item{means}{if \code{!is.null(y)} mean of y in non-empty bins}
\item{devs}{if \code{!is.null(y)} standard deviations of y in non-empty bins}
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\note{ This function has been adapted from the code of function binning in package sm.
}
\seealso{See Also as \code{\link{aws.irreg}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
