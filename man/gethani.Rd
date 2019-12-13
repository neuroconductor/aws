\name{auxiliary}
\alias{gethani}
\alias{getvofh}
\alias{sofmchi}
\alias{residualSpatialCorr}
\alias{residualVariance}
\title{
 Auxiliary functions (for internal use)
}
\description{
Function \code{gethani} determines a bandwidth that leads to, for the specified kernel, a variance reduction for a non-adaptive kernel estimate by a factor of \code{value}.
\code{getvofh} calculates the sum of location weights for a given bandwidth vector and kernel.
\code{sofmchi} precomputes the variance of a non-central chi distribution 
with \code{2*L} degrees of freedom as a function of the noncentrality parameter for an interval \code{c(0,to)}.
Functions \code{residualVariance} and \code{residualSpatialCorr} are
used in package \code{fmri} to calculate variances and spatial correlations from residual objects. 
}
\usage{
gethani(x, y, lkern, value, wght, eps = 0.01)
getvofh(bw, lkern, wght)
sofmchi(L, to = 50, delta = 0.01)
residualVariance(residuals, mask, resscale = 1, compact = FALSE)
residualSpatialCorr(residuals, mask, lags = c(5, 5, 3), compact = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
lower bound of search interval
}
  \item{y}{
upper bound of search interval
}
  \item{lkern}{
code for location kernel 
}
  \item{value}{
target sum of location weights
}
  \item{wght}{
relative size of voxel dimensions \code{c(0,0)} for 1D and \code{c(w1,0)}
for 2D problems.
}
  \item{eps}{
attempted precision for bandwidth search
}
\item{bw}{
vector of bandwidths, length equal to 1,2 or 3 depending on the dimensionality of the problem.
}
\item{L}{number of effective coils, \code{2*L} is the degree of freedom of the non-central chi distribution.

}
\item{to}{
upper interval bound.
}
\item{delta}{
discretization width.
}
\item{residuals}{
array of residuals, if\code{compact} only containing voxel with mask, otherwise for complete data cubes.
}
\item{mask}{
mask of active voxel (e.g. brain masks)
}
\item{resscale}{
scale for residuals (residuals may be scaled for optimal integer*2 storage)
}
\item{compact}{
logical, determines if only information for voxel within mask or full for full data cubes is given.
}
\item{lags}{
positive integer vector of length 3, maximum lags for spatial correlations
for each coordinate direction to be computed
}






}
\details{
These are auxiliary functions not to be used by the user. They are only exported to be available for internal use in packages \code{fmri}, \code{dti}, \code{qMRI} and \code{adimpro}.
}
\value{
\code{gethani} returns a vector of bandwidths,
\code{getvofh} returns the variance reduction that would be obtained with a kernel estimate employing the specified kernel and bandwidth,
\code{sofmchi} returns a list with, e.g., components \code{ncp} and \code{s2} containing vectors of noncentralityparameter values and corresponding variances, respectively, for the specified noncentral Chi distribution, 
\code{residualVariance} returns a vector (compact==TRUE) or array(compact==FALSE) of voxelwise residual variances,
\code{residualSpatialCorr} returns an array of dimension \code{lags} 
containing spatial correlations.

}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\note{
These functions are for internal use only. They are only exported to be available in other packages.
}

\keyword{ misc }% use one of  RShowDoc("KEYWORDS")

