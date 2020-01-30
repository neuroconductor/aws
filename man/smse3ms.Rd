\name{smse3ms}
\alias{smse3}
\alias{smse3ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptive smoothing in orientation space SE(3)
}
\description{
The functions perform adaptive weights smoothing for data in orientation space SE(3),
e.g. diffusion weighted MR data,
with spatial coordinates given by voxel location within a mask and spherical information given
by gradient direction. Observations can belong to different shells characterized by b-value \code{bv}.
The data provided should only refer to voxel within mask.
}
\usage{
smse3ms(sb, s0, bv, grad, kstar, lambda, kappa0, mask, sigma,
    ns0 = 1, ws0 = 1, vext = NULL, ncoils = 1, verbose = FALSE, usemaxni = TRUE)
smse3(sb, s0, bv, grad, mask, sigma, kstar, lambda, kappa0,
    ns0 = 1, vext = NULL, vred = 4, ncoils = 1, model = 0, dist = 1,
    verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sb}{
  2D array of diffion weighted data, first dimension refers to index ov voxel within the mask,
  second dimension to the number diffusion weighted images.
}
  \item{s0}{
  vector of length \code{sum(mask)} containing values within mask of an average non-diffusion-weigthed image.
}
  \item{bv}{
   vector of b-values.
}
  \item{grad}{
   matrix of gradient directions with \code{dim(grad)[1]==3}.
}
  \item{kstar}{
  number of steps in adaptive weights smoothing.
}
  \item{lambda}{
  Scale parameter in adaptation
}
  \item{kappa0}{
determines amount of smoothing on the sphere. Larger values correspond to stronger smoothing
on the sphere. If \code{kappa0=NULL} a value is that corresponds to a variace reduction
with factor \code{vred} on the sphere.
}
  \item{mask}{
3D image defining a mask (logical)
}
  \item{sigma}{
  Error standard deviation. Assumed to be known and homogeneous in the current implementation.
  A reasonable estimate may be defined
  as the modal value of standard deviations obtained using method \code{getsdofsb}.
}
  \item{ns0}{
  Actual number of non-diffusion-weigthed images used to obtain \code{s0} by averaging.
}
  \item{ns0}{
  Weight for non-diffusion-weigthed images in statistical penalty.
}
  \item{vext}{
  Voxel extensions.
}
  \item{ncoils}{
  Effective number of receiver coils (in case of e.g. GRAPPA reconstructions),
  should be 1 in case of SENSE reconstructions. \code{2*ncoils} is the number of degrees of freedom of
  the intensity distribution used.
}
  \item{verbose}{
If \code{verbose=TRUE} additional reports are given.
}
  \item{usemaxni}{
If \code{"usemaxni==TRUE"} a strikter penalization is used.
}
\item{vred}{
Used if \code{kappa0=NULL} to specify the variance reduction on the
sphere when suggesting a value of \code{kappa0}.
}
\item{model}{
Determines which quantities are smoothed. Possible values are
\code{"Chi"} for observed values (assumed to be distributed as noncentral
Chi with \code{2*ncoils} degrees of freedom), \code{"Chi2"} for squares of
observed values (assumed to be distributed as noncentral
Chi-squared with \code{2*ncoils} degrees of freedom). \code{"Gapprox"}
and \code{"Gapprox2"} use a Gaussian approximation for the noncentral
Chi distribution to smooth ovserved and squared values, respectively.
}
\item{dist}{
Distance in SE3. Reasonable values are 1 (default, see Becker et.al. 2012),
2 ( a slight modification of 1: with k6^2 instead of abs(k6)) and
3 (using a 'naive' distance on the sphere)
}
}
\value{
The functions return lists with main results in components
\code{th} and \code{th0} containing the smoothed data.
}
\references{
Joerg Polzehl, Karsten Tabelow (2019). Magnetic Resonance Brain Imaging:
Modeling and Data Analysis Using R. Springer, Use R! series. Doi:10.1007/978-3-030-29184-6.

S. Becker, K. Tabelow, H.U. Voss, A. Anwander, R. Heidemann, J. Polzehl.
Position-orientation adaptive smoothing of diffusion weighted magnetic resonance
data (POAS). Medical Image Analysis, 2012, 16, 1142-1155. DOI:10.1016/j.media.2012.05.007.

S. Becker, K. Tabelow, S. Mohammadi, N. Weiskopf, J. Polzehl. Adaptive smoothing
of multi-shell diffusion-weighted magnetic resonance data by msPOAS.
Neuroimage, 2014, 95, 90-105. DOI:10.1016/j.neuroimage.2014.03.053.
}

\author{
J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\note{
These functions are intended to be used internally in package \code{dti} only.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% use one of  RShowDoc("KEYWORDS")
