\name{awsLocalSigma}
\alias{awsLocalSigma}
\alias{awslinsd}
\alias{estimateSigmaCompl}
\alias{AFLocalSigma}
\alias{estGlobalSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D variance estimation
}
\description{
Functions for 3D variance estimation
}
\usage{
awsLocalSigma(y, steps, mask, ncoils, vext = c(1, 1), lambda = 5,
    minni = 2, hsig = 5, sigma = NULL, family = c("NCchi", "Gauss"),
    verbose = NULL, trace = FALSE, u = NULL)
awslinsd(y, hmax = NULL, hpre = NULL, h0 = NULL, mask = NULL,
    ladjust = 1, wghts = NULL, varprop = 0.1, A0, A1)
AFLocalSigma(y, ncoils, level = NULL, mask = NULL, h = 2, hadj = 1,
    vext = c(1, 1))
estGlobalSigma(y, mask = NULL, ncoils = 1, steps = 16, vext = c(1, 1),
    lambda = 20, hinit = 2, hadj = 1, q = 0.25, level = NULL,
    sequence = FALSE, method = c("awsVar", "awsMAD", "AFmodevn",
                "AFmodem1chi", "AFbkm2chi", "AFbkm1chi"))
estimateSigmaCompl(magnitude, phase, mask, kstar = 20, kmin = 8, hsig = 5,
        lambda = 12, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    3D array of image intensities.
}
  \item{steps}{
number of steps in adapive weights smoothing, used to reveal the unerlying
   mean structure.

}
  \item{mask}{
restrict computations to voxel in mask, if \code{is.null(mask)} all voxel are used.
   In function \code{estGlobalSigma} mask should refer to background for \code{method \%in\% c("modem1chi","bkm2chi","bkm1chi")} and to voxel within the head for 
   \code{method=="modevn"}. 
}
  \item{ncoils}{
effective number of coils, or equivalently number of effective degrees of freedom of non-central chi distribution divided by 2.
}
  \item{vext}{
  voxel extentions or relative voxel extensions
}
  \item{lambda}{
   scale parameter in adaptive weights smoothing
}
  \item{minni}{
   minimal bandwidth for calculating local variance estimates 
}
  \item{hsig}{
  bandwwidth for median filter
}
  \item{sigma}{
   optional initial global variance estimate
}
  \item{family}{
  type of distribution, either noncentral Chi ("NCchi") or Gaussian ("Gauss")
}
  \item{verbose}{
   if \code{verbose==TRUE} density plots 
   and quantiles of local estimates of \code{sigma} are provided.
}
  \item{trace}{
   if \code{trace==TRUE} intermediate results for each step are
   returned in component tergs for all voxel in mask.
}
  \item{u}{if \code{verbose==TRUE} an array of noncentrality paramters for
      comparisons. Internal use for tests only
}
\item{hmax}{
   maximal bandwidth
}
\item{hpre}{
   minimal bandwidth
}
\item{h0}{
   bandwidth vector characterizing to spatial correlation as correlation induced
   by convolution with a Gaussian kernel
}
\item{ladjust}{
   correction factor for lambda
}
\item{wghts}{
      relative voxel extensions
}
\item{varprop}{
    defines a lower bound for the estimated variance as \code{varprop*mean(sigma2hat}
}
\item{A0}{
    select voxel with \code{A0 < theta < A1} to estimate parameters of the variance model
}
\item{A1}{
    select voxel with \code{A0 < theta < A1} to estimate parameters of the variance model
}
\item{level}{
    threshold for mask definition
}
\item{h}{
    bandwidth for local variance estimates.
}
\item{hinit}{
    minimal bandwidth for local variance estimates with \code{method="awsxxx"}.
}
\item{hadj}{
    bandwidth for mode estimation
}
\item{q}{
    Quantile for interquantile estimate of standard deviation
}
\item{sequence}{
    logical, return sequence of estimated variances for iterative methods.
}
\item{method}{
    determines variance estimation method
}
\item{magnitude}{
    magnitude of complex 3D image
}
\item{phase}{
    phase of complex 3D image
}
  \item{kstar}{
number of steps in adapive weights smoothing, used to reveal the unerlying
   mean structure.

}
\item{kmin}{
    iteration to start adaptation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
all functions return lists with variance estimates in component \code{sigma}
}
\references{
K. Tabelow, H.U. Voss, J. Polzehl, Local estimation of the noise level in MRI using structural adaptation, Medical Image Analysis, 20 (2015), pp. 76--86.

}
\author{
J\"org Polzehl \email{polzehl@wias-berlin.de}
}

% R documentation directory.
\keyword{ smooth }% use one of  RShowDoc("KEYWORDS")
\keyword{ utilities }% __ONLY ONE__ keyword per line
