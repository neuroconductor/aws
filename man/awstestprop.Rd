\name{awstestprop}
\alias{awstestprop}
\alias{pawstestprop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagation condition for adaptive weights smoothing
}
\description{
The function enables testing of the propagation condition in order to select
appropriate values for the parameter \code{lambda} in function \code{aws}.
}
\usage{
awstestprop(dy, hmax, theta = 1, family = "Gaussian", lkern = "Triangle",
            aws = TRUE, memory = FALSE, shape = 2, homogeneous=TRUE, varadapt=FALSE,
            ladjust = 1, spmin=0.25, seed = 1, minlevel=1e-6, maxz=25, diffz=.5,
            maxni=FALSE, verbose=FALSE)
pawstestprop(dy, hmax, theta = 1, family = "Gaussian", lkern = "Triangle",
             aws = TRUE, patchsize=1, shape = 2,
             ladjust = 1, spmin = 0.25, seed = 1, minlevel = 1e-6,
             maxz = 25, diffz = .5, maxni = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dy}{ Dimension of grid used in 1D, 2D or 3D. May also be specified as an array of values.
  In this case data are generated with parameters \code{dy-mean(dy)+theta} and the propagation condition
  is testet as if \code{theta} is the true parameter. This can be used to study properties for a
  slighty misspecified structural assumption.
}
  \item{hmax}{
    Maximum bandwidth.
}
  \item{theta}{
    Parameter determining the distribution in case of
        \code{family \%in\% c("Poisson","Bernoulli")}
}
  \item{family}{
\code{family} specifies the probability distribution. Default is \code{family="Gaussian"}, also implemented
   are "Bernoulli", "Poisson", "Exponential", "Volatility", "Variance" and "NCchi". \code{family="Volatility"} specifies a Gaussian distribution with
   expectation 0 and unknown variance. \code{family="Volatility"} specifies that \code{p*y/theta} is distributed as \eqn{\chi^2} with \code{p=shape}
   degrees of freedom. \code{family="NCchi"} uses a noncentral Chi distribution with \code{p=shape} degrees of freedom and noncentrality parameter \code{theta}.
}
  \item{lkern}{
character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian"}
  \item{aws}{
logical: if TRUE structural adaptation (AWS) is used. }
  \item{patchsize}{
patchsize in case of paws. }
  \item{memory}{
logical: if TRUE stagewise aggregation is used as an additional
   adaptation scheme.}
  \item{shape}{
Allows to specify an additional shape parameter for certain family models. Currently only used for family="Variance", that is \eqn{\chi}-Square distributed observations
  with \code{shape} degrees of freedom.}
  \item{homogeneous}{
    if \code{homgeneous==FALSE} and \code{family==Gaussian} then create heterogeneous variances according to
    a chi-squared distribution with number of degrees of freedom given by \code{sphere}}
  \item{varadapt}{
    if \code{varadapt==TRUE} use inverse of variance reduction instead of sum of weights in definition of statistical penalty.}
  \item{ladjust}{
Factor to increase the default value of lambda }
  \item{spmin}{Determines the form (size of the plateau) in the adaptation kernel.
    Not to be changed by the user.}
  \item{seed}{Seed value for random generator.}
  \item{minlevel}{Minimum exceedence probability to use in contour plots.}
  \item{maxz}{Maximum of z-scale in plots.}
  \item{diffz}{Gridlength in z}
  \item{maxni}{If TRUE use \eqn{max_{l<=k}(N_i^{(l)}} instead of \eqn{(N_i^{(k)}} in the definition of the statistical penalty.}
  \item{verbose}{If TRUE provide additional information.}
}
\details{
Estimates exceedence probabilities

Results for intermediate steps are provided as contour plots. For a good choice of lambda
(ladjust) the contours up to probabilities of \code{1e-5} should be vertical.
}
\value{
A list with components
\item{h}{Sequence of bandwidths used}
\item{z}{\code{seq(0,30,.5)}, the quantiles exceedence probabilities refer to}
\item{prob}{the matrix of exceedence probabilities, columns corresponding to \code{h}}
\item{probna}{the matrix of exceedence probabilities for corresponding nonadaptive estimates, columns corresponding to \code{h}}
}
\references{
Becker (2013)
%% ~put references to the literature/web site here ~
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
\code{\link{aws}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ nonparametric}
