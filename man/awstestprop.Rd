\name{awstestprop}
\alias{awstestprop}
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
        aws = TRUE, memory = FALSE, shape = 2, ladjust = 1, seed = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dy}{ Dimension of grid used in 1D, 2D or 3D
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
   are "Bernoulli", "Poisson", "Exponential", "Volatility" and "Variance". \code{family="Volatility"} specifies a Gaussian distribution with
   expectation 0 and unknown variance. \code{family="Volatility"} specifies that \code{p*y/theta} is distributed as \eqn{\chi^2} with \code{p=shape}
   degrees of freedom.
}
  \item{lkern}{
character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian"}
  \item{aws}{
logical: if TRUE structural adaptation (AWS) is used. }
  \item{memory}{
logical: if TRUE stagewise aggregation is used as an additional
   adaptation scheme.}
  \item{shape}{
Allows to specify an additional shape parameter for certain family models. Currently only used for family="Variance", that is \eqn{\chi}-Square distributed observations
  with \code{shape} degrees of freedom.}
  \item{ladjust}{
Factor to increase the default value of lambda }
  \item{seed}{
Seed value for random generator.
}
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
\item{pofalpha}{the matrix of exceedence probabilities with respect to a grid of \code{alpha}
      for the direct propagation condition}
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