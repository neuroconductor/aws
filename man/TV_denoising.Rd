\name{TV_denoising}
\alias{TV_denoising}
\alias{TGV_denoising}
\alias{TV_denoising_colour}
\alias{TGV_denoising_colour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   TV/TGV denoising of image data
}
\description{
   Total variation and total generalized variation are classical energy minimizing methods
   for image denoising.
}
\usage{
TV_denoising(datanoisy, alpha,
        iter = 1000, eps = 1e-06, epssup = 1e-04, scale = 1)
TGV_denoising(datanoisy, alpha, beta,
        iter = 1000, eps = 1e-06, epssup = 1e-04, scale = 1)
TV_denoising_colour(datanoisy, alpha,
        iter = 1000, eps = 1e-06, epssup = 1e-04, scale = 1)
TGV_denoising_colour(datanoisy, alpha, beta,
        iter = 1000, eps = 1e-06, epssup = 1e-04, scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datanoisy}{
    noisy 2D image data matrix. In case of \code{TV_denoising_colour}
    and \code{TGV_denoising_colour} and array with third dimension refering to
    RGB channels.
}
  \item{alpha}{
    TV regularization parameter.
}
  \item{beta}{
    additional TGV regularization parameter.
}
  \item{iter}{
    max. number of iterations
}
  \item{eps}{
    requested accuracy for mean image correction
}
  \item{epssup}{
    requested accuracy for max (over pixel) image correction
}
  \item{scale}{
    image scale
}
}
\details{
   Reimplementation of original matlab code by Kostas Papafitsoros (WIAS).
}
\value{
   TV/TGV reconstructed image data (2D or 3D array)
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{http://www.wias-berlin.de/people/polzehl/}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% use one of  RShowDoc("KEYWORDS")
