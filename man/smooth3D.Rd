\name{smooth3D}
\alias{smooth3D}
\alias{medianFilter3D}
\alias{aws3Dmask}
\alias{aws3Dmaskfull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary 3D smoothing routines
}
\description{
\code{smooth3D} and \code{medianFilter3D} are auxiliary functions for non-adaptive
smoothing of 3D image data using kernel or median smoothing. Both function restrict
to sub-areas determined by a mask.  The functions are used in packages \code{dti} and
\code{qMRI}.

Functions \code{aws3Dmask} and \code{aws3Dmaskfull} perform adaptive weights smoothing on
statistical parametric maps in fMRI. Variability of results is determined from smoothed
(using the same weighting schemes) residuals in order to correctly account for spatial
correlation. These functions are intended to be used internally in package \code{fmri}.
They have been moved here because they share significant parts of the
openMP parallelized Fortran code underlying function \code{aws}.
}
\usage{
smooth3D(y, lkern = "Gaussian", weighted = FALSE, sigma2 = NULL, mask = NULL,
         h = NULL, wghts = NULL)
medianFilter3D(y, h = 10, mask = NULL)
aws3Dmask(y, mask, lambda, hmax, res = NULL, sigma2 = NULL, lkern = "Gaussian",
   skern = "Plateau", weighted = TRUE, u = NULL, wghts = NULL,
   h0 = c(0, 0, 0), testprop = FALSE)
aws3Dmaskfull(y, mask, lambda, hmax, res = NULL, sigma2 = NULL, lkern = "Gaussian",
       skern = "Plateau", weighted = TRUE, u = NULL, wghts = NULL,
       testprop = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
3D array of data in case of functions \code{smooth3D} and \code{medianFilter3D}. For \code{aws3Dmask*}
with \code{!is.null(mask)} a vector of length \code{sum(mask)} containing only data values within the specified mask.
}
\item{lkern}{character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian".
 The default "Triangle" is equivalent to using an Epanechnikov kernel, "Quadratic" and  "Cubic" refer to a Bi-weight and Tri-weight
  kernel, see Fan and Gijbels (1996). "Gaussian" is a truncated (compact support) Gaussian kernel.
  This is included for comparisons only and should be avoided due to its large computational costs.
}
\item{weighted}{
logical: use inverse variances as weights.
}
\item{sigma2}{
  \code{sigma2} allows to specify the variance of data entries.
}
\item{mask}{
optional logical mask, same dimensionality as \code{y}
}
\item{h}{
 bandwidth to use. In case of \code{lkern="Gaussian"} this is in FWHM (full width half maximum) units.
 Value refers to first voxel dimension.
}
\item{wghts}{
   voxel dimensions. Defaults to \code{c(1,1,1)}
}
\item{lambda}{
kritical scale parameter in hypothesis testing (adaptive weights smoothing)
}
\item{hmax}{
maximum bandwidth for adaptive weights smoothing
}
\item{res}{
array of residuals with dimension \code{c(nres,sum(mask))}.
}
\item{skern}{ \code{skern} specifies the kernel for the statistical
penalty. Defaults to "Plateau", the alternatives are "Triangle" and "Exp".
"Plateau" specifies a kernel that is equal to 1 in the interval (0,.3),
decays linearly in (.3,1) and is 0 for arguments larger than 1.
\code{lkern="Plateau"} and \code{lkern="Triangle"} allow for much faster computation (saves up
to 50\% CPU-time). \code{lkern="Plateau"} produces a less random weighting scheme.
}
\item{u}{
For test purposes in simulations: noisless 3D data.
}
\item{h0}{
Vector of 3 bandwidths corresponding to a Gaussian kernel that would produce a
comparable spatial correlation by convoluting iid data.
}
\item{testprop}{
logical: test the validity of a propagation condition for the specified value of \code{lambda}.
}
}
\value{
Functions \code{smooth3D} and \code{medianFilter3D} return a 3D array.
Functions \code{awsmask*} return a list with smoothed values of \code{y} in
component \code{theta} and smoothed residuals in component \code{res}.
}

\author{ Joerg Polzehl \email{polzehl@wias-berlin.de},
  Karsten Tabelow \email{tabelow@wias-berlin.de}  }
}
\note{
Functions \code{awsmask*} are used intenally in package \code{fmri}. They refer to the
situation, typical for fMRI, where the data are spatially correlated and this correlation
can be accessed using residuals with respect to a model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% use one of  RShowDoc("KEYWORDS")
\keyword{ misc }% __ONLY ONE__ keyword per line
