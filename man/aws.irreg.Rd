\name{aws.irreg}
\alias{aws.irreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ local constant AWS for irregular (1D/2D) design }
\description{
The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient Gaussian models on a 1D or 2D irregulat design. 
The function allows for a paramertic (polynomial) mean-variance dependence. 
}
\usage{
aws.irreg(y, x, hmax = NULL, aws=TRUE, memory=FALSE, varmodel = "Constant", lkern = "Triangle", aggkern = "Uniform", sigma2 = NULL, nbins = 100, hpre = NULL, henv = NULL, ladjust =1, varprop = 0.1, graph = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The observed response vector (length n) }
  \item{x}{ Design matrix, dimension n x d, \code{d \%in\% 1:2} }
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Unit is binwidth in the first dimension.}
  \item{aws}{ logical: if TRUE structural adaptation (AWS) is used. }
  \item{memory}{ logical: if TRUE stagewise aggregation is used as an additional
   adaptation scheme. }
  \item{varmodel}{ determines the model that relates variance to mean. Either
"Constant", "Linear" or "Quadratic". }
  \item{lkern}{ character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian"}
  \item{aggkern}{ character: kernel used in stagewise aggregation, either "Triangle" or "Uniform" }
  \item{sigma2}{\code{sigma2} allows to specify the variance in case of \code{varmodel="Constant"}, estimated if not given.  }
  \item{nbins}{ numer of bins, can be NULL, a positive integer or a vector of positive integers (length d) }
  \item{hpre}{ smoothing bandwidth for initial variance estimate }
  \item{henv}{ radius of balls around each observed design point where
       estimates will be calculated }
  \item{ladjust}{ factor to increase the default value of lambda }
  \item{varprop}{ exclude the largest 100*varprop\% squared residuals when estimating the error variance }
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
}
\details{
  Data are first binned (1D/2D), then aws is performed on all datapoints within 
  distance <= henv of nonempty bins.  
}
\value{
returns anobject of class \code{aws} with slots
                        \item{y = "numeric"}{y}
                        \item{dy = "numeric"}{dim(y)}
                        \item{x = "numeric"}{x}
                        \item{ni = "integer"}{number of observations per bin}
                        \item{mask = "logical"}{bins where parameters have been estimated}
                        \item{theta = "numeric"}{Estimates of regression function, \code{length: length(y)}}
                        \item{mae = "numeric"}{numeric(0)}
                        \item{var = "numeric"}{approx. variance of the estimates of the regression function. Please note that this does not reflect variability due to randomness of weights.}
                        \item{xmin = "numeric"}{vector of minimal x-values (bins)}
                        \item{xmax = "numeric"}{vector of maximal x-values (bins)}
                        \item{wghts = "numeric"}{relative binwidths}
                        \item{degree = "integer"}{0}
                        \item{hmax  = "numeric"}{effective hmax}
                        \item{sigma2  = "numeric"}{provided or estimated error variance}
                        \item{scorr = "numeric"}{0}
                        \item{family = "character"}{"Gaussian"}
                        \item{shape = "numeric"}{numeric(0)}
                        \item{lkern  = "integer"}{integer code for lkern,
                        1="Plateau", 2="Triangle", 3="Quadratic", 4="Cubic", 5="Gaussian"}
                        \item{lambda = "numeric"}{effective value of lambda}
                        \item{ladjust = "numeric"}{effective value of ladjust}
                        \item{aws = "logical"}{aws}
                        \item{memory = "logical"}{memory}
                        \item{homogen = "logical"}{FALSE}
                        \item{earlystop = "logical"}{FALSE}
                        \item{varmodel = "character"}{varmodel}
                        \item{vcoef = "numeric"}{estimated coefficients in variance model}
                        \item{call = "function"}{the arguments of the call to \code{aws}}
}
\references{ Joerg Polzehl, Vladimir Spokoiny, in V. Chen, C.; Haerdle, W. and Unwin, A. (ed.) Handbook of Data Visualization Structural adaptive smoothing by propagation-separation methods Springer-Verlag, 2008, 471-492
 }
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ See also  \code{\link{lpaws}}, \code{link{awsdata}}, \code{\link{lpaws}} }
\examples{
require(aws)
# 1D local constant smoothing
\dontrun{demo(irreg_ex1)}
# 2D local constant smoothing
\dontrun{demo(irreg_ex2)}
}
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
