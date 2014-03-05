\name{lpaws}
\alias{lpaws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Local polynomial smoothing by AWS}
\description{
The function allows for structural adaptive smoothing using a local 
polynomial (degree <=2) structural assumption. Response variables are assumed to be observed on a 1 or 2 dimensional regular grid. 
}
\usage{
lpaws(y, degree = 1, hmax = NULL, aws = TRUE, memory = FALSE, lkern = "Triangle", 
      homogen = TRUE, earlystop = TRUE, aggkern = "Uniform", sigma2 = NULL, 
      hw = NULL, ladjust = 1, u = NULL, graph = FALSE, demo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Response, either a vector (1D) or matrix (2D). The corresponding design is assumed to be a regular grid in 1D or 2D, respectively. }
  \item{degree}{Polynomial degree of the local model }
  \item{hmax}{ maximal bandwidth }
  \item{aws}{ logical: if TRUE structural adaptation (AWS) is used. }
  \item{memory}{ logical: if TRUE stagewise aggregation is used as an additional
   adaptation scheme. }
  \item{lkern}{ character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian"}
  \item{homogen}{ logical: if TRUE the function tries to determine regions 
   where weights can be fixed to 1. This may increase speed.}
  \item{earlystop}{ logical: if TRUE the function tries to determine points 
   where the homogeneous region is unlikely to change in further steps. This may increase speed.}
  \item{aggkern}{ character: kernel used in stagewise aggregation, either "Triangle" or "Uniform" }
  \item{sigma2}{ Error variance, the value is estimated if not provided. }
  \item{hw}{ Regularisation bandwidth, used to prevent from unidentifiability of
   local estimates for small bandwidths. }
  \item{ladjust}{ factor to increase the default value of lambda }
  \item{u}{ a "true" value of the regression function, may be provided to 
   report risks at each iteration. This can be used to test the propagation condition with \code{u=0} }
  \item{graph}{ logical: If TRUE intermediate results are illustrated graphically. May significantly slow down the computations in 2D. Please
   avoid using the default \code{X11()} on systems build with \code{cairo}, use
   \code{X11(type="Xlib")} instead (faster by a factor of 30). }
  \item{demo}{ logical: if TRUE wait after each iteration }
}
\value{ returns anobject of class \code{aws} with slots
                        \item{y = "numeric"}{y}
                        \item{dy = "numeric"}{dim(y)}
                        \item{x = "numeric"}{numeric(0)}
                        \item{ni = "integer"}{integer(0)}
                        \item{mask = "logical"}{logical(0)}
                        \item{theta = "numeric"}{Estimates of regression function and derivatives, \code{length: length(y)*(degree+1)}}
                        \item{mae = "numeric"}{Mean absolute error for each iteration step if u was specified,  numeric(0) else}
                        \item{var = "numeric"}{approx. variance of the estimates of the regression function. Please note that this does not reflect variability due to randomness of weights.}
                        \item{xmin = "numeric"}{numeric(0)}
                        \item{xmax = "numeric"}{numeric(0)}
                        \item{wghts = "numeric"}{numeric(0)}
                        \item{degree = "integer"}{degree}
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
                        \item{homogen = "logical"}{homogen}
                        \item{earlystop = "logical"}{eralustop}
                        \item{varmodel = "character"}{"Constant"}
                        \item{vcoef = "numeric"}{numeric(0)}
                        \item{call = "function"}{the arguments of the call to \code{lpaws}}
}
\references{Joerg Polzehl, Vladimir Spokoiny, in V. Chen, C.; Haerdle, W. and Unwin, A. (ed.) Handbook of Data Visualization Structural adaptive smoothing by propagation-separation methods Springer-Verlag, 2008, 471-492
}
\author{Joerg Polzehl \email{polzehl@wias-berlin.de} }
\note{If you specify \code{graph=TRUE} for 2D problems 
   avoid using the default \code{X11()} on systems build with \code{cairo}, use
   \code{X11(type="Xlib")} instead (faster by a factor of 30). 
}
\seealso{\code{link{awsdata}},\code{\link{aws}}, \code{\link{aws.irreg}}}
\examples{
library(aws)
# 1D local polynomial smoothing
\dontrun{demo(lpaws_ex1)}
# 2D local polynomial smoothing
\dontrun{demo(lpaws_ex2)}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
