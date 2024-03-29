\name{aws}
\alias{aws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AWS for local constant models on a grid}
\description{The function implements the propagation separation approach to
nonparametric smoothing (formerly introduced as Adaptive weights smoothing)
for varying coefficient likelihood models on a 1D, 2D or 3D grid. For "Gaussian"
models, i.e. regression with additive "Gaussian" errors, a homoskedastic
or heteroskedastic model is used depending on the content of \code{sigma2}
}
\usage{
aws(y,hmax=NULL, mask=NULL, aws=TRUE, memory=FALSE, family="Gaussian",
                lkern="Triangle", aggkern="Uniform",
                sigma2=NULL, shape=NULL, scorr=0, spmin=0.25,
		            ladjust=1,wghts=NULL,u=NULL,graph=FALSE,demo=FALSE,
                testprop=FALSE,maxni=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{array \code{y} containing the observe response (image intensity) data. \code{dim(y)} determines the dimensionality and extend of the grid design.}
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=250, 12, 5}
  for 1D, 2D, 3D images, respectively.
  In case of \code{lkern="Gaussian"} the bandwidth is assumed to be given in full width half maximum (FWHM) units, i.e., \code{0.42466} times gridsize.   }
  \item{aws}{ logical: if TRUE structural adaptation (AWS) is used. }
	\item{mask}{ optional logical mask, same dimensionality as \code{y}}
  \item{memory}{ logical: if TRUE stagewise aggregation is used as an additional
   adaptation scheme. }
  \item{family}{\code{family} specifies the probability distribution. Default is \code{family="Gaussian"}, also implemented
   are "Bernoulli", "Poisson", "Exponential", "Volatility", "Variance" and "NCchi". \code{family="Volatility"} specifies a Gaussian distribution with
   expectation 0 and unknown variance. \code{family="Volatility"} specifies that \code{p*y/theta} is distributed as \eqn{\chi^2} with \code{p=shape}
   degrees of freedom. \code{family="NCchi"} uses a noncentral Chi distribution with \code{p=shape} degrees of freedom and noncentrality parameter \code{theta}}
  \item{lkern}{ character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian".
   The default "Triangle" is equivalent to using an Epanechnikov kernel, "Quadratic" and  "Cubic" refer to a Bi-weight and Tri-weight
    kernel, see Fan and Gijbels (1996). "Gaussian" is a truncated (compact support) Gaussian kernel.
    This is included for comparisons only and should be avoided due to its large computational costs.}
  \item{aggkern}{ character: kernel used in stagewise aggregation, either "Triangle" or "Uniform" }
  \item{sigma2}{\code{sigma2} allows to specify the variance in case of \code{family="Gaussian"}. Not used if \code{family!="Gaussian"}.
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==length(y)} then \code{sigma2}
   is assumed to contain the pointwise variance of \code{y} and a heteroscedastic variance model is used.}
  \item{shape}{Allows to specify an additional shape parameter for certain family models. Currently only used for family="Variance", that is \eqn{\chi}-Square distributed observations
  with \code{shape} degrees of freedom. }
  \item{scorr}{The vector \code{scorr} allows to specify a first order correlations of the noise for each coordinate direction,
    defaults to 0 (no correlation). }
  \item{spmin}{Determines the form (size of the plateau) in the adaptation kernel.
    Not to be changed by the user.}
  \item{ladjust}{ factor to increase the default value of lambda }
  \item{wghts}{ \code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space. }
  \item{u}{ a "true" value of the regression function, may be provided to
   report risks at each iteration. This can be used to test the propagation condition with \code{u=0} }
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
  \item{demo}{ If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{testprop}{If set this provides diagnostics for testing the propagation condition. The values of \code{y} should correspond to the specified
   family and a global model. }
  \item{maxni}{If TRUE use \eqn{max_{l<=k}(N_i^{(l)}} instead of \eqn{(N_i^{(k)}} in the definition of the statistical penalty.}
}
\details{The function implements the propagation separation approach to
nonparametric smoothing (formerly introduced as Adaptive weights smoothing)
for varying coefficient likelihood models  on a 1D, 2D or 3D grid. For "Gaussian"
models, i.e. regression with additive "Gaussian" errors, a homoskedastic
or heteroskedastic model is used depending on the content of \code{sigma2}.
\code{aws==FALSE} provides the stagewise aggregation procedure from Belomestny and Spokoiny (2004).
\code{memory==FALSE} provides Adaptive weights smoothing without control by stagewise aggregation.

The essential parameter in the procedure is a critical value \code{lambda}. This parameter has an
   interpretation as a significance level of a test for equivalence of two local
   parameter estimates. Optimal values mainly depend on the choosen \code{family}.
   Values set internally are choosen to fulfil a propagation condition, i.e. in case of a
   constant (global) parameter value and large \code{hmax} the procedure
   provides, with a high probability, the global (parametric) estimate.
   More formally we require the parameter \code{lambda}
   to be specified such that
   \eqn{\bf{E} |\hat{\theta}^k - \theta| \le (1+\alpha) \bf{E} |\tilde{\theta}^k - \theta|}
   where \eqn{\hat{\theta}^k} is the aws-estimate in step \code{k} and \eqn{\tilde{\theta}^k}
   is corresponding nonadaptive estimate using the same bandwidth (\code{lambda=Inf}).
   The value of lambda can be adjusted by specifying the factor \code{ladjust}. Values \code{ladjust>1} lead to an less effective adaptation while \code{ladjust<<1} may lead to random segmentation
   of, with respect to a constant model, homogeneous regions.

   The numerical complexity of the procedure is mainly determined by \code{hmax}. The number
   of iterations is approximately \code{Const*d*log(hmax)/log(1.25)} with \code{d} being the dimension
   of \code{y} and the constant depending on the kernel \code{lkern}. Comlexity in each iteration step is \code{Const*hakt*n} with \code{hakt}
   being the actual bandwith in the iteration step and \code{n} the number of design points.
   \code{hmax} determines the maximal possible variance reduction.

}
\value{
returns anobject of class \code{aws} with slots
                        \item{y = "numeric"}{y}
                        \item{dy = "numeric"}{dim(y)}
                        \item{x = "numeric"}{numeric(0)}
                        \item{ni = "integer"}{integer(0)}
                        \item{mask = "logical"}{logical(0)}
                        \item{theta = "numeric"}{Estimates of regression function, \code{length: length(y)}}
                        \item{mae = "numeric"}{Mean absolute error for each iteration step if u was specified,  numeric(0) else}
                        \item{var = "numeric"}{approx. variance of the estimates of the regression function. Please note that this does not reflect variability due to randomness of weights.}
                        \item{xmin = "numeric"}{numeric(0)}
                        \item{xmax = "numeric"}{numeric(0)}
                        \item{wghts = "numeric"}{numeric(0), ratio of distances \code{wghts[-1]/wghts[1]}}
                        \item{degree = "integer"}{0}
                        \item{hmax  = "numeric"}{effective hmax}
                        \item{sigma2  = "numeric"}{provided or estimated error variance}
                        \item{scorr = "numeric"}{scorr}
                        \item{family = "character"}{family}
                        \item{shape = "numeric"}{shape}
                        \item{lkern  = "integer"}{integer code for lkern,
                        1="Plateau", 2="Triangle", 3="Quadratic", 4="Cubic", 5="Gaussian"}
                        \item{lambda = "numeric"}{effective value of lambda}
                        \item{ladjust = "numeric"}{effective value of ladjust}
                        \item{aws = "logical"}{aws}
                        \item{memory = "logical"}{memory}
                        \item{homogen = "logical"}{homogen}
                        \item{earlystop = "logical"}{FALSE}
                        \item{varmodel = "character"}{"Constant"}
                        \item{vcoef = "numeric"}{numeric(0)}
                        \item{call = "function"}{the arguments of the call to \code{aws}}
}
\references{J. Polzehl, K. Tabelow (2019). Magnetic Resonance Brain Imaging:
             Modeling and Data Analysis Using R. Springer, Use R! series. Appendix A.
             Doi:10.1007/978-3-030-29184-6.

             J. Polzehl, K. Papafitsoros, K. Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
             Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06.

             J. Polzehl, V. Spokoiny,  Adaptive Weights Smoothing with applications to image restoration,
						 J. R. Stat. Soc. Ser. B Stat. Methodol.  62 ,  (2000)  , pp. 335--354. DOI:10.1111/1467-9868.00235.

             J. Polzehl, V. Spokoiny,  Propagation-separation approach for local likelihood estimation,
						 Probab. Theory Related Fields  135 (3),  (2006)  , pp. 335--362. DOI:10.1007/s00440-005-0464-1.
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{http://www.wias-berlin.de/people/polzehl/}}
\note{
use \code{setCores='number of threads'} to enable parallel execution.
}

\seealso{See also \code{\link{paws}}, \code{\link{lpaws}}, \code{\link{vaws}},\code{link{awsdata}}, \code{\link{aws.irreg}}, \code{\link{aws.gaussian}}}
\examples{
require(aws)
# 1D local constant smoothing
\dontrun{demo(aws_ex1)}
\dontrun{demo(aws_ex2)}
# 2D local constant smoothing
\dontrun{demo(aws_ex3)}
}
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
