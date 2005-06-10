\name{caws}
\alias{caws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Propagation-Separation approach - likelihood models}
\description{This function implements the propagation-separation
approach to adaptive estimation for exponential models with local 
constant expectation (likelihood approach). 
Function \code{psaws} provides a simpler interface to this function.
}
\usage{
caws(y, x = NULL, qlambda = NULL, lkern = "Triangle", family = "Gaussian", sigma2 = NULL, shape = NULL, hinit = NULL, hincr = NULL, hmax = NULL, heta = NULL, tau1 = NULL, tau2 = NULL, eta0 = NULL, NN = FALSE, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed values at location \code{x}. 
           In case of \code{x==NULL} (second parameter) \code{y} is assumed to be 
           observed on a one, two or three-dimensional grid. The dimension of 
           \code{y} determines if one, two or three-dimensional AWS is used.}
  \item{x}{ \code{x} is either \code{NULL}, in this case \code{y} is assumed 
           to be observed on a grid, or is a matrix, with rows corresponding to variables, 
           containing the design points where \code{y} is observed. }
  \item{qlambda}{ \code{qlambda} determines the scale parameter \code{lambda}
           for the stochastic penalty. The scaling parameter in the stochastic 
           penalty \code{lambda} is choosen as the \code{qlambda}-quantile
           of a Chi-square-distribution with number of parameters in the polynomial 
           model as degrees of freedom. If \code{qlambda==NULL} a standard value 
           depending on \code{model} and \code{symmetric} is choosen. 
	   \code{qlambda} == 1 eliminates the stochastic penalty resulting in a 
	   stagewise aggregation procedure.}
  \item{lkern}{ \code{lkern} determines the location kernel to be used. Options 
           are \code{"Uniform"}, \code{"Triangle"}, \code{"Quadratic"}, 
           \code{"Cubic"} and \code{"Exponential"}. Default is \code{"Triangle"}. 
           The Kernel operates on the squared distance, so \code{"Triangle"}
           corresponds to the use of an Epanechnikov kernel in kernel smoothing. 
           \code{"Exponential"} requires larger values of \code{hmax} and 
           therefore more iterations to reach comparable results. }
  \item{family}{ \code{family} determines the distribution type of \code{y}. 
           Currently implemented exponential families are \code{"Poisson"}, \code{"Bernoulli"},
           \code{"Gaussian"}, \code{"Exponential"}, \code{"Weibull"}, \code{"Volatility"}
           (Estimation of the scale parameter of a Gaussian distribution).
           Defaults to \code{"Poisson"}.}
  \item{sigma2}{\code{sigma2} is used to specify the variance with \code{family=="Gaussian".
                Variance will be estimated if \code{sigma2=NULL}. }
  \item{shape}{used for additional parameters of the specified distribution if needed.
               Exception: use \code{sigma2} for variance with \code{family=="Gaussian"}
  \item{hinit}{ \code{hinit} Initial bandwidth for the location penalty. 
               Appropriate value is choosen in case of \code{hinit==NULL}  }
  \item{hincr}{  \code{hincr} \code{hincr^(1/d)}, with \code{d} the 
           dimensionality of the design, is used as a factor to increase the 
           bandwidth between iterations. Defauts to \code{hincr==1.25}}
  \item{hmax}{ \code{hmax} Maximal bandwidth to be used. Determines the 
           number of iterations and is used as the stopping rule. }
  \item{heta}{ \code{heta} Minimal bandwidth for using memory control. Can be used
              to switch memory control off.}
  \item{tau1}{ \code{tau1} and \code{tau2} determine the scaling parameter for 
                 the penalty in memory control. Default values are provided that 
		 fulfil the propagation condition for default setting of \code{qlambda}.}
  \item{tau2}{ \code{tau2} and \code{tau2}determine the scaling parameter for 
                 the penalty in memory control. Default values are provided that 
		 fulfil the propagation condition for default setting of \code{qlambda}.}
  \item{eta0}{ \code{eta0} minimal memory value. }
  \item{NN}{ If \code{NN==TRUE} use nearest neighbor-rules instead of 
           distances in the location term. }
  \item{u}{ \code{u} used to supply values of the true regression function 
           for test purposes to calculate  Mean Squared Error (MSE) and 
           Mean Absolute Error (MAE) }
  \item{graph}{ \code{graph} if \code{TRUE} results are displayed after each 
           iteration step. }
  \item{demo}{ \code{demo} if \code{TRUE} after each iteration step results 
           are displayed and the process waits for user interaction. }
  \item{wghts}{ Specifies wghts for distance evaluation on a bi- or trivariate grid.
                Allows for anisotropic local neighborhoods. If \code{wghts=NULL}
                isotropic neighborhoods are used.}
}
\details{
  This function implements the propagation-separation approach 
  to adaptive estimation for exponential models with local 
  constant expectation (likelihood approach), see Polzehl and Spokoiny (2004a).
  The procedure is very similar to the adaptive weights approach 
  in Polzehl and Spokoiny (2002). The main difference lies in an additional
  control step that allows for "optimal" rates of estimation in case of 
  piecewise smooth functions.

  Function \code{psaws} provides a simpler interface to this function.
  
  The following exponential families (specified by parameter \code{family}) are
  currently implemented:
   \describe{
   \item{Binary response}{\code{family="Bernoulli"}}
   \item{Poisson regression}{\code{family="Poisson"} This family allows e.g. for 
   density estimation or for the analysis of poisson count data on a grid.}
   \item{Exponential regression}{\code{family="Exponential"} Applications of this family
   include e.g. test for constant (over time) failure rates and estimation of tail indicies.}
   \item{Gaussian regression}{\code{family="Gaussian"} This is the local constant homoskedastic
    regression model with additive sub-gaussian errors. }
   \item{Weibull regression}{\code{family="Weibull"} Applications in reliability analysis.}
   \item{Volatility model}{\code{family="Volatility"} Implements a local constant variance
   model for "Gaussian" data with expectation 0.}
   }
   
   The essential parameter in the procedure is \code{qlambda}. This parameter has an 
   interpretation as a significance level of a test for equivalence of two local
   parameter estimates. Optimal values mainly depend on the choosen \code{family}. 
   Default values provided are choosen to fulfil the propagation, i.e. in case of a 
   constant (global) parameter value and large \code{hmax} the procedure should 
   provide, with a high probability, the global (parametric) estimate.
   The optimal values only slightly depend on the model parameters, i.e. the
   default parameters should work in most situations. Larger values of \code{qlambda}
   may lead to oversmoothing, small values of \code{qlambda} lead to a random segmentation
   of homogeneous regions. 
   
   The numerical complexity of the procedure is mainly determined by \code{hmax}. The number
   of iterations is \code{d*log(hmax/hinit)/log(hincr)} with \code{d} being the dimension 
   of \code{y}. Comlexity in each iteration step is \code{Const*hakt*n} with \code{hakt}
   being the actual bandwith in the iteration step and \code{n} the number of design points.
   \code{hmax} determines the maximal possible variance reduction.
}
\value{ object of class \code{laws}.\code{family} with components
  \item{theta }{Parameter estimates}
  \item{ni}{sum of weights in each design point}
  \item{y}{values provided in \code{y}}
  \item{x}{values provided in \code{x}}
  \item{call}{actual function call}

}
\references{ ~put references to the literature/web site here ~ }
\item{ }{Polzehl, J. and Spokoiny, V. (2004b). \emph{Spatially adaptive regression estimation: Propagation-separation approach}, WIAS-Preprint 998.}
\item{ }{Polzehl, J. and Spokoiny, V. (2004a). \emph{Propagation-Separation Approach for Local Likelihood Estimation}, 
Manuscript, see webside.}
\item{ }{Polzehl, J. and Spokoiny, V. (2003). \emph{Varying coefficient 
regression modeling by adaptive weights smoothing}, WIAS-Preprint 818.}
\item{ }{Polzehl, J. and Spokoiny, V. (2002). \emph{Local likelihood 
modelling by adaptive weights smoothing}, WIAS-Preprint 787.}
\item{ }{ Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp. 335-354}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html}}

\seealso{ SEE ALSO \code{\link{aws}}, \code{\link{laws}}, \code{\link{psaws}} }

\examples{
#
#   Gaussian regression
#
u0<-
###
###    Artificial PET data
###
x <- 1:128
f12 <- function(x,y){
.5+2*((1.5*(x-64)^2+(y-64)^2<3025)) +
((x-64)^2+(y-104)^2<16)-1.5*((x-92)^2+(y-84)^2<25)+
((x-78)^2+(y-84)^2<25)-1.5*((x-64)^2+(y-84)^2<36)+
((x-50)^2+(y-84)^2<36)-1.5*((x-36)^2+(y-84)^2<25)+
((x-92)^2+(y-64)^2<25)-1.5*((x-78)^2+(y-64)^2<16)+
((x-64)^2+(y-64)^2<16)-1.5*((x-50)^2+(y-64)^2<25)+
((x-36)^2+(y-64)^2<25)-1.5*((x-92)^2+(y-44)^2<36)+
((x-78)^2+(y-44)^2<36)-1.5*((x-64)^2+(y-44)^2<25)+
((x-50)^2+(y-44)^2<25)-1.5*((x-36)^2+(y-44)^2<16)+
((x-64)^2+(y-24)^2<16)
}
u0 <- 5*outer(x,x,"f12")
set.seed(1)
y <- matrix(rpois(u0,u0),128,128)
# use larger hmax for good results
yhat <- caws(y,family="Poisson",hmax=20,graph=TRUE,u=u0)$theta
par(mfrow=c(1,3),mar=c(3,3,3,.5),mgp=c(2,1,0))
image(y,col=gray((0:255)/255))
title("Observed image")
image(yhat,col=gray((0:255)/255))
title("AWS-Reconstruction")
image(u0,col=gray((0:255)/255))
title("True image")
rm(u0,y,yhat,x)
###
###   VOLATITILTY ESTIMATION
###
###   artificial example
###
sigma <- c(rep(1,125),rep(2,125),rep(.5,125),rep(1,125))
set.seed(1)
x <- rnorm(sigma,0,sigma)
tmp <- caws(x,family="Volatility",u=sigma^2,graph=TRUE,hmax=1000)
par(mfrow=c(1,1),mar=c(3,3,3,1))
plot(abs(x),col=3,xlab="time t",ylab=expression(abs( R[t] )))
lines(tmp$theta,col=1,lwd=2)
lines(sigma,col=2,lwd=2,lty=3)
legend(350,5.5,c("PS-AWS","true sigma"),
                 lwd=c(2,2),lty=c(1,3),col=c(1,2))
title(expression(paste("Estimates of  ",sigma(t))))
rm(tmp,sigma,x)
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
