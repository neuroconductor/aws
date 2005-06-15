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
aws(y, qlambda = NULL, qtau = NULL, family = "Gaussian", lkern = "Triangle", aggkern = "Uniform", sigma2 = NULL, shape = NULL, hinit = NULL, hincr = NULL, hmax = NULL, heta = NULL, eta0 = NULL, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed response data. \code{dim(y)} determines the dimensionality and extend of the grid design.}
  \item{qlambda}{ \code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,1)}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  \code{qlambda} can be selected for a given model, i.e. specification of \code{family}, as the smallest value that fulfils a propagation condition.
  This means that in a parametric, in this case constant, model the resulting estimate is, for a large \code{hmax} almost parametric. 
  The default values of \code{qlambda} for the different specifications of \code{family} are selected to fulfil this condition. Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate.}
  \item{qtau}{Stagewise Aggregation, see Belomestny and Spokoiny (2004) is used as an adaptive control step if \code{qtau<1}. 
   \code{qtau} determines the scaling parameter \code{\tau} in the Stagewise Aggregation algorithm in the same way as \code{qlambda} does
   for AWS. Default values are again selected, for the probability distribution specified with the \code{family} argument, by a propagation condition.}
  \item{family}{\code{family} specifies the probability distribution. Default is \code{family="Gaussian"}, also implemented
   are "Bernoulli", "Poisson", "Exponential", "Weibull" and "Volatility". \code{family="Volatility"} specifies a Gaussian distribution with
   expectation 0 and unknown variance}
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
   \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle"}
  \item{sigma2}{\code{sigma2} allows to specify the variance in case of \code{family="Gaussian"}. Not used if \code{family!="Gaussian"}.
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==length(y)} then \code{sigma2}
   is assumed to contain the pointwise variance of \code{y} and a heteroscedastic variance model is used.}
  \item{shape}{\code{shape} allows to specify additional parameters if \code{family!="Gaussian"}, e.g. for \code{family="Weibull"} }
  \item{hinit}{\code{hinit} specifies the initial bandwidth. Defaults to \code{hinit=1} }
  \item{hincr}{\code{hincr} specifies the factor used to increase the size of local neigborhoods after each iteration. The bandwidth is increased by
   a factor \code{hincr^(1/dd)} with \code{dd} specifying the dimensionality of the grid. Defaults to \code{hincr=1.25}  }
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=250, 12, 5} for \code{dd=1, 2, 3}, respectively.}
  \item{heta}{ \code{heta} specifies the minimal bandwidth to use with stagewise aggregation.  }
  \item{eta0}{ \code{eta0} minimal mixing coefficient in stagewise aggregation, defaults to 0 in case of \code{family="Gaussian"} and values
   larger than 0 otherwise to avoid problems at the boundary of the parameter space.}
  \item{u}{\code{u} allows to specify the true parameter. This is only used to test the algorithm and to select the smoothing parameters 
    \code{qlambda} and \code{qtau} by a propagation condition. If \code{u} is specified MSE and MAE of the estimates are 
    printed for each iteration step.}
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
  \item{demo}{ If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{wghts}{ \code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space. }
}
\details{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient likelihood models  on a 1D, 2D or 3D grid. For "Gaussian"
models, i.e. regression with additive "Gaussian" errors, a homoskedastic 
or heteroskedastic model is used depending on the content of \code{sigma2}.
\code{qlambda>=1} provides the stagewise aggregation procedure from Belomestny and Spokoiny (2004).
\code{qtau>=1} provides Adaptive weights smoothing without control by stagewise aggregation.  

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
\value{
  \item{theta }{Contains the parameter estimates \code{thetahat} on the grid. Dimension coincides with \code{dim(y)}.}
  \item{ni }{Contains the sum of weights, i.e. \code{trace(W_i)}, in all grid points \code{i}.}
  \item{var }{Contains an estimate of \code{var(thetahat)}, in all grid points \code{i}, conditional on the weights in the last
  iteration step being independent form \code{y}. In this sense it only is an approximation of the true variance.}
  \item{y }{The observed data.}
  \item{call }{The arguments of the function call.}
}
\references{ 
\item{ }{Polzehl, J. and Spokoiny, V. (2004a). \emph{Propagation-Separation Approach for Local Likelihood Estimation}, 
WIAS-Preprint 1000.}
\item{ }{Belomestny, D. and Spokoiny, V. (2004a). \emph{Local likelihood modeling via stagewise aggregation}, 
Manuscript, see webside.}
\item{ }{Polzehl, J. and Spokoiny, V. (2004b). \emph{Spatially adaptive regression estimation: Propagation-separation approach}, WIAS-Preprint 998.}
\item{ }{Polzehl, J. and Spokoiny, V. (2003). \emph{Varying coefficient 
regression modeling by adaptive weights smoothing}, WIAS-Preprint 818.}
\item{ }{Polzehl, J. and Spokoiny, V. (2002). \emph{Local likelihood 
modelling by adaptive weights smoothing}, WIAS-Preprint 787.}
\item{ }{ Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp. 335-354}}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html}}
##\note{}
\seealso{See Also  \code{\link{paws}}, \code{\link{vaws}}}
\examples{
require(aws)
#  
#  Local constant image from Polzehl and Spokoiny (2000)
#
     xy <- rbind(rep(0:255,256),rep(0:255,rep(256,256)))
     indw <- c(1:12,29:48,73:100,133:168,209:256)
     w0 <- matrix(rep(FALSE,256*256),ncol=256)
     w0[indw,] <- TRUE
     w0[,indw] <- !w0[,indw]
     w0 <- w0-.5
     
     w0[((xy[1,]-129)^2+(xy[2,]-129)^2)<=10000&((xy[1,]-129)^2+(xy[2,]-129)^2)>=4900] <- 0
     w0[abs(xy[1,]-xy[2,])<=20&((xy[1,]-129)^2+(xy[2,]-129)^2)<4900] <- 0
     w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625] <- 0
     
     w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625&xy[2,]>27&xy[2,]<31] <- -.5
     
     w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625&xy[1,]>223&xy[1,]<227] <- .5
     w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625] <- 0
     
     w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625&xy[1,]>27&xy[1,]<31] <- -.5
     
     w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625&xy[2,]>223&xy[2,]<227] <- .5
     w0[((xy[2,]-225)^2+(xy[1,]-225)^2)+1*(xy[2,]-225)*(xy[1,]-225)<=400] <- 0
     w0[((xy[2,]-30)^2+(xy[1,]-30)^2)<=256] <- 0
     rm(xy,indw)
     set.seed(1)
     sigma <- .5
     y <- w0+rnorm(w0,0,sigma)
#   increase hmax for better results
     tmp <- aws(y,hmax=5,graph=TRUE)
     par(mfrow=c(1,3))
     image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
     image(tmp$theta,zlim=range(y),col=gray((0:255)/255),xaxt="n",yaxt="n")
     image(w0,zlim=range(y),col=gray((0:255)/255),xaxt="n",yaxt="n")
     rm(y,w0,tmp,sigma)
###
###    Artificial Poisson data
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
yhat <- aws(y,family="Poisson",hmax=10,graph=TRUE,u=u0,qtau=1)$theta
par(mfrow=c(1,3),mar=c(3,3,3,.5),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed image")
image(yhat,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("AWS-Reconstruction")
image(u0,col=gray((0:255)/255),xaxt="n",yaxt="n")
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
tmp <- aws(x,family="Volatility",u=sigma^2,qtau=1,graph=TRUE,hmax=1000)
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
