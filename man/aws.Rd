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
aws(y,hmax=NULL,qlambda=NULL,qtau=NULL,family="Gaussian",
                sigma2=NULL,scorr=0,shape=NULL,wghts=NULL,graph=FALSE,demo=FALSE,
		lkern="Triangle",aggkern="Uniform",
		spmin=0,homogen=TRUE,lseq=NULL,u=NULL,testprop=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed response data. \code{dim(y)} determines the dimensionality and extend of the grid design.}
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=250, 12, 5} for \code{dd=1, 2, 3}, respectively.}
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
   are "Bernoulli", "Poisson", "Exponential", "Volatility" and "Variance". \code{family="Volatility"} specifies a Gaussian distribution with
   expectation 0 and unknown variance. \code{family="Volatility"} specifies that \code{p*y/theta} is distributed as \code{\Chi^2} with \code{p=shape}
   degrees of freedom.}
  \item{sigma2}{\code{sigma2} allows to specify the variance in case of \code{family="Gaussian"}. Not used if \code{family!="Gaussian"}.
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==length(y)} then \code{sigma2}
   is assumed to contain the pointwise variance of \code{y} and a heteroscedastic variance model is used.}
  \item{scorr}{The vector \code{scorr} allows to specify a first order correlations of the noise for each coordinate direction,
    defaults to 0 (no correlation). }
  \item{shape}{Allows to specify an additional shape parameter for certain family models. Currently only used for family="Variance", that is \code{\Chi}-Square distributed observations
  with \code{shape} degrees of freedom. }
  \item{wghts}{ \code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space. }
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
  \item{demo}{ If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
  \item{aggkern}{ \code{aggkern} specifies the kernel \eqn{K_{agg}} used in the adaptive control step if \code{qtau<1} Defaults to "Uniform", the alternative choice is
   "Triangle"}
  \item{spmin}{Specifies the size of a plateau of the stochastic kernel.}
  \item{homogen}{If \code{homogen==TRUE} the algorithm tries to identify, in each point and each iteration, a sphere of maximal radius, such that
     all statistical penalties are less then \code{spmin} (the value specifying the plateau of the stochastic kernel). The statistical penalty
     is set to 0 for all points \code{j} within the sphere centered in point \code{i} for the following iteration.  }
  \item{lseq}{ \code{lseq} allows to increase the value of the scaling parameter \code{\lambda} for the first \code{length(lseq)} iteration steps by
  the factor specified in \code{lseq}. Defaults to NULL. In this case a default sequence is used that fulfils the propagation condition. }
  \item{u}{\code{u} allows to specify the true parameter. This is only used to test the algorithm and to select the smoothing parameters 
    \code{qlambda} and \code{qtau} by a propagation condition. If \code{u} is specified MSE and MAE of the estimates are 
    printed for each iteration step.
  \item{testprop}{If set this provides diagnostics for testing the propagation condition. The values of \code{y} should correspond to the specified
   family and a global model. }
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
   More formally we require the parameter \code{qlambda} and eventually \code{lseq}
   to be specified such that
   \code{\bf{E} |\hat{\theta}^k - \theta| \le (1+\alpha) \bf{E} |\tilde{\theta}^k - \theta|}
   where \code{\hat{\theta}^k} is the aws-estimate in step \code{k} and \code{\tilde{\theta}^k}
   is corresponding nonadaptive estimate using the same bandwidth (\code{qlambda=1}).
   Default values are selected to fulfil this condition for \code{\alpha=0.1}.
   
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
  \item{hmax }{Bandwidth used in the last iteration.}
  \item{mae }{If \code{u} is specified this parameter contains the MAE observed for all iteration steps. This information
  is useful to see if a parameter specification fulfils the propagation condition, see example.}
  \item{lseq }{Faktor used to scale \code{\lambda} for each iteration step.}
  \item{call }{The arguments of the function call.}
}
\references{ 
\item{ }{Polzehl, J. and Spokoiny, V. (2006). \emph{Propagation-Separation Approach for Local Likelihood Estimation}, 
Probab. Theory & Relat. Fields, in print.}
\item{ }{Belomestny, D. and Spokoiny, V. (2004a). \emph{Local likelihood modeling via stagewise aggregation}, 
WIAS-Preprint 1000.}
\item{ }{Polzehl, J. and Spokoiny, V. (2004b). \emph{Spatially adaptive regression estimation: Propagation-separation approach}, WIAS-Preprint 998.}
\item{ }{ Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp. 335-354}}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html}}
##\note{}
\seealso{See Also  \code{\link{lpaws}}}
\examples{
require(aws)
#
#    test of propagation condition for "Gaussian" case
#
    set.seed(1); y<-rnorm(10000)
    ttt<-aws(y,hmax=100,u=0,qlambda=1,qtau=1) # non-adaptive estimates
    ttt.975 <- aws(y,hmax=100,u=0,qlambda=.975,qtau=1) # aws
    ttt.975$mae/ttt$mae-1
#   shows that our propagation condition is fulfilled with alpha=.1
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
