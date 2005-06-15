\name{vaws}
\alias{vaws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AWS for multiple local constant models on a grid  }
\description{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient "Gaussian" models on a 1D, 2D or 3D grid. 
The observations in each gridpoint are assumed to consist of a vector of
independent components. The underlying structure, i.e. location of discontinuities,
is assumed to be the same in all components.
A homoskedastic 
or heteroskedastic model is used depending on the content of \code{sigma2}. In the latter case
all components are assumed to have the same structure in the variance.
}
\usage{
vaws(y, qlambda = NULL, qtau = NULL, lkern = "Triangle", aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, heta = NULL, eta0 = 0, u = NULL, graph = FALSE, wghts = NULL, vwghts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the vector-valued observed response data. The first dimension determines the number of components of \code{y} in each design point while \code{dim(y)[-1]} gives the dimensionality and extend of the grid design. }
  \item{qlambda}{\code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,dim(y)[1])}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  \code{qlambda} can be selected as the smallest value that fulfils a propagation condition.
  This means that in a parametric, in this case constant, model the resulting estimate is, for a large \code{hmax} almost parametric. 
  The default values of \code{qlambda} is selected to fulfil this condition. Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate. }
  \item{qtau}{Stagewise Aggregation, see Belomestny and Spokoiny (2004) is used as an adaptive control step if \code{qtau<1}. 
   \code{qtau} determines the scaling parameter \code{\tau} in the Stagewise Aggregation algorithm in the same way as \code{qlambda} does
   for AWS. Default values are again selected, for the probability distribution specified with the \code{family} argument, by a propagation condition.}
  \item{lkern}{\code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
  \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle"}
  \item{sigma2}{\code{sigma2} allows to specify the variance of the observations up to a specific constant factor for each 
      component. This factor is provided with \code{vwghts}.
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==prod(dim(y)[-1])} then \code{sigma2}
   is assumed to contain the pointwise variance of \code{y} and a heteroscedastic variance model is used.}
  \item{hinit}{\code{hinit} specifies the initial bandwidth. Defaults to \code{hinit=1} }
  \item{hincr}{\code{hincr} specifies the factor used to increase the size of local neigborhoods after each iteration. The bandwidth is increased by
   a factor \code{hincr^(1/dd)} with \code{dd} specifying the dimensionality of the grid. Defaults to \code{hincr=1.25} }
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=250, 12, 5} for \code{dd=1, 2, 3}, respectively.}
  \item{heta}{ \code{heta} specifies the minimal bandwidth to use with stagewise aggregation.  }
  \item{u}{\code{u} allows to specify the true parameter. This is only used to test the algorithm and to select the smoothing parameters 
    \code{qlambda} and \code{qtau} by a propagation condition. If \code{u} is specified MSE and MAE of the estimates are 
    printed for each iteration step. }
  \item{graph}{If  \code{graph=TRUE} intermediate results for the first component are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
  \item{wghts}{ \code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space.  }
  \item{vwghts}{ \code{vwghts} specifies a specific constant factor for the variance of each 
      component. defaults to \code{vwghts=rep(1,dim(y)[1])} }
}
\details{
  AWS for multiple local constant models on a grid  }
\description{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient "Gaussian" models on a 1D, 2D or 3D grid. 
The observations in each gridpoint are assumed to consist of a vector of
independent components. The underlying structure, i.e. location of discontinuities,
is assumed to be the same in all components.
A homoskedastic 
or heteroskedastic model is used depending on the content of \code{sigma2}. In the latter case
all components are assumed to have the same structure in the variance.

A typical application of this function occurs in case of registered 2D or 3D images.

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
   
\code{qlambda>=1} provides the stagewise aggregation procedure from Belomestny and Spokoiny (2004).
\code{qtau>=1} provides Adaptive weights smoothing without control by stagewise aggregation.  
   
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
~put references to the literature/web site here ~ }
\seealso{See Also  \code{\link{aws}}}
\examples{
require(aws)
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
     y <- rbind(rnorm(w0,w0,sigma),rnorm(w0,1-w0,sigma/sqrt(2)))
     u <- rbind(as.vector(w0),as.vector(1-w0))
     dim(y) <- dim(u) <- c(2,dim(w0))
#   increase hmax for better results  
#      vwghts is used to specify that the second component has lower variance 
     tmp <- vaws(y,hmax=5,graph=TRUE,u=u,vwghts=c(1,.5))
     par(mfrow=c(2,2))
     image(y[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
     image(y[2,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
     image(tmp$theta[1,,],zlim=range(y[1,,]),col=gray((0:255)/255),xaxt="n",yaxt="n")
     image(tmp$theta[2,,],zlim=range(y[2,,]),col=gray((0:255)/255),xaxt="n",yaxt="n")
     rm(y,w0,tmp,sigma)

}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
