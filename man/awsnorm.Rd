\name{awsnorm}
\alias{awsnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AWS for local constant heteroskedastic regression}
\description{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient Gaussian models with unknown mean and variance on a 1D, 2D or 3D grid.}
\usage{
awsnorm(y, qlambda = NULL, qtau = NULL, lkern = "Triangle", aggkern = "Uniform", hinit = NULL, hincr = NULL, hmax = NULL, heta = NULL, eta0 = NULL, u = NULL, graph = FALSE, wghts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{\code{y} contains the observed response data. \code{dim(y)} determines the dimensionality and extend of the grid design.}
  \item{qlambda}{ \code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,1)}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  \code{qlambda} can be selected  as the smallest value that fulfils a propagation condition.
  This means that in a parametric, in this case constant in mean and variance, model the resulting estimate is, for a large \code{hmax} almost parametric. 
  The default value of \code{qlambda} is selected to fulfil this condition. Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate.}
  \item{qtau}{Stagewise Aggregation, see Belomestny and Spokoiny (2004) is used as an adaptive control step if \code{qtau<1}. 
   \code{qtau} determines the scaling parameter \code{\tau} in the Stagewise Aggregation algorithm in the same way as \code{qlambda} does
   for AWS. Default values are again selected by a propagation condition.}
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
   \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle"}
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
\details{
The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient Gaussian models with unknown mean and variance on a 1D, 2D or 3D grid. 
\code{qtau>=1} provides Adaptive weights smoothing without control by stagewise aggregation.  

The essential parameter in the procedure is \code{qlambda}. This parameter has an 
   interpretation as a significance level of a test for equivalence of two local
   parameter estimates. 
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
  \item{mu }{Contains the estimated mean}
  \item{sigma }{Contains the estimated variance}
  \item{bi }{Contains the sum of weights, i.e. \code{trace(W_i)}, in all grid points \code{i}.}
  \item{args}{Contains the call of the function}
}
\references{ }
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
\seealso{See Also  \code{\link{aws}}, \code{\link{vawsnorm}}}

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, qlambda = NULL, qtau = NULL, lkern = "Triangle", 
    aggkern = "Uniform", hinit = NULL, hincr = NULL, hmax = NULL, 
    heta = NULL, eta0 = NULL, graph = FALSE, u = NULL, wghts = NULL, 
    eps = 5.1e-06) 
{
    KLnorm <- function(mu1, mu2, sig1, sig2) {
        log(sig1/sig2) - 1 + (sig2 + (mu1 - mu2)^2)/sig1
    }
    updtheta <- function(zobj, tobj, cpar, aggkern) {
        heta <- cpar$heta
        eta0 <- cpar$eta0
        tau1 <- cpar$tau1
        kstar <- cpar$kstar
        hakt <- zobj$hakt
        tau <- tau1 * (2 + max(kstar - log(hakt), 0))
        hakt <- zobj$hakt
        bi0 <- zobj$bi0
        bi <- zobj$bi
        n <- length(bi)
        munew <- zobj$ami/bi
        sigmanew <- zobj$asi/bi - munew * munew
        sigmanew <- sigmanew * bi/(bi - 1) + cpar$eps
        sigma <- tobj$sigma
        sigmanew[is.na(sigmanew)] <- sigma[is.na(sigmanew)]
        mu <- tobj$mu
        sigmanew[tobj$fix] <- sigma[tobj$fix]
        munew[tobj$fix] <- mu[tobj$fix]
        if (hakt > heta) {
            eta <- switch(aggkern, Uniform = (1 - eta0) * as.numeric(bi0/tau * 
                KLnorm((1 - eta0) * munew + eta0 * mu, mu, (1 - 
                  eta0) * sigmanew + eta0 * sigma, sigma) > 1) + 
                eta0, Triangle = (1 - eta0) * pmin(1, bi0/tau * 
                KLnorm((1 - eta0) * munew + eta0 * mu, mu, (1 - 
                  eta0) * sigmanew + eta0 * sigma, sigma)) + 
                eta0)
        }
        else {
            eta <- rep(eta0, n)
        }
        eta[tobj$fix] <- 1
        bi <- (1 - eta) * bi + eta * tobj$bi
        mu <- (1 - eta) * munew + eta * mu
        sigma <- (1 - eta) * sigmanew + eta * sigma
        list(mu = mu, sigma = sigma, bi = bi, eta = eta, fix = (eta == 
            1))
    }
    args <- match.call()
    spmax <- 5
    dy <- dim(y)
    yyt <- y * y
    if (is.null(qlambda)) 
        qlambda <- 0.96
    if (is.null(qtau)) 
        if (qlambda == 1) {
            if (is.null(dy)) 
                qtau <- 0.8
            else qtau <- 0.25
        }
        else qtau <- 0.92
    d1 <- 2
    if (qtau < 1) 
        tau1 <- qchisq(qtau, d1)
    else tau1 <- 1e+50
    if (aggkern == "Triangle") 
        tau1 <- 2.5 * tau1
    if (is.null(eta0)) 
        eta0 <- 0
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, 2)
    if (qlambda < 1) 
        lambda <- 2 * qchisq(qlambda, 2)
    else lambda <- 1e+50
    if (is.null(dy)) {
        form = "uni"
        n1 <- n <- length(y)
        n2 <- n3 <- 1
        ddim <- 1
        kstar <- log(100)
    }
    if (length(dy) == 2) {
        form = "bi"
        n1 <- dy[1]
        n2 <- dy[2]
        n3 <- 1
        n <- n1 * n2
        ddim <- 2
        kstar <- log(15)
        hincr <- sqrt(hincr)
    }
    if (length(dy) == 3) {
        form = "tri"
        n1 <- dy[1]
        n2 <- dy[2]
        n3 <- dy[3]
        n <- n1 * n2 * n3
        n1 <- dy[1]
        n2 <- dy[2]
        n3 <- dy[3]
        n <- n1 * n2 * n3
        ddim <- 3
        kstar <- log(5)
        hincr <- hincr^(1/3)
    }
    if (length(dy) > 3) 
        return("AWS for more than 3 dimensional grids is not implemented")
    if (is.null(wghts)) 
        wghts <- c(1, 1, 1)
    if (is.null(hinit) || hinit < 10^(1/ddim)) 
        hinit <- 10^(1/ddim)
    if (is.null(hincr) || hincr <= 1) 
        hincr <- 1.25
    if (is.null(heta)) 
        heta <- max(4, hinit + 1)
    if (is.null(hmax)) {
        if (is.null(dim(y))) 
            hmax <- 250
        if (length(dim(y)) == 2) 
            hmax <- 12
        if (length(dim(y)) == 3) 
            hmax <- 5
    }
    cpar <- list(heta = heta, tau1 = tau1, eta0 = eta0, eps = eps, 
        kstar = kstar)
    hinit <- hinit/wghts[1]
    hmax <- hmax/wghts[1]
    wghts <- (wghts[2:3]/wghts[1])
    tobj <- list(bi = rep(1, n), mu = y, sigma = yyt, fix = rep(FALSE, 
        n))
    zobj <- list(asi = yyt, ami = y, bi0 = rep(1, n))
    bi0old <- rep(1, n)
    hakt <- hinit
    lambda0 <- lambda
    if (hinit > 1) 
        lambda0 <- 1e+50
    while (hakt <= hmax) {
        zobj <- .Fortran("cawsnorm", as.double(y), as.double(yyt), 
            as.logical(tobj$fix), as.integer(n1), as.integer(n2), 
            as.integer(n3), hakt = as.double(hakt), as.double(lambda0), 
            as.double(tobj$mu), as.double(tobj$sigma), bi = as.double(tobj$bi), 
            bi0 = as.double(zobj$bi0), ami = as.double(zobj$ami), 
            asi = as.double(zobj$asi), as.integer(lkern), as.double(spmax), 
            as.double(wghts), PACKAGE = "aws")[c("bi", "bi0", 
            "ami", "asi", "hakt")]
        if (hakt > n1/2) 
            zobj$bi0 <- hincr * biold
        biold <- zobj$bi0
        dim(zobj$ami) <- dy
        dim(zobj$asi) <- dy
        dim(zobj$bi) <- dy
        dim(zobj$bi0) <- dy
        tobj <- updtheta(zobj, tobj, cpar, aggkern)
        dim(tobj$mu) <- dy
        dim(tobj$sigma) <- dy
        dim(tobj$bi) <- dy
        dim(tobj$eta) <- dy
        if (graph) {
            if (ddim == 1) {
                par(mfrow = c(1, 3), mar = c(1, 1, 3, 0.25), 
                  mgp = c(2, 1, 0))
                plot(y)
                lines(tobj$mu, col = 2)
                title("Observed data and estimated mean")
                plot(sqrt(tobj$sigma), col = 2, type = "l")
                title("Estimated standard deviation")
                plot(tobj$bi, type = "l")
                title(paste("Sum of weights  min=", signif(min(tobj$bi), 
                  3), " max=", signif(max(tobj$bi), 3)))
            }
            if (ddim == 2) {
                par(mfrow = c(2, 2), mar = c(1, 1, 3, 0.25), 
                  mgp = c(2, 1, 0))
                image(y, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title("Observed Image")
                image(tobj$mu, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                image(sqrt(tobj$sigma), col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Estimated SD  h=", signif(hakt, 
                  3), " min=", signif(sqrt(min(tobj$sigma)), 
                  3), " max=", signif(sqrt(max(tobj$sigma)), 
                  3)))
                image(tobj$bi, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title(paste("Sum of weights  min=", signif(min(tobj$bi), 
                  3), " max=", signif(max(tobj$bi), 3)))
            }
            if (ddim == 3) {
                n3h <- n3%/%2
                par(mfrow = c(2, 2), mar = c(1, 1, 3, 0.25), 
                  mgp = c(2, 1, 0))
                image(y[, , n3h], col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title("Observed Image")
                image(tobj$mu[, , n3h], col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                image(sqrt(tobj$sigma[, , n3h]), col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Estimated SD  h=", signif(hakt, 
                  3), " min=", signif(sqrt(min(tobj$sigma[, , 
                  n3h])), 3), " max=", signif(sqrt(max(tobj$sigma[, 
                  , n3h])), 3)))
                image(tobj$bi[, , n3h], col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Sum of weights  min=", signif(min(tobj$bi[, 
                  , n3h]), 3), " max=", signif(max(tobj$bi[, 
                  , n3h]), 3)))
            }
        }
        if (!is.null(u)) 
            cat("bandwidth: ", signif(hakt, 3), "   MSE: ", mean((tobj$mu - 
                u)^2), "   MAE: ", mean(abs(tobj$mu - u)), "\n")
        hakt <- hakt * hincr
        lambda0 <- lambda
        gc()
    }
    list(mu = tobj$mu, sigma = tobj$sigma, bi = tobj$bi, args = args)
  }
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
