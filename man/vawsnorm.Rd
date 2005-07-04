\name{vawsnorm}
\alias{vawsnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ AWS for local constant multivariate heteroskedastic regression }
\description{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for multivariate varying coefficient Gaussian models with unknown mean and variance on a 1D grid.
}
\usage{
vawsnorm(y, qlambda = NULL, qtau = NULL, model = "full", lkern = "Triangle", aggkern = "Uniform", hinit = NULL, hincr = NULL, hmax = NULL, heta = NULL, eta0 = NULL, graph = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{\code{y} contains the observed response data. \code{dim(y)[1]} 
  is the number of components of the response, \code{dim(y)[2]} 
  the number of observations.}
  \item{qlambda}{ \code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,1)}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is (for \code{qtau=1} a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  \code{qlambda} can be selected  as the smallest value that fulfils a propagation condition. 
  This means that in a parametric, in this case constant in mean and variance, model the resulting estimate is, for a large \code{hmax} almost parametric. 
  The default values of \code{qlambda}  are selected to fulfil this condition. Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate.}
  \item{qtau}{Stagewise Aggregation, see Belomestny and Spokoiny (2004) is used as an adaptive control step if \code{qtau<1}. 
   \code{qtau} determines the scaling parameter \code{\tau} in the Stagewise Aggregation algorithm in the same way as \code{qlambda} does
   for AWS. Default values are again selected by a propagation condition.}
  \item{model}{ \code{model} specifies the model for the variance. 
  Alternatives are \code{model="full"} for a general varying coefficient covariance structure 
  and \code{model="corr"} for varying coefficient variances and a 
  common (global) correlation structure.}
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
   \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle"}
  \item{hinit}{\code{hinit} specifies the initial bandwidth. Defaults to \code{hinit=1} }
  \item{hincr}{\code{hincr} specifies the factor used to increase the size of local neigborhoods after each iteration. The bandwidth is increased by
   a factor \code{hincr}. Defaults to \code{hincr=1.25}  }
  \item{hmax}{\code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=100 \code{dim(y)[1]}}
  \item{heta}{ \code{heta} specifies the minimal bandwidth to use with stagewise aggregation.  }
  \item{eta0}{ \code{eta0} minimal mixing coefficient in stagewise aggregation, defaults to 0 in case of \code{family="Gaussian"} and values
   larger than 0 otherwise to avoid problems at the boundary of the parameter space.}
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
}
\details{
The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for multivariate varying coefficient Gaussian models with unknown mean and variance on a 1D grid. 
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
\seealso{See Also  \code{\link{aws}}, \code{\link{awsnorm}}, \code{\link{vawsvola}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, qlambda = NULL, qtau = NULL, model = "full", lkern = "Triangle", 
    aggkern = "Uniform", hinit = NULL, hincr = NULL, hmax = NULL, 
    heta = NULL, eta0 = NULL, graph = FALSE) 
{
    KLnorm <- function(mu1, mu2, sig1, sig2) {
        ds <- dim(sig1)[1]
        n <- dim(sig1)[2]
        d <- dim(mu1)[1]
        .Fortran("klmnormn", as.double(mu1), as.double(mu2), 
            as.double(sig1), as.double(sig2), as.integer(d), 
            as.integer(ds), as.integer(n), double(d * d), double(d * 
                d), double(d * d), kld = double(n), PACKAGE = "aws")$kld
    }
    updtheta <- function(zobj, tobj, cpar, aggkern) {
        heta <- cpar$heta
        eta0 <- cpar$eta0
        tau1 <- cpar$tau1
        kstar <- cpar$kstar
        d <- cpar$d
        hakt <- zobj$hakt
        tau <- tau1 * (2 + max(kstar - log(hakt), 0))
        hakt <- zobj$hakt
        bi0 <- zobj$bi0
        bi <- zobj$bi
        n <- length(bi)
        yyt <- tobj$yyt
        sigma <- tobj$sigma
        mu <- tobj$mu
        munew <- t(t(zobj$ami)/bi)
        mumut <- matrix(0, ds, n)
        m <- 1
        for (i in 1:d) for (j in 1:i) {
            mumut[m, ] <- munew[i, ] * munew[j, ]
            m <- m + 1
        }
        sigmanew <- t(t(zobj$asi)/bi) - mumut
        sigmanew <- t(t(sigmanew) * bi/(bi - d))
        m <- 1
        for (i in 1:d) for (j in 1:i) {
            sigmanew[m, bi < (d + 1)] <- sigma[m, bi < (d + 1)]
            m <- m + 1
        }
        ind <- (1:d) * ((1:d) + 1)/2
        for (i in ind) {
            msi <- mean(sigmanew[i, ])
            sigmanew[i, sigmanew[i, ] < 1e-05 * msi] <- 1e-05 * 
                msi
        }
        if (cpar$model != "full") {
            sn <- sqrt(sigmanew[ind, ])
            m <- 1
            for (i in 1:d) for (j in 1:i) {
                if (i != j) 
                  sigmanew[m, ] <- mean((yyt[m, ] - mumut[m, 
                    ])/(sn[i, ] * sn[j, ])) * sn[i, ] * sn[j, 
                    ]
                m <- m + 1
            }
        }
        sigmanew[, tobj$fix] <- sigma[, tobj$fix]
        munew[, tobj$fix] <- mu[, tobj$fix]
        if (hakt > heta) {
            eta <- switch(aggkern, Uniform = , (1 - eta0) * as.numeric(bi0/tau * 
                KLnorm((1 - eta0) * munew + eta0 * mu, mu, (1 - 
                  eta0) * sigmanew + eta0 * sigma, sigma) > 1) + 
                eta0, Triangle = (1 - eta0) * pmin(1, bi0/tau * 
                KLnorm((1 - eta0) * munew + eta0 * mu, mu, (1 - 
                  eta0) * sigmanew + eta0 * sigma, sigma)) + 
                eta0)
        }
        else {
            eta <- rep(0, n)
        }
        eta[tobj$fix] <- 1
        bi <- (1 - eta) * bi + eta * tobj$bi
        for (i in 1:d) {
            mu[i, ] <- (1 - eta) * munew[i, ] + eta * mu[i, ]
        }
        for (i in 1:dim(sigma)[1]) {
            sigma[i, ] <- (1 - eta) * sigmanew[i, ] + eta * sigma[i, 
                ]
        }
        list(mu = mu, sigma = sigma, bi = bi, eta = eta, fix = (eta == 
            1), yyt = yyt)
    }
    args <- match.call()
    spmax <- 5
    d <- dim(y)[1]
    n <- dim(y)[2]
    ds <- d * (d + 1)/2
    yyt <- matrix(0, ds, n)
    m <- 1
    for (i in 1:d) for (j in 1:i) {
        yyt[m, ] <- y[i, ] * y[j, ]
        m <- m + 1
    }
    if (is.null(qlambda)) 
        qlambda <- 0.98
    if (is.null(qtau)) 
        if (qlambda == 1) 
            qtau <- 0.85
        else qtau <- 0.96
    if (model == "full") 
        d1 <- (d + 1)/2 + 1
    else d1 <- 2
    if (qtau < 1) 
        tau1 <- qchisq(qtau, 2)
    else tau1 <- 1e+50
    if (aggkern == "Triangle") 
        tau1 <- 2.5 * tau1
    if (is.null(eta0)) 
        eta0 <- 0.25
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, 2)
    if (qlambda < 1) 
        lambda <- 2 * qchisq(qlambda, d1 * d)
    else lambda <- 1e+50
    if (is.null(hinit) || hinit < max(10, 2 * d)) 
        hinit <- max(10, 2 * d)
    if (is.null(hincr) || hincr <= 1) 
        hincr <- 1.25
    if (is.null(hmax)) 
        hmax <- 100 * d
    if (is.null(heta)) 
        heta <- max(2 * (d1 + d), hinit + 1)
    cpar <- list(heta = heta, tau1 = tau1, eta0 = eta0, model = model, 
        kstar = log(d * 100), d = d)
    tobj <- list(yyt = yyt, bi = rep(1, n), mu = y, sigma = yyt, 
        fix = rep(FALSE, n))
    zobj <- list(asi = yyt, ami = y, bi0 = rep(1, n))
    bi0old <- rep(1, n)
    hakt <- hinit
    lambda0 <- lambda
    if (hinit > 1) 
        lambda0 <- 1e+50
    while (hakt <= hmax) {
        zobj <- .Fortran("cvawsnun", as.double(y), as.double(yyt), 
            as.logical(tobj$fix), as.integer(n), as.integer(d), 
            as.integer(ds), hakt = as.double(hakt), as.double(lambda0), 
            as.double(tobj$mu), as.double(tobj$sigma), bi = as.double(tobj$bi), 
            bi0 = as.double(zobj$bi0), ami = as.double(zobj$ami), 
            asi = as.double(zobj$asi), as.integer(lkern), as.double(spmax), 
            double(ds), double(d), double(d * d), double(d * 
                d), double(d * d), PACKAGE = "aws")[c("bi", "bi0", 
            "ami", "asi", "hakt")]
        if (hakt > n/2) 
            zobj$bi0 <- hincr * biold
        dim(zobj$ami) <- c(d, n)
        dim(zobj$asi) <- c(ds, n)
        biold <- zobj$bi0
        cat("Now update for h=", hakt, "\n")
        tobj <- updtheta(zobj, tobj, cpar, aggkern)
        if (graph) {
            par(mfrow = c(1, 3), mar = c(3, 3, 3, 0.2), mgp = c(2, 
                1, 0))
            plot(y[1, ], ylim = range(y[1, ], tobj$mu[1, ]), 
                col = 3)
            lines(tobj$mu[1, ], lwd = 2)
            title(paste("Mean  h=", signif(hakt, 3)))
            plot(yyt[1, ] - tobj$mu[1, ]^2, ylim = range(yyt[1, 
                ] - tobj$mu[1, ]^2, tobj$sigma[1, ]), col = 3)
            lines(tobj$sigma[1, ], lwd = 2)
            title(paste("Variance  h=", signif(hakt, 3)))
            plot(tobj$bi, type = "l", ylim = range(0, tobj$bi))
            lines(tobj$eta * max(tobj$bi), col = 2)
            title("Sum of weights and eta")
        }
        hakt <- hakt * hincr
        lambda0 <- lambda
        gc()
    }
    list(mu = tobj$mu, sigma = tobj$sigma, bi = tobj$bi, args = args)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
