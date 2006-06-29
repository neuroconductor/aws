\name{glmaws}
\alias{glmaws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Varying coefficient generalized linear models}
\description{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for varying coefficient generalized linear models.
}
\usage{
glmaws(y, degree = 1, family = "Gaussian", qlambda = NULL, heta = NULL, qtau = NULL, lkern = "Triangle", skern="Exponential",aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, hw = NULL, lseq = NULL, iter = 1000, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL, spmax = NULL, eps = 1e-08, showwghts = FALSE, conf = FALSE, usevar = TRUE, kstar = NULL, tau2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Observed data on a 1D or 2D grid (1D currently implemented}
  \item{degree}{ Degree of polynomial basis. \code{degree=1,2} are currently implemented. For polynomial degree 
  0 use \code{aws}. }
  \item{family}{ This argument determines the family of distributions for \code{y}. Currently implemented are
  \code{family="Poisson"}, \code{family="Gaussian"}, \code{family="Bernoulli"} and \code{family="Exponential"}.
  The canonical link function is used in each case.}
  \item{qlambda}{ \code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,degree+1)}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  \code{qlambda} can be selected for a given model, i.e. specification of \code{family}, as the smallest value that fulfils a propagation condition.
  This means that in a parametric, in this case generalized linear, model the resulting estimate is, for a large \code{hmax} almost parametric. 
  The default values of \code{qlambda} for the different specifications of \code{family} are selected to fulfil this condition. Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate. Defaults currently only correct for \code{family="Poisson"}.}
  \item{heta}{ \code{heta} specifies the minimal bandwidth to use with stagewise aggregation.  }
  \item{qtau}{ Stagewise Aggregation, see Belomestny and Spokoiny (2004) is used as an adaptive control step if \code{qtau<1}. 
   \code{qtau} determines the scaling parameter \code{\tau} in the Stagewise Aggregation algorithm in the same way as \code{qlambda} does
   for AWS. Default values are again selected, for the probability distribution specified with the \code{family} argument, by a propagation condition.}
  \item{lkern}{  \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
  \item{skern}{ \code{skern} specifies the kernel \eqn{K_s}. Defaults to "Exp", the altenative choice is "Triangle".}
   \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle"}
  \item{sigma2}{\code{sigma2} allows to specify the variance in case of \code{family="Gaussian"}. Not used if \code{family!="Gaussian"}.
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==length(y)} then \code{sigma2}
   is assumed to contain the pointwise variance of \code{y} and a heteroscedastic variance model is used. }
  \item{hinit}{\code{hinit} specifies the initial bandwidth. Defaults to \code{hinit=1}. }
  \item{hincr}{\code{hincr} specifies the factor used to increase the size of local neigborhoods after each iteration. The bandwidth is increased by
   a factor \code{hincr^(1/dd)} with \code{dd} specifying the dimensionality of the grid. Defaults to \code{hincr=1.25}  }
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth.}
  \item{hw}{Determines a minimal regularity of weights. Weights are spread if the number of observations
   with positive weights is small and the parameters may therefore numerically non-identifiable. Larger values of \code{hw} lead to more regularisation. }
  \item{lseq}{\code{lseq} allows to increase the value of the scaling parameter \code{\lambda} for the first \code{length(lseq)} iteration steps by
  the factor specified in \code{lseq}. Defaults to NULL. In this case a default sequence is used that fulfils the propagation condition. }
  \item{iter}{ Maximal number of iteration to numerically solve the estimating equation for a given weighting scheme.
   Should be large (default 1000). }
  \item{u}{ \code{u} allows to specify the true parameter. This is only used to test the algorithm and to select the smoothing parameters 
    \code{qlambda} and \code{qtau} by a propagation condition. If \code{u} is specified MSE and MAE of the estimates are 
    printed for each iteration step.}
  \item{graph}{ If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
  \item{demo}{ If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{wghts}{\code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space. }
  \item{spmax}{ ~~Describe \code{spmax} here~~ }
  \item{eps}{ ~~Describe \code{eps} here~~ }
  \item{showwghts}{}
  \item{conf}{ Compute confidence information  }
  \item{usevar}{ Use variance instead of B_i }
  \item{kstar}{ Stagewise aggregation parameter. take defaults ... }
  \item{tau2}{ Stagewise aggregation parameter. take defaults ... }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list of class \code{glmaws} with components
  \item{yhat }{Estimate of E(y|x)}
  \item{theta }{Estimated parameters }
  \item{confint }{approximative 0.95-Confidence intervals for yhat, conditional on the
     weighting scheme. This ignores the data-dependence of the weights.}
  \item{y}{The data}
  \item{mae}{if \code{!is.null(u)} this containes the Mean Absolute Error. Needed when testing the Propagation Condition. Otherwise NULL.}
  \item{lseq}{Sequence of correction factors used for lambda.}
  \item{call}{arguments of the call to glmaws}}  
}
\references{ ~put references to the literature/web site here ~ }
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html}}
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, degree = 1, family = "Gaussian", qlambda = NULL, 
    heta = NULL, qtau = NULL, lkern = "Triangle", aggkern = "Uniform", 
    sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, hw = NULL, 
    lseq = NULL, iter = 50, u = NULL, graph = FALSE, demo = FALSE, 
    wghts = NULL, spmax = 5, eps = 1e-08, showwghts = FALSE, 
    conf = FALSE, usevar = TRUE, kstar = NULL, tau2 = NULL) 
{
    Pardist <- function(mcode, Bi0, dtheta) {
        dp1 <- dim(dtheta)[1]
        dp2 <- dim(Bi0)[1]
        if (mcode == 1) {
            dist <- 0
            for (i in 1:dp1) for (j in 1:dp1) dist <- dist + 
                dtheta[i, ] * Bi0[i + j - 1, ] * dtheta[j, ]
        }
        if (mcode == 2) {
            ind <- matrix(c(1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 
                9, 3, 5, 6, 8, 9, 10, 4, 7, 8, 11, 12, 13, 5, 
                8, 9, 12, 13, 14, 6, 9, 10, 13, 14, 15), 6, 6)[1:dp1, 
                1:dp1]
            dist <- 0
            for (i in 1:dp1) for (j in 1:dp1) dist <- dist + 
                dtheta[i, , ] * Bi0[ind[i, j], , ] * dtheta[j, 
                  , ]
        }
        if (mcode == 3) {
            for (i in 1:dp1) dist <- dist + dtheta[i, ] * Bi0[i + 
                (i - 1) * i/2, ] * dtheta[i, ]
            for (j in 2:dp1) for (i in 1:(j - 1)) dist <- dist + 
                2 * dtheta[i, ] * Bi0[i + (j - 1) * j/2, ] * 
                  dtheta[j, ]
        }
        dist
    }
    updtheta <- function(zobj, tobj, cpar, aggkern, hakt) {
        regdiff <- function(x1, x2) {
            dx <- dim(x1)
            x <- x1 - x2
            dim(x) <- dx
            x
        }
        eta0 <- 0
        heta <- cpar$heta
        tau1 <- cpar$tau1
        tau2 <- cpar$tau2
        kstar <- cpar$kstar
        hakt <- zobj$hakt
        tau <- 2 * (tau1 + tau2 * max(kstar - log(hakt), 0))
        mcode <- cpar$mcode
        mfamily <- cpar$mfamily
        bi0 <- zobj$bi0
        bi <- zobj$bi
        bi2 <- zobj$bi2
        theta <- tobj$theta
        thetanew <- zobj$theta
        dp1 <- dim(theta)[1]
        n <- dim(theta)[2]
        if (mfamily == 2) {
            failed <- (1:n)[(hakt^(0:(dp1 - 1))) \%*\% abs(thetanew) > 
                100]
            failed <- failed[!tobj$fix[failed]]
            tobj$fix[failed] <- TRUE
            if (length(failed) > 0) 
                cat("Old estimate kept in", failed, "\n")
        }
        if (mfamily == 3) {
            failed <- (1:n)[(hakt^(0:(dp1 - 1))) \%*\% abs(thetanew) > 
                100]
            failed <- failed[!tobj$fix[failed]]
            tobj$fix[failed] <- TRUE
            if (length(failed) > 0) 
                cat("Old estimate kept in", failed, "\n")
        }
        dd <- dim(theta)
        if (hakt > heta) {
            eta <- switch(aggkern, Uniform = (1 - eta0) * as.numeric(Pardist(mcode, 
                bi0, regdiff(thetanew, theta)) > tau/2.5) + eta0, 
                Triangle = (1 - eta0) * pmin(1, Pardist(mcode, 
                  bi0, regdiff(thetanew, theta))/tau) + eta0)
            if (length(dd) > 2) 
                dim(eta) <- dd[-1]
        }
        else {
            eta <- rep(eta0, prod(dim(bi)[-1]))
            dim(eta) <- dim(bi)[-1]
        }
        eta[tobj$fix] <- 1
        thetanew[, tobj$fix] <- 0
        bi[, tobj$fix] <- 0
        bi2[, tobj$fix] <- 0
        dp1 <- dim(zobj$theta)[1]
        dp2 <- dim(bi)[1]
        etadd <- outer(rep(1, dp1), eta)
        theta <- (1 - etadd) * thetanew + etadd * theta
        etadd <- outer(rep(1, dp2), eta)
        bi <- (1 - etadd) * bi + etadd * tobj$bi
        bi2 <- (1 - etadd) * bi2 + etadd * tobj$bi2
        list(theta = theta, bi = bi, bi2 = bi2, eta = eta, fix = (eta == 
            1))
    }
    if (degree == 0) 
        return("use aws for local constant models")
    if (degree > 3) 
        return("no defaults for parameters available")
    mfamily <- switch(family, Gaussian = 1, Poisson = 2, Bernoulli = 3, 
        Exponential = 4, 1)
    args <- match.call()
    if (is.null(qlambda)) {
        if (is.null(dim(y))) {
            qlambda <- switch(family, Gaussian = switch(degree, 
                0.65, 0.966), Poisson = switch(degree, 0.35, 
                0.5), Bernoulli = switch(degree, 0.985, 0.999), 
                Exponential = switch(degree, 0.65, 0.999), Volatility = switch(degree, 
                  0.65, 0.999))
        }
        else qlambda <- switch(degree, 0.65, 0.92)
    }
    if (qlambda < 0.6) 
        warning("Inappropriate value of qlambda")
    if (is.null(dim(y))) {
        if (is.null(lseq)) 
            lseq <- switch(family, Gaussian = 1.3, Bernoulli = 1, 
                Exponential = 1.7, Poisson = 1, 1.3)
    }
    if (demo && !graph) 
        graph <- TRUE
    dy <- dim(y)
    if (is.null(hinit) || hinit <= 0) 
        hinit <- 1
    if (is.null(dy)) {
        form <- "uni"
        ddim <- 1
        n <- length(y)
        dp1 <- degree + 1
    }
    if (length(dy) == 2) {
        form <- "bi"
        ddim <- 2
        if (is.null(wghts)) 
            wghts <- c(1, 1)
        hinit <- hinit/wghts[1]
        hmax <- hmax/wghts[1]
        wghts <- (wghts[2]/wghts[1])
        n1 <- dy[1]
        n2 <- dy[2]
        n <- n1 * n2
        if (degree > 2) 
            return("bivariate aws on a grid is not implemented for degree>2")
        dp1 <- switch(degree + 1, 1, 3, 6)
    }
    if (length(dy) > 2) 
        return("polynomial AWS for more than 2 dimensional grids is not implemented")
    mae <- NULL
    if (is.null(hincr)) 
        hincr <- 1.25^(1/ddim)
    if (is.null(sigma2)) {
        if (family == "Gaussian") {
            sigma2 <- IQRdiff(y)^2
            cat("sigma^2=", sigma2, "\n")
            sigma2 <- sigma2
        }
        else sigma2 <- 1
    }
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, 2)
    if (is.null(dy)) {
        if (is.null(heta)) 
            heta <- switch(degree, 10, 25, 100)
        if (is.null(qtau)) 
            tau <- switch(degree, 100, 500, 4000)
        else if (qtau < 1) 
            tau <- qchisq(qtau, dp1)
        else {
            tau <- 1e+40
            heta <- 1e+40
        }
        if (is.null(tau2)) 
            tau2 <- (degree + 1)^2 * tau
        if (is.null(kstar)) 
            kstar <- log(switch(degree, 150, 300, 600))
    }
    else {
        if (is.null(heta)) 
            heta <- switch(degree, 3, 4)
        tau <- switch(degree, 4, 12)
        tau2 <- (degree + 1)^2 * tau
        kstar <- log(switch(degree, 15, 30))
    }
    if (qlambda >= 1) 
        lamakt <- 1e+50
    else lamakt <- 2 * qchisq(qlambda, dp1) * sigma2
    cat("Value of lambda", signif(lamakt, 3), "\n")
    eta0 <- 0
    cpar <- list(heta = heta, eta0 = eta0, tau1 = tau * sigma2, 
        tau2 = tau2 * sigma2, kstar = kstar, mfamily = mfamily)
    steps <- as.integer(log(hmax/hinit)/log(hincr) + 1)
    if (is.null(lseq)) 
        lseq <- 1
    if (length(lseq) < steps) 
        lseq <- c(lseq, rep(1, steps - length(lseq)))
    lseq <- lseq[1:steps]
    k <- 1
    if (form == "uni" && degree > 0) {
        confint <- NULL
        cpar$mcode <- mcode <- 1
        dp1 <- degree + 1
        dp2 <- degree + dp1
        dp3 <- (dp1 + 1) * dp1/2
        if (is.null(hinit) || hinit < 1) 
            hinit <- 1
        cb <- matrix(0, dp1, dp1)
        for (i in (1:dp1)) cb[i:dp1, i] <- choose((i:dp1) - 1, 
            i - 1)
        hakt <- hinit
        tobj <- list(bi = matrix(0, dp2, n), bi2 = matrix(0, 
            dp2, n), theta = matrix(0, dp1, n), fix = rep(FALSE, 
            n))
        theta <- matrix(0, dp1, n)
        tobj$theta[1, ] <- theta[1, ] <- switch(mfamily, Gaussian = mean(y), 
            Poisson = log(mean(y)), Bernoulli = log(mean(y)/(1 - 
                mean(y))), Exponential = 1/mean(y))
        biold <- biold2 <- matrix(0, dp2, n)
        bii <- matrix(0, dp3, n)
        zobj <- list(bi0 = biold, bi02 = biold, ni = rep(hinit, 
            n))
        lamakt0 <- 1e+50
        if (is.null(hw)) 
            hw <- switch(ddim, switch(degree, 5, 10), degree + 
                0.1)
        else hw <- max(hw, degree + 0.1)
        while (hakt <= hmax) {
            if (showwghts) {
                zobj <- .Fortran("glawsunw", as.integer(n), as.integer(dp1), 
                  as.integer(dp2), as.integer(dp3), as.double(y), 
                  fix = as.logical(tobj$fix), as.integer(mfamily), 
                  as.double(tobj$theta), theta = as.double(tobj$theta), 
                  as.double(bii), bi = as.double(tobj$bi), bi2 = as.double(tobj$bi2), 
                  bi0 = as.double(zobj$bi0), bi02 = as.double(zobj$bi02), 
                  ai = double(dp1), ni = as.double(zobj$ni), 
                  as.double(lamakt0), hakt = as.double(hakt), 
                  as.double(hw), as.integer(lkern), as.double(cb), 
                  double(dp1 * dp1), double(dp1), double(dp2), 
                  wghts = double(n * n), double(n), as.double(spmax), 
                  as.integer(iter), double(n), PACKAGE = "aws")[c("fix", 
                  "theta", "bi", "bi0", "bi2", "bi02", "ni", 
                  "hakt", "wghts")]
            }
            else {
                zobj <- .Fortran("glawsuni", as.integer(n), as.integer(dp1), 
                  as.integer(dp2), as.integer(dp3), as.double(y), 
                  fix = as.logical(tobj$fix), as.integer(mfamily), 
                  as.double(tobj$theta), theta = as.double(tobj$theta), 
                  as.double(bii), bi = as.double(tobj$bi), bi2 = as.double(tobj$bi2), 
                  bi0 = as.double(zobj$bi0), bi02 = as.double(zobj$bi02), 
                  ai = double(dp1), ni = as.double(zobj$ni), 
                  as.double(lamakt0), hakt = as.double(hakt), 
                  as.double(hw), as.integer(lkern), as.double(cb), 
                  double(dp1 * dp1), double(dp1), double(dp2), 
                  double(n), double(n), as.double(spmax), as.integer(iter), 
                  double(n), PACKAGE = "aws")[c("fix", "theta", 
                  "bi", "bi0", "bi2", "bi02", "ni", "hakt")]
            }
            dim(zobj$theta) <- c(dp1, n)
            dim(zobj$bi) <- c(dp2, n)
            dim(zobj$bi0) <- c(dp2, n)
            dim(zobj$bi2) <- c(dp2, n)
            dim(zobj$bi02) <- c(dp2, n)
            failed <- (1:n)[tobj$fix != zobj$fix]
            zobj$bi0[, failed] <- biold[, failed]
            zobj$bi02[, failed] <- biold2[, failed]
            if (length(failed) > 0) 
                cat("No convergence in", failed, "\n")
            tobj$fix[failed] <- TRUE
            gc()
            if (hakt > n/2) {
                zobj$bi0 <- hincr * biold
                zobj$bi02 <- hincr * biold2
            }
            biold <- zobj$bi0
            biold2 <- zobj$bi02
            tobj <- updtheta(zobj, tobj, cpar, aggkern, hakt)
            if (usevar) {
                bii <- matrix(.Fortran("bibi2ibi", as.integer(n), 
                  as.integer(dp1), as.integer(dp2), as.integer(dp3), 
                  as.double(tobj$bi), as.double(tobj$bi2), bii = double(dp3 * 
                    n), double(dp1 * dp1), double(dp1 * dp1), 
                  PACKAGE = "aws")$bii, dp3, n)
            }
            else {
                bii <- tobj$bi[c(1:3, 3:5, 4:7)[1:dp3], ]
            }
            gc()
            if (conf) {
                confint <- matrix(.Fortran("confuni", as.integer(n), 
                  as.integer(dp1), as.integer(dp2), as.double(tobj$theta[1, 
                    ]), as.double(tobj$bi), as.double(tobj$bi2), 
                  conf = double(2 * n), double(dp1 * dp1), PACKAGE = "aws")$conf, 
                  2, n)
            }
            if (graph) {
                if (showwghts) 
                  par(mfrow = c(1, 3), mar = c(3, 3, 3, 1))
                else par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
                plot(y, ylim = range(y), col = 3)
                if (!is.null(u)) 
                  lines(u, col = 2)
                lines(switch(mfamily, Gaussian = tobj$theta[1, 
                  ], Poisson = exp(tobj$theta[1, ]), Bernoulli = exp(tobj$theta[1, 
                  ])/(1 + exp(tobj$theta[1, ])), Exponential = 1/tobj$theta[1, 
                  ]), lwd = 2)
                if (conf) {
                  lines(switch(mfamily, Gaussian = confint[1, 
                    ], Poisson = exp(confint[1, ]), Bernoulli = exp(confint[1, 
                    ])/(1 + exp(confint[1, ])), Exponential = 1/confint[1, 
                    ]), lwd = 1, col = 4)
                  lines(switch(mfamily, Gaussian = confint[2, 
                    ], Poisson = exp(confint[2, ]), Bernoulli = exp(confint[2, 
                    ])/(1 + exp(confint[2, ])), Exponential = 1/confint[2, 
                    ]), lwd = 1, col = 4)
                }
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                plot(zobj$ni, type = "l", ylim = c(0, max(zobj$ni)))
                lines(tobj$eta * max(zobj$ni), col = 2)
                title("Sum of weights and eta")
                if (showwghts) 
                  image(matrix(zobj$wghts, n, n), zlim = c(0, 
                    1), col = gray((0:255)/255))
            }
            if (!is.null(u)) {
                mtheta <- switch(mfamily, Gaussian = tobj$theta[1, 
                  ], Poisson = exp(tobj$theta[1, ]), Bernoulli = exp(tobj$theta[1, 
                  ])/(1 + exp(tobj$theta[1, ])), Exponential = 1/tobj$theta[1, 
                  ])
                cat("bandwidth: ", signif(hakt, 3), "eta==1", 
                  sum(tobj$eta == 1), "   MSE: ", mean((mtheta - 
                    u)^2), "   MAE: ", mean(abs(mtheta - u)), 
                  "\n")
                mae <- c(mae, signif(mean(abs(mtheta - u)), 3))
            }
            if (demo) 
                readline("Press return")
            hakt <- hakt * hincr
            lamakt0 <- lamakt * lseq[k]
            k <- k + 1
            gc()
        }
    }
    if (form == "bi" && p > 0) {
        cpar$mcode <- 2
        dp1 <- switch(p + 1, 1, 3, 6)
        dp2 <- switch(p + 1, 1, 6, 15)
        dpm <- 1
        ind <- matrix(c(1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, 3, 
            5, 6, 8, 9, 10, 4, 7, 8, 11, 12, 13, 5, 8, 9, 12, 
            13, 14, 6, 9, 10, 13, 14, 15), 6, 6)[1:dp1, 1:dp1]
        if (is.null(hinit) || hinit < p + 0.1) 
            hinit <- p + 0.1
        if (demo) 
            readline("Press return")
        hakt <- hinit
        tobj <- list(bi = array(0, c(dp2, n1, n2)), theta = array(0, 
            c(dp1, n1, n2)), fix = matrix(FALSE, n1, n2))
        biold <- array(0, c(dp2, n1, n2))
        zobj <- list(ai = array(0, c(dp1, n1, n2)), bi0 = biold)
        lamakt0 <- 1e+50
        while (hakt <= hmax) {
            zobj <- .Fortran("cpawsbi", as.integer(n1), as.integer(n2), 
                as.integer(dp1), as.integer(dp2), as.double(y), 
                as.logical(tobj$fix), as.double(tobj$theta), 
                bi = as.double(tobj$bi), bi0 = as.double(zobj$bi0), 
                ai = as.double(zobj$ai), as.double(lamakt0), 
                hakt = as.double(hakt), as.integer(lkern), double(dp1 * 
                  dp1), double(dp1), double(dp2), double(dp2), 
                double(dp2), double(dp1), as.integer(ind), as.double(wghts), 
                as.double(spmax), PACKAGE = "aws")[c("ai", "bi", 
                "bi0", "hakt")]
            gc()
            dim(zobj$ai) <- c(dp1, n1, n2)
            dim(zobj$bi) <- c(dp2, n1, n2)
            dim(zobj$bi0) <- c(dp2, n1, n2)
            if (hakt > min(n1, n2)/2) 
                zobj$bi0 <- hincr * hincr * biold
            biold <- zobj$bi0
            tobj <- updtheta(zobj, tobj, cpar)
            gc()
            lamakt0 <- lamakt
            if (graph) {
                par(mfrow = c(2, 2), mar = c(1, 1, 3, 0.25), 
                  mgp = c(2, 1, 0))
                image(y, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title("Observed Image")
                image(tobj$theta[1, , ], col = gray((0:255)/255), 
                  zlim = range(y), xaxt = "n", yaxt = "n")
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                image(tobj$bi[1, , ], col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title("Sum of weights")
                image(matrix(tobj$eta, n1, n2), col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n", zlim = c(0, 1))
                title("eta")
            }
            if (!is.null(u)) 
                cat("bandwidth: ", signif(hakt, 3), " eta==1: ", 
                  sum(tobj$eta == 1), "   MSE: ", mean((tobj$theta[1, 
                    , ] - u)^2), "   MAE: ", mean(abs(tobj$theta[1, 
                    , ] - u)), "\n")
            hakt <- hakt * hincr
            gc()
        }
    }
    z <- list(theta = tobj$theta, confint = confint, y = y, x = x, 
        mae = mae, lseq = c(0, lseq[-steps]), call = args)
    class(z) <- "aws"
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
