\name{cpaws}
\alias{cpaws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cpaws(y, x = NULL, p = 1, sigma2 = NULL, qlambda = NULL, heta = NULL, tau1 = NULL, tau2 = NULL, eta0 = 0, lkern = "Triangle", hinit = NULL, hincr = NULL, hmax = NULL, NN = FALSE, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL, spmax = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{sigma2}{ ~~Describe \code{sigma2} here~~ }
  \item{qlambda}{ ~~Describe \code{qlambda} here~~ }
  \item{heta}{ ~~Describe \code{heta} here~~ }
  \item{tau1}{ ~~Describe \code{tau1} here~~ }
  \item{tau2}{ ~~Describe \code{tau2} here~~ }
  \item{eta0}{ ~~Describe \code{eta0} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{hinit}{ ~~Describe \code{hinit} here~~ }
  \item{hincr}{ ~~Describe \code{hincr} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{NN}{ ~~Describe \code{NN} here~~ }
  \item{u}{ ~~Describe \code{u} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{demo}{ ~~Describe \code{demo} here~~ }
  \item{wghts}{ ~~Describe \code{wghts} here~~ }
  \item{spmax}{ ~~Describe \code{spmax} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x = NULL, p = 1, sigma2 = NULL, qlambda = NULL, 
    heta = NULL, tau1 = NULL, tau2 = NULL, eta0 = 0, lkern = "Triangle", 
    hinit = NULL, hincr = NULL, hmax = NULL, NN = FALSE, u = NULL, 
    graph = FALSE, demo = FALSE, wghts = NULL, spmax = 5) 
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
    gettheta <- function(mcode, ai, bi) {
        if (mcode == 1) {
            n <- dim(ai)[2]
            dp1 <- dim(ai)[1]
            dp2 <- dim(bi)[1]
            theta <- matrix(.Fortran("mpawsuni", as.integer(n), 
                as.integer(dp1), as.integer(dp2), as.double(ai), 
                as.double(bi), theta = double(dp1 * n), double(dp1 * 
                  dp1), PACKAGE = "awsnew")$theta, dp1, n)
        }
        if (mcode == 2) {
            n1 <- dim(ai)[2]
            n2 <- dim(ai)[3]
            n <- n1 * n2
            dp1 <- dim(ai)[1]
            dp2 <- dim(bi)[1]
            ind <- matrix(c(1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 
                9, 3, 5, 6, 8, 9, 10, 4, 7, 8, 11, 12, 13, 5, 
                8, 9, 12, 13, 14, 6, 9, 10, 13, 14, 15), 6, 6)[1:dp1, 
                1:dp1]
            theta <- array(.Fortran("mpawsbi", as.integer(n), 
                as.integer(dp1), as.integer(dp2), as.double(ai), 
                as.double(bi), theta = double(dp1 * n), double(dp1 * 
                  dp1), as.integer(ind), PACKAGE = "awsnew")$theta, 
                c(dp1, n1, n2))
        }
        if (mcode == 3) {
            n <- dim(ai)[2]
            dp1 <- dim(ai)[1]
            dp2 <- dim(bi)[1]
            theta <- matrix(.Fortran("mpawsmul", as.integer(n), 
                as.integer(dp1), as.integer(dp2), as.double(ai), 
                as.double(bi), double(dp1 * n), double(dp1 * 
                  dp1), PACKAGE = "awsnew")[[6]], dp1, n)
        }
        theta
    }
    updtheta <- function(zobj, tobj, cpar) {
        heta <- cpar$heta
        eta0 <- cpar$eta0
        tau1 <- cpar$tau1
        tau2 <- cpar$tau2
        kstar <- cpar$kstar
        hakt <- zobj$hakt
        tau <- 2 * (tau1 + tau2 * max(kstar - log(hakt), 
            0))
        mcode <- cpar$mcode
        bi0 <- zobj$bi0
        bi <- zobj$bi
        thetanew <- gettheta(mcode, zobj$ai, bi)
        theta <- tobj$theta
        if (hakt > heta) {
            eta <- (1 - eta0) * pmin(1, Pardist(mcode, bi0, thetanew - 
                theta)/tau) + eta0
        }
        else {
            eta <- rep(eta0, prod(dim(bi)[-1]))
            dim(eta) <- dim(bi)[-1]
        }
        eta[tobj$fix] <- 1
        dp1 <- dim(zobj$ai)[1]
        dp2 <- dim(bi)[1]
        if (mcode == 2) {
            for (i in 1:dp1) theta[i, , ] <- (1 - eta) * thetanew[i, 
                , ] + eta * theta[i, , ]
            for (i in 1:dp2) bi[i, , ] <- (1 - eta) * bi[i, , 
                ] + eta * tobj$bi[i, , ]
        }
        else {
            for (i in 1:dp1) theta[i, ] <- (1 - eta) * thetanew[i, 
                ] + eta * theta[i, ]
            for (i in 1:dp2) bi[i, ] <- (1 - eta) * bi[i, ] + 
                eta * tobj$bi[i, ]
        }
        list(theta = theta, bi = bi, eta = eta, fix = (eta == 
            1))
    }
    if (p == 0) 
        return("use caws for local constant models")
    if (p > 3) 
        return("no defaults for parameters available")
    args <- match.call()
    if (is.null(qlambda)) {
        if (is.null(dim(y))) 
            qlambda <- switch(p, 0.65, 0.966, 0.966)
        else qlambda <- switch(p, 0.65, 0.92)
    }
    if (qlambda < 0.6) 
        return("Inappropriate value of qlambda")
    if (demo && !graph) 
        graph <- TRUE
    gridded <- is.null(x)
    dy <- dim(y)
    if (gridded) {
        if (is.null(hinit) || hinit <= 0) 
            hinit <- p + 1
        if (is.null(dy)) {
            form <- "uni"
            ddim <- 1
            n <- length(y)
            dp1 <- p + 1
        }
        if (length(dy) == 2) {
            form <- "bi"
            ddim <- 2
            if (is.null(wghts)) 
                wghts <- c(1, 1)
            hinit <- hinit/wghts[1]
            hmax <- hmax/wghts[1]
            wghts <- (wghts[2]/wghts[1])^2
            n1 <- dy[1]
            n2 <- dy[2]
            n <- n1 * n2
            if (p > 2) 
                return("bivariate aws on a grid is not implemented for p>2")
            dp1 <- switch(p + 1, 1, 3, 6)
        }
        if (length(dy) > 2) 
            return("polynomial AWS for more than 2 dimensional grids is not implemented")
    }
    else {
        dx <- dim(x)
        ddim <- 1
        if (is.null(dx)) {
            form <- "uni"
            n <- length(x)
            if (n != length(y)) 
                return("incompatible lengths of x and y")
            ox <- order(x)
            x <- x[ox]
            y <- y[ox]
            dp1 <- p + 1
        }
        else {
            px <- dx[1]
            n <- dx[2]
            if (p > 1) {
                p <- 1
                cat("p is set to 1, the maximal polynomial degree implemented")
            }
            form <- "multi"
            if (n != length(y)) 
                return("incompatible dimensions of x and y")
            if (is.null(wghts) || length(wghts) != px) 
                wghts <- rep(1, px)
            dp1 <- 1 + p * px
            if (NN) {
                ihmax <- trunc(hmax)
                if (ihmax > n) 
                  ihmax <- n
                neighbors <- matrix(0, ihmax, n)
                for (i in 1:n) {
                  if (px == 1) 
                    adist <- (x - x[i])^2
                  else adist <- wghts %*% ((x - x[, i])^2)
                  neighbors[, i] <- order(adist)[1:ihmax]
                }
            }
            else {
                ihmax <- n
                ddim <- px
                neighbors <- distmat <- matrix(0, n, n)
                for (i in 1:n) {
                  if (px == 1) 
                    adist <- (x - x[i])^2
                  else adist <- wghts %*% ((x - x[, i])^2)
                  od <- order(adist)
                  distmat[, i] <- adist[od]
                  neighbors[, i] <- od
                }
                gc()
                distmat <- sqrt(distmat)
                maxdist <- apply(distmat, 1, max)
                ihmax <- sum(maxdist <= hmax)
                distmat <- distmat[1:ihmax, ]
                neighbors <- neighbors[1:ihmax, ]
                maxdist <- maxdist[1:ihmax]
                gc()
            }
        }
        if (length(y) != n) 
            return("incompatible dimensions of x and y")
    }
    if (is.null(hincr)) 
        hincr <- 1.25^(1/ddim)
    if (is.null(sigma2)) {
        sigma2 <- IQRdiff(y)^2
        cat("sigma^2=", sigma2, "\n")
    }
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, 2)
    if (is.null(dy)) {
        if (is.null(heta)) 
            heta <- switch(p, 10, 25, 100)
        if (is.null(tau1)) 
            tau1 <- switch(p, 100, 500, 4000)
        if (is.null(tau2)) 
            tau2 <- switch(p, 500, 3000, 80000)
        kstar <- log(switch(p, 150, 300, 600))
    }
    else {
        if (is.null(heta)) 
            heta <- switch(p, 3, 4)
        if (is.null(tau1)) 
            tau1 <- switch(p, 4, 12)
        if (is.null(tau2)) 
            tau2 <- switch(p, 30, 50)
        kstar <- log(switch(p, 15, 30))
    }
    if (qlambda >= 1) 
        lamakt <- 1e+50
    else lamakt <- 2 * qchisq(qlambda, dp1) * sigma2
    cpar <- list(heta = heta, eta0 = eta0, tau1 = tau1 * 
        sigma2, tau2 = tau2 * sigma2, kstar = kstar)
    if (form == "uni" && p > 0) {
        cpar$mcode <- 1
        if (gridded) 
            x <- 1:length(y)
        dp1 <- p + 1
        dp2 <- p + dp1
        dxp <- max(diff(x, p)) + 0.1
        if (is.null(hinit) || hinit < dxp) 
            hinit <- dxp
        cb <- matrix(0, dp1, dp1)
        for (i in (1:dp1)) cb[i:dp1, i] <- choose((i:dp1) - 1, 
            i - 1)
        hakt <- hinit
        tobj <- list(bi = matrix(0, dp2, n), theta = matrix(0, 
            dp1, n), fix = rep(FALSE, n))
        bi0old <- matrix(0, dp2, n)
        zobj <- list(bi0 = bi0old, ai = matrix(0, dp1, n))
        lamakt0 <- 1e+50
        while (hakt <= hmax) {
            zobj <- .Fortran("cpawsuni", as.integer(n), as.integer(dp1), 
                as.integer(dp2), as.double(x), as.double(y), 
                as.logical(tobj$fix), as.double(tobj$theta), 
                as.double(tobj$bi), bi = as.double(tobj$bi), 
                bi0 = as.double(zobj$bi0), ai = as.double(zobj$ai), 
                as.double(lamakt0), hakt = as.double(hakt), as.integer(lkern), 
                as.double(cb), double(dp1 * dp1), double(dp1 * 
                  dp1), double(dp1 * dp1), double(dp1), double(dp1), 
                double(dp2), double(dp1), as.double(spmax), PACKAGE = "awsnew")[c("ai", 
                "bi", "bi0", "hakt")]
            dim(zobj$ai) <- c(dp1, n)
            dim(zobj$bi) <- c(dp2, n)
            dim(zobj$bi0) <- c(dp2, n)
            if (hakt > n/2) 
                zobj$bi0 <- hincr * biold
            biold <- zobj$bi0
            tobj <- updtheta(zobj, tobj, cpar)
            lamakt0 <- lamakt
            if (graph) {
                par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
                plot(x, y, ylim = range(y, tobj$theta[1, ]), 
                  col = 3)
                if (!is.null(u)) 
                  lines(x, u, col = 2)
                lines(x, tobj$theta[1, ], lwd = 2)
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                plot(x, tobj$bi[1, ], type = "l", ylim = c(0, 
                  max(tobj$bi[1, ])))
                lines(x, tobj$eta * max(tobj$bi[1, ]), col = 2)
                title("Sum of weights and eta")
            }
            if (!is.null(u)) 
                cat("bandwidth: ", signif(hakt, 3), "eta==1", 
                  sum(tobj$eta == 1), "   MSE: ", mean((tobj$theta[1, 
                    ] - u)^2), "   MAE: ", mean(abs(tobj$theta[1, 
                    ] - u)), "\n")
            if (demo) 
                readline("Press return")
            hakt <- hakt * hincr
            gc()
        }
    }
    if (gridded && form == "bi" && p > 0) {
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
                as.double(tobj$bi), bi = as.double(tobj$bi), 
                bi0 = as.double(zobj$bi0), ai = as.double(zobj$ai), 
                as.double(lamakt0), hakt = as.double(hakt), as.integer(lkern), 
                double(dp1 * dp1), double(dp1 * dp1), double(dp1 * 
                  dp1), double(dp1), double(dp1), double(dp2), 
                double(dp2), double(dp2), double(dp1), as.integer(ind), 
                as.double(wghts), as.double(spmax), PACKAGE = "awsnew")[c("ai", 
                "bi", "bi0", "hakt")]
            dim(zobj$ai) <- c(dp1, n1, n2)
            dim(zobj$bi) <- c(dp2, n1, n2)
            dim(zobj$bi0) <- c(dp2, n1, n2)
            if (hakt > min(n1, n2)/2) 
                zobj$bi0 <- hincr * hincr * biold
            biold <- zobj$bi0
            tobj <- updtheta(zobj, tobj, cpar)
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
    if (form == "multi") {
        cpar$mcode <- 3
        dp1 <- 1 + p * px
        dp2 <- dp1 * (dp1 + 1)/2
        bi <- matrix(0, dp2, n)
        theta <- ai <- matrix(0, dp1, n)
        if (NN) {
            if (is.null(hinit) || hinit < (p + 0.1)) 
                hinit <- p + 0.1
            info <- 1
            while (info > 0) {
                ihinit <- trunc(hinit)
                z <- .Fortran("ipawsmnn", as.integer(n), as.integer(px), 
                  as.integer(dp1), as.integer(dp2), as.double(x), 
                  as.double(y), as.integer(neighbors[1:ihinit, 
                    ]), as.integer(ihinit), as.double(hinit), 
                  bi = as.double(bi), ai = as.double(ai), theta = as.double(theta), 
                  as.integer(lkern), double(dp1 * dp1), double(dp1), 
                  info = as.integer(info), PACKAGE = "awsnew")[c("ai", 
                  "bi", "theta", "info")]
                info <- z$info
                hinit <- hinit + 1
            }
            theta <- matrix(z$theta, dp1, n)
            bi <- bi0 <- matrix(z$bi, dp2, n)
            ai <- z$ai
            if (!is.null(u)) 
                cat("bandwidth: ", signif(hinit, 3), "   MSE: ", 
                  mean((theta[1, ] - u)^2), "   MAE: ", mean(abs(theta[1, 
                    ] - u)), "\n")
            hakt <- hinit * hincr
            while (hakt <= hmax) {
                ihakt <- min(ihmax, trunc(hakt))
                z <- .Fortran("cpawsmnn", as.integer(n), as.integer(px), 
                  as.integer(dp1), as.integer(dp2), as.double(x), 
                  as.double(y), as.integer(neighbors[1:ihakt, 
                    ]), as.integer(ihakt), as.double(theta), 
                  as.double(bi), bi = as.double(bi), bi0 = as.double(bi0), 
                  as.double(ai), ai = as.double(ai), as.double(lamakt), 
                  as.double(hakt), as.integer(lkern), double(dp1 * 
                    dp1), double(dp1 * dp1), double(dp1 * dp1), 
                  double(dp1), double(dp1), double(dp1), as.double(spmax), 
                  PACKAGE = "awsnew")[c("ai", "bi", "bi0")]
                ai <- matrix((1 - eta) * z$ai + eta * ai, dp1, 
                  n)
                bi <- matrix((1 - eta) * z$bi + eta * bi, dp2, 
                  n)
                bi0 <- matrix((1 - eta) * z$bi0 + eta * bi0, 
                  dp2, n)
                theta <- gettheta(mcode, ai, bi)
                if (!is.null(u)) 
                  cat("bandwidth: ", signif(hakt, 3), "   MSE: ", 
                    mean((theta[1, ] - u)^2), "   MAE: ", mean(abs(theta[1, 
                      ] - u)), "\n")
                hakt <- hakt * hincr
                gc()
            }
        }
        else {
            dpd <- dp1 + 1
            if (is.null(hinit)) 
                hinit <- maxdist[dpd]
            info <- 1
            while (info > 0) {
                if (hinit <= maxdist[dpd]) 
                  hinit <- maxdist[dpd]
                ihinit <- sum(maxdist <= hinit)
                z <- .Fortran("ipawsmul", as.integer(n), as.integer(px), 
                  as.integer(dp1), as.integer(dp2), as.double(x), 
                  as.double(y), as.integer(neighbors[1:ihinit, 
                    ]), as.double(distmat[1:ihinit, ]), as.integer(ihinit), 
                  as.double(hinit), bi = as.double(bi), ai = as.double(ai), 
                  theta = as.double(theta), as.integer(lkern), 
                  double(dp1 * dp1), double(dp1), info = integer(1), 
                  PACKAGE = "awsnew")[c("ai", "bi", "theta", 
                  "info")]
                info <- z$info
                dpd <- dpd + 1
            }
            theta <- matrix(z$theta, dp1, n)
            bi <- bi0 <- matrix(z$bi, dp2, n)
            ai <- z$ai
            if (!is.null(u)) 
                cat("bandwidth: ", signif(hinit, 3), "   MSE: ", 
                  mean((theta[1, ] - u)^2), "   MAE: ", mean(abs(theta[1, 
                    ] - u)), "\n")
            hakt <- hinit * hincr
            while (hakt <= hmax) {
                ihakt <- sum(maxdist <= hakt)
                z <- .Fortran("cpawsmul", as.integer(n), as.integer(px), 
                  as.integer(dp1), as.integer(dp2), as.double(x), 
                  as.double(y), as.integer(neighbors[1:ihakt, 
                    ]), as.double(distmat[1:ihakt, ]), as.integer(ihakt), 
                  as.double(theta), as.double(bi), bi = as.double(bi), 
                  bi0 = as.double(bi0), as.double(ai), ai = as.double(ai), 
                  as.double(lamakt), as.double(hakt), as.integer(lkern), 
                  double(dp1 * dp1), double(dp1 * dp1), double(dp1 * 
                    dp1), double(dp1), double(dp1), double(dp1), 
                  as.double(spmax), PACKAGE = "awsnew")[c("ai", 
                  "bi", "bi0")]
                ai <- matrix((1 - eta) * z$ai + eta * ai, dp1, 
                  n)
                bi <- matrix((1 - eta) * z$bi + eta * bi, dp2, 
                  n)
                bi0 <- matrix((1 - eta) * z$bi0 + eta * bi0, 
                  dp2, n)
                theta <- gettheta(mcode, ai, bi)
                if (!is.null(u)) 
                  cat("bandwidth: ", signif(hakt, 3), "   MSE: ", 
                    mean((theta[1, ] - u)^2), "   MAE: ", mean(abs(theta[1, 
                      ] - u)), "\n")
                hakt <- hakt * hincr
                gc()
            }
        }
    }
    z <- list(theta = tobj$theta, y = y, x = x, call = args)
    class(z) <- "aws"
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
