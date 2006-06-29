\name{aws.gaussian}
\alias{aws.gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aws.gaussian(y, hmax = NULL, hpre = NULL, qlambda = NULL, qtau = NULL, varmodel = "Constant", varpar = NULL, varprop = 0.1, scorr = 0, wghts = NULL, graph = FALSE, demo = FALSE, lkern = "Triangle", skern = "Triangle", aggkern = "Uniform", spmin = 0, spmax = 5, lseq = NULL, u = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{hpre}{ ~~Describe \code{hpre} here~~ }
  \item{qlambda}{ ~~Describe \code{qlambda} here~~ }
  \item{qtau}{ ~~Describe \code{qtau} here~~ }
  \item{varmodel}{ ~~Describe \code{varmodel} here~~ }
  \item{varpar}{ ~~Describe \code{varpar} here~~ }
  \item{varprop}{ ~~Describe \code{varprop} here~~ }
  \item{scorr}{ ~~Describe \code{scorr} here~~ }
  \item{wghts}{ ~~Describe \code{wghts} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{demo}{ ~~Describe \code{demo} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{skern}{ ~~Describe \code{skern} here~~ }
  \item{aggkern}{ ~~Describe \code{aggkern} here~~ }
  \item{spmin}{ ~~Describe \code{spmin} here~~ }
  \item{spmax}{ ~~Describe \code{spmax} here~~ }
  \item{lseq}{ ~~Describe \code{lseq} here~~ }
  \item{u}{ ~~Describe \code{u} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
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
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, hmax = NULL, hpre = NULL, qlambda = NULL, qtau = NULL, 
    varmodel = "Constant", varpar = NULL, varprop = 0.1, scorr = 0, 
    wghts = NULL, graph = FALSE, demo = FALSE, lkern = "Triangle", 
    skern = "Triangle", aggkern = "Uniform", spmin = 0, spmax = 5, 
    lseq = NULL, u = NULL) 
{
    args <- match.call()
    dy <- dim(y)
    if (length(dy) > 3) 
        stop("AWS for more than 3 dimensional grids is not implemented")
    if (!(varmodel \%in\% c("Constant", "Linear", "Quadratic"))) 
        stop("Model for variance not implemented")
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, Gaussian = 5, 2)
    skern <- switch(skern, Exp = 1, Triangle = 2, 2)
    cpar <- setawsdefaults(dy, mean(y), "Gaussian", skern, aggkern, 
        qlambda, qtau, lseq, hmax, 1, spmax)
    lambda <- cpar$lambda
    hmax <- cpar$hmax
    lseq <- cpar$lseq
    shape <- cpar$shape
    d <- cpar$d
    cpar$heta <- 20^(1/d)
    hinit <- cpar$hinit
    hincr <- cpar$hincr
    spmax <- cpar$spmax
    n <- length(y)
    zfamily <- awsgfamily(y, scorr, lambda, cpar)
    cpar <- zfamily$cpar
    lambda <- zfamily$lambda
    sigma2 <- zfamily$sigma2
    h0 <- zfamily$h0
    y <- zfamily$y
    rm(zfamily)
    if (lkern == 5) {
        hmax <- hmax * 0.42445 * 4
        hinit <- 0.42445 * 4
    }
    if (demo && !graph) 
        graph <- TRUE
    n <- length(y)
    n1 <- switch(d, n, dy[1], dy[1])
    n2 <- switch(d, 1, dy[2], dy[2])
    n3 <- switch(d, 1, 1, dy[3])
    if (is.null(wghts)) 
        wghts <- c(1, 1, 1)
    hinit <- hinit/wghts[1]
    hmax <- hmax/wghts[1]
    wghts <- (wghts[2:3]/wghts[1])
    tobj <- list(bi = rep(1, n), bi2 = rep(1, n), theta = y/shape, 
        fix = rep(FALSE, n))
    zobj <- list(ai = y, bi0 = rep(1, n))
    biold <- rep(1, n)
    vred <- rep(1, n)
    mae <- NULL
    hakt <- hinit * hincr
    hakt0 <- hinit * hincr
    lambda0 <- lambda
    lambda0 <- 1e+50
    if (is.null(hpre)) 
        hpre <- 20^(1/d)
    dlw <- (2 * trunc(hpre/c(1, wghts)) + 1)[1:d]
    hobj <- .Fortran("caws", as.double(y), as.logical(tobj$fix), 
        as.integer(n1), as.integer(n2), as.integer(n3), hakt = as.double(hpre), 
        as.double(1e+40), as.double(tobj$theta), bi = as.double(tobj$bi), 
        bi2 = double(n), bi0 = as.double(zobj$bi0), ai = as.double(zobj$ai), 
        as.integer(cpar$mcode), as.integer(lkern), as.integer(skern), 
        as.double(spmin), as.double(spmax), double(prod(dlw)), 
        as.double(wghts), PACKAGE = "aws", DUP = FALSE)[c("bi", 
        "ai")]
    hobj$theta <- hobj$ai/hobj$bi
    dim(hobj$theta) <- dim(hobj$bi) <- dy
    steps <- as.integer(log(hmax/hinit)/log(hincr))
    cat("Progress:")
    for (k in 1:steps) {
        dlw <- (2 * trunc(hakt/c(1, wghts)) + 1)[1:d]
        if (scorr[1] >= 0.1) 
            lambda0 <- lambda0 * Spatialvar.gauss(hakt0/0.42445/4, 
                h0, d)/Spatialvar.gauss(hakt0/0.42445/4, 1e-05, 
                d)
        hakt0 <- hakt
        zobj <- .Fortran("cgaws", as.double(y), as.logical(tobj$fix), 
            as.double(sigma2), as.integer(n1), as.integer(n2), 
            as.integer(n3), hakt = as.double(hakt), as.double(lambda0), 
            as.double(tobj$theta), bi = as.double(tobj$bi), bi2 = double(n), 
            bi0 = as.double(zobj$bi0), gi = double(n), vred = double(n), 
            ai = as.double(zobj$ai), as.integer(cpar$mcode), 
            as.integer(lkern), as.integer(skern), as.double(spmin), 
            as.double(spmax), double(prod(dlw)), as.double(wghts), 
            PACKAGE = "aws", DUP = FALSE)[c("bi", "bi0", "bi2", 
            "vred", "ai", "gi", "hakt")]
        vred[!tobj$fix] <- zobj$vred[!tobj$fix]
        dim(zobj$ai) <- dy
        if (hakt > n1/2) 
            zobj$bi0 <- hincr^d * biold
        biold <- zobj$bi0
        tobj <- updtheta(zobj, tobj, cpar)
        tobj$gi <- zobj$gi
        dim(tobj$theta) <- dy
        dim(tobj$bi) <- dy
        dim(tobj$eta) <- dy
        if (graph) {
            if (d == 1) {
                oldpar <- par(mfrow = c(1, 2), mar = c(3, 3, 
                  3, 0.2), mgp = c(2, 1, 0))
                plot(y, ylim = range(y, tobj$theta), col = 3)
                if (!is.null(u)) 
                  lines(u, col = 2)
                lines(tobj$theta, lwd = 2)
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                plot(tobj$bi, type = "l", ylim = range(0, tobj$bi))
                lines(tobj$eta * max(tobj$bi), col = 2)
                title("Sum of weights and eta")
            }
            if (d == 2) {
                oldpar <- par(mfrow = c(2, 2), mar = c(1, 1, 
                  3, 0.25), mgp = c(2, 1, 0))
                image(y, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title(paste("Observed Image  min=", signif(min(y), 
                  3), " max=", signif(max(y), 3)))
                zlim <- quantile(tobj$theta, c(0.001, 0.999))
                image(array(pmax(pmin(tobj$theta, zlim[2]), zlim[1]), 
                  dy), col = gray((0:255)/255), xaxt = "n", yaxt = "n")
                title(paste("Reconstruction  h=", signif(hakt, 
                  3), " min=", signif(min(tobj$theta), 3), " max=", 
                  signif(max(tobj$theta), 3)))
                image(tobj$bi, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title(paste("Sum of weights: min=", signif(min(tobj$bi), 
                  3), " mean=", signif(mean(tobj$bi), 3), " max=", 
                  signif(max(tobj$bi), 3)))
                image(tobj$eta, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n", zlim = c(0, 1))
                title("eta")
            }
            if (d == 3) {
                oldpar <- par(mfrow = c(2, 2), mar = c(1, 1, 
                  3, 0.25), mgp = c(2, 1, 0))
                image(y[, , n3\%/\%2 + 1], col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Observed Image  min=", signif(min(y), 
                  3), " max=", signif(max(y), 3)))
                zlim <- quantile(tobj$theta, c(0.001, 0.999))
                image(array(pmax(pmin(tobj$theta[, , n3\%/\%2 + 
                  1], zlim[2]), zlim[1]), dy[-3]), col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Reconstruction  h=", signif(hakt, 
                  3), " min=", signif(min(tobj$theta), 3), " max=", 
                  signif(max(tobj$theta), 3)))
                image(tobj$bi[, , n3\%/\%2 + 1], col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n")
                title(paste("Sum of weights: min=", signif(min(tobj$bi), 
                  3), " mean=", signif(mean(tobj$bi), 3), " max=", 
                  signif(max(tobj$bi), 3)))
                image(tobj$eta[, , n3\%/\%2 + 1], col = gray((0:255)/255), 
                  xaxt = "n", yaxt = "n", zlim = c(0, 1))
                title("eta")
            }
            par(oldpar)
        }
        if (!is.null(u)) {
            cat("bandwidth: ", signif(hakt, 3), "eta==1", sum(tobj$eta == 
                1), "   MSE: ", signif(mean((tobj$theta - u)^2), 
                3), "   MAE: ", signif(mean(abs(tobj$theta - 
                u)), 3), " mean(bi)=", signif(mean(tobj$bi), 
                3), "\n")
            mae <- c(mae, signif(mean(abs(tobj$theta - u)), 3))
        }
        if (demo) 
            readline("Press return")
        vobj <- awsgsigma2(y, hobj, tobj, varmodel, varpar, varprop, 
            h0)
        sigma2 <- vobj$sigma2inv
        coef <- vobj$coef
        rm(vobj)
        hakt <- hakt * hincr
        x <- 1.25^(k - 1)
        scorrfactor <- x/(3^d * prod(scorr) * prod(h0) + x)
        lambda0 <- lambda * lseq[k] * scorrfactor
        cat(paste(signif(sum(hincr^(2 * (1:k)))/sum(hincr^(2 * 
            (1:steps))) * 100, 2), "\% ", sep = ""))
        gc()
    }
    cat("\n")
    if (length(sigma2) == n) {
        vartheta <- tobj$bi2/tobj$bi^2
    }
    else {
        vartheta <- sigma2 * tobj$bi2/tobj$bi^2
        vred <- tobj$bi2/tobj$bi^2
    }
    vartheta <- vartheta/Spatialvar.gauss(hakt/0.42445/4, h0 + 
        1e-05, d) * Spatialvar.gauss(hakt/0.42445/4, 1e-05, d)
    z <- list(theta = tobj$theta, sigma2 = 1/sigma2, ni = tobj$bi, 
        var = vartheta, vred = vred, y = y, varcoef = coef, hmax = hakt/hincr, 
        mae = mae, lseq = c(0, lseq), call = args)
    class(z) <- "aws.gaussian"
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
