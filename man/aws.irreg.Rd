\name{aws.irreg}
\alias{aws.irreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aws.irreg(y, x, d = 2, hmax = NULL, hpre = NULL, qlambda = NULL, qtau = NULL, varmodel = "Constant", sigma2 = NULL, varprop = 0.1, graph = FALSE, lkern = "Triangle", skern = "Triangle", aggkern = "Uniform", spmin = 0, spmax = 5, lseq = NULL, nbins = 100, henv = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{hpre}{ ~~Describe \code{hpre} here~~ }
  \item{qlambda}{ ~~Describe \code{qlambda} here~~ }
  \item{qtau}{ ~~Describe \code{qtau} here~~ }
  \item{varmodel}{ ~~Describe \code{varmodel} here~~ }
  \item{sigma2}{ ~~Describe \code{sigma2} here~~ }
  \item{varprop}{ ~~Describe \code{varprop} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{skern}{ ~~Describe \code{skern} here~~ }
  \item{aggkern}{ ~~Describe \code{aggkern} here~~ }
  \item{spmin}{ ~~Describe \code{spmin} here~~ }
  \item{spmax}{ ~~Describe \code{spmax} here~~ }
  \item{lseq}{ ~~Describe \code{lseq} here~~ }
  \item{nbins}{ ~~Describe \code{nbins} here~~ }
  \item{henv}{ ~~Describe \code{henv} here~~ }
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
function (y, x, d = 2, hmax = NULL, hpre = NULL, qlambda = NULL, 
    qtau = NULL, varmodel = "Constant", sigma2 = NULL, varprop = 0.1, 
    graph = FALSE, lkern = "Triangle", skern = "Triangle", aggkern = "Uniform", 
    spmin = 0, spmax = 5, lseq = NULL, nbins = 100, henv = NULL) 
{
    args <- match.call()
    n <- length(y)
    dx <- dim(x)
    if (!(d \%in\% 1:2)) 
        stop("this version is for 1D and 2D only")
    if ((d == 1 && length(x) != length(y)) || (d == 2 && (is.null(dx) || 
        dx[1] != n || dx[2] != 2))) 
        stop("incorrect size of x")
    if (!(varmodel \%in\% c("Constant", "Linear", "Quadratic"))) 
        stop("Model for variance not implemented")
    require(sm)
    zbins <- binning(x, y, nbins = n^(1/d)/2)
    given.var <- !is.null(sigma2)
    if (!given.var) {
        sigma20 <- mean(zbins$devs[zbins$x.freq > 1]/(zbins$x.freq[zbins$x.freq > 
            1] - 1))
        cat("Preliminary variance estimate:", sigma20, "\n")
    }
    else {
        coef <- sigma2[1]
    }
    zbins <- binning(x, y, nbins = nbins)
    ni <- t(zbins$table.freq)
    mask <- ni > 0
    if (!is.null(henv)) 
        mask <- .Fortran("mask", as.logical(mask), mask = as.logical(mask), 
            as.integer(nbins), as.integer(switch(d, 1, nbins)), 
            as.integer(max(0, henv)), PACKAGE = "aws", DUP = FALSE)$mask
    yy <- rep(mean(y), length(mask))
    dim(yy) <- dim(mask) <- dim(ni)
    yy[ni > 0] <- zbins$means
    nn <- length(yy)
    if (given.var) {
        if (length(sigma2) != nn) 
            sigma2 <- rep(sigma2[1], nn)
        sigma2 <- 1/sigma2
    }
    else {
        sigma2 <- 1/rep(sigma20, nn)
    }
    if (d == 2) 
        wghts <- diff(range(x[, 1]))/diff(range(x[, 2]))
    else wghts <- 1
    if (d == 2) 
        dy <- dim(yy) <- dim(sigma2) <- c(nbins, nbins)
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, Gaussian = 5, 2)
    skern <- switch(skern, Exp = 1, Triangle = 2, 2)
    cpar <- setawsdefaults(dim(yy), mean(y), "Gaussian", skern, 
        aggkern, qlambda, qtau, lseq, hmax, 1, spmax)
    lambda <- 2 * cpar$lambda
    cpar$tau1 <- cpar$tau1 * 2
    cpar$tau2 <- cpar$tau2 * 2
    hmax <- cpar$hmax
    lseq <- cpar$lseq
    shape <- cpar$shape
    cpar$heta <- 20^(1/d)
    hinit <- cpar$hinit
    hincr <- cpar$hincr
    spmax <- cpar$spmax
    cpar$heta <- 1e+10
    if (lkern == 5) {
        hmax <- hmax * 0.42445 * 4
        hinit <- 0.42445 * 4
    }
    n1 <- nbins
    n2 <- switch(d, 1, nbins)
    dy <- switch(d, NULL, c(n1, n2))
    tobj <- list(bi = ni, bi2 = ni^2, theta = yy/shape, fix = rep(FALSE, 
        nn))
    zobj <- list(ai = yy, bi0 = rep(1, nn))
    biold <- ni
    vred <- ni
    hakt <- hinit * hincr
    hakt0 <- hinit * hincr
    lambda0 <- lambda
    lambda0 <- 1e+50
    if (is.null(hpre)) 
        hpre <- (20 * nn/n)^(1/d)
    dlw <- (2 * trunc(hpre/c(1, wghts)) + 1)[1:d]
    hobj <- .Fortran("cawsmask", as.double(yy), as.logical(ni > 
        0), as.integer(ni), as.logical(tobj$fix), as.integer(n1), 
        as.integer(n2), hakt = as.double(hpre), as.double(1e+40), 
        as.double(tobj$theta), bi = as.double(tobj$bi), bi2 = double(nn), 
        bi0 = as.double(zobj$bi0), ai = as.double(zobj$ai), as.integer(cpar$mcode), 
        as.integer(lkern), as.integer(skern), as.double(spmin), 
        as.double(spmax), double(prod(dlw)), as.double(wghts), 
        PACKAGE = "aws", DUP = FALSE)[c("bi", "ai")]
    hobj$theta <- hobj$ai/hobj$bi
    hobj$theta[ni == 0] <- mean(hobj$theta[ni > 0])
    dim(hobj$theta) <- dim(hobj$bi) <- dy
    steps <- as.integer(log(hmax/hinit)/log(hincr))
    cat("Progress:")
    for (k in 1:steps) {
        dlw <- (2 * trunc(hakt/c(1, wghts)) + 1)[1:d]
        hakt0 <- hakt
        zobj <- .Fortran("cgawsmas", as.double(yy), as.logical(mask), 
            as.integer(ni), as.logical(tobj$fix), as.double(sigma2), 
            as.integer(n1), as.integer(n2), hakt = as.double(hakt), 
            as.double(lambda0), as.double(tobj$theta), bi = as.double(tobj$bi), 
            bi2 = double(nn), bi0 = as.double(zobj$bi0), vred = double(nn), 
            ai = as.double(zobj$ai), as.integer(cpar$mcode), 
            as.integer(lkern), as.integer(skern), as.double(spmin), 
            as.double(spmax), double(prod(dlw)), as.double(wghts), 
            PACKAGE = "aws", DUP = FALSE)[c("bi", "bi0", "bi2", 
            "vred", "ai", "hakt")]
        vred[!tobj$fix] <- zobj$vred[!tobj$fix]
        dim(zobj$ai) <- dy
        if (hakt > n1/2) 
            zobj$bi0 <- hincr^d * biold
        biold <- zobj$bi0
        tobj <- updtheta(zobj, tobj, cpar)
        tobj$vred <- vred
        tobj$theta[tobj$bi == 0] <- mean(tobj$theta[ni > 0])
        dim(tobj$vred) <- dy
        dim(tobj$theta) <- dy
        dim(tobj$bi) <- dy
        dim(tobj$eta) <- dy
        if (graph) {
            if (d == 1) {
                oldpar <- par(mfrow = c(1, 2), mar = c(3, 3, 
                  3, 0.2), mgp = c(2, 1, 0))
                plot((1:nn)[ni > 0], yy[ni > 0], ylim = range(yy, 
                  tobj$theta[mask]), col = 3)
                points((1:nn)[mask & ni == 0], yy[mask & ni == 
                  0], col = 4)
                lines((1:nn)[mask], tobj$theta[mask], lwd = 2)
                title(paste("Reconstruction  h=", signif(hakt, 
                  3)))
                plot((1:nn), tobj$bi, type = "l", ylim = range(0, 
                  tobj$bi))
                points((1:nn)[ni > 0], max(tobj$bi)/max(ni) * 
                  ni[ni > 0], col = 3)
                points((1:nn)[mask], rep(0, sum(mask)), col = 4)
                title("Sum of weights, ni and mask")
            }
            if (d == 2) {
                oldpar <- par(mfrow = c(2, 2), mar = c(1, 1, 
                  3, 0.25), mgp = c(2, 1, 0))
                image(yy, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title(paste("Observed Image  min=", signif(min(yy[mask]), 
                  3), " max=", signif(max(yy[mask]), 3)))
                zlim <- quantile(tobj$theta, c(0.001, 0.999))
                image(array(pmax(pmin(tobj$theta, zlim[2]), zlim[1]), 
                  dy), col = gray((0:255)/255), xaxt = "n", yaxt = "n")
                title(paste("Reconstruction  h=", signif(hakt, 
                  3), " min=", signif(min(tobj$theta[mask]), 
                  3), " max=", signif(max(tobj$theta[mask]), 
                  3)))
                image(tobj$bi, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title(paste("Sum of weights: min=", signif(min(tobj$bi[mask]), 
                  3), " mean=", signif(mean(tobj$bi[mask]), 3), 
                  " max=", signif(max(tobj$bi), 3)))
                image(mask, col = gray((0:255)/255), xaxt = "n", 
                  yaxt = "n")
                title("mask")
            }
            par(oldpar)
        }
        if (!given.var) {
            sigma2 <- awsisigma2(yy, hobj, tobj, ni, sigma20, 
                varmodel, varprop)
        }
        hakt <- hakt * hincr
        x <- 1.25^(k - 1)
        lambda0 <- lambda * lseq[k]
        cat(paste(signif(sum(hincr^(2 * (1:k)))/sum(hincr^(2 * 
            (1:steps))) * 100, 2), "\% ", sep = ""))
        gc()
    }
    cat("\n")
    if (length(sigma2) == nn) {
        vartheta <- tobj$bi2/tobj$bi^2
    }
    else {
        vartheta <- sigma2 * tobj$bi2/tobj$bi^2
        vred <- tobj$bi2/tobj$bi^2
    }
    z <- list(theta = tobj$theta, sigma2 = 1/sigma2, bi = tobj$bi, 
        var = vartheta, vred = vred, y = yy, ni = ni, varcoef = coef, 
        hmax = hakt/hincr, lseq = c(0, lseq), call = args, zbins = zbins, 
        x = x)
    class(z) <- "aws.gaussian"
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
