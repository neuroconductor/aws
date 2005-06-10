\name{cawsdens}
\alias{cawsdens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cawsdens(y, ngrid = NULL, nempty = NULL, qlambda = NULL, eta0 = NULL, lkern = "Triangle", fu = NULL, hinit = 1, hincr = 1.25, hmax = NULL, heta = NULL, delta1 = NULL, delta2 = NULL, graph = FALSE, demo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{ngrid}{ ~~Describe \code{ngrid} here~~ }
  \item{nempty}{ ~~Describe \code{nempty} here~~ }
  \item{qlambda}{ ~~Describe \code{qlambda} here~~ }
  \item{eta0}{ ~~Describe \code{eta0} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{fu}{ ~~Describe \code{fu} here~~ }
  \item{hinit}{ ~~Describe \code{hinit} here~~ }
  \item{hincr}{ ~~Describe \code{hincr} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{heta}{ ~~Describe \code{heta} here~~ }
  \item{delta1}{ ~~Describe \code{delta1} here~~ }
  \item{delta2}{ ~~Describe \code{delta2} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{demo}{ ~~Describe \code{demo} here~~ }
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
function (y, ngrid = NULL, nempty = NULL, qlambda = NULL, eta0 = NULL, 
    lkern = "Triangle", fu = NULL, hinit = 1, hincr = 1.25, hmax = NULL, 
    heta = NULL, delta1 = NULL, delta2 = NULL, graph = FALSE, 
    demo = FALSE) 
{
    args <- match.call()
    dy <- dim(y)
    u <- NULL
    if (is.null(dy)) {
        n <- length(y)
        if (is.null(ngrid)) 
            ngrid <- as.integer(2 * n)
        if (is.null(nempty)) 
            nempty <- trunc(0.1 * ngrid)
        dd <- 1
        ry <- range(y)
        dry <- diff(ry)
        ry[1] <- ry[1] - dry/ngrid * nempty
        ry[2] <- ry[2] + dry/ngrid * nempty
        dry <- dry * (1 + 2/ngrid * nempty)
        ngrid <- ngrid + 2 * nempty
        bin <- numeric(ngrid)
        ind <- trunc((y - ry[1])/dry * ngrid) + 1
        bin[as.integer(levels(factor(ind)))] <- table(ind)
        xgrid <- seq(ry[1] + dry/2/ngrid, ry[2] + dry/2/ngrid, 
            length = ngrid)
        if (!is.null(fu)) 
            u <- fu(xgrid)
    }
    else {
        dd <- dy[1]
        n <- dy[2]
        if (is.null(ngrid)) 
            ngrid <- as.integer(2 * n^(1/dd))
        if (is.null(nempty)) 
            nempty <- trunc(0.1 * ngrid)
        ry <- apply(y, 1, range)
        dry <- ry[2, ] - ry[1, ]
        ry[1, ] <- ry[1, ] - dry/ngrid * nempty
        ry[2, ] <- ry[2, ] + dry/ngrid * nempty
        dry <- dry * (1 + 2/ngrid * nempty)
        if (length(ngrid) == 1) 
            ngrid <- rep(ngrid, dd)
        if (length(ngrid) != dd) 
            return("incompatible length of ngrid")
        ngrid[1] <- ngrid[1] + 2 * nempty
        ngrid[2] <- ngrid[2] + 2 * nempty
        if (dd > 2) 
            ngrid[3] <- ngrid[3] + 2 * nempty
        if (dd > 3) 
            return("not implemented for more than three dimensions")
        if (dd == 2) 
            bin <- matrix(0, ngrid[1], ngrid[2])
        if (dd > 2) 
            bin <- array(0, ngrid)
        ind <- matrix(0, n, dd)
        for (i in 1:dd) ind[, i] <- trunc((y[i, ] - ry[1, i])/dry[i] * 
            ngrid[i]) + 1
        if (dd == 2) {
            bin[as.integer(levels(factor(ind[, 1]))), as.integer(levels(factor(ind[, 
                2])))] <- table(ind[, 1], ind[, 2])
            xgrid <- list(seq(ry[1, 1] + dry[1]/2/ngrid[1], ry[2, 
                1] + dry[1]/2/ngrid[1], length = ngrid[1]), seq(ry[1, 
                2] + dry[2]/2/ngrid[2], ry[2, 2] + dry[2]/2/ngrid[2], 
                length = ngrid[2]))
            if (!is.null(fu)) 
                u <- fu(xgrid[[1]], xgrid[[2]])
        }
        else {
            bin[as.integer(levels(factor(ind[, 1]))), as.integer(levels(factor(ind[, 
                2]))), as.integer(levels(factor(ind[, 3])))] <- table(ind[, 
                1], ind[, 2], ind[, 3])
            xgrid <- list(seq(ry[1, 1] + dry[1]/2/ngrid[1], ry[2, 
                1] + dry[1]/2/ngrid[1], length = ngrid[1]), seq(ry[1, 
                2] + dry[2]/2/ngrid[2], ry[2, 2] + dry[2]/2/ngrid[2], 
                length = ngrid[2]), seq(ry[1, 3] + dry[3]/2/ngrid[3], 
                ry[2, 3] + dry[3]/2/ngrid[3], length = ngrid[3]))
            if (!is.null(fu)) 
                u <- fu(xgrid[[1]], xgrid[[2]], xgrid[[3]])
        }
    }
    if (!is.null(fu)) 
        u <- u * n * prod(dry/ngrid)
    dens <- caws(bin, qlambda = qlambda, family = "Poisson", 
        eta0 = eta0, lkern = lkern, hinit = hinit, hincr = hincr, 
        hmax = hmax, graph = graph, demo = demo, u = u, heta = heta, 
        delta1 = delta1, delta2 = delta2)$theta
    dens <- dens/sum(dens)/prod(dry/ngrid)
    z <- list(bin = bin, dens = dens, xgrid = xgrid, call = args)
    class(z) <- "psaws.density"
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
