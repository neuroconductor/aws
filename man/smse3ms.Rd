\name{smse3ms}
\alias{smse3}
\alias{smse3ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptive smoothing in orientation space SE(3)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smse3ms(sb, s0, bv, grad, kstar, lambda, kappa0, mask, sigma,
    ns0 = 1, vext = NULL, ncoils = 1, verbose = FALSE, usemaxni = TRUE)
smse3(sb, s0, bv, grad, mask, sigma, kstar, lambda, kappa0,
    ns0 = 1, vext = NULL, vred = 4, ncoils = 1, model = 0, dist = 1,
    verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sb}{
%%     ~~Describe \code{sb} here~~
}
  \item{s0}{
%%     ~~Describe \code{s0} here~~
}
  \item{bv}{
%%     ~~Describe \code{bv} here~~
}
  \item{grad}{
%%     ~~Describe \code{grad} here~~
}
  \item{kstar}{
%%     ~~Describe \code{kstar} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{kappa0}{
%%     ~~Describe \code{kappa0} here~~
}
  \item{mask}{
%%     ~~Describe \code{mask} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{ns0}{
%%     ~~Describe \code{ns0} here~~
}
  \item{vext}{
%%     ~~Describe \code{vext} here~~
}
  \item{ncoils}{
%%     ~~Describe \code{ncoils} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{usemaxni}{
%%     ~~Describe \code{usemaxni} here~~
}
\item{vred}{
%%     ~~Describe \code{vred} here~~
}
\item{model}{
%%     ~~Describe \code{eps} here~~
}
\item{dist}{
%%     ~~Describe \code{eps} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sb, s0, bv, grad, kstar, lambda, kappa0, mask, sigma,
    ns0 = 1, vext = NULL, ncoils = 1, verbose = FALSE, usemaxni = TRUE)
{
    ngrad <- length(bv)
    ddim <- dim(mask)
    msstructure <- getnext3g(grad, bv)
    bv <- msstructure$bv
    nshell <- as.integer(msstructure$nbv)
    ubv <- msstructure$ubv
    nvoxel <- sum(mask)
    nbv <- length(bv)
    if (dim(sb)[1] != nvoxel || length(s0) != nvoxel || dim(sb)[2] !=
        nbv) {
        stop("smse3ms - sb and s0 should only contain data within mask")
    }
    dsigma <- dim(sigma)
    if (is.null(dsigma)) {
        s0 <- s0/sigma
        sb <- sb/sigma
    }
    else if (dsigma[2] == nshell + 1) {
        s0 <- s0/sigma[, 1]
        for (shnr in 1:nshell) {
            indbv <- (1:length(bv))[bv == ubv[shnr]]
            sb[, indbv] <- sb[, indbv]/sigma[, shnr + 1]
        }
    }
    else if (dsigma[2] == nbv + 1) {
        s0 <- s0/sigma[, 1]
        for (ibv in 1:nbv) {
            sb[, ibv] <- sb[, ibv]/sigma[, ibv + 1]
        }
    }
    else {
        warning("incompatible number of sigma images, needs to be 1, nshell+1 or nbv+1\n")
    }
    position <- array(0, ddim)
    position[mask] <- 1:nvoxel
    if (is.null(vext))
        vext <- c(1, 1)
    varstats <- sofmchi(ncoils)
    z <- list(th = array(1, c(nvoxel, nbv)), th0 = rep(1, nvoxel),
        ni = array(1, c(nvoxel, nbv)), ni0 = rep(1, nvoxel))
    gradstats <- getkappasmsh3(grad, msstructure)
    hseq <- gethseqfullse3msh(kstar, gradstats, kappa0, vext = vext)
    nind <- as.integer(hseq$n * 1.25)
    if (usemaxni) {
        ni <- array(1, c(nvoxel, nbv))
        ni0 <- rep(1, nvoxel)
    }
    prt0 <- Sys.time()
    cat("adaptive smoothing in SE3, kstar=", kstar, if (verbose)
        "\n"
    else " ")
    kinit <- if (lambda < 1e+10)
        0
    else kstar
    mc.cores <- setCores(, reprt = FALSE)
    gc()
    for (k in kinit:kstar) {
        hakt <- hseq$h[, k + 1]
        t0 <- Sys.time()
        thnimsh <- interpolatesphere1(z$th, z$th0, z$ni, z$ni0,
            msstructure, mask)
        rm(z)
        gc()
        t1 <- Sys.time()
        param <- lkfullse3msh(hakt, kappa0/hakt, gradstats, vext,
            nind)
        hakt0 <- mean(hakt)
        param0 <- lkfulls0(hakt0, vext, nind)
        vs2 <- varstats$s2[findInterval(thnimsh$mstheta, varstats$mu,
            all.inside = TRUE)]/2
        vs02 <- varstats$s2[findInterval(thnimsh$msth0, varstats$mu,
            all.inside = TRUE)]/2
        t2 <- Sys.time()
        z <- .Fortran(C_adsmse3s, as.double(sb), as.double(s0),
            as.double(thnimsh$mstheta), as.double(thnimsh$msni),
            as.double(thnimsh$msth0), as.double(thnimsh$msni0),
            as.double(vs2), as.double(vs02), as.integer(position),
            as.integer(nvoxel), as.integer(nshell + 1), as.integer(ddim[1]),
            as.integer(ddim[2]), as.integer(ddim[3]), as.integer(ngrad),
            as.double(lambda), as.double(ns0), as.integer(param$ind),
            as.double(param$w), as.integer(param$n), as.integer(param0$ind),
            as.double(param0$w), as.integer(param0$n), th = double(nvoxel *
                nbv), ni = double(nvoxel * nbv), th0 = double(nvoxel),
            ni0 = double(nvoxel), double(ngrad * mc.cores), double(ngrad *
                mc.cores), double((nshell + 1) * mc.cores), double((nshell +
                1) * mc.cores), double((nshell + 1) * mc.cores))[c("ni",
            "th", "ni0", "th0")]
        t3 <- Sys.time()
        gc()
        rm(thnimsh, vs2, vs02)
        gc()
        if (usemaxni) {
            ni <- z$ni <- if (usemaxni)
                pmax(ni, z$ni)
            ni0 <- z$ni0 <- if (usemaxni)
                pmax(ni0, z$ni0)
        }
        dim(z$ni) <- dim(z$th) <- c(nvoxel, nbv)
        if (is.null(dsigma)) {
            z$th0 <- z$th0 * sigma/sqrt(ns0)
            z$th <- z$th * sigma
        }
        else if (dsigma[2] == nshell + 1) {
            for (shnr in 1:nshell) {
                z$th0 <- z$th0 * sigma[, 1]/sqrt(ns0)
                indbv <- (1:length(bv))[bv == ubv[shnr]]
                z$th[, indbv] <- z$th[, indbv] * sigma[, shnr +
                  1]
            }
        }
        else if (dsigma[2] == nbv + 1) {
            for (ibv in 1:nbv) {
                z$th0 <- z$th0 * sigma[, 1]/sqrt(ns0)
                z$th[, ibv] <- z$th[, ibv] * sigma[, ibv + 1]
            }
        }
        if (verbose) {
            cat("k:", k, "h_k:", signif(max(hakt), 3), " quartiles of ni",
                signif(quantile(z$ni), 3), "mean of ni", signif(mean(z$ni),
                  3), "\n              quartiles of ni0", signif(quantile(z$ni0),
                  3), "mean of ni0", signif(mean(z$ni0), 3),
                " time elapsed:", format(difftime(Sys.time(),
                  prt0), digits = 3), "\n")
            cat("interpolation:", format(difftime(t1, t0), digits = 3),
                "param:", format(difftime(t2, t1), digits = 3),
                "smoothing:", format(difftime(t3, t2), digits = 3),
                "\n")
        }
        else {
            cat(".")
        }
    }
    list(th = z$th, th0 = z$th0, ni = z$ni, ni0 = z$ni0, hseq = hseq,
        kappa0 = kappa0, lambda = lambda)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
