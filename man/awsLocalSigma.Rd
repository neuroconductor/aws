\name{awsLocalSigma}
\alias{awsLocalSigma}
\alias{awslinsd}
\alias{estimateSigmaCompl}
\alias{AFLocalSigma}
\alias{estGlobalSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D variance estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
awsLocalSigma(y, steps, mask, ncoils, vext = c(1, 1), lambda = 5,
    minni = 2, hsig = 5, sigma = NULL, family = c("NCchi", "Gauss"),
    verbose = NULL, trace = FALSE, u = NULL)
awslinsd(y, hmax = NULL, hpre = NULL, h0 = NULL, mask = NULL,
    ladjust = 1, wghts = NULL, varprop = 0.1, A0, A1)
AFLocalSigma(y, ncoils, level = NULL, mask = NULL, h = 2, hadj = 1,
    vext = c(1, 1))
estGlobalSigma(y, mask = NULL, ncoils = 1, steps = 16, vext = c(1, 1),
    lambda = 20, h0 = 2, hadj = 1, q = 0.25, qni = 0.8, level = NULL,
    sequence = FALSE, method = c("awsVar", "awsMAD", "AFmodevn",
                "AFmodem1chi", "AFbkm2chi", "AFbkm1chi"))
estimateSigmaCompl(magnitude, phase, mask, kstar = 20, kmin = 8, hsig = 5,
        lambda = 12, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{steps}{
%%     ~~Describe \code{steps} here~~
}
  \item{mask}{
%%     ~~Describe \code{mask} here~~
}
  \item{ncoils}{
%%     ~~Describe \code{ncoils} here~~
}
  \item{vext}{
%%     ~~Describe \code{vext} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{minni}{
%%     ~~Describe \code{minni} here~~
}
  \item{hsig}{
%%     ~~Describe \code{hsig} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
}
\item{hmax}{
%%     ~~Describe \code{eps} here~~
}
\item{hpre}{
%%     ~~Describe \code{eps} here~~
}
\item{h0}{
%%     ~~Describe \code{eps} here~~
}
\item{ladjust}{
%%     ~~Describe \code{eps} here~~
}
\item{wghts}{
%%     ~~Describe \code{eps} here~~
}
\item{varprop}{
%%     ~~Describe \code{eps} here~~
}
\item{A0}{
%%     ~~Describe \code{eps} here~~
}
\item{A1}{
%%     ~~Describe \code{eps} here~~
}
\item{level}{
%%     ~~Describe \code{eps} here~~
}
\item{h}{
%%     ~~Describe \code{eps} here~~
}
\item{hadj}{
%%     ~~Describe \code{eps} here~~
}
\item{q}{
%%     ~~Describe \code{eps} here~~
}
\item{qni}{
%%     ~~Describe \code{eps} here~~
}
\item{sequence}{
%%     ~~Describe \code{eps} here~~
}
\item{method}{
%%     ~~Describe \code{eps} here~~
}
\item{magnitude}{
%%     ~~Describe \code{eps} here~~
}
\item{phase}{
%%     ~~Describe \code{eps} here~~
}
\item{kstar}{
%%     ~~Describe \code{eps} here~~
}
\item{kmin}{
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
function (y, steps, mask, ncoils, vext = c(1, 1), lambda = 5,
    minni = 2, hsig = 5, sigma = NULL, family = c("NCchi", "Gauss"),
    verbose = NULL, trace = FALSE, u = NULL)
{
    if (trace)
        tergs <- array(0, c(steps, 4, sum(mask)))
    else tergs <- NULL
    family <- match.arg(family, c("NCchi", "Gauss", "Gaussian"))
    if (family == "NCchi") {
        varstats <- sofmchi(ncoils)
        th <- seq(0, 30, 0.01)
        z <- fncchiv(th, varstats)
        minz <- min(z)
        th <- th[z > min(z)]
        minth <- min(th)
        z <- z[z > min(z)]
        nz <- nls(z ~ (a * th + b * th * th + c * th * th * th +
            d)/(a * th + b * th * th + c * th * th * th + 1 +
            d), data = list(th = th, z = z), start = list(a = 1,
            b = 1, c = 1, d = 1))
        vpar <- c(minth, minz, coef(nz))
    }
    if (length(vext) == 3)
        vext <- vext[2:3]/vext[1]
    ddim <- dim(y)
    n <- prod(ddim)
    if (length(ddim) != 3)
        stop("first argument should be a 3-dimentional array")
    if (is.null(mask))
        mask <- array(TRUE, ddim)
    if (length(mask) != n)
        stop("dimensions of data array and mask should coincide")
    if (is.null(sigma)) {
        sigma <- IQQdiff(y, mask, 0.25, verbose = verbose)
        if ("NCchi" == family) {
            sigma <- iniSigmaNCchi(y, mask, 0.25, ncoils, sigma)
            sigma <- iniSigmaNCchi(y, mask, 0.25, ncoils, sigma)
            sigma <- iniSigmaNCchi(y, mask, 0.25, ncoils, sigma)
        }
    }
    if (verbose) {
        mslice <- (ddim[3] + 1)/2
        ymslice <- y[, , mslice]
        ymslice[!mask[, , mslice]] <- 0
        if (!is.null(u) && "NCchi" == family) {
            par(mfrow = c(2, 4), mar = c(3, 3, 3, 1), mgp = c(2,
                1, 0))
        }
        else {
            par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), mgp = c(2,
                1, 0))
        }
    }
    else {
        cat("step")
    }
    th <- y
    ksi <- array(y^2, ddim)
    ni <- array(1, ddim)
    sigma <- array(sigma, ddim)
    mc.cores <- setCores(, reprt = FALSE)
    parammd <- getparam3d(hsig, vext)
    if (family == "NCchi") {
        nfb <- 200 * ncoils
        x <- 1:nfb
        flb <- x + log(besselI(x, ncoils - 1, TRUE))
    }
    for (i in 1:steps) {
        h <- 1.25^((i - 1)/3)
        param <- getparam3d(h, vext)
        nw <- length(param$w)
        if (family == "NCchi") {
            z <- .Fortran(C_awslchi2, as.double(y), as.double(ksi),
                ni = as.double(ni), as.double(sigma), as.double(vpar),
                as.double(ncoils), as.integer(mask), as.integer(ddim[1]),
                as.integer(ddim[2]), as.integer(ddim[3]), as.integer(param$ind),
                as.double(param$w), as.integer(nw), as.double(minni),
                double(nw * mc.cores), double(nw * mc.cores),
                as.double(lambda), as.integer(mc.cores), as.integer(floor(ncoils)),
                double(floor(ncoils) * mc.cores), th = double(n),
                sigman = double(n), ksi = double(n), as.double(flb),
                as.integer(nfb))[c("ni", "ksi", "th", "sigman")]
            thchi <- z$th
            ksi <- z$ksi
            thchi[!mask] <- 0
        }
        else {
            z <- .Fortran(C_awslgaus, as.double(y), as.double(th),
                ni = as.double(ni), as.double(sigma), as.integer(mask),
                as.integer(ddim[1]), as.integer(ddim[2]), as.integer(ddim[3]),
                as.integer(param$ind), as.double(param$w), as.integer(nw),
                as.double(minni), as.double(lambda), th = double(n),
                sigman = double(n))[c("ni", "th", "sigman")]
        }
        th <- array(z$th, ddim)
        ni <- array(z$ni, ddim)
        z$sigman[z$sigman == 0] <- median(z$sigman[z$sigman >
            0])
        nmask <- sum(mask)
        if (verbose)
            cat("local estimation in step ", i, " completed",
                format(Sys.time()), "\n")
        if (any(ni[mask] > minni)) {
            z$sigman[mask] <- z$sigman[mask] * (1 + runif(nmask,
                -1e-04, 1e-04))
            nwmd <- length(parammd$w)
            sigma <- .Fortran(C_mediansm, as.double(z$sigman),
                as.integer(mask), as.integer(ddim[1]), as.integer(ddim[2]),
                as.integer(ddim[3]), as.integer(parammd$ind),
                as.integer(nwmd), double(nwmd * mc.cores), as.integer(mc.cores),
                sigman = double(n))$sigman
        }
        dim(sigma) <- ddim
        mask[sigma == 0] <- FALSE
        if (verbose)
            cat("local median smoother in step ", i, " completed",
                format(Sys.time()), "\n")
        if (verbose) {
            meds <- median(sigma[mask])
            means <- mean(sigma[mask])
            image(ymslice, col = grey(0:255/255))
            title(paste("S  max=", signif(max(y[mask]), 3), " median=",
                signif(median(y[mask]), 3)))
            image(th[, , mslice], col = grey(0:255/255))
            title(paste("E(S)  max=", signif(max(th[mask]), 3),
                " median=", signif(median(th[mask]), 3)))
            image(sigma[, , mslice], col = grey(0:255/255), zlim = c(0,
                max(sigma[mask])))
            title(paste("sigma max=", signif(max(sigma[mask]),
                3), " median=", signif(meds, 3)))
            image(ni[, , mslice], col = grey(0:255/255))
            title(paste("Ni    max=", signif(max(ni[mask]), 3),
                " median=", signif(median(ni[mask]), 3)))
            plot(density(sigma[mask]), main = "density of sigma")
            plot(density(ni[mask]), main = "density of Ni")
            cat("mean sigma", means, "median sigma", meds, "sd sigma",
                sd(sigma[mask]), "\n")
            if (!is.null(u) && "NCchi" == family) {
                thchims <- fncchir(th/sigma, varstats) * sigma
                thchims[!mask] <- u[!mask]
                image(abs(thchims - u)[, , mslice], col = grey(0:255/255))
                title("abs Error in thchims")
                plot(density((thchims - u)[mask]), main = "density of thchims-u")
                cat("MAE(th)", mean(abs(thchi - u)[mask]), "RMSE(th)",
                  sqrt(mean((thchi - u)[mask]^2)), "MAE(thms)",
                  mean(abs(thchims - u)[mask]), "RMSE(thms)",
                  sqrt(mean((thchims - u)[mask]^2)), "\n")
            }
        }
        else {
            cat(" ", i)
        }
        if (trace) {
            tergs[i, 1, ] <- ni[mask]
            tergs[i, 2, ] <- th[mask]
            tergs[i, 3, ] <- z$sigma[mask]
            tergs[i, 4, ] <- sigma[mask]
        }
    }
    if (!verbose)
        cat("\n")
    if (family == "NCchi") {
        thchi <- fncchir(th/sigma, varstats) * sigma
        thchi[!mask] <- 0
    }
    else {
        thchi <- NULL
    }
    invisible(list(sigma = sigma, sigmal = array(z$sigman, ddim),
        theta = th, thchi = thchi, ni = ni, tergs = tergs, mask = mask))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
