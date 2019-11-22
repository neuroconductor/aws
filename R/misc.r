#########################################################################
#
#   functions to handle the  noncentral chi case (mcode=6)
#
#########################################################################
sofmchi <- function(L) {
  minlev <- sqrt(2) * gamma(L + .5) / gamma(L)
  x <- seq(0, 50, .01)
  mu <-
    sqrt(pi / 2) * gamma(L + 1 / 2) / gamma(1.5) / gamma(L) * hyperg_1F1(-0.5, L,-x ^
                                                                           2 / 2, give = FALSE, strict = TRUE)
  s2 <- 2 * L + x ^ 2 - mu ^ 2
  s <- sqrt(s2)
  ## return list containing values of noncentrality parameter (ncp),
  ## mean (mu), standard deviation (sd) and variance (s2) to be used
  ## in variance modeling
  list(
    ncp = x,
    mu = mu,
    s = s,
    s2 = s2,
    minlev = minlev,
    L = L
  )
}

fncchiv <- function(mu, varstats) {
  mu <- pmax(varstats$minlev, mu)
  ind <-
    findInterval(mu,
                 varstats$mu,
                 rightmost.closed = FALSE,
                 all.inside = FALSE)
  varstats$s2[ind]
}
#########################################################################
#
#   binning in 1D -- 3D (adapted from binning function in package sm
#
#########################################################################
binning <- function (x, y, nbins, xrange = NULL) {
  dx <- dim(x)
  if (is.null(dx))
    d <- 1
  else
    d <- dx[2]
  if (d > 3) {
    warning("Binning only implemented in 1D, 2D and 3D")
    return(NULL)
  }
  if (length(nbins) < d || any(nbins < 2)) {
    warning("Invalid values for nbins")
    return(NULL)
  }
  if (!is.null(y) && length(y) * d != length(x)) {
    warning("Dimensions of design matrix incompatible with length of response vector")
    return(NULL)
  }
  if (is.null(xrange)) {
    xrange <- if (d == 1)
      range(x)
    else
      apply(x, 2, range)
  } else {
    if ((d == 1 &&
         length(xrange) != 2) || (d > 1 && any(dim(xrange) != c(2, d)))) {
      warning("Dimensions of xrange incorrect ")
      return(NULL)
    }
    xrange <-
      if (d == 1)
        range(x, xrange)
    else
      apply(rbind(x, xrange), 2, range)
  }
  xnames <- if (d > 1)
    dimnames(x)[[2]]
  else
    names(x)
  breaks.x1 <- seq(xrange[1], xrange[2], length = nbins[1] + 1)
  if (d > 1)
    breaks.x2 <- seq(xrange[1, 2], xrange[2, 2], length = nbins[2] + 1)
  if (d > 2)
    breaks.x3 <- seq(xrange[1, 3], xrange[2, 3], length = nbins[3] + 1)
  f1 <- cut(if (d == 1)
    x
    else
      x[, 1], breaks = breaks.x1)
  if (d > 1)
    f2 <- cut(x[, 2], breaks = breaks.x2)
  if (d > 2)
    f3 <- cut(x[, 3], breaks = breaks.x3)
  freq <- switch(d, table(f1), table(f1, f2), table(f1, f2, f3))
  dimnames(freq) <- NULL
  midpoints.x1 <-
    (breaks.x1[-1] + breaks.x1[-(nbins[1] + 1)]) / 2
  if (d > 1)
    midpoints.x2 <- (breaks.x2[-1] + breaks.x2[-(nbins[2] + 1)]) / 2
  if (d > 2)
    midpoints.x3 <- (breaks.x3[-1] + breaks.x3[-(nbins[3] + 1)]) / 2
  z1 <- midpoints.x1
  if (d > 1)
    z2 <- midpoints.x2
  if (d > 2)
    z3 <- midpoints.x3
  X <- switch(d, z1,
              cbind(rep(z1, length(z2)),
                    rep(z2, rep(
                      length(z1), length(z2)
                    ))),
              cbind(rep(z1, length(z2) * length(z3)),
                    rep(z2, rep(
                      length(z1) * length(z3), length(z2)
                    )),
                    rep(z3, rep(
                      length(z1) * length(z2), length(z3)
                    ))))
  X.f <- as.vector(freq)
  id <- (X.f > 0)
  if (d > 1)
    X <- X[id,]
  else
    X <- X[id]
  if (d > 1)
    dimnames(X) <- list(NULL, xnames)
  else
    names(X) <- xnames
  X.f <- X.f[id]
  result <- list(
    x = X,
    x.freq = X.f,
    midpoints.x1 = midpoints.x1,
    midpoints.x2 = if (d > 1)
      midpoints.x2
    else
      NULL,
    midpoints.x3 = if (d > 2)
      midpoints.x3
    else
      NULL,
    breaks.x1 = breaks.x1,
    breaks.x2 = if (d > 1)
      breaks.x2
    else
      NULL,
    breaks.x3 = if (d > 2)
      breaks.x3
    else
      NULL,
    table.freq = freq
  )
  if (!is.null(y) && !all(is.na(y))) {
    result$means <- as.numeric(tapply(y, switch(
      d, list(f1),
      list(f1, f2), list(f1, f2, f3)
    ),
    mean))[id]
    result$devs <- as.numeric(tapply(y, switch(
      d, list(f1),
      list(f1, f2), list(f1, f2, f3)
    ),
    function(x)
      sum((x - mean(
        x
      )) ^ 2)))[id]
  }
  result
}
Varcor.gauss <- function(h) {
  #
  #   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
  #
  #   in case of colored noise that was produced by smoothing with lkern and bandwidth h
  #
  h <- pmax(h / 2.3548, 1e-5)
  ih <- trunc(4 * h) + 1
  dx <- 2 * ih + 1
  d <- length(h)
  penl <- dnorm(((-ih[1]):ih[1]) / h[1])
  if (d == 2)
    penl <- outer(penl, dnorm(((-ih[2]):ih[2]) / h[2]), "*")
  if (d == 3)
    penl <-
    outer(penl, outer(dnorm(((
      -ih[2]
    ):ih[2]) / h[2]), dnorm(((
      -ih[3]
    ):ih[3]) / h[3]), "*"), "*")
  2 * sum(penl) ^ 2 / sum(diff(penl) ^ 2)
}


Spatialvar.gauss <- function(h, h0, d) {
  #
  #   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in
  #
  #   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
  #
  #   Spatialvariance(lkern,h,h0,d)/Spatialvariance(lkern,h,1e-5,d) gives the
  #   a factor for lambda to be used with bandwidth h
  #
  h0 <- max(1e-5, h0)
  h <- h / 2.3548
  if (length(h) == 1)
    h <- rep(h, d)
  ih <- trunc(4 * h)
  ih <- pmax(1, ih)
  dx <- 2 * ih + 1
  penl <- dnorm(((-ih[1]):ih[1]) / h[1])
  if (d == 2)
    penl <-
    outer(dnorm(((-ih[1]):ih[1]) / h[1]), dnorm(((-ih[2]):ih[2]) / h[2]), "*")
  if (d == 3)
    penl <-
    outer(dnorm(((-ih[1]):ih[1]) / h[1]), outer(dnorm(((
      -ih[2]
    ):ih[2]) / h[2]), dnorm(((
      -ih[3]
    ):ih[3]) / h[3]), "*"), "*")
  dim(penl) <- dx
  h0 <- h0 / 2.3548
  if (length(h0) == 1)
    h0 <- rep(h0, d)
  ih <- trunc(4 * h0)
  ih <- pmax(1, ih)
  dx0 <- 2 * ih + 1
  x <- ((-ih[1]):ih[1]) / h0[1]
  penl0 <- dnorm(((-ih[1]):ih[1]) / h0[1])
  if (d == 2)
    penl0 <-
    outer(dnorm(((-ih[1]):ih[1]) / h0[1]), dnorm(((-ih[2]):ih[2]) / h0[2]), "*")
  if (d == 3)
    penl0 <-
    outer(dnorm(((-ih[1]):ih[1]) / h0[1]), outer(dnorm(((
      -ih[2]
    ):ih[2]) / h0[2]), dnorm(((
      -ih[3]
    ):ih[3]) / h0[3]), "*"), "*")
  dim(penl0) <- dx0
  penl0 <- penl0 / sum(penl0)
  dz <- dx + dx0 - 1
  z <- array(0, dz)
  if (d == 1) {
    for (i1 in 1:dx0) {
      ind1 <- c(0:(i1 - 1), (dz - dx0 + i1):dz + 1)
      ind1 <- ind1[ind1 <= dz][-1]
      z[-ind1] <- z[-ind1] + penl * penl0[i1]
    }
  } else if (d == 2) {
    for (i1 in 1:dx0[1])
      for (i2 in 1:dx0[2]) {
        ind1 <- c(0:(i1 - 1), (dz[1] - dx0[1] + i1):dz[1] + 1)
        ind1 <- ind1[ind1 <= dz[1]][-1]
        ind2 <- c(0:(i2 - 1), (dz[2] - dx0[2] + i2):dz[2] + 1)
        ind2 <- ind2[ind2 <= dz[2]][-1]
        z[-ind1, -ind2] <- z[-ind1, -ind2] + penl * penl0[i1, i2]
      }
  } else if (d == 3) {
    for (i1 in 1:dx0[1])
      for (i2 in 1:dx0[2])
        for (i3 in 1:dx0[3]) {
          ind1 <- c(0:(i1 - 1), (dz[1] - dx0[1] + i1):dz[1] + 1)
          ind1 <- ind1[ind1 <= dz[1]][-1]
          ind2 <- c(0:(i2 - 1), (dz[2] - dx0[2] + i2):dz[2] + 1)
          ind2 <- ind2[ind2 <= dz[2]][-1]
          ind3 <- c(0:(i3 - 1), (dz[3] - dx0[3] + i3):dz[3] + 1)
          ind3 <- ind3[ind3 <= dz[3]][-1]
          z[-ind1, -ind2, -ind3] <- z[-ind1, -ind2, -ind3] + penl * penl0[i1, i2, i3]
        }
  }
  sum(z ^ 2) / sum(z) ^ 2
}

geth.gauss <- function(corr, step = 1.01) {
  #   get the   bandwidth for lkern corresponding to a given correlation
  #
  #  keep it simple result does not depend on d
  #
  if (corr < 0.1) {
    h <- 1e-5
  } else {
    h <- .8
    z <- 0
    while (z < corr) {
      h <- h * step
      z <- get.corr.gauss(h, interv = 2)
    }
    h <- h / step
  }
  h
}


get.corr.gauss <- function(h, interv = 1) {
  #
  #   Calculates the correlation of
  #   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
  #   Result does not depend on d for "Gaussian" kernel !!
  h <- h / 2.3548 * interv
  ih <- trunc(4 * h + 2 * interv - 1)
  dx <- 2 * ih + 1
  penl <- dnorm(((-ih):ih) / h)
  sum(penl[-(1:interv)] * penl[-((dx - interv + 1):dx)]) / sum(penl ^
                                                                 2)
}

residualVariance <- function(residuals, mask, resscale=1, compact=FALSE){
   nt <- dim(residuals)[1]
   nvoxel <- sum(mask)
   if(!compact){
      ddim <- dim(mask)
      dim(residuals) <- c(nt,prod(ddim))
      residuals <- residuals[,mask]
   }
   z <- .Fortran(C_ivar,as.double(residuals),
                        as.double(resscale),
                        as.integer(nvoxel),
                        as.integer(nt),
                        var = double(nvoxel))$var
   if(compact){
      resvar <- z
   } else {
      resvar <- array(0,ddim)
      resvar[mask] <- z
   }
   resvar
}

sweepMean <- function()

residualSpatialCorr <- function(residuals, mask, lags=c(5,5,3), compact=FALSE){
   nt <- dim(residuals)[1]
   ddim <- dim(mask)
   if(compact){
#  for the current code we need to expand
      res <- array(0,c(nt,prod(ddim)))
      res[,mask] <- residuals
      residuals <- res
      rm(res)
   }
   corr <- .Fortran(C_imcorr, as.double(residuals), as.integer(mask),
                    as.integer(ddim[1]), as.integer(ddim[2]), as.integer(ddim[3]),
                    as.integer(nt), scorr = double(prod(lags)), as.integer(lags[1]),
                    as.integer(lags[2]), as.integer(lags[3]))$scorr
   dim(corr) <- lags
   corr  
}
