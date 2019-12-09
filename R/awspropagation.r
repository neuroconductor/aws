awstestprop <- function(dy,
                        hmax,
                        theta = 1,
                        family = "Gaussian",
                        lkern = "Triangle",
                        aws = TRUE,
                        memory = FALSE,
                        shape = 2,
                        homogeneous = TRUE,
                        varadapt = FALSE,
                        ladjust = 1,
                        spmin = 0.25,
                        seed = 1,
                        minlevel = 1e-6,
                        maxz = 25,
                        diffz = .5,
                        maxni = FALSE,
                        verbose = FALSE) {
  if (length(dy) > 3) {
    cat(
      "maximum array dimension is 3\n contents of first argument will be interpreted as array of
      parameters\n"
    )
    nnn <- length(dy)
  } else {
    nnn <- prod(dy)
  }
  if (minlevel < 5 / nnn) {
    minlevel <- 5 / nnn
    cat("minlevel reset to",
        minlevel,
        "due to insufficient size of test sample\n")
  }
  set.seed(seed)
  par(
    mfrow = c(1, 1),
    mar = c(3, 3, 3, 3),
    mgp = c(2, 1, 0)
  )
  if (length(dy) <= 3) {
    y <- array(switch(
      family,
      "Gaussian" = rnorm(nnn),
      "Poisson" = rpois(nnn, theta),
      "Exponential" = rexp(nnn, 1),
      "Bernoulli" = rbinom(nnn, 1, theta),
      "Volatility" = rnorm(nnn),
      "Variance" = rchisq(nnn, shape) / shape,
      "NCchi" = sqrt(rchisq(nnn, shape, theta ^ 2))
    ), dy)
  } else {
    ddy <- if (!is.null(dim(dy)))
      dim(dy)
    else
      length(dy)
    y <- array(switch(
      family,
      "Gaussian" = rnorm(nnn, dy - mean(dy)),
      "Poisson" = rpois(nnn, dy - mean(dy) + theta),
      "Exponential" = rexp(nnn, dy - mean(dy) + 1),
      "Bernoulli" = rbinom(nnn, 1, dy - mean(dy) + theta),
      "Volatility" = rnorm(nnn, dy - mean(dy)),
      "Variance" = rchisq(nnn, dy - mean(dy) + shape) /
        (dy - mean(dy) + shape),
      "NCchi" = sqrt(rchisq(nnn, shape, (
        dy - mean(dy) + theta
      ) ^ 2))
    ), ddy)
    dy <- ddy
  }
  cat("minlevel ",
      minlevel,
      "due to insufficient size of test sample\n")
  if (family == "NCchi") {
    varstats <-
      sofmchi(shape / 2) # precompute table of mean, sd and var for
    #
    #   NCchi for noncentral chi with shape=degrees of freedom and theta =NCP
    #
  }
  z <- seq(0, maxz, diffz)
  nz <- length(z)
  elevel <- trunc(log(1e-6, 10))
  levels <- as.vector(outer(c(.5, .2, .1), 10 ^ (-0:elevel)))
  levels <- levels[levels >= minlevel]
  wghts <- switch(length(dy), c(0, 0), c(1, 0), c(1, 1))
  cpar <-
    setawsdefaults(dy,
                   mean(y),
                   family,
                   lkern,
                   "Uniform",
                   aws,
                   memory,
                   ladjust,
                   hmax,
                   shape,
                   wghts)
  lambda <- cpar$lambda
  hmax <- cpar$hmax
  shape <- cpar$shape
  d <- cpar$d
  n <- length(y)
  if (!homogeneous & family == "Gaussian") {
    sigma2 <- array(rchisq(prod(dy), shape) / shape, dy)
  } else
    sigma2 <- 1
  zfamily <- awsfamily(family, y, sigma2, shape, 0, lambda, cpar)
  cpar <- zfamily$cpar
  lambda <- zfamily$lambda
  sigma2 <- zfamily$sigma2
  h0 <- zfamily$h0
  y <- zfamily$y
  lkern <- cpar$lkern
  rm(zfamily)
  n1 <- switch(d, n, dy[1], dy[1])
  n2 <- switch(d, 1, dy[2], dy[2])
  n3 <- switch(d, 1, 1, dy[3])
  maxvol <- cpar$maxvol
  # k <- cpar$k
  k <- 1
  kstar <- cpar$kstar
  h <- numeric(kstar)
  if (k > 1)
    h[1:(k - 1)] <- 1 + (0:(k - 2)) * .001
  exceedence  <-
    exceedencena  <-
    matrix(0, nz, kstar) # this is used to store exceedence probabilities for adaptive and nonadaptive estimates
  zobj <- zobj0 <- list(ai = y, bi = rep(1, n))
  bi <- rep(1, n)
  yhat <- y / shape
  lambda0 <- 1e50
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  #
  #  get initial conditions for a comparison
  #
  if (family == "Bernoulli")
    y0 <- (10 * y + 1) / 12
  if (family == "Poisson")
    y0 <- y + .1
  #
  #  this corresponds to the regularization used to avoid Inf distances
  #
  kldistnorm1 <- function(th1, y, df) {
    L <- df / 2
    m1 <-
      sqrt(pi / 2) * gamma(L + 1 / 2) / gamma(1.5) / gamma(L) * hyperg_1F1(-0.5, L, -th1 ^
                                                                             2 / 2, give = FALSE, strict = TRUE)
    (m1 - y) ^ 2 / 2 / (2 * L + th1 ^ 2 - m1 ^ 2)
  }
  KLdist0 <- switch(
    family,
    "Gaussian" = y ^ 2 / 2,
    "Poisson" = (theta - y0 + y0 * (log(y0) - log(theta))),
    "Exponential" = (log(y) - 1 + 1 / y),
    "Bernoulli" = (y0 * log(y0 / theta) +
                     (1 - y0) * log((1 - y0) / (1 - theta))),
    "Volatility" = (log(y) - 1 + 1 / y) / 2,
    "Variance" = shape / 2 * (log(y) - 1 + 1 / y),
    "NCchi" = kldistnorm1(theta, y, shape)
  )
  exceedence0 <- .Fortran(C_exceed,
    as.double(KLdist0),
    as.integer(length(KLdist0)),
    as.double(z),
    as.integer(nz),
    exprob = double(nz)
  )$exprob
  #
  #  now iterate
  #
  t0 <- Sys.time()
  cat("using lambda=", lambda, "\n")
  while (k <= kstar) {
    t1 <- Sys.time()
    hakt0 <-
      gethani(1, 1.25 * hmax, lkern, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, lkern, 1.25 ^ k, wghts, 1e-4)
    cat("step", k, "hakt", hakt, "\n")
    if (lkern == 5) {
      #  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
      hakt <- hakt * 0.42445 * 4
    }
    h[k] <- hakt
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)[1:d]
    #
    #   get nonadaptive estimate
    #
    if (!homogeneous & family == "Gaussian") {
      zobj0 <- .Fortran(C_chaws,
        as.double(y),
        as.double(sigma2),
        as.integer(1:n),#position
        as.integer(n1),
        as.integer(n2),
        as.integer(n3),
        hakt = as.double(hakt),
        as.double(1e40),
        double(n),
        bi = double(n),
        bi2 = double(n),
        double(n),
        ai = double(n),
        as.integer(cpar$mcode),
        as.integer(lkern),
        as.double(spmin),
        double(prod(dlw)),
        as.double(wghts)
      )[c("bi", "bi2", "ai", "hakt")]
    } else {
      zobj0 <- .Fortran(C_caws,
        as.double(y),
        as.integer(1:n),# position for full mask
        as.integer(n1),
        as.integer(n2),
        as.integer(n3),
        hakt = as.double(hakt),
        as.double(1e40),
        double(n),
        bi = double(n),
        bi2 = double(n),
        double(n),
        ai = double(n),
        as.integer(cpar$mcode),
        as.integer(lkern),
        as.double(spmin),
        double(prod(dlw)),
        as.double(wghts)
      )[c("bi", "bi2", "ai", "hakt")]
    }
    if (family %in% c("Bernoulli", "Poisson"))
      zobj0 <- regularize(zobj0, family)
    yhat0 <- zobj0$ai / zobj0$bi
    dim(yhat0) <- dy
    bi <- zobj0$bi
    ih <- as.integer(hakt)
    ind1 <- (ih + 1):(dy[1] - ih)
    if (length(dy) > 1)
      ind2 <- (ih + 1):(dy[2] - ih)
    if (length(dy) > 2)
      ind3 <- (ih + 1):(dy[3] - ih)
    yhat0 <-
      switch(length(dy), yhat0[ind1], yhat0[ind1, ind2], yhat0[ind1, ind2, ind3])
    ni <- max(zobj0$bi)
    KLdist0 <- switch(
      family,
      "Gaussian" = yhat0 ^ 2 / 2,
      "Poisson" = (theta - yhat0 + yhat0 * (log(yhat0) -
                                              log(theta))),
      "Exponential" = (log(yhat0) - 1 + 1 / yhat0),
      "Bernoulli" = (yhat0 * log(yhat0 / theta) +
                       (1 - yhat0) * log((1 - yhat0) / (1 -
                                                          theta))),
      "Volatility" = (log(yhat0) - 1 + 1 / yhat0) / 2,
      "Variance" = shape / 2 * (log(yhat0) - 1 + 1 /
                                  yhat0),
      "NCchi" = kldistnorm1(theta, yhat0, shape)
    )
    exceedencena[, k] <- .Fortran(C_exceed,
      as.double(KLdist0),
      as.integer(length(KLdist0)),
      as.double(z / ni),
      as.integer(nz),
      exprob = double(nz)
    )$exprob
    #
    #   get adaptive estimate
    #
    if (!homogeneous & family == "Gaussian") {
      zobj <- .Fortran(C_chaws,
        as.double(y),
        as.double(sigma2),
        as.integer(1:n),# position for full mask
        as.integer(n1),
        as.integer(n2),
        as.integer(n3),
        hakt = as.double(hakt),
        as.double(lambda0),
        as.double(yhat),
        bi = as.double(bi),
        bi2 = double(n),
        double(n),
        ai = double(n),
        as.integer(cpar$mcode),
        as.integer(lkern),
        as.double(spmin),
        double(prod(dlw)),
        as.double(wghts)
      )[c("bi", "bi2", "ai", "hakt")]
    } else {
      if (cpar$mcode != 6) {
        zobj <- .Fortran(C_caws,
          as.double(y),
          as.integer(1:n),# position for full mask
          as.integer(n1),
          as.integer(n2),
          as.integer(n3),
          hakt = as.double(hakt),
          as.double(lambda0),
          as.double(yhat),
          bi = as.double(bi),
          bi2 = double(n),
          double(n),
          ai = double(n),
          as.integer(cpar$mcode),
          as.integer(lkern),
          as.double(spmin),
          double(prod(dlw)),
          as.double(wghts)
        )[c("bi", "bi2", "ai", "hakt")]
      } else {
        zobj <- .Fortran(C_caws6,
          as.double(y),
          as.integer(1:n),# position for full mask
          as.integer(n1),
          as.integer(n2),
          as.integer(n3),
          hakt = as.double(hakt),
          as.double(lambda0),
          as.double(yhat),
          as.double(fncchiv(yhat, varstats) / 2),
          bi = as.double(bi),
          bi2 = double(n),
          double(n),
          ai = double(n),
          as.integer(lkern),
          as.double(spmin),
          double(prod(dlw)),
          as.double(wghts)
        )[c("bi", "bi2", "ai", "hakt")]
      }
    }
    if (family %in% c("Bernoulli", "Poisson"))
      zobj <- regularize(zobj, family)
    dim(zobj$ai) <- dy
    yhat <- zobj$ai / zobj$bi
    dim(yhat) <- dy
    if (varadapt)
      bi <- bi ^ 2 / zobj$bi2
    if (maxni)
      bi <- pmax(bi, zobj$bi)
    else
      bi <- zobj$bi
    lambda0 <- lambda
    yhat0 <-
      switch(length(dy), yhat[ind1], yhat[ind1, ind2], yhat[ind1, ind2, ind3])
    KLdist1 <- switch(
      family,
      "Gaussian" = yhat0 ^ 2 / 2,
      "Poisson" = (theta - yhat0 + yhat0 * (log(yhat0) -
                                              log(theta))),
      "Exponential" = (log(yhat0) - 1 + 1 / yhat0),
      "Bernoulli" = (yhat0 * log(yhat0 / theta) +
                       (1 - yhat0) * log((1 - yhat0) / (1 -
                                                          theta))),
      "Volatility" = (log(yhat0) - 1 + 1 / yhat0) / 2,
      "Variance" = shape / 2 * (log(yhat0) - 1 + 1 /
                                  yhat0),
      "NCchi" = kldistnorm1(theta, yhat0, shape)
    )
    exceedence[, k] <- .Fortran(C_exceed,
      as.double(KLdist1),
      as.integer(length(KLdist1)),
      as.double(z / ni),
      as.integer(nz),
      exprob = double(nz)
    )$exprob

    contour(
      z,
      0:k,
      cbind(exceedence0, exceedence[, 1:k]),
      levels = levels,
      ylab = "step",
      xlab = "z",
      main = paste(family, length(dy), "-dim. ladj=", ladjust, " Exceed. Prob.")
    )
    contour(
      z,
      0:k,
      cbind(exceedence0, exceedencena[, 1:k]),
      levels = levels,
      ylab = "step",
      xlab = "z",
      add = TRUE,
      col = 2,
      lty = 3
    )
    yaxp <- par("yaxp")
    at <-
      unique(as.integer(seq(yaxp[1], yaxp[2], length = yaxp[3] + 1)))
    at <- at[at > 0 & at <= k]
    axis(4, at = at, labels = as.character(signif(h[at], 3)))
    mtext("bandwidth", 4, 1.8)
    if (max(total) > 0) {
      cat("Progress:", signif(total[k], 2) * 100, "% . ", sep = "")
    }
    t2 <- Sys.time()
    tpar <-
      if (family %in% c("Bernoulli", "Poisson"))
        paste("theta=", theta)
    else
      ""
    cat(
      family,
      "(dim:",
      length(dy),
      tpar,
      ") ni=",
      ni,
      " Time: Step",
      format(signif(difftime(t2, t1), 3)),
      "Total",
      format(signif(difftime(t2, t0), 3)),
      "\n"
    )
    k <- k + 1
    gc()
  }
  if (family %in% c("Bernoulli", "Poisson"))
    y <- y0

  list(
    h = h,
    z = z,
    prob = exceedence,
    probna = exceedencena,
    y = if (verbose)
      y
    else
      NULL ,
    theta = if (verbose)
      yhat
    else
      NULL,
    levels = levels,
    family = family
  )
}

pawstestprop <- function(dy,
                        hmax,
                        theta = 1,
                        family = "Gaussian",
                        lkern = "Triangle",
                        aws = TRUE,
                        patchsize=1,
                        shape = 2,
                        ladjust = 1,
                        spmin = 0.25,
                        seed = 1,
                        minlevel = 1e-6,
                        maxz = 25,
                        diffz = .5,
                        maxni = FALSE,
                        verbose = FALSE) {
  varadapt <- FALSE
  patchadapt <- "max"
  if (length(dy) > 3) {
    cat(
      "maximum array dimension is 3\n contents of first argument will be interpreted as array of
      parameters\n"
    )
    nnn <- length(dy)
  } else {
    nnn <- prod(dy)
  }
  if (minlevel < 5 / nnn) {
    minlevel <- 5 / nnn
    cat("minlevel reset to",
        minlevel,
        "due to insufficient size of test sample\n")
  }
  set.seed(seed)
  par(
    mfrow = c(1, 1),
    mar = c(3, 3, 3, 3),
    mgp = c(2, 1, 0)
  )
  if (length(dy) <= 3) {
    y <- array(switch(
      family,
      "Gaussian" = rnorm(nnn),
      "Poisson" = rpois(nnn, theta),
      "Exponential" = rexp(nnn, 1),
      "Bernoulli" = rbinom(nnn, 1, theta),
      "Volatility" = rnorm(nnn),
      "Variance" = rchisq(nnn, shape) / shape,
      "NCchi" = sqrt(rchisq(nnn, shape, theta ^ 2))
    ), dy)
  } else {
    ddy <- if (!is.null(dim(dy)))
      dim(dy)
    else
      length(dy)
    y <- array(switch(
      family,
      "Gaussian" = rnorm(nnn, dy - mean(dy)),
      "Poisson" = rpois(nnn, dy - mean(dy) + theta),
      "Exponential" = rexp(nnn, dy - mean(dy) + 1),
      "Bernoulli" = rbinom(nnn, 1, dy - mean(dy) + theta),
      "Volatility" = rnorm(nnn, dy - mean(dy)),
      "Variance" = rchisq(nnn, dy - mean(dy) + shape) /
        (dy - mean(dy) + shape),
      "NCchi" = sqrt(rchisq(nnn, shape, (
        dy - mean(dy) + theta
      ) ^ 2))
    ), ddy)
    dy <- ddy
  }
  cat("minlevel ",
      minlevel,
      "due to insufficient size of test sample\n")
  if (family == "NCchi") {
    varstats <-
      sofmchi(shape / 2) # precompute table of mean, sd and var for
    #
    #   NCchi for noncentral chi with shape=degrees of freedom and theta =NCP
    #
  }
  z <- seq(0, maxz, diffz)
  nz <- length(z)
  elevel <- trunc(log(1e-6, 10))
  levels <- as.vector(outer(c(.5, .2, .1), 10 ^ (-0:elevel)))
  levels <- levels[levels >= minlevel]
  wghts <- switch(length(dy), c(0, 0), c(1, 0), c(1, 1))
  cpar <-
    setawsdefaults(dy,
                   mean(y),
                   family,
                   lkern,
                   "Uniform",
                   aws,
                   FALSE,
                   ladjust,
                   hmax,
                   shape,
                   wghts)
  lambda <- cpar$lambda
  # additional adjustments for taking the maximum of s_{ij} over patches
  # adjusted using simulations such that for homogeneous structures
  # loss by adaptation < 1% in MAE and <.1 in PSNR
  if(length(dy)==1){
     ladjust <- switch(patchsize,.55,.55,.5,.5)*ladjust
  } else if(length(dy)==2){
     ladjust <- switch(patchsize,1.1,1.1,1.2)*ladjust
} else {
     ladjust <- switch(patchsize,1.44,1.8)*ladjust
}
  hmax <- cpar$hmax
  shape <- cpar$shape
  d <- cpar$d
  n <- length(y)
  sigma2 <- 1
  zfamily <- awsfamily(family, y, sigma2, shape, 0, lambda, cpar)
  cpar <- zfamily$cpar
  lambda <- zfamily$lambda
  sigma2 <- zfamily$sigma2
  h0 <- zfamily$h0
  y <- zfamily$y
  lkern <- cpar$lkern
  rm(zfamily)
  n1 <- switch(d, n, dy[1], dy[1])
  n2 <- switch(d, 1, dy[2], dy[2])
  n3 <- switch(d, 1, 1, dy[3])
  maxvol <- cpar$maxvol
  # k <- cpar$k
  k <- 1
  kstar <- cpar$kstar
  h <- numeric(kstar)
  if (k > 1)
    h[1:(k - 1)] <- 1 + (0:(k - 2)) * .001
  exceedence  <-
    exceedencena  <-
    matrix(0, nz, kstar) # this is used to store exceedence probabilities for adaptive and nonadaptive estimates
  zobj <- zobj0 <- list(ai = y, bi = rep(1, n))
  bi <- rep(1, n)
  yhat <- y / shape
  lambda0 <- 1e50
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  #
  #  get initial conditions for a comparison
  #
  if (family == "Bernoulli")
    y0 <- (10 * y + 1) / 12
  if (family == "Poisson")
    y0 <- y + .1
  #
  #  this corresponds to the regularization used to avoid Inf distances
  #
  kldistnorm1 <- function(th1, y, df) {
    L <- df / 2
    m1 <-
      sqrt(pi / 2) * gamma(L + 1 / 2) / gamma(1.5) / gamma(L) * hyperg_1F1(-0.5, L, -th1 ^
                                                                             2 / 2, give = FALSE, strict = TRUE)
    (m1 - y) ^ 2 / 2 / (2 * L + th1 ^ 2 - m1 ^ 2)
  }
  KLdist0 <- switch(
    family,
    "Gaussian" = y ^ 2 / 2,
    "Poisson" = (theta - y0 + y0 * (log(y0) - log(theta))),
    "Exponential" = (log(y) - 1 + 1 / y),
    "Bernoulli" = (y0 * log(y0 / theta) +
                     (1 - y0) * log((1 - y0) / (1 - theta))),
    "Volatility" = (log(y) - 1 + 1 / y) / 2,
    "Variance" = shape / 2 * (log(y) - 1 + 1 / y),
    "NCchi" = kldistnorm1(theta, y, shape)
  )
  exceedence0 <- .Fortran(C_exceed,
    as.double(KLdist0),
    as.integer(length(KLdist0)),
    as.double(z),
    as.integer(nz),
    exprob = double(nz)
  )$exprob
  #
  #  define patch area
  #
  patchsize <- min(patchsize,5-length(dy))
  np1 <- if (patchsize > 0)
    2 * patchsize + 1
  else
    7
  ## patchsize == 0 includes immediate neighbors only
  np2 <- if (n2 > 1)
    2 * patchsize + 1
  else
    1
  np3 <- if (n3 > 1)
    2 * patchsize + 1
  else
    1
  t0 <- Sys.time()
  cat("using lambda=", lambda, "\n")
  mask <- array(TRUE,dy)
  #
  #  now iterate
  #
  while (k <= kstar) {
    t1 <- Sys.time()
    hakt0 <-
      gethani(1, 1.25 * hmax, lkern, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, lkern, 1.25 ^ k, wghts, 1e-4)
    cat("step", k, "hakt", hakt, "\n")
    if (lkern == 5) {
      #  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
      hakt <- hakt * 0.42445 * 4
    }
    h[k] <- hakt
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)[1:d]
    #
    #  define interior mask
    #
    ih <- trunc(h)+patchsize
    if(nnn==1) mask[c(1:ih,n1+1-1:ih)] <- FALSE
    if(nnn==2) mask[c(1:ih,n1+1-1:ih),
                    c(1:ih,n2+1-1:ih)] <- FALSE
    if(nnn==3) mask[c(1:ih,n1+1-1:ih),
                    c(1:ih,n2+1-1:ih),
                    c(1:ih,n3+1-1:ih)] <- FALSE
    #
    #   get nonadaptive estimate
    #
    zobj0 <- .Fortran(C_caws,
      as.double(y),
      as.integer(1:n),# position for full mask
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(1e40),
      double(n),
      bi = double(n),
      bi2 = double(n),
      double(n),
      ai = double(n),
      as.integer(cpar$mcode),
      as.integer(lkern),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts)
    )[c("bi", "bi2", "ai", "hakt")]
    if (family %in% c("Bernoulli", "Poisson"))
      zobj0 <- regularize(zobj0, family)
    yhat0 <- zobj0$ai / zobj0$bi
    dim(yhat0) <- dy
    bi <- zobj0$bi
    ih <- as.integer(hakt)
    ind1 <- (ih + 1):(dy[1] - ih)
    if (length(dy) > 1)
      ind2 <- (ih + 1):(dy[2] - ih)
    if (length(dy) > 2)
      ind3 <- (ih + 1):(dy[3] - ih)
    yhat0 <-
      switch(length(dy), yhat0[ind1], yhat0[ind1, ind2], yhat0[ind1, ind2, ind3])
    ni <- max(zobj0$bi)
    KLdist0 <- switch(
      family,
      "Gaussian" = yhat0 ^ 2 / 2,
      "Poisson" = (theta - yhat0 + yhat0 * (log(yhat0) -
                                              log(theta))),
      "Exponential" = (log(yhat0) - 1 + 1 / yhat0),
      "Bernoulli" = (yhat0 * log(yhat0 / theta) +
                       (1 - yhat0) * log((1 - yhat0) / (1 -
                                                          theta))),
      "Volatility" = (log(yhat0) - 1 + 1 / yhat0) / 2,
      "Variance" = shape / 2 * (log(yhat0) - 1 + 1 /
                                  yhat0),
      "NCchi" = kldistnorm1(theta, yhat0, shape)
    )
    exceedencena[, k] <- .Fortran(C_exceedm,
      as.double(KLdist0),
      as.integer(length(KLdist0)),
      as.double(z / ni),
      as.integer(nz),
      exprob = double(nz),
      as.integer(mask)
    )$exprob
    #
    #   get adaptive estimate
    #

        zobj <- .Fortran(C_pcaws,
          as.double(y),
          as.double(1:n),# full mask
          as.integer(n1),
          as.integer(n2),
          as.integer(n3),
          hakt = as.double(hakt),
          as.double(lambda0),
          as.double(yhat),
          as.double(zobj$bi),
          bi2 = double(n),
          bi = double(n), #biout
          theta = double(n),
          as.integer(cpar$mcode),
          as.integer(lkern),
          as.double(spmin),
          double(prod(dlw)),
          as.double(wghts),
          as.integer(patchsize))[c("bi", "bi2", "theta", "hakt")]
    if (family %in% c("Bernoulli", "Poisson"))
      zobj <- regularize(zobj, family)
    yhat <- zobj$theta
    dim(yhat) <- dy
    if (varadapt)
      bi <- bi ^ 2 / zobj$bi2
    if (maxni)
      bi <- pmax(bi, zobj$bi)
    else
      bi <- zobj$bi
    lambda0 <- lambda
    yhat0 <-
      switch(length(dy), yhat[ind1], yhat[ind1, ind2], yhat[ind1, ind2, ind3])
    KLdist1 <- switch(
      family,
      "Gaussian" = yhat0 ^ 2 / 2,
      "Poisson" = (theta - yhat0 + yhat0 * (log(yhat0) -
                                              log(theta))),
      "Exponential" = (log(yhat0) - 1 + 1 / yhat0),
      "Bernoulli" = (yhat0 * log(yhat0 / theta) +
                       (1 - yhat0) * log((1 - yhat0) / (1 -
                                                          theta))),
      "Volatility" = (log(yhat0) - 1 + 1 / yhat0) / 2,
      "Variance" = shape / 2 * (log(yhat0) - 1 + 1 /
                                  yhat0),
      "NCchi" = kldistnorm1(theta, yhat0, shape)
    )
    exceedence[, k] <- .Fortran(C_exceedm,
      as.double(KLdist1),
      as.integer(length(KLdist1)),
      as.double(z / ni),
      as.integer(nz),
      exprob = double(nz),
      as.integer(mask)
    )$exprob

    contour(
      z,
      0:k,
      cbind(exceedence0, exceedence[, 1:k]),
      levels = levels,
      ylab = "step",
      xlab = "z",
      main = paste(family, length(dy), "-dim. ladj=", ladjust, " Exceed. Prob.")
    )
    contour(
      z,
      0:k,
      cbind(exceedence0, exceedencena[, 1:k]),
      levels = levels,
      ylab = "step",
      xlab = "z",
      add = TRUE,
      col = 2,
      lty = 3
    )
    yaxp <- par("yaxp")
    at <-
      unique(as.integer(seq(yaxp[1], yaxp[2], length = yaxp[3] + 1)))
    at <- at[at > 0 & at <= k]
    axis(4, at = at, labels = as.character(signif(h[at], 3)))
    mtext("bandwidth", 4, 1.8)
    if (max(total) > 0) {
      cat("Progress:", signif(total[k], 2) * 100, "% . ", sep = "")
    }
    t2 <- Sys.time()
    tpar <-
      if (family %in% c("Bernoulli", "Poisson"))
        paste("theta=", theta)
    else
      ""
    cat(
      family,
      "(dim:",
      length(dy),
      tpar,
      ") ni=",
      ni,
      " mean(bi/ni)=",
      mean(bi[mask])/ni,
      " Time: Step",
      format(signif(difftime(t2, t1), 3)),
      "Total",
      format(signif(difftime(t2, t0), 3)),
      "\n"
    )
    k <- k + 1
    gc()
  }
  if (family %in% c("Bernoulli", "Poisson"))
    y <- y0

  list(
    h = h,
    z = z,
    prob = exceedence,
    probna = exceedencena,
    y = if (verbose)
      y
    else
      NULL ,
    theta = if (verbose)
      yhat
    else
      NULL,
    levels = levels,
    family = family
  )
}
