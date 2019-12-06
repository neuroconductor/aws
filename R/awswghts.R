#
#   these routines are not public (in namespace)
#   they have been used to compute weighting schemes for illustrations
#   in the JSS 2019 article
#
awsweights <- function(awsobj, spmin = 0.25, inx = NULL) {
  if (awsobj@degree != 0 || awsobj@varmodel != "Constant" ||
      any(awsobj@scorr != 0))
    stop("Adjustment for correlations is not implemented")
  ##  if is.null(inx)  the complete weight configuration will be
  ##  computed. This may be huge and exceed memory
  ##
  ##
  dy <- awsobj@dy
  n1 <- dy[1]
  ldy <- length(dy)
  if (is.null(ldy))
    ldy <- 1
  if (ldy > 1)
    n2 <- dy[2]
  else
    n2 <- 1
  if (ldy == 3)
    n3 <- dy[3]
  else
    n3 <- 1
  n <- n1 * n2 * n3
  hakt <- awsobj@hmax * 1.25 ^ (1 / ldy)
  ## bandwidth for an additional step of aws
  lambda0 <- awsobj@lambda
  ## this reflects lambda*sigma^2 from aws
  yhat <- awsobj@theta
  bi <- awsobj@ni
  mcode <- switch(
    awsobj@family,
    Gaussian = 1,
    Bernoulli = 2,
    Poisson = 3,
    Exponential = 4,
    Volatility = 4,
    Variance = 5
  )
  lkern <- awsobj@lkern
  hakt <- rep(hakt, ldy)
  dlw <- (2 * trunc(hakt) + 1)
  if (!is.null(inx)) {
    dinx <- dim(inx)
    if (is.null(dinx) & n2 == 1) {
      ## distinguish between univariate problems and multiple points of interest
      ## and 2D/3D problems with single point of interest
      dinx <- c(1, length(inx))
      dim(inx) <- dinx
    }
    if (is.null(dinx)) {
      anzx <- 1
      linx <- length(inx)
      ix <- inx[1]
      iy <- if (linx > 1)
        inx[2]
      else
        1
      iz <- if (linx > 2)
        inx[3]
      else
        1
    } else {
      linx <- dinx[1]
      anzx <- dinx[2]
      ix <- inx[1,]
      iy <- if (linx > 1)
        inx[2,]
      else
        rep(1, anzx)
      iz <- if (linx > 2)
        inx[3,]
      else
        rep(1, anzx)
    }
    zobj <- .Fortran(C_cawsw1,
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      as.integer(ix),
      as.integer(iy),
      as.integer(iz),
      as.integer(anzx),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(yhat),
      bi = as.double(bi),
      as.integer(mcode),
      as.integer(lkern),
      as.double(spmin),
      double(prod(dlw)),
      wghts = double(n * anzx)
    )$wghts
    dim(zobj) <- if (anzx == 1)
      dy
    else
      c(dy, anzx)
  } else{
    if (n > 128 ^ 2)
      stop("Weight scheme would use more than 2 GB memory, please specify locations in inx ")
    zobj <- .Fortran(C_cawsw,
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(yhat),
      bi = as.double(bi),
      as.integer(mcode),
      as.integer(lkern),
      as.double(spmin),
      double(prod(dlw)),
      wghts = double(n * n)
    )$wghts
    dim(zobj) <- c(dy, dy)
  }
  zobj
}

pawswghts <- function(awsobj,
           patchsize,
           position)
  {
    #
    #   patch based version (patches of patchsize neighbors in each direction)
    #   compute weighting scheme in position corresponding to ths situation
    #   after the last iteration step of patch based aws
    #   Gaussian model only
    #
    hmax <- awsobj@hmax
    theta <- awsobj@theta
    bi <- awsobj@ni
    sigma2 <- awsobj@sigma2
    lambda <- awsobj@lambda
    lkern <- awsobj@lkern
    spmin <- .25
    wghts <- awsobj@wghts
    mcode <- 1

    n <- length(theta)
    dy <- dim(theta)
    if (is.null(dy))
    dy <- length(theta)
    d <- length(dy)
    if (d > 3)
      stop("AWS for more than 3 dimensional grids is not implemented")
    #
    #   family dependent transformations that depend on the value of family
    #
    n1 <- switch(d, n, dy[1], dy[1])
    n2 <- switch(d, 1, dy[2], dy[2])
    n3 <- switch(d, 1, 1, dy[3])
    i1 <- position[1]
    i2 <- if(n2>1) position[2] else 1
    i3 <- if(n2>1) position[3] else 1
    #
    #    Initialize  for the iteration
    #
    hakt <- 1.25^(1/d)*hmax
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)[1:d]
    np1 <- if (patchsize > 0) 2 * patchsize + 1 else 1
      ## patchsize == 0 includes immediate neighbors only
    np2 <- if (n2 > 1) 2 * patchsize + 1 else 1
    np3 <- if (n3 > 1) 2 * patchsize + 1 else 1
      # all other cases
    wi <- .Fortran(C_pawswght,
          as.integer(n1),
          as.integer(n2),
          as.integer(n3),
          as.integer(i1),
          as.integer(i2),
          as.integer(i3),
          hakt = as.double(hakt),
          as.double(lambda),
          as.double(theta),
          as.double(bi),
          as.integer(mcode),
          as.integer(lkern),
          as.double(spmin),
          double(prod(dlw)),
          as.double(wghts),
          as.integer(patchsize),
          w=double(n1*n2*n3))[["w"]]
      dim(wi) <- dim(theta)
      wi
    }
