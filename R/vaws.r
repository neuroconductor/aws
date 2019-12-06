vaws <- function(y,
                 kstar = 16,
                 sigma2 = 1,
                 mask = NULL,
                 scorr = 0,
                 spmin = 0.25,
                 ladjust = 1,
                 wghts = NULL,
                 u = NULL,
                 maxni = FALSE) {
  args <- match.call()
  dy <- dim(y)
  nvec <- dy[1]
  dy <- dy[-1]
  d <- length(dy)
  if (length(dy) > 3)
    stop("Vector AWS for more than 3 dimensional grids is not implemented")
  lambda <-
    2 * sigma2 * ladjust * switch(d,
                                  qchisq(pchisq(14.6, 1), nvec),
                                  ## 1D
                                  qchisq(pchisq(9.72, 1), nvec),
                                  ## 2D
                                  qchisq(pchisq(8.82, 1), nvec))## 3D
  if (is.null(wghts))
    wghts <- c(1, 1, 1)
  wghts <-
    switch(length(dy), c(0, 0), c(wghts[1] / wghts[2], 0), wghts[1] / wghts[2:3])
  n1 <- switch(d, dy, dy[1], dy[1])
  n2 <- switch(d, 1, dy[2], dy[2])
  n3 <- switch(d, 1, 1, dy[3])
  n <- n1 * n2 * n3
  if (is.null(mask)) {
      mask <- array(TRUE, dy)
  }
  dmask <- dim(mask)
  if(is.null(dmask)) dmask <- length(mask)
  nvoxel <- sum(mask)
  position <- array(0,dmask)
  position[mask] <- 1:nvoxel
  dim(y) <- c(nvec,n)
  y <- y[,mask]
  if(!is.null(u)){
     dim(u) <- c(nvec,n)
     u <- u[,mask]
  }
  # reduce to voxel in mask
  h0 <- 0
  if (any(scorr > 0)) {
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0))
      h0[i] <- geth.gauss(scorr[i])
    if (length(h0) < d)
      h0 <- rep(h0[1], d)
    cat("Corresponding bandwiths for specified correlation:",
        h0,
        "\n")
  }
  hseq <- 1
  zobj <- list(bi = rep(1, voxel), theta = y)
  bi <- zobj$bi
  cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  k <- 1
  hmax <- 1.25 ^ (kstar / d)
  lambda0 <- lambda
  mae <- NULL
  while (k <= kstar) {
    hakt0 <- gethani(1, 1.25 * hmax, 2, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4)
    cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
    hseq <- c(hseq, hakt)
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)[1:d]
    if (scorr[1] >= 0.1)
      lambda0 <-
      lambda0 * Spatialvar.gauss(hakt0 / 0.42445 / 4, h0, d) / Spatialvar.gauss(hakt0 /
                                                                                  0.42445 / 4, 1e-5, d)
    zobj <- .Fortran(C_vaws,
      as.double(y),
      as.integer(pos),
      as.integer(nvec),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      bi = as.double(zobj$bi),
      vred = double(nvoxel),
      theta = double(nvec * nvoxel),
      as.integer(mc.cores),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts),
      double(nvec * mc.cores)
    )[c("bi", "theta", "hakt","vred")]
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    if (!is.null(u)) {
      cat(
        "bandwidth: ",
        signif(hakt, 3),
        "   MSE: ",
        signif(mean((zobj$theta - u) ^ 2), 3),
        "   MAE: ",
        signif(mean(abs(zobj$theta - u)), 3),
        " mean(bi)=",
        signif(mean(zobj$bi), 3),
        "\n"
      )
      mae <- c(mae, signif(mean(abs(zobj$theta - u)), 3))
    }
    x <- 1.25 ^ k
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (max(total) > 0) {
      cat(signif(total[k], 2) * 100, "% . ", sep = "")
    }
    k <- k + 1
    gc()
  }
  cat("\n")
  y0 <- theta <- array(0,c(nvec,n))
  y0[,mask] <- y
  theta[,mask] <- zobj$theta
  dim(y0) <- dim(theta) <- c(nvec,dy)
  bi <- vred <- array(dy)
  bi[mask] <- zobj$bi
  vred[mask] <- zobj$vred
  rm(zobj)
  awsobj(
    y0,
    theta,
    sigma2 * vred,
    hakt,
    sigma2,
    lkern = 1L,
    lambda,
    ladjust,
    aws = TRUE,
    memory = FALSE,
    args,
    hseq = hseq,
    homogen = FALSE,
    earlystop = FALSE,
    family = "Gaussian",
    wghts = wghts,
    mae = mae,
    psnr = NULL,
    ni = bi
  )
}
vawscov <- function(y,
                      kstar = 16,
                      invcov = NULL,
                      mask = NULL,
                      scorr = 0,
                      spmin = 0.25,
                      ladjust = 1,
                      wghts = NULL,
                      u = NULL,
                      maxni = FALSE) {
  ##
  ##  this is the version with full size invcov (triangular storage)
  ##  needed for MPM
  ##
  args <- match.call()
  dy <- dim(y)
  nvec <- dy[1]
  dy <- dy[-1]
  d <- length(dy)
  if (length(dy) > 3)
    stop("Vector AWS for more than 3 dimensional grids is not implemented")
  lambda <- 2 * ladjust * switch(d,
                                 qchisq(pchisq(14.6, 1), nvec),
                                 ## 1D
                                 qchisq(pchisq(9.72, 1), nvec),
                                 ## 2D
                                 qchisq(pchisq(8.82, 1), nvec))## 3D
  if (is.null(wghts))
    wghts <- c(1, 1, 1)
  wghts <-
    switch(length(dy), c(0, 0), c(wghts[1] / wghts[2], 0), wghts[1] / wghts[2:3])
  n1 <- switch(d, dy, dy[1], dy[1])
  n2 <- switch(d, 1, dy[2], dy[2])
  n3 <- switch(d, 1, 1, dy[3])
  n <- n1 * n2 * n3
  if (is.null(mask)) {
      mask <- array(TRUE, dy)
  }
  dmask <- dim(mask)
  if(is.null(dmask)) dmask <- length(mask)
  nvoxel <- sum(mask)
  position <- array(0,dmask)
  position[mask] <- 1:nvoxel
  dim(y) <- c(nvec,n)
  y <- y[,mask]
  nvd <- nvec * (nvec + 1) / 2
  dim(invcov) <- c(nvd,n)
  invcov <- invcov[,mask]
  if(!is.null(u)){
     dim(u) <- c(nvec,n)
     u <- u[,mask]
  }
  h0 <- 0
  if (any(scorr > 0)) {
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0))
      h0[i] <- geth.gauss(scorr[i])
    if (length(h0) < d)
      h0 <- rep(h0[1], d)
    cat("Corresponding bandwiths for specified correlation:",
        h0,
        "\n")
  }
  hseq <- 1
  zobj <- list(bi = rep(1, nvoxel), theta = y)
  bi <- zobj$bi
  cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  k <- 1
  hmax <- 1.25 ^ (kstar / d)
  lambda0 <- lambda
  mae <- NULL
  while (k <= kstar) {
    hakt0 <- gethani(1, 1.25 * hmax, 2, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4)
    cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
    hseq <- c(hseq, hakt)
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)[1:d]
    if (scorr[1] >= 0.1)
      lambda0 <-
      lambda0 * Spatialvar.gauss(hakt0 / 0.42445 / 4, h0, d) / Spatialvar.gauss(hakt0 /
                                                                                  0.42445 / 4, 1e-5, d)
    zobj <- .Fortran(C_vaws2cov,
      as.double(y),
      as.integer(mask),
      as.integer(nvec),
      as.integer(nvec * (nvec + 1) / 2),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      bi = as.double(zobj$bi),
      vred = double(nvoxel),
      theta = double(nvec * nvoxel),
      as.double(invcov),
      as.integer(mc.cores),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts),
      double(nvec * mc.cores),
      double(nvec * mc.cores),
      double(nvec * (nvec + 1) / 2 * mc.cores)
    )[c("bi", "theta", "hakt", "vred")]
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    if (!is.null(u)) {
      cat(
        "bandwidth: ",
        signif(hakt, 3),
        "   MSE: ",
        signif(mean((zobj$theta - u) ^ 2), 3),
        "   MAE: ",
        signif(mean(abs(zobj$theta - u)), 3),
        " mean(bi)=",
        signif(mean(zobj$bi), 3),
        "\n"
      )
      mae <- c(mae, signif(mean(abs(zobj$theta - u)), 3))
    }
    x <- 1.25 ^ k
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (max(total) > 0) {
      cat(signif(total[k], 2) * 100, "% . ", sep = "")
    }
    k <- k + 1
    gc()
  }
  cat("\n")
  theta0 <- array(0,c(nvec,n))
  dim(theta) <- c(nvec,n)
  theta[,mask] <- zobj$theta
  dim(theta) <- c(nvec,dy)
  bi <- vred <- array(dy)
  bi[mask] <- zobj$bi
  vred[mask] <- zobj$vred
  rm(zobj)
  y0  <- array(0,c(nvec,n))
  y0[,mask] <- y
  dim(y0) <- c(nvec,dy)
  rm(y)
  invcov0 <- array(0,c(nvd,n))
  invcov0[,mask] <- invcov
  dim(invcov0) <- c(nvd,dy)
  rm(invcov)
  awsobj(
    y0,
    theta,
    sweep(invcov, 2:(d + 1), vred, "*"),
    hakt,
    invcov0,
    lkern = 1L,
    lambda,
    ladjust,
    aws = TRUE,
    memory = FALSE,
    args,
    hseq = hseq,
    homogen = FALSE,
    earlystop = FALSE,
    family = "Gaussian",
    wghts = wghts,
    mae = mae,
    psnr = NULL,
    ni = bi
  )
}
