vpaws <- function(y,
                  kstar = 16,
                  sigma2 = 1,
                  mask = NULL,
                  scorr = 0,
                  spmin = 0.25,
                  ladjust = 1,
                  wghts = NULL,
                  u = NULL,
                  maxni = FALSE,
                  patchsize = 1) {
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
  if (is.null(mask))
    mask <- rep(TRUE, n)
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
  zobj <- list(bi = rep(1, n), theta = y)
  bi <- zobj$bi
  cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  np1 <- 2 * patchsize + 1
  np2 <- if (n2 > 1)
    2 * patchsize + 1
  else
    1
  np3 <- if (n3 > 1)
    2 * patchsize + 1
  else
    1
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
    zobj <- .Fortran(C_pvaws,
      as.double(y),
      as.logical(mask),
      as.integer(nvec),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      bi = as.double(zobj$bi),
      theta = double(nvec * n),
      as.integer(mc.cores),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts),
      double(nvec * mc.cores),
      as.integer(np1),
      as.integer(np2),
      as.integer(np3),
      double(nvec * np1 * np2 * np3 * mc.cores),
      double(np1 * np2 * np3 * mc.cores)
    )[c("bi", "theta", "hakt")]
    dim(zobj$theta) <- c(nvec, dy)
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    dim(zobj$bi) <- dy
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
  #   list(y=y,theta=zobj$theta,hakt=hakt,sigma2=sigma2,lambda=lambda,
  #        ladjust=ladjust,args=args,wghts=wghts,mae=mae,ni=zobj$bi)
  awsobj(
    y,
    zobj$theta,
    sigma2 / zobj$bi,
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
    ni = zobj$bi
  )
}

vpawscov <- function(y,
                      kstar = 16,
                      invcov = NULL,
                      mask = NULL,
                      scorr = 0,
                      spmin = 0.25,
                      ladjust = 1,
                      wghts = NULL,
                      u = NULL,
                      maxni = FALSE,
                      patchsize = 1) {
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
  if (is.null(mask))
    mask <- rep(TRUE, n)
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
  zobj <- list(bi = rep(1, n), theta = y)
  bi <- zobj$bi
  cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  np1 <- 2 * patchsize + 1
  np2 <- if (n2 > 1)
    2 * patchsize + 1
  else
    1
  np3 <- if (n3 > 1)
    2 * patchsize + 1
  else
    1
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
    zobj <- .Fortran(C_pvaws2,
      as.double(y),
      as.logical(mask),
      as.integer(nvec),
      as.integer(nvec * (nvec + 1) / 2),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      bi = as.double(zobj$bi),
      theta = double(nvec * n),
      as.double(invcov),
      as.integer(mc.cores),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts),
      double(nvec * mc.cores),
      as.integer(np1),
      as.integer(np2),
      as.integer(np3),
      double(nvec * np1 * np2 * np3 * mc.cores),
      double(nvec * (nvec + 1) / 2 * np1 * np2 * np3 *
               mc.cores),
      double(np1 * np2 * np3 * mc.cores)
    )[c("bi", "theta", "hakt")]
    dim(zobj$theta) <- c(nvec, dy)
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    dim(zobj$bi) <- dy
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
  #   list(y=y,theta=zobj$theta,hakt=hakt,invcov=invcov,lambda=lambda,
  #        ladjust=ladjust,args=args,wghts=wghts,mae=mae,ni=zobj$bi)
  awsobj(
    y,
    zobj$theta,
    sweep(invcov, 2:(d + 1), zobj$bi, "*"),
    hakt,
    invcov,
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
    ni = zobj$bi
  )
}
