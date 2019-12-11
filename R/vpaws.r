vpaws <- function(y,
                  kstar = 16,
                  sigma2 = 1,
                  invcov = NULL,
                  mask = NULL,
                  scorr = 0,
                  spmin = 0.25,
                  ladjust = 1,
                  wghts = NULL,
                  u = NULL,
                  patchsize = 1) {
  args <- match.call()
  dy <- dim(y)
  nvec <- dy[1]
  n <- prod(dy[-1])
  if(is.null(mask)){
     dy <- dy[-1]
   } else {
     dy <- dim(mask)
     if(is.null(dy)) dy <- length(mask)
   }
   d <- length(dy)

   if (d > 3)
    stop("Vector AWS for more than 3 dimensional grids is not implemented")
  #
  #   set appropriate defaults
  #
  if(is.null(mask)) mask <- array(TRUE,dy)
  nvoxel <- sum(mask)
  position <- array(0,dy)
  position[mask] <- 1:nvoxel
  condensed <- n==nvoxel & n<prod(dy)
  if(!condensed & prod(dy)!=n) stop("incompatible mask and data")
  if(!condensed & !is.null(u)){
     if(length(dim(u))==d){
        u <- u[mask]
     } else if(length(dim(u))==(d+1)&dim(u)[1]==nvec){
        dim(u) <- c(nvec,prod(dy))
        u <- u[,mask]
     }
  }
  if(!is.null(invcov)) sigma2 <- 1
# set sigma2 to an uninformative value if invcov is given
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
  if(!condensed){
     dim(y) <- c(nvec,prod(dy))
     y <- y[,mask]
   }
   if(!is.null(invcov)){
     dinvcov <- dim(invcov)
     nvd <- dinvcov[1]
     if(nvd!=(nvec+1)*nvec/2) stop(paste("First dimension of invcov should be",
                     nvd,"(dense storage)"))
     if(condensed){
        if(dinvcov[2]!=nvoxel) stop("invcov should contain information for all voxel in mask")
     } else {
        if(prod(dinvcov[-1])!=n) stop("invcov should contain information for all voxel")
        dim(invcov) <- c(nvd,n)
        invcov <- invcov[,mask]
     }
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
      lambda0 * Spatialvar.gauss(hakt0 / 0.42445 / 4, h0, d) /
                Spatialvar.gauss(hakt0 / 0.42445 / 4, 1e-5, d)
    if(is.null(invcov)){
    zobj <- .Fortran(C_pvaws,
      as.double(y),
      as.integer(position),
      as.integer(nvec),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      as.double(zobj$bi),
      bi = double(nvoxel),
      theta = double(nvec * nvoxel),
      as.integer(mc.cores),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts),
      double(nvec * mc.cores),
      as.integer(np1),
      as.integer(np2),
      as.integer(np3))[c("bi", "theta", "hakt")]
    } else {
      zobj <- .Fortran(C_pvaws2,
        as.double(y),
        as.integer(position),
        as.integer(nvec),
        as.integer(nvd),
        as.integer(n1),
        as.integer(n2),
        as.integer(n3),
        hakt = as.double(hakt),
        as.double(lambda0),
        as.double(zobj$theta),
        as.double(zobj$bi),
        bi = double(nvoxel), #binn
        theta = double(nvec * nvoxel),
        as.double(invcov),
        as.integer(mc.cores),
        as.double(spmin),
        double(prod(dlw)),
        as.double(wghts),
        double(nvec * mc.cores),
        as.integer(np1),
        as.integer(np2),
        as.integer(np3))[c("bi", "theta", "hakt")]
}
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
  sigma2new <- if(is.null(invcov)) zobj$bi/sigma2 else
     sweep(invcov, 2, zobj$bi, "*")
# this is the new inverse covariance
  if(!condensed){
#     expand results
     y0 <- y
     y <- array(0,c(nvec,prod(dy)))
     y[,mask] <- y0
     rm(y0)
     dim(y) <- c(nvec,dy)
     theta <- array(0,c(nvec,prod(dy)))
     theta[,mask] <- zobj$theta
     dim(theta) <- c(nvec,dy)
     if(!is.null(invcov)){
       sigma20 <- sigma2new
       sigma2new <- array(0,c(nvd,prod(dy)))
       sigma2new[,mask] <- sigma20
       rm(sigma20)
       dim(sigma2new) <- c(nvd,dy)
     }
     ni <- array(0,dy)
     ni[mask] <- zobj$bi
} else {
#     set dimensions
     dim(y) <- c(nvec,nvoxel)
     theta <- array(zobj$theta,c(nvec,nvoxel))
     ni <- array(zobj$bi,nvoxel)
}
  awsobj(
    y,
    theta,
    sigma2new,
    hakt,
    sigma2,
    lkern = 1L,
    lambda,
    ladjust,
    aws = TRUE,
    nvec = as.integer(nvec),
    memory = FALSE,
    args,
    hseq = hseq,
    homogen = FALSE,
    mask = mask,
    earlystop = FALSE,
    family = "Gaussian",
    wghts = wghts,
    mae = mae,
    psnr = NULL,
    ni = ni
  )
}
