#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and
#    Volatility models
#
#    emphazises on the propagation-separation approach
#
#    Copyright (C) 2006 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#
#     default parameters:  see function setawsdefaults
#
aws.segment <- function(y,
                        level,
                        delta = 0,
                        hmax = NULL,
                        hpre = NULL,
                        mask = NULL,
                        varmodel = "Constant",
                        lkern = "Triangle",
                        scorr = 0,
                        ladjust = 1,
                        wghts = NULL,
                        u = NULL,
                        varprop = .1,
                        ext = 0,
                        graph = FALSE,
                        demo = FALSE,
                        fov = NULL)
{
  setawsthresh <- function(d, kstar, ladjust, ext) {
    ladjust <- min(ladjust, switch(d, 1.68, 2.37, 2.44))
    #  for ladjust>=switch(d,1.68,2.37,2.44) use nonadaptive thresholds
    switch(
      d,
      1.65  + 0.1952 * log(kstar) - 0.0473 * ladjust - 0.6771 * ext,
      1.729 + 0.2831 * log(kstar) - 0.2715 * ladjust - 0.4576 * ext,
      1.696 + 0.4010 * log(kstar) - 0.2975 * ladjust - 0.4273 * ext
    )
  }
  #    first check arguments and initialize
  #
  args <- match.call()
  n <- length(y)
  dy <- dim(y)
  if (length(dy) > 3)
    stop("AWS for more than 3 dimensional grids is not implemented")
  if (!(varmodel %in% c("Constant", "Linear", "Quadratic")))
    stop("Model for variance not implemented")
  #
  #   check for segmentation parameters
  #
  if (is.null(level) ||
      is.null(delta) || level + delta > max(y) || level - delta < min(y)) {
    stop(
      paste(
        "Inproper specifications for level ",
        level,
        " or delta ",
        delta,
        "\n Values specified outside range of y"
      )
    )
  }
  #
  #   set appropriate defaults
  #
  if (is.null(dy)) {
    d <- 1
  } else {
    d <- length(dy)
  }
  if (is.null(wghts))
    wghts <- c(1, 1, 1)
  wghts <-
    switch(d, c(0, 0), c(wghts[1] / wghts[2], 0), wghts[1] / wghts[2:3])
    if (is.null(mask)) {
      if (is.null(dy))
        mask <- rep(TRUE, length(y))
      else
        mask <- array(TRUE, dy)
    } else {
  ## these things need full data cubes
      u <- NULL
      graph <- FALSE
      demo <- FALSE
    }
    dmask <- dim(mask)
    if(is.null(dmask)) dmask <- length(mask)
    nvoxel <- sum(mask)
    position <- array(0,dmask)
    position[mask] <- 1:nvoxel
    # reduce to voxel in mask
    y <- y[mask]
  cpar <-
    setawsdefaults(dy,
                   mean(y),
                   "Gaussian",
                   lkern,
                   "Uniform",
                   TRUE,
                   FALSE,
                   ladjust,
                   hmax,
                   1,
                   wghts)
  lkern <- cpar$lkern
  lambda <- 1.25 * cpar$lambda # Gaussian case
  maxvol <- cpar$maxvol
  k <- cpar$k
  kstar <- cpar$kstar
  cpar$tau1 <- cpar$tau1 * 2
  cpar$tau2 <- cpar$tau2 * 2
  hmax <- cpar$hmax
  shape <- cpar$shape
  if (lkern == 5) {
    #  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax * 0.42445 * 4
  }
  #
  #    set threshold
  #
  thresh <- setawsthresh(d, kstar, ladjust, ext)
  beta <- switch(d, 1.06, 1.42, 1.33)
  #  optimised for sample sizes
  #     d=1: n=c(2000,4000,8000)
  #     d=2: n=c(256^2,512^2,1024^2)=c(65536,262144,1048576)
  #     d=3: n=c(32^3,64^2*32,128^2*32)=c(32768,131072,524288)
  #
  #   family dependent transformations
  #
  zfamily <- awsgfamily(y, scorr, d)
  sigma2 <- zfamily$sigma2
  varest <- 1 / sigma2
  h0 <- zfamily$h0
  rm(zfamily)
  if (lkern == 5) {
    #  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax * 0.42445 * 4
    hinit <- 0.42445 * 4
  }
  if (demo && !graph)
    graph <- TRUE
  # now check which procedure is appropriate
  ##  this is the version on a grid
  n1 <- switch(d, n, dy[1], dy[1])
  n2 <- switch(d, 1, dy[2], dy[2])
  n3 <- switch(d, 1, 1, dy[3])
  if (is.null(fov))
    fov <- nvoxel
  interval <-
    if (delta > 0)
      paste("Central interval: (", level - delta, ",", level + delta, ")")
  else
    paste("Level: ", level)
  cat(
    "Running segmentation algorithm with following parameters:\n",
    interval,
    "\n",
    "Critical value: ",
    thresh,
    "  Extension: ",
    ext,
    "  Field of view: ",
    fov,
    "\n"
  )
  #
  #    Initialize  for the iteration
  #
  fix <- rep(FALSE, nvoxel)
  zobj <-
    list(
      ai = y,
      bi0 = rep(1, nvoxel),
      bi2 = rep(1, nvoxel),
      bi = rep(1, nvoxel),
      theta = y / shape,
      fix = fix
    )
  segment <- rep(0,nvoxel)
  gi <- gi2 <- rep(1,nvoxel)
  mae <- NULL
  lambda0 <-
    1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
  #
  #   produce a presmoothed estimate to stabilize variance estimates
  #
  if (is.null(hpre))
    hpre <- 20 ^ (1 / d)
  dlw <- (2 * trunc(hpre / c(1, wghts)) + 1)[1:d]
  hobj <- .Fortran(C_caws,
    as.double(y),
    as.integer(position),
    as.integer(n1),
    as.integer(n2),
    as.integer(n3),
    as.double(hpre),
    as.double(1e40),
    as.double(zobj$theta),
    bi = double(nvoxel),
    double(nvoxel),
    as.double(zobj$bi0),
    ai = double(nvoxel),
    as.integer(cpar$mcode),
    as.integer(lkern),
    as.double(0.25),
    double(prod(dlw)),
    as.double(wghts)
  )[c("bi", "ai")]
  hobj$theta <- hobj$ai / hobj$bi
  #
  #   iteratate until maximal bandwidth is reached
  #
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  while (k <= kstar) {
    hakt0 <- gethani(1, 10, lkern, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 10, lkern, 1.25 ^ k, wghts, 1e-4)
    if (lkern == 5) {
      #  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
      hakt <- hakt * 0.42445 * 4
    }
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)[1:d]
    if (scorr[1] >= 0.1)
      lambda0 <-
      lambda0 * Spatialvar.gauss(hakt0 / 0.42445 / 4, h0, d) / Spatialvar.gauss(hakt0 /
                                                                                  0.42445 / 4, 1e-5, d)
    # Correction for spatial correlation depends on h^{(k)}
    hakt0 <- hakt
    # heteroskedastic Gaussian case
    zobj <- .Fortran(C_segment,
      as.double(y),
      as.integer(position),
      fix=as.integer(fix),
      as.double(level),
      as.double(delta),
      as.double(sigma2),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      bi = as.double(zobj$bi),
      bi2 = double(nvoxel),
      bi0 = as.double(zobj$bi0),
      gi = double(nvoxel),
      gi2 = double(nvoxel),
      theta = as.double(zobj$theta),
      as.integer(lkern),
      as.double(0.25),
      double(prod(dlw)),
      as.double(wghts),
      as.integer(segment),
      # previous segmentation array
      segment = as.integer(segment),
      # new array for segment (-1,0,1)
      as.double(beta),
      as.double(thresh),
      as.double(ext),
      as.double(fov),
      varest = as.double(varest)
    )[c("fix",
        "bi",
        "bi0",
        "bi2",
        "gi2",
        "segment",
        "theta",
        "gi",
        "hakt",
        "varest")]
    gi[!fix] <- zobj$gi[!fix]
    gi2[!fix] <- zobj$gi2[!fix]
    if (hakt > n1 / 2)
      zobj$bi0 <- rep(max(zobj$bi), nvoxel)
    segment <- zobj$segment
    varest <- zobj$varest
    fix <- as.logical(zobj$fix)
    if (graph) {
        dim(zobj$theta) <-
        dim(gi) <- dim(segment) <- dim(zobj$bi) <- dmask
      #
      #     Display intermediate results if graph == TRUE
      #
      if (d == 1) {
        oldpar <- par(
          mfrow = c(1, 3),
          mar = c(3, 3, 3, .2),
          mgp = c(2, 1, 0)
        )
        plot(y, ylim = range(y, zobj$theta), col = 3)
        if (!is.null(u))
          lines(u, col = 2)
        lines(zobj$theta, lwd = 2)
        title(paste("Reconstruction  h=", signif(hakt, 3)))
        plot(
          segment,
          type = "l",
          main = "Segmentation result",
          ylim = c(-1, 1)
        )
        plot(zobj$bi, type = "l", ylim = range(0, zobj$bi))
        title("Sum of weights")
      }
      if (d == 2) {
        oldpar <- par(
          mfrow = c(2, 2),
          mar = c(1, 1, 3, .25),
          mgp = c(2, 1, 0)
        )
        image(array(y,dy),
              col = grey((0:255) / 255),
              xaxt = "n",
              yaxt = "n")
        title(paste(
          "Observed Image  min=",
          signif(min(y), 3),
          " max=",
          signif(max(y), 3)
        ))
        image(
          array(zobj$theta,dy),
          col = grey((0:255) / 255),
          xaxt = "n",
          yaxt = "n"
        )
        title(paste(
          "Reconstruction  h=",
          signif(hakt, 3),
          " min=",
          signif(min(zobj$theta), 3),
          " max=",
          signif(max(zobj$theta), 3)
        ))
        image(
          array(zobj$gi,dy),
          col = grey((0:255) / 255),
          xaxt = "n",
          yaxt = "n"
        )
        title(paste(
          "Sum of weights: min=",
          signif(min(zobj$gi), 3),
          " mean=",
          signif(mean(zobj$gi), 3),
          " max=",
          signif(max(zobj$gi), 3)
        ))
        image(
          array(segment,dy),
          col = grey((0:255) / 255),
          xaxt = "n",
          yaxt = "n",
          zlim = c(-1, 1)
        )
        title("Segmentation result")
      }
      if (d == 3) {
        oldpar <- par(
          mfrow = c(2, 2),
          mar = c(1, 1, 3, .25),
          mgp = c(2, 1, 0)
        )
        image(array(y,dy)[, , n3 %/% 2 + 1],
              col = grey((0:255) / 255),
              xaxt = "n",
              yaxt = "n")
        title(paste(
          "Observed Image  min=",
          signif(min(y), 3),
          " max=",
          signif(max(y), 3)
        ))
        image(
          array(zobj$theta,dy)[, , n3 %/% 2 + 1],
          col = grey((0:255) / 255),
          xaxt = "n",
          yaxt = "n"
        )
        title(paste(
          "Reconstruction  h=",
          signif(hakt, 3),
          " min=",
          signif(min(zobj$theta), 3),
          " max=",
          signif(max(zobj$theta), 3)
        ))
        image(
          array(zobj$bi,dy)[, , n3 %/% 2 + 1],
          col = grey((0:255) / 255),
          xaxt = "n",
          yaxt = "n"
        )
        title(paste(
          "Sum of weights: min=",
          signif(min(zobj$bi), 3),
          " mean=",
          signif(mean(zobj$bi), 3),
          " max=",
          signif(max(zobj$bi), 3)
        ))
        image(
          array(segment,dy)[, , n3 %/% 2 + 1],
          col = grey((0:255) / 255),
          xaxt = "n",
          yaxt = "n",
          zlim = c(-1, 1)
        )
        title("Segmentation result")
      }
      par(oldpar)
    }
    #
    #    Calculate MAE and MSE if true parameters are given in u
    #    this is for demonstration and testing for propagation (parameter adjustments)
    #    only.
    #
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
    if (demo)
      readline("Press return")
    #
    #   Prepare for next iteration
    #
    #
    #   Create new variance estimate
    #
    if (sum(zobj$fix) < nvoxel / 4) {
      vobj <-
        awsgsigma2(y, hobj, zobj[c("theta", "gi", "gi2")], varmodel, varprop)
      sigma2 <- vobj$sigma2inv
      coef <- vobj$coef
      rm(vobj)
    }
    x <- 1.25 ^ (k - 1)
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (max(total) > 0) {
      cat(
        "step:",
        k,
        "  hakt=",
        hakt,
        "  progress:",
        signif(total[k], 2) * 100,
        "\n",
        sep = ""
      )
    }
    k <- k + 1
    gc()
  }
  cat("\n")
  ###
  ###            end iterations now prepare results
  ###
  ###   component var contains an estimate of Var(zobj$theta)
  ###
  vartheta <- y0 <- theta <- segment <- sigma0 <- array(0, dmask)
  vartheta[mask] <- zobj$bi2 / zobj$bi ^ 2
  vartheta <-
    vartheta / Spatialvar.gauss(hakt / 0.42445 / 4, h0 + 1e-5, d) * Spatialvar.gauss(hakt /
                                                    0.42445 / 4, 1e-5, d)
  sigma0[mask] <- 1/sigma2
  y0[mask] <- y
  theta[mask] <- zobj$theta
  segment[mask] <- zobj$segment
  awssegmentobj(
    y0,
    zobj$theta,
    segment,
    vartheta,
    level,
    delta,
    hakt,
    sigma0,
    lkern,
    lambda,
    ladjust,
    TRUE,
    FALSE,
    args,
    homogen = FALSE,
    earlystop = FALSE,
    family = "Gaussian",
    wghts = wghts,
    scorr = scorr,
    varmodel = varmodel,
    vcoef = coef,
    mae = mae
  )
}
