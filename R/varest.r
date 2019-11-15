estimateSigmaCompl <- function(magnitude,phase,mask,kstar=20,kmin=8,hsig=5,lambda=12,verbose=TRUE){
  ## kmin = 10 corresponds to an initial bandwidth of 1.47 giving positive weight to direct neighbors and
  ## 2D diagonal neigbors
  args <- sys.call(-1)
  sdim <- dim(mask)
#  if(!is.numeric(magnitude)){
#    if (verbose) cat("reading Magnitude file ... ")
#    R <- readNIfTI(magnitude, reorient = FALSE)
#  } else {
    R <- magnitude
#  }
#  if(!is.numeric(phase)){
#    if (verbose) cat("reading Phase file ... ")
#    Ph <- readNIfTI(phase, reorient = FALSE)
#  } else {
    Ph <- phase
#  }
  ComplImg <- array(0,c(2,sdim))
  ComplImg[1,,,] <- R*cos(Ph)
  ComplImg[2,,,] <- R*sin(Ph)
  ## find the number of usable cores
  mc.cores <- setCores(, reprt = FALSE)
  ##
  ##  start smoothing and variance estimation
  ##
  n <- prod(sdim)
  lambda0 <- 1e40
  sigma2 <- array(1e10,sdim)
  # just inilitialize with something large, first step is nonadaptive due to lambda0
  k <- kmin
  hmax <- 1.25^(kstar/3)
  ## preparations for median smoothing
  nwmd <- (2*as.integer(hsig)+1)^3
  parammd <- .Fortran(C_paramw3,
                      as.double(hsig),
                      as.double(c(1,1)),
                      ind=integer(3*nwmd),
                      w=double(nwmd),
                      n=as.integer(nwmd))[c("ind","w","n")]
  nwmd <- parammd$n
  parammd$ind <- parammd$ind[1:(3*nwmd)]
  dim(parammd$ind) <- c(3,nwmd)

  if (verbose) pb <- txtProgressBar(min = 0, max = kstar-kmin+1, style = 3)
  bi <- array(1,sdim)
  zobj <- list(theta=ComplImg, bi=bi)
  if (verbose) {
    mae <- NULL
    protocol <- matrix("", kstar-kmin+1, 1, dimnames = list(paste("step", kmin:kstar), "protocol"))
  }
  while (k <= kstar) {
    ## determine the actual bandwidth for this step
    hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, c(1,1), 1e-4)

    ## we need the (approx.) size of the weigthing scheme array
    dlw <- (2*trunc(hakt/c(1, 1, 1))+1)[1:3]

    ## perform the actual adaptive smoothing
    zobj <- .Fortran(C_vaws2cov,
                     as.double(ComplImg),
                     as.logical(mask),
                     as.integer(2),
                     as.integer(sdim[1]),
                     as.integer(sdim[2]),
                     as.integer(sdim[3]),
                     hakt = as.double(hakt),
                     as.double(lambda0),
                     as.double(zobj$theta),
                     as.double(sigma2),
                     bi = as.double(zobj$bi),
                     theta = double(2*n),
                     sigma2 = double(n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(c(1,1)),
                     double(2 * mc.cores))[c("bi", "theta", "hakt","sigma2")]
    ##
    ##  now get local median variance estimates
    ##
    dim(zobj$sigma2) <- sdim
    sigma2 <- .Fortran(C_mediansm,
                       as.double(zobj$sigma2),
                       as.logical(mask),
                       as.integer(sdim[1]),
                       as.integer(sdim[2]),
                       as.integer(sdim[3]),
                       as.integer(parammd$ind),
                       as.integer(nwmd),
                       double(nwmd*mc.cores), # work(nw,nthreds)
                       as.integer(mc.cores),
                       sigma2n = double(n))$sigma2n/0.6931
    # sigma2n containes sum of 2 independent squared residuals
    # 0.6931 approximates  median \chi_2 /2
    # needed to get correct results
    ## use maximum ni
    bi <- zobj$bi <- pmax(bi, zobj$bi)

    ## some verbose stuff
    if (verbose) {
      protocol[k-kmin+1,1] <- paste("bandwidth: ", signif(hakt, 3),
                                    "sigma: mean: ", signif(sqrt(mean(sigma2[mask])),3),
                                    "median: ", signif(sqrt(median(sigma2[mask])),3),
                                    "sd: ", signif(sd(sqrt(sigma2[mask])),3),
                                    "median(bi):", signif(median(zobj$bi[mask]),3),
                                    "max(bi):", signif(max(zobj$bi[mask]),3))
      setTxtProgressBar(pb, k-kmin+1)
    }

    ## go for next iteration
    k <- k+1
    lambda0 <- lambda
    gc()
  }
  dim(zobj$theta) <- c(2,sdim)
  # return estimated parameters of rician distribution
  list(sigma=array(sqrt(sigma2),sdim),
       theta=array(sqrt(zobj$theta[1,,,]^2+zobj$theta[2,,,]^2),sdim),
       sigmal=array(sqrt(zobj$sigma2),sdim),mask=mask,
       protocol=protocol,args=args)
}

medianFilter3D <- function(sigma2, hsig=10, mask=NULL){
  sdim <- dim(sigma2)
  n <- prod(sdim)
  if(length(sdim)!=3) stop("obj needs to be of class 'array' (3D) or 'sigmaEstSENSE'")
  if(is.null(mask)) mask <- array(TRUE,sdim)
  if(any(dim(mask)!=sdim)) stop("dimensions do not coinside")
  nwmd <- (2*as.integer(hsig)+1)^3
  parammd <- .Fortran(C_paramw3,
                      as.double(hsig),
                      as.double(c(1,1)),
                      ind=integer(3*nwmd),
                      w=double(nwmd),
                      n=as.integer(nwmd))[c("ind","w","n")]
  nwmd <- parammd$n
  parammd$ind <- parammd$ind[1:(3*nwmd)]
  dim(parammd$ind) <- c(3,nwmd)
  mc.cores <- setCores(, reprt = FALSE)
  sigma2hat <- .Fortran(C_mediansm,
                     as.double(sigma2),
                     as.logical(mask),
                     as.integer(sdim[1]),
                     as.integer(sdim[2]),
                     as.integer(sdim[3]),
                     as.integer(parammd$ind),
                     as.integer(nwmd),
                     double(nwmd*mc.cores), # work(nw,nthreds)
                     as.integer(mc.cores),
                     sigma2n = double(n))$sigma2n/0.6931
  dim(sigma2hat) <- sdim
  sigma2hat
}
