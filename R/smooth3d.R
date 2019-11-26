smooth3D <- function(y,lkern="Gaussian",weighted=FALSE,sigma2=NULL,mask=NULL,h=NULL,
                     wghts=NULL) {
  #
  #  3D nonadaptive smoothing using a mask
  #  y is assumed to only contain voxel within mask
  #
  d <- 3
  dmask <- dim(mask)
  nvoxel <- sum(mask)
  position <- array(0,dmask)
  position[mask] <- 1:nvoxel
  dy <- dim(y)
  if(is.null(dy)){
     ly <- length(y)
     nv <- 1
   } else {
     ly <- dy[1]
     nv <- dy[2]
   }
  if(ly!=nvoxel){
      stop("smooth3d: y should have length equal to sum(mask)")
  }
  if(is.null(sigma2)) {
      weighted <- FALSE
  } else {
    if(length(sigma2)!=nvoxel) weighted <- FALSE
    sigma2 <- 1/sigma2
  }
  if (is.null(h)) h <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
  }
  if (is.null(wghts)) wghts <- c(1,1,1)
  if(is.null(mask)) mask <- array(TRUE,dy[1:3])
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
  dlw <- (2*trunc(hmax/c(1,wghts))+1)[1:d]
  ysmooth <- .Fortran(C_smooth3d,
                     as.double(y),
                     as.double(sigma2),
                     as.integer(position),
                     as.integer(weighted),
                     as.integer(nvoxel),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(nv),
                     hakt=as.double(h),
                     thnew=double(nvoxel*dv),
                     as.integer(lkern),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(dv))$thnew
array(ysmooth,dy)
}

aws3Dmask <- function(y, mask, lambda, hmax, res=NULL, sigma2=NULL,
                      lkern="Gaussian", skern="Plateau", weighted=TRUE,
                      u=NULL, wghts=NULL, h0=c(0,0,0),
                      testprop=FALSE){
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  uses sum of weights and correction factor C(h,g) in statistical penalty
#
  dmask <- dim(mask)
  n1 <- dmask[1]
  n2 <- dmask[2]
  n3 <- dmask[3]
  nvoxel <- sum(mask)
  position <- array(0,dmask)
  position[mask] <- 1:nvoxel
  if (length(dy)!=nvoxel) {
    stop("aws3Dmask: y needs to have length sum(mask)")
  }
  # set the code for the kernel (used in lkern) and set lambda
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  skern <- switch(skern,
                  Triangle=2,
                  Plateau=1,
                  Exponential=3,
                  2)
  if(skern%in%c(1,2)) {
    # to have similar preformance compared to skern="Exp"
    lambda <- 4/3*lambda
    if(skern==1) spmin <- .3
    spmax <- 1
  } else {
      spmin <- 0
      spmax <- 4
  }
  # set hinit and hincr if not provided
  hinit <- 1

  # define hmax
  if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
    hinit <- min(hinit,hmax)
  }
  if (qlambda == 1) hinit <- hmax

  # define hincr
  hincr <- 1.25
  hincr <- hincr^(1/3)

  # determine corresponding bandwidth for specified correlation
  if(is.null(h0)) h0 <- rep(0,3)

  # estimate variance in the gaussian case if necessary
  if (is.null(sigma2)) {
    sigma2 <- IQRdiff(as.vector(y))^2
    if (any(h0>0)) sigma2<-sigma2*Varcor.gauss(h0)*Spatialvar.gauss(h0,1e-5,3)
  }
  if (length(sigma2)==1) sigma2 <- rep(sigma2,nvoxel)
  # deal with homoskedastic Gaussian case by extending sigma2
  if (length(sigma2)!=nvoxel) stop("sigma2 does not have length 1 or same length as y")
  lambda <- lambda*2
  sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes

  # Initialize  list for bi and theta
  if (is.null(wghts)) wghts <- c(1,1,1)
  hinit <- hinit/wghts[1]
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
  tobj <- list(bi= rep(1,nvoxel))
  theta <- y
  maxvol <- aws::getvofh(hmax,lkern,wghts)
  kstar <- as.integer(log(maxvol)/log(1.25))
  steps <- kstar+1

  if(lambda < 1e10) k <- 1 else k <- kstar
  hakt <- hinit
  hakt0 <- hinit
  lambda0 <- lambda
  if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
  scorr <- numeric(3)
  if(h0[1]>0) scorr[1] <-  get.corr.gauss(h0[1],2)
  if(h0[2]>0) scorr[2] <-  get.corr.gauss(h0[2],2)
  if(h0[3]>0) scorr[3] <-  get.corr.gauss(h0[3],2)
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  if(testprop) {
    if(is.null(u)) u <- 0
    propagation <- NULL
  }
  if(!is.null(u)) mae <- NULL
##
##  determine number of cores to use
##
  mc.cores <- setCores(,reprt=FALSE)

  # run single steps to display intermediate results
  while (k<=kstar) {
      hakt0 <- aws::gethani(1,3,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- aws::gethani(1,3,lkern,1.25^k,wghts,1e-4)
      hakt.oscale <- if(lkern==3) bw2fwhm(hakt/4) else hakt
      cat("step",k,"bandwidth",signif(hakt.oscale,3)," ")
    dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:3]
#  need bandwidth in voxel for Spaialvar.gauss, h0 is in voxel
    if (any(h0>0)) lambda0 <- lambda0 * Spatialvar.gauss(bw2fwhm(hakt0)/4/c(1,wghts),h0,3)/
      Spatialvar.gauss(h0,1e-5,3)/Spatialvar.gauss(bw2fwhm(hakt0)/4/c(1,wghts),1e-5,3)
        # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)}  all bandwidth-arguments in FWHM
    hakt0 <- hakt
    theta0 <- theta
    bi0 <- tobj$bi
    #
    #   need these values to compute variances after the last iteration
    #
    tobj <- .Fortran(C_chaws2,as.double(y),
                     as.double(sigma2),
                     as.integer(position),
                     as.integer(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(theta0),
                     bi=as.double(bi0),
                     thnew=double(nvoxel),
                     as.integer(lkern),
                     as.integer(skern),
                     as.double(spmin),
                     as.double(spmax),
                     double(prod(dlw)),
                     as.double(wghts))[c("bi","thnew","hakt")]
    gc()
    theta <- tobj$thnew
    dim(tobj$bi) <- dy
    if(testprop) {
      pobj <- .Fortran(C_chaws2,as.double(y),
                       as.double(sigma2),
                       as.integer(position),
                       as.integer(weighted),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(1e50),
                       as.double(theta0),
                       bi=as.double(bi0),
                       thnew=double(nvoxel),
                       as.integer(lkern),
		                   as.integer(skern),
	                     as.double(spmin),
		                   as.double(spmax),
		                   double(prod(dlw)),
		                   as.double(wghts))[c("bi","thnew")]
      gc()
      propagation <- c(propagation,sum(abs(theta-pobj$thnew))/sum(abs(pobj$thnew-u)))
      cat("Propagation with alpha=",max(propagation),"\n")
      cat("alpha values:","\n")
      print(rbind(hakt.oscale,signif(propagation[-1],3)))
    }
    if (!is.null(u)) {
      cat("bandwidth: ",signif(hakt.oscale,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
          signif(mean((theta-u)^2),3),"   MAE: ",signif(mean(abs(theta-u)),3),
          " mean(bi)=",signif(mean(tobj$bi),3),"\n")
      mae <- c(mae,signif(mean(abs(theta-u)),3))
    } else if (max(total) >0) {
      cat(signif(total[k],2)*100,"%                 \r",sep="")
     }
    k <- k+1
#  adjust lambda for the high intrinsic correlation between  neighboring estimates
    c1 <- (prod(h0+1))^(1/3)
    c1 <- 2.7214286 - 3.9476190*c1 + 1.6928571*c1*c1 - 0.1666667*c1*c1*c1
    x <- (prod(1.25^(k-1)/c(1,wghts)))^(1/3)
    scorrfactor <- (c1+x)/(c1*prod(h0+1)+x)
    lambda0 <- lambda*scorrfactor
    gc()
  }

  #   Now compute variance of theta and variance reduction factor (with respect to the spatially uncorrelated situation
  if(!is.null(residuals)){
    nres <- dim(residuals)[1]
    vartheta0 <- residualVariance(residuals, mask, resscale=1)
    residuals <- .Fortran(C_ihaws2,
                     as.double(residuals),
                     as.double(sigma2),
                     as.integer(position),
                     as.integer(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(nres),
                     hakt=as.double(tobj$hakt),
                     as.double(lambda0),
                     as.double(theta0),
                     as.integer(mc.cores),
                     as.double(bi0),
                     resnew=double(nvoxel*nres),
                     as.integer(lkern),
                     as.integer(skern),
                     as.double(spmin),
                     as.double(spmax),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(nres*mc.cores))$resnew
    dim(residuals) <- c(nres,nvoxel)
    vartheta <- residualVariance(residuals, mask, resscale=1)
    vred <- vartheta/vartheta0
    vartheta <- vred/sigma2
    lags <- pmin(c(5,5,3),dmask-1)
    scorr <- residualSpatialCorr(residuals,mask,lags)
  }
  list(theta=theta, ni=tobj$bi, var=vartheta, vred=vred, mae=mae,
       alpha=propagation, hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),
       scorr=scorr, res=residuals, mask=mask)
}
