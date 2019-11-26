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
