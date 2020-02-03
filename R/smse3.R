smse3 <- function(sb, s0, bv, grad, mask, sigma, kstar, lambda, kappa0,
                  ns0=1, vext = NULL, vred = 4,
                  ncoils = 1,
                  model = 0,
                  dist = 1,
                  verbose = FALSE){
  # dist determines distance on sphere (can take 1:3), see getkappas
  #  data need to be scaled by sigma
  if(is.null(kappa0)){
    #  select kappa based on variance reduction on the sphere
    if(is.null(vred)||!is.numeric(vred)||vred<1) {
      stop("snse3 You need to specify either kappa0 or vred")
    }
    kappa0 <- suggestkappa(grad,vred,dist)$kappa
  }
  if(model>=2) varstats <- sofmchi(ncoils)
  #
  #  model=0  uses approx of noncentral Chi and smoothes Y
  #  model=1  uses approx of noncentral Chi2 and smoothes Y^2
  #  model=2  uses Gaussian approx of noncentral Chi2 and smoothes Y^2
  #  model=3  uses Gaussian approx of noncentral Chi and smoothes Y
  multishell <- sd(bv) > mean(bv)/50
  ngrad <- length(bv)
    if(multishell) {
    msstructure <- getnext3g(grad,bv)
    model <- 3
    nshell <- msstructure$nbv
  }
  ddim <- dim(mask)
  nvoxel <- sum(mask)
  nbv <- length(bv)
  if(dim(sb)[1]!=nvoxel || length(s0)!=nvoxel || dim(sb)[2]!=nbv){
     stop("smse3 - sb and s0 should only contain data within mask")
  }
  # define position of voxel in mask within the 3D cube
      position <- array(0,ddim)
      position[mask] <- 1:nvoxel
  # set relation between voxel lengths
      if(is.null(vext)) vext <- c(1,1)
  ##
  ##  rescale
  ##
  sb <- sb/sigma
  s0 <- s0/sigma
  if(model==1){
    #
    #   use squared values for Chi^2
    #
    sb <- sb^2
    s0 <- s0^2
    sigma <- sigma^2
  }
  th0 <- s0
  ni0 <- rep(1,nvoxel)
  if(multishell){
    gradstats <- getkappasmsh(grad, msstructure, dist=dist)
    hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
  } else {
    gradstats <- getkappas(grad, dist=dist)
    hseq <- gethseqfullse3(kstar,gradstats,kappa0,vext=vext)
    hseq$h <- cbind(rep(1,ngrad),hseq$h)
  }
  nind <- as.integer(hseq$n*1.25)#just to avoid n being to small due to rounding
  hseq <- hseq$h
  # make it nonrestrictive for the first step
  ni <- rep(1,nvoxel)
  minlevel <- if(model==1) 2*ncoils else sqrt(2)*gamma(ncoils+.5)/gamma(ncoils)
  minlevel0 <- if(model==1) 2*ns0*ncoils else sqrt(2)*gamma(ns0*ncoils+.5)/gamma(ns0*ncoils)
  z <- list(th=array(minlevel,c(nvoxel,nbv)), ni = ni)
  th0 <- rep(minlevel0,nvoxel)
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
  kinit <- if(lambda<1e10) 0 else kstar
  mc.cores <- setCores(,reprt=FALSE)
  for(k in kinit:kstar){
    gc()
    hakt <- hseq[,k+1]
    if(multishell){
      thmsh <- interpolatesphere(z$th,msstructure)
      param <- lkfullse3msh(hakt,kappa0/hakt,gradstats,vext,nind)
        z <- .Fortran(C_adsmse3m,
                      as.double(sb),#y
                      as.double(thmsh),#th
                      ni=as.double(z$ni),#ni
                      as.double(fncchiv(thmsh,varstats)),#sthi
                      as.integer(position),#mask
                      as.integer(nvoxel),
                      as.integer(nshell),# number of shells
                      as.integer(ddim[1]),#n1
                      as.integer(ddim[2]),#n2
                      as.integer(ddim[3]),#n3
                      as.integer(ngrad),#ngrad
                      as.double(lambda),#lambda
                      as.integer(mc.cores),#ncores
                      as.integer(param$ind),#ind
                      as.double(param$w),#w
                      as.integer(param$nind),#n
                      th=double(nvoxel*ngrad),#thn
                      double(ngrad*mc.cores),#sw
                      double(ngrad*mc.cores),#swy
                      double(nshell*mc.cores),#si
                      double(nshell*mc.cores))[c("ni","th")]
        dim(z$th) <- dim(z$ni) <- c(nvoxel,ngrad)
        gc()
    } else {
      param <- lkfullse3(hakt,kappa0/hakt,gradstats,vext,nind)
        z <- .Fortran(C_adsmse3p,
                      as.double(sb),
                      as.double(z$th),
                      ni=as.double(z$ni/if(model>=2) fncchiv(z$th,varstats) else 1),
                      as.integer(position),#mask
                      as.integer(nvoxel),
                      as.integer(ddim[1]),
                      as.integer(ddim[2]),
                      as.integer(ddim[3]),
                      as.integer(ngrad),
                      as.double(lambda),
                      as.integer(ncoils),
                      as.integer(mc.cores),
                      as.integer(param$ind),
                      as.double(param$w),
                      as.integer(param$n),
                      th=double(nvoxel*ngrad),
                      double(nvoxel*ngrad),#ldf (to precompute lgamma)
                      double(ngrad*mc.cores),
                      double(ngrad*mc.cores),
                      as.integer(model))[c("ni","th")]
        gc()
    }
    if(verbose){
      dim(z$ni)  <- c(prod(ddim),ngrad)
      cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni),3),
          "mean of ni",signif(mean(z$ni),3),
          " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
    } else {
      cat(".")
    }
    param <- reduceparam(param)
    z0 <- .Fortran(C_asmse30p,
                   as.double(s0),
                   as.double(th0),
                   ni=as.double(ni0/if(model==2) fncchiv(th0,varstats) else 1),
                   as.integer(position),#mask
                   as.integer(nvoxel),
                   as.integer(ddim[1]),
                   as.integer(ddim[2]),
                   as.integer(ddim[3]),
                   as.double(lambda),
                   as.integer(ncoils*ns0),
                   as.integer(param$ind),
                   as.double(param$w),
                   as.integer(param$n),
                   as.integer(param$starts),
                   as.integer(param$nstarts),
                   th0=double(nvoxel),
                   double(nvoxel),#ldf (to precompute lgamma)
                   double(param$nstarts),#swi
                   as.integer(model))[c("ni","th0")]
    th0 <- z0$th0
    ni0 <- z0$ni
    rm(z0)
    gc()
    if(verbose){
      cat("End smoothing s0: quartiles of ni",signif(quantile(ni0[mask]),3),
          "mean of ni",signif(mean(ni0[mask]),3),
          " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
    }
  }
  list(th=z$th*sigma, th0=th0*sigma, ni=z$ni, ni0=ni0, hseq=hseq, kappa0=kappa0, lambda=lambda)
}

smse3ms <- function(sb, s0, bv, grad, kstar, lambda, kappa0,
                  mask, sigma, ns0=1, ws0=1,
                  vext = NULL,
                  ncoils = 1,
                  verbose = FALSE,
                  usemaxni = TRUE){
#
#   determine structure of the space
#
    ngrad <- length(bv)
    ddim <- dim(mask)
    msstructure <- getnext3g(grad, bv)
    bv <- msstructure$bv
    nshell <- as.integer(msstructure$nbv)
    ubv <- msstructure$ubv
    nvoxel <- sum(mask)
    nbv <- length(bv)
# check dimensions
    if(dim(sb)[1]!=nvoxel || length(s0)!=nvoxel || dim(sb)[2]!=nbv){
       stop("smse3ms - sb and s0 should only contain data within mask")
    }
#
#  rescale so that we have Chi-distributed values
#
    dsigma <- dim(sigma)
    if(is.null(dsigma)){
       s0 <- s0/sigma
       sb <- sb/sigma
    } else if(dsigma[2]==nshell+1){
      s0 <- s0/sigma[,1]
      for (shnr in 1:nshell) {
        indbv <- (1:length(bv))[bv == ubv[shnr]]
        sb[, indbv] <- sb[, indbv]/ sigma[, shnr+1]
      }
    } else if(dsigma[2]==nbv+1){
      s0 <- s0/sigma[,1]
      for (ibv in 1:nbv) {
        sb[, ibv] <- sb[, ibv]/ sigma[, ibv+1]
      }
    } else {
      warning("incompatible number of sigma images, needs to be 1, nshell+1 or nbv+1\n")
    }
# define position of voxel in mask within the 3D cube
    position <- array(0,ddim)
    position[mask] <- 1:nvoxel
# set relation between voxel lengths
    if(is.null(vext)) vext <- c(1,1)
# approximate the noise distribution
    varstats <- sofmchi(ncoils)
# only keep information within mask
    z <- list(th = array(1,c(nvoxel,nbv)), th0 = rep(1,nvoxel),
              ni = array(1,c(nvoxel,nbv)), ni0 = rep(1,nvoxel))
# generate sequence od bandwidths
    gradstats <- getkappasmsh3(grad, msstructure)
    hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
    nind <- as.integer(hseq$n*1.25)
    if(usemaxni){
      ni <- array(1,c(nvoxel,nbv))
      ni0 <- rep(1,nvoxel)
    }
    prt0 <- Sys.time()
    cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
    kinit <- if(lambda<1e10) 0 else kstar
    mc.cores <- setCores(,reprt=FALSE)
    gc()
    for(k in kinit:kstar){
      hakt <- hseq$h[,k+1]
      t0 <- Sys.time()
      thnimsh <- interpolatesphere1(z$th,z$th0,z$ni,z$ni0,msstructure)
      rm(z)
      gc()
      t1 <- Sys.time()
      param <- lkfullse3msh(hakt,kappa0/hakt,gradstats,vext,nind)
      hakt0 <- mean(hakt)
      param0 <- lkfulls0(hakt0,vext,nind)
# calls to base::findInterval
      vs2 <- varstats$s2[findInterval(thnimsh$mstheta, varstats$mu, all.inside = TRUE)]/2
      vs02 <- varstats$s2[findInterval(thnimsh$msth0, varstats$mu, all.inside = TRUE)]/2
      t2 <- Sys.time()
### need to rewrite adsmse3s to only use data in mask
      z <- .Fortran(C_adsmse3s,
                    as.double(sb),#y
                    as.double(s0),#y0
                    as.double(thnimsh$mstheta),#th
                    as.double(thnimsh$msni),#ni/si^2
                    as.double(thnimsh$msth0),#th0
                    as.double(thnimsh$msni0),#ni0/si^2
                    as.double(vs2),#var/2
                    as.double(vs02),#var/2 for s0
                    as.integer(position),
                    as.integer(nvoxel),
                    as.integer(nshell+1),#ns number of shells
                    as.integer(ddim[1]),#n1
                    as.integer(ddim[2]),#n2
                    as.integer(ddim[3]),#n3
                    as.integer(ngrad),#ngrad
                    as.double(lambda),#lambda
                    as.double(ws0),# wghts0 rel. weight for s0 image
                    as.integer(param$ind),#ind
                    as.double(param$w),#w
                    as.integer(param$n),#n
                    as.integer(param0$ind),#ind0
                    as.double(param0$w),#w0
                    as.integer(param0$n),#n0
                    th=double(nvoxel*nbv),#thn
                    ni=double(nvoxel*nbv),#nin
                    th0=double(nvoxel),#th0n
                    ni0=double(nvoxel),#ni0n
                    double(ngrad*mc.cores),#sw
                    double(ngrad*mc.cores),#swy
                    double((nshell+1)*mc.cores),#thi
                    double((nshell+1)*mc.cores),#nii
                    double((nshell+1)*mc.cores))[c("ni","th","ni0","th0")]
      t3 <- Sys.time()
      gc()
      rm(thnimsh,vs2,vs02)
      gc()
      if(usemaxni){
        ni <- z$ni <- if(usemaxni) pmax(ni,z$ni)
        ni0 <- z$ni0 <- if(usemaxni) pmax(ni0,z$ni0)
      }
      dim(z$ni) <- dim(z$th) <- c(nvoxel,nbv)
      if(verbose){
        cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni),3),
            "mean of ni",signif(mean(z$ni),3),
            "\n              quartiles of ni0",signif(quantile(z$ni0),3),
            "mean of ni0",signif(mean(z$ni0),3),
            " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
        cat("interpolation:",format(difftime(t1,t0),digits=3),
            "param:",format(difftime(t2,t1),digits=3),
            "smoothing:",format(difftime(t3,t2),digits=3),"\n")
      } else {
        cat(".")
      }
    }
    #
    #   now rescale results with sigma
    #
    #  s0 was scaled differently
          if(is.null(dsigma)){
            z$th0 <- z$th0*sigma/sqrt(ns0)
            z$th <- z$th*sigma
          } else if(dsigma[2]==nshell+1){
            for (shnr in 1:nshell) {
              z$th0 <- z$th0*sigma[,1]/sqrt(ns0)
              indbv <- (1:length(bv))[bv == ubv[shnr]]
              z$th[, indbv] <- z$th[, indbv]* sigma[, shnr+1]
            }
          }
          else if(dsigma[2]==nbv+1){
            for (ibv in 1:nbv) {
              z$th0 <- z$th0*sigma[,1]/sqrt(ns0)
              z$th[, ibv] <- z$th[, ibv]* sigma[, ibv+1]
            }
          }
    list(th=z$th, th0=z$th0, ni=z$ni, ni0=z$ni0, hseq=hseq, kappa0=kappa0, lambda=lambda)
  }
