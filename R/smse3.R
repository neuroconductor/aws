smse3 <- function(sb, s0, bv, grad, sigma, kstar, lambda, kappa0,
                  mask,
                  vext = NULL,
                  ncoils = 1,
                  ws0 = 1,
                  level = NULL,
                  verbose = FALSE,
                  usemaxni = TRUE){ ### should arguments only be voxel in mask ??
#
#   determine structure of the space
#
    msstructure <- getnext3g(grad, bvalues)
    bv <- msstructure$bv
    nshell <- as.integer(msstructure$nbv)
    ubv <- msstructure$ubv
    ddim0 <- dim(mask)
    nvoxel <- sum(mask)
    nbv <- length(bv)
# check dimensions
    if(dim(sb)[1]!=nvoxel || length(s0)!=nvoxel || dim(sb)[2]!=nbv){
       stop("aws::smse3 - sb and s0 should only contain data within mask")
    }
# define position of voxel in mask within the 3D cube
    position <- array(0,ddim0)
    position[mask] <- 1:nvoxel
# set relation between voxel lengths
    if(is.null(vext)) vext <- c(1,1)
# approximate the noise distribution
    varstats <- sofmchi(ncoils)
# only keep information within mask
    z <- list(th = array(1,c(nvoxel,nbv)), th0 = rep(1,nvoxel),
              ni = array(1,c(nvoxel,nbv)), ni0 = rep(1,nvoxel))
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
      thnimsh <- interpolatesphere1(z$th,z$th0,z$ni,z$ni0,msstructure,mask)
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
                    as.integer(nshell+1),#ns number of shells
                    as.integer(ddim0[1]),#n1
                    as.integer(ddim0[2]),#n2
                    as.integer(ddim0[3]),#n3
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
      if(memrelease) rm(thnimsh,vs2,vs02)
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
    ngrad <- ngrad+1
    #
    #  one s0 image only
    #
    si <- array(object@si,c(ddim,ngrad))
    #
    #  back to original scale
    #
    if (length(sigma) > 1) {
      if(length(dim(sigma)) == 3) {
        si[xind, yind, zind, 1] <-  z$th0/sqrt(ns0)*sigma
      } else {
        si[xind, yind, zind, 1] <-  z$th0/sqrt(ns0)*sigma[, , , 1]
      }
    } else {
      si[xind, yind, zind, 1] <-  z$th0/sqrt(ns0)*sigma
    }
    #  go back to original s0 scale
    #  for the DWI we need to scale back differently
    if (sigmacase == 1) {
      si[xind, yind, zind, -1] <- z$th*sigma
    } else if (sigmacase == 2) {
      si[xind, yind, zind, -1] <- sweep(z$th, 1:3, sigma, "*")
    } else if (sigmacase == 3) { # now sigma is an array with (identical(dim(sigma), ddim) == TRUE)
      for (shnr in 1:nshell) {
        indth <- (1:length(bv))[bv == ubv[shnr]]
        indsi <- indth+1
        si[xind, yind, zind, indsi] <- sweep(z$th[, , , indth], 1:3, sigma[, , , shnr+1], "*")
        }
    } else {
        si[xind, yind, zind, -1] <- z$th*sigma[,,,-1]
    }

                  }
