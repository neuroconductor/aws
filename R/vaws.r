vaws <- function(y,kstar=16,
                 sigma2=1,mask=NULL,scorr=0,spmin=0.25,
                 ladjust=1,wghts=NULL,u=NULL,
                 maxni=FALSE){
   args <- match.call()
   dy<-dim(y)
   nvec <- dy[1]
   dy <- dy[-1]
   d <- length(dy)
   if(length(dy)>3) stop("Vector AWS for more than 3 dimensional grids is not implemented")
   lambda <- 2*sigma2*ladjust*switch(d,qchisq(pchisq(14.6,1),nvec),## 1D
                      qchisq(pchisq(9.72,1),nvec),## 2D
                      qchisq(pchisq(8.82,1),nvec))## 3D
   if(is.null(wghts)) wghts <- c(1,1,1)
   wghts <- switch(length(dy),c(0,0),c(wghts[1]/wghts[2],0),wghts[1]/wghts[2:3])
   n1 <- switch(d,n,dy[1],dy[1])
   n2 <- switch(d,1,dy[2],dy[2])
   n3 <- switch(d,1,1,dy[3])
   n <- n1*n2*n3
   if(is.null(mask)) mask <- rep(TRUE,n)
   h0 <- 0
   if(any(scorr>0)) {
         h0<-numeric(length(scorr))
         for(i in 1:length(h0))
         h0[i]<-geth.gauss(scorr[i])
         if(length(h0)<d) h0<-rep(h0[1],d)
         cat("Corresponding bandwiths for specified correlation:",h0,"\n")
   }
   zobj<-list(bi= rep(1,n), theta= y)
   bi <- zobj$bi
   cat("Progress:")
   total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
   mc.cores <- setCores(,reprt=FALSE)
   k <- 1
   hmax <- 1.25^(kstar/d)
   lambda0 <- lambda
   mae <- NULL
   while (k<=kstar) {
      hakt0 <- gethani(1,1.25*hmax,2,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,1.25*hmax,2,1.25^k,wghts,1e-4)
      cat("step",k,"hakt",hakt,"time",format(Sys.time()),"\n")
      dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
      if(scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
      zobj <- .Fortran("vaws",as.double(y),
                       as.logical(mask),
                       as.integer(nvec),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi=as.double(zobj$bi),
                       theta=double(nvec*n),
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec*mc.cores),
                       PACKAGE="aws")[c("bi","theta","hakt")]
      dim(zobj$theta)<-c(nvec,dy)
      if(maxni) bi <- zobj$bi <- pmax(bi,zobj$bi)
      dim(zobj$bi)<-dy
      if(!is.null(u)) {
      cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    signif(mean((zobj$theta-u)^2),3),"   MAE: ",
		    signif(mean(abs(zobj$theta-u)),3)," mean(bi)=",
		    signif(mean(zobj$bi),3),"\n")
                    mae<-c(mae,signif(mean(abs(zobj$theta-u)),3))
		    }
      x<-1.25^k
      scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
      lambda0<-lambda*scorrfactor
      if (max(total) >0) {
          cat(signif(total[k],2)*100,"% . ",sep="")
       }
      k <- k+1
      gc()
   }
   cat("\n")
   list(y=y,theta=zobj$theta,hakt=hakt,sigma2=sigma2,lambda=lambda,
        ladjust=ladjust,args=args,wghts=wghts,mae=mae,ni=zobj$bi)
}
          
vsegm <- function(theta,bi,mask,level){
   dth <- dim(theta)
   dy <- dth[-1]
   nv <- dth[1]
   n1 <- dy[1]
   n2 <- dy[2]
   n3 <- dy[3]
   n <- n1*n2*n3
   segm <- .Fortran("vsegmen0",
                 as.double(theta),
                 as.double(bi),
                 as.logical(mask),
                 as.double(level),
                 double(prod(2*dy+1)),
                 as.integer(nv),
                 as.integer(2*n1+1),
                 as.integer(2*n2+1),
                 as.integer(2*n3+1),
                 integer(n),
                 segm=integer(n),
                 PACKAGE="aws")$segm
   dim(segm) <- dy
   segm          
}

fillsegm <- function(segm,mask,slevel){
##
##  replace segment number>slevel with a segment number 
##  from an adjacent region (in mask and segm<=slevel)
##
   ds <- dim(segm)
   n1 <- ds[1]
   n2 <- ds[2]
   n3 <- ds[3]
   n <- n1*n2*n3
   ntbl <- sum(segm[mask]>slevel)
   z <- .Fortran("fillsegm",
                 nseg=as.integer(segm),
                 as.integer(n1),
                 as.integer(n2),
                 as.integer(n3),
                 as.logical(mask),
                 integer(3*ntbl),
                 as.integer(ntbl),
                 as.integer(slevel),
                 PACKAGE="aws")$nseg
    array(z,ds)
}

vdetrend <- function(theta,mask,h){
##
##  remove spatial trends using a median filter
##  trend removal will be everywhere as long as there are
##  active (mask==TRUE) voxel within a ball of radius h 
##
   dth <- dim(theta)
   dy <- dth[-1]
   nv <- dth[1]
   n1 <- dy[1]
   n2 <- dy[2]
   n3 <- dy[3]
   n <- n1*n2*n3
   newtheta <- theta
   mc.cores <- setCores(,reprt=FALSE)
   nwmd <- (2*as.integer(h)+1)^3
   parammd <- .Fortran("paramw3",
                      as.double(h),
                      as.double(c(1,1)),
                      ind=integer(3*nwmd),
                      w=double(nwmd),
                      n=as.integer(nwmd),
                      PACKAGE = "aws")[c("ind","n")]
   nwmd <- parammd$n
   parammd$ind <- parammd$ind[1:(3*nwmd)]
   dim(parammd$ind) <- c(3,nwmd)
   nind <- (2*h+1)^3
   for(k in 1:nv){
      z <- .Fortran("medsm1",
                    as.double(theta[k,,,]),
                    as.logical(mask),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.integer(parammd$ind),
                    as.integer(nwmd),
                    double(2*nwmd*mc.cores),
                    as.integer(mc.cores),
                    thnew=double(n),
                    PACKAGE="aws")$thnew
      newtheta[k,,,] <- z
   }
   theta-newtheta
}
