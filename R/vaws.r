vaws <- function(y,kstar=16,homogen=TRUE,
                 sigma2=1,scorr=0,spmin=0.25,
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
   h0 <- 0
   if(any(scorr>0)) {
         h0<-numeric(length(scorr))
         for(i in 1:length(h0))
         h0[i]<-geth.gauss(scorr[i])
         if(length(h0)<d) h0<-rep(h0[1],d)
         cat("Corresponding bandwiths for specified correlation:",h0,"\n")
   }
   zobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y, bi0= rep(1,n))
   bi <- zobj$bi
   hhom <- rep(1,n)
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
      cat("step",k,"hakt",hakt,"\n")
      dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
      if(scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
      zobj <- .Fortran("vaws",as.double(y),
                       as.integer(nvec),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi=as.double(zobj$bi),
                       bi2=double(n),
                       bi0=double(n),
                       theta=double(nvec*n),
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec*mc.cores),
                       PACKAGE="aws")[c("bi","bi0","bi2","theta","hakt","hhom")]
      dim(zobj$theta)<-c(nvec,dy)
      if(maxni) bi <- zobj$bi <- pmax(bi,zobj$bi)
      dim(zobj$bi)<-dy
      if(homogen) hhom <- zobj$hhom
      if(!is.null(u)) {
      cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    signif(mean((zobj$theta-u)^2),3),"   MAE: ",
		    signif(mean(abs(zobj$theta-u)),3)," mean(bi)=",
		    signif(mean(zobj$bi),3),"mean hhom",signif(mean(hhom),3),"\n")
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
                 