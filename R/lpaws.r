#
#    R - functions  for  Adaptive Weights Smoothing (AWS)
#    in (generalized) local polynomial regression models in 1D and 2D         
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
##############################################################################
#
#   Local polynomal AWS (Gaussian case on a grid) max. polynomial degree 2  
#
##############################################################################
lpaws <- function(y,degree=1,hmax=NULL,qlambda=NULL,qtau=NULL,lkern="Triangle",
                  homogen=TRUE,earlystop=TRUE,
                  aggkern="Uniform",sigma2=NULL,hinit=NULL,hw=NULL,
                  lseq=NULL,u=NULL,graph=FALSE,demo=FALSE,spmin=.25)
{ 
#
#          Auxilary functions
#
Pardist <- function(d,Bi0,dtheta){
#  local polynomial uni  mcode=1
#  local polynomial bi   mcode=2
   dp1 <- switch(d,dim(dtheta)[2],dim(dtheta)[3])
   dp2 <- switch(d,dim(Bi0)[2],dim(Bi0)[3])
   if(d==1){
      dist <- 0
      for(i in 1:dp1) for(j in 1:dp1) dist <- dist+dtheta[,i]*Bi0[,i+j-1]*dtheta[,j]
   }
   if(d==2){
         ind <- matrix(c(1, 2, 3, 4, 5, 6,
                         2, 4, 5, 7, 8, 9,
                         3, 5, 6, 8, 9,10,
                         4, 7, 8,11,12,13,
                         5, 8, 9,12,13,14,
                         6, 9,10,13,14,15),6,6)[1:dp1,1:dp1,drop=FALSE]
                dist <- 0
                for(i in 1:dp1) for(j in 1:dp1) dist <- dist+dtheta[,,i]*Bi0[,,ind[i,j]]*dtheta[,,j]
   }
   dist
   }
#
#     Compute theta
#
gettheta <- function(d,ai,bi){
if(d==1){
   n <- dim(ai)[1]
   dp1 <- dim(ai)[2]
   dp2 <- dim(bi)[2]
   ind <- matrix(c(1, 2, 3,
                   2, 3, 4,
                   3, 4, 5),3,3)[1:dp1,1:dp1]
}  else {
   n1 <- dim(ai)[1]
   n2 <- dim(ai)[2]
   n <- n1*n2
   dp1 <- dim(ai)[3]
   dp2 <- dim(bi)[3]
   ind <- matrix(c(1, 2, 3, 4, 5, 6,
                   2, 4, 5, 7, 8, 9,
                   3, 5, 6, 8, 9,10,
                   4, 7, 8,11,12,13,
                   5, 8, 9,12,13,14,
                   6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
} 
   theta <- .Fortran("mpaws2",
                  as.integer(n),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=double(dp1*n),
                  double(dp1*dp1),
                  as.integer(ind),PACKAGE="aws")$theta
   dim(theta) <- switch(d,c(n,dp1),c(n1,n2,dp1))
   theta
}
#
#  update theta (memory step)
#
updtheta <- function(d,zobj,fix,cpar,aggkern,bikm2,bi2km2,theta){
    heta <- cpar$heta
    tau1 <- cpar$tau1
    tau2 <- cpar$tau2
    kstar <- cpar$kstar
    hakt <- zobj$hakt
    tau <- 2*(tau1+tau2*max(kstar-log(hakt),0))
    bi <- zobj$bi
    bi2 <- zobj$bi2
    thetanew <- gettheta(d,zobj$ai,bi)
    n<-length(fix)
    if(any(fix)) thetanew[array(fix,dim(thetanew))] <- theta[rep(fix,dp1)]
    if(max(bikm2)<1) heta <- 1e40
#  we don't have initializations for bikm2 and theta 
    if (hakt>heta) {
	eta <- rep(pmin(1,Pardist(d,bikm2,thetanew-theta)/tau),dp2)
    } else {
        eta <- rep(0,n*dp2)
    }
    if(any(fix)) eta[rep(fix,dp2)] <- 1
    if(any(eta>0)){
       bi <- (1-eta)*bi + eta * bikm2
       bi2 <- (1-eta)*bi2 + eta * bi2km2
       eta <- array(eta[1:n],dim(theta))
       theta <- (1-eta)*thetanew + eta * theta
       eta <- eta[1:n]
       if(d==2) dim(eta) <- dim(theta)[1:2]
    } else {
      theta <- thetanew
    }
    eta <- eta[1:n]
    if(d==2) dim(eta) <- dim(theta)[1:2]
    list(theta=theta,bi=bi,bi2=bi2,eta=eta,fix=as.vector(eta==1))
  }
#
#          Main function body
#
#    first check arguments and initialize                                 
#
args <- match.call()
mae <- NULL
#
#     set approriate defaults
#
dy <- dim(y)
if(is.null(dy)) d<-1 else d<-length(dy) 
if(d>2) return(warning("Local polynomial PS is only implemented for 1 and 2 dimensional grids"))
if(!(degree %in% 0:2)) return(warning("Local polynomial PS is only implemented for degrees up to 2"))
if(d==1){
dp1 <-  degree+1
dp2 <- degree+dp1
n <- length(y)
} else {
n1 <- dy[1]
n2 <- dy[2]
n <- n1*n2
dp1 <- switch( degree+1,1,3,6)
dp2 <- switch( degree+1,1,6,15)
}
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,
	            Gaussian=5,2)
if(is.null(qlambda)) qlambda <-switch(d,
                                     switch( degree+1,.98,.92,.92),
                                     switch( degree+1,.96,.75,.92))
if(earlystop) nfix <- switch(d,
                             switch(degree+1,2,5,20),
                             switch(degree+1,2,3,4)) else nfix <- n
if(is.null(spmin)) spmin <- switch(degree+1,.3,.25,.25)
if(spmin<0|spmin>1) spmin <- switch(degree+1,.3,.25,.25)
  if (qlambda>=1) {
    # thats stagewise aggregation with kernel specified by aggkern
    if(is.null(qtau)) qtau <- switch(d,
                                     switch( degree+1,.6,.8,.94),
                                     switch( degree+1,.55,.45,.45))
    if (qtau==1) {
        tau1 <- 1e50 
	hinit <- heta <- hmax
    } else {
	tau1 <- qchisq(qtau,degree+1)
	hinit <- heta <- degree+1.1
    }
    if (aggkern=="Triangle") tau1 <- 2.5*tau1
    tau2 <- tau1/2
  } else {
    if (is.null(qtau)) qtau<-switch(degree+1,1,.68,.86)
    if (qtau>=1) {
      tau1 <- 1e50 
      heta <- 1e40
    } else {
      tau1 <- qchisq(qtau,degree+1)^(2/d)
      heta <- 1.25^((degree+2)/d)
    }
    if (aggkern=="Triangle") tau1 <- 2.5*tau1
    tau2 <- tau1/2
  }
  cat("tau1",tau1,"\n")
  if (qlambda<1) lambda <- qchisq(qlambda,dp1) else lambda <- 1e50
   lambda<-lambda*2
cat("Value of lambda",lambda,"\n")
if(is.null(hinit)) hinit <- 1 
hincr <- 1.25^(1/d)
if (is.null(hmax)) hmax <- switch(d,100,12)
kstar <- switch(d,log(250*dp1),log(15*dp1))
cpar <- list(heta=heta,tau1=tau1,tau2=tau2,dy=dy,kstar=kstar)
    if(is.null(sigma2)) {
       sigma2 <- IQRdiff(as.vector(y))^2
       cat("Estimated error variance",signif(sigma2,3),"\n")
       }
    if (length(sigma2)==1) {
      #   homoskedastic Gaussian case
      lambda <- lambda*sigma2*2 
      cpar$tau1 <- cpar$tau1*sigma2*2 
      cpar$tau2 <- cpar$tau2*sigma2*2 
    } else if (length(sigma2)!=n) {
      cpar$tau1 <- cpar$tau1*sigma2*2 
      cpar$tau2 <- cpar$tau2*sigma2*2 
      lambda <- lambda*2 
    } else {
      #   heteroskedastic Gaussian case
      if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as img")
      lambda <- lambda*2 
      cpar$tau1 <- cpar$tau1*2 
      cpar$tau2 <- cpar$tau2*2 
      sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
    }
#  tobj <- list(bi= rep(1,n*dp2), bi2= rep(1,n*dp2), theta= rep(0,n*dp1), fix=rep(FALSE,n))
  bi0old <- rep(1,n*dp2)
  bi <- array(rep(1,n*dp2),c(switch(d,n,dy),dp2))
  bi2 <- array(rep(0,n*dp2),c(switch(d,n,dy),dp2))
  bikm1 <- array(rep(1,n*dp2),c(switch(d,n,dy),dp2))
  bi2km1 <- array(rep(0,n*dp2),c(switch(d,n,dy),dp2))
  theta <- thetakm1 <- array(rep(0,n*dp1),c(switch(d,n,dy),dp1))
  hhom <- rep(1,n)
  fix <- rep(FALSE,n)
  ind <- switch(d,
                matrix(c(1, 2, 3,
                         2, 3, 4, 
                         3, 4, 5),3,3)[1:dp1,1:dp1],
                matrix(c(1, 2, 3, 4, 5, 6,
                         2, 4, 5, 7, 8, 9,
                         3, 5, 6, 8, 9,10,
                         4, 7, 8,11,12,13,
                         5, 8, 9,12,13,14,
                         6, 9,10,13,14,15),6,6)[1:dp1,1:dp1])
  if(is.null(hw)) hw<-switch(d,degree+1.1,degree+.1) else hw<-max(hw,degree+.1)
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)
  if (is.null(lseq)) lseq <- 1
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  k <- 1
  hakt0 <- hakt <- hinit
  lambda0 <- 1e50
  #
  #   run single steps to display intermediate results
  #
  while (hakt<=hmax) {
    twohp1<-2*trunc(hakt)+1
    twohhwp1<-2*trunc(hakt+hw)+1
    if (length(sigma2)==n) {
      # heteroskedastic Gaussian case
      hakt0 <- hakt
      zobj <- switch(d,
                     .Fortran("awsph1",
		       as.double(y),
                       as.double(sigma2),
                       as.logical(fix),
                       as.integer(nfix),
                       as.integer(n),
                       as.integer(degree),
		       as.double(hw),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(theta),
                       bi=as.double(bi),
                       bi2=double(n*dp2),
                       bi0=double(n*dp2),
                       ai=double(n*dp1),
                       as.integer(lkern),
                       as.double(spmin),
                       double(twohp1),# array for location weights
                       double(twohp1),# array for general weights
                       double(twohhwp1),# array for smoothed location weights
                       double(twohhwp1),# array for smoothed general weights
                       as.integer(ind),
                       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt")],
                     .Fortran("awsph2",
		       as.double(y),
                       as.double(sigma2),
                       as.logical(fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(degree),
		       as.double(hw),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(theta),
                       bi=as.double(bi),
                       bi2=double(n*dp2),
                       bi0=double(n*dp2),
                       ai=double(n*dp1),
                       as.integer(lkern),
                       as.double(spmin),
                       double(twohp1*twohp1),# array for location weights
                       double(twohp1*twohp1),# array for general weights
                       double(twohhwp1*twohhwp1),# array for smoothed location weights
                       double(twohhwp1*twohhwp1),# array for smoothed general weights
                       as.integer(ind),
                       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt","hhom","fix")])
    } else {
      # all other cases
      zobj <- switch(d,
                     .Fortran("awsp1",
		       as.double(y),
                       fix=as.logical(fix),
                       as.integer(nfix),
                       as.integer(n),
                       as.integer(degree),
		       as.double(hw),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(theta),
                       bi=as.double(bi),
                       bi2=double(n*dp2),
                       bi0=double(n*dp2),
                       ai=double(n*dp1),
                       as.integer(lkern),
                       as.double(spmin),
                       double(twohp1),# array for location weights
                       double(twohp1),# array for general weights
                       double(twohhwp1),# array for smoothed location weights
                       double(twohhwp1),# array for smoothed general weights
                       as.integer(ind),
                       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt","hhom","fix")],
                     .Fortran("awsp2",
		       as.double(y),
                       as.logical(fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(degree),
		       as.double(hw),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(theta),
                       bi=as.double(bi),
                       bi2=double(n*dp2),
                       bi0=double(n*dp2),
                       ai=double(n*dp1),
                       as.integer(lkern),
                       as.double(spmin),
                       double(twohp1*twohp1),# array for location weights
                       double(twohp1*twohp1),# array for general weights
                       double(twohhwp1*twohhwp1),# array for smoothed location weights
                       double(twohhwp1*twohhwp1),# array for smoothed general weights
                       as.integer(ind),
                       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt")])
    }
    gc()
    dim(zobj$ai) <- c(switch(d,n,dy),dp1)
    if (hakt>n^(1/d)/2) zobj$bi0 <- hincr^d*biold
    biold <- zobj$bi0
    dim(zobj$bi0)<-c(switch(d,n,dy),dp2)
    if(!homogen) zobj$hhom<-rep(1,n)
    tobj <- updtheta(d,zobj,fix,cpar,aggkern,bikm1,bi2km1,thetakm1)
    if(!is.null(zobj$hhom)) hhom <- zobj$hhom else hhom<-rep(1,n)
    fix <- tobj$fix
    fix[zobj$fix] <- TRUE
    rm(zobj)
    gc()
    dim(tobj$theta) <- c(switch(d,n,dy),dp1)
    dim(tobj$bi) <- c(switch(d,n,dy),dp2)
    dim(tobj$bi2) <- c(switch(d,n,dy),dp2)
    dim(tobj$eta) <- switch(d,NULL,dy)
       dim(bikm1)<-dim(bi2km1)<-dim(bi)<-dim(bi2) <- c(n,dp2)
       dim(thetakm1)<-dim(theta) <- c(n,dp1)
       bikm1[!fix,] <- bi[!fix,]
       bi2km1[!fix,] <- bi2[!fix,]
       thetakm1[!fix,] <- theta[!fix,]
       dim(bikm1)<-dim(bi2km1)<-dim(bi)<-dim(bi2) <- c(switch(d,n,dy),dp2)
       dim(thetakm1)<-dim(theta) <- c(switch(d,n,dy),dp1)
    bi <- tobj$bi
    bi2 <- tobj$bi2
    theta <- tobj$theta
    eta <- tobj$eta
    rm(tobj)
    gc()
    if (graph) {
      if(d==1){
      oldpar<-par(mfrow=c(1,2),mar=c(3,3,3,.25),mgp=c(2,1,0))
      plot(y)
      lines(theta[,1],col=2)
      title("Observed data and estimate")
      plot(bi[,1],type="l",ylim=c(0,max(bi[,1])))
      lines(eta*max(bi[,1]),col=2)
      lines(hhom/max(hhom)*max(bi[,1]),col=3)
      title(paste("hakt=",signif(hakt,3),"bi,eta and hhom"))
      } else {
      oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
      image(y,xaxt="n",yaxt="n",col=gray((0:255)/255))
      title("Observed Image")
      image(theta[,,1],xaxt="n",yaxt="n",col=gray((0:255)/255))
      title(paste("Reconstruction  h=",signif(hakt,3)," Range ",signif(min(theta[,,1]),3),"-",signif(max(theta[,,1]),3)))
      image(bi[,,1],xaxt="n",yaxt="n",col=gray((0:255)/255))
      title(paste("Sum of weights: min=",signif(min(bi[,,1]),3)," mean=",signif(mean(bi[,,1]),3)," max=",signif(max(bi[,,1]),3)))
      image(eta,xaxt="n",yaxt="n",col=gray((0:255)/255))
      title(paste("eta   max=",signif(max(eta),3)))
    }
    par(oldpar)
    }
    if (!is.null(u)) {
      th <- switch(d,theta[,1],theta[,,1])
       cat("bandwidth: ",signif(hakt,3),"eta==1",sum(eta==1),"   MSE: ",
          signif(mean((th-u)^2),3),"   MAE: ",signif(mean(abs(th-u)),3),"mean hhom",signif(mean(hhom),3),"\n")
      mae<-c(mae,signif(mean(abs(th-u)),3))
    }
    if (demo) readline("Press return")
    hakt <- hakt*hincr
    lambda0<-lambda*lseq[k]
    k<-k+1
    gc()
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 
  ###   component var contains an estimate of Var(theta) if and aggkern="Uniform", or if qtau1=1 
  ###   
  if (length(sigma2)==n) {
    # heteroskedastic Gaussian case 
    vartheta <- bi2/bi^2
  } else {
    vartheta <- sigma2[1]*bi2/bi^2
  }
  z <- list(theta=theta,
            ni=switch(d,bi[,1],bi[,,1]),
	    sigma2=sigma2,
            var=vartheta,
            hmax=hakt/hincr,
            mae=mae,
            lseq=c(0,lseq[-steps]),
            call=args)
  class(z) <- "awslp.gaussian"
  z
}




