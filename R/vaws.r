#
#    R - function  vaws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for vector valued local constant Gaussian models                              #
# 
#
#    Copyright (C) 2002 Weierstrass-Institut fuer
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
#     default parameters:  
#       
#                        heta=2   qtau=.95   qlambda=.96
#
#
vaws <- function(y,qlambda=NULL,qtau=NULL,lkern="Triangle",aggkern="Uniform",sigma2=NULL,
                 hinit=NULL,hincr=NULL,hmax=NULL,heta=NULL,
                 eta0=0,u=NULL,graph=FALSE,wghts=NULL,vwghts=NULL)
{
#
#
IQRdiff <- function(y) IQR(diff(y))/1.908
#
#
updtheta<-function(zobj,tobj,cpar){
heta<-cpar$heta
eta0<-cpar$eta0
tau1<-cpar$tau1
tau2<-cpar$tau2
kstar<-cpar$kstar
vw<-cpar$vw
dd<-cpar$dd
hakt<-zobj$hakt
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
bi0<-zobj$bi0
bi<-zobj$bi
bi2<-zobj$bi2
n<-cpar$n
thetanew<-zobj$ai/outer(rep(1,dd),bi)
theta<-tobj$theta
if(hakt>heta) {
eta<-vw%*%matrix((thetanew-theta)^2,dd,n)/tau
dim(eta)<-dim(bi0)
eta<-bi0*eta
eta<-(1-eta0)*(1-pmax(0,1-eta))+eta0
} else {
eta <- rep(eta0,length(bi))
}
dim(eta)<-dim(bi)
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
bi2 <- (1-eta)*bi2 + eta * tobj$bi2
etadd<-outer(rep(1,dd),eta)
theta <- (1-etadd)*thetanew + etadd * theta
list(theta=theta,bi=bi,bi2=bi2,eta=eta,fix=(eta==1))
}
#
#    first check arguments and initialize
#
args <- match.call()
spmax <- 5
dy<-dim(y)
if(is.null(dy)) return("dim(y) should exist for multivariate data")
dd<-dy[1]
if(is.null(heta)) {
if(length(dy)==2) heta<-2
if(length(dy)==3) heta<-2
if(length(dy)==4) heta<-2
}
if(is.null(qtau)) {
if(length(dy)==2) qtau<-.95 
if(length(dy)==3) qtau<-.95
if(length(dy)==4) qtau<-.95 
}
if(qtau<1) tau1<-2*qchisq(qtau,1) else {
tau1<-1e50
heta<- 1e50
}  
if(aggkern=="Triangle") tau1<-2.5*tau1
tau2<-tau1/2
if(is.null(hmax)){
if(length(dy)==2) hmax<-250    # uses a maximum of about 500 points
if(length(dy)==3) hmax<-12   # uses a maximum of about 450 points
if(length(dy)==4) hmax<-5    # uses a maximum of about 520 points
}
cpar<-list(heta=heta,tau1=tau1,tau2=tau2,eta0=eta0,vw=vwghts,dd=dd)
#
#          if not provided set default value for qlambda 
#
if(is.null(qlambda)) qlambda <- .98
# qlambda>=1  eliminates the use of the stochastic penalty
#
#  check if family is implemented and set the code for the family (used in kldist) 
#
if(is.null(cpar$eta0)) cpar$eta0<-0
#
#    set the code for the kernel (used in lkern) and set lambda
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
#
#      get lambda as quantile of appropriate chisq,
#                rescale to be consistent with the paper in  lambda
#
if(is.null(vwghts)) vwghts<-rep(1,dd) 
if(qlambda<1) lambda <- 2*qchisq(qlambda,dd)*sum(vwghts)/dd else lambda <- 1e50
#
#   the factor sum(vwghts)/dd  is just a guess how to correct lambda if vwghts are not equal 
#   resulting value may be to small ... .
#
vwghts<-1/vwghts 
# vwghts originally contains a componentwise factor for the variance, we need the inverse of this as a weight
#
#     now set hinit and hincr if not provided
#
if(is.null(hinit)||hinit<1) hinit <- 1
if(is.null(hincr)||hincr<=1) hincr <-1.25
#
#   estimate variance in the gaussian case if necessary
#
    if(is.null(sigma2)) {
        sigma2 <- mean(apply(y,1,IQRdiff)^2*vwghts)
	cat("Estimated variance: ", signif(sigma2,4),"\n")
	}
    if(length(sigma2)==1){
#   homoskedastic Gaussian case
    lambda <- lambda*sigma2 
    cpar$tau1 <- cpar$tau1*sigma2 
    cpar$tau2 <- cpar$tau2*sigma2
    } else {
#   heteroskedastic Gaussian case
    if((dd*length(sigma2))!=length(y)) {
        cat("sigma2 does not have length 1 or same length as y")
	return("sigma2 does not have length 1 or same length as y")
	}
    sigma2<-1/sigma2 #  taking the invers yields simpler formulaes 
    }
dy <- dim(y)
if(length(dy)==2) {
   form <- "uni"
   ddim  <- 1
   n1 <- n <- dy[2]
   n2 <- n3 <- 1
   cpar$n<-n
   cpar$kstar<-log(100)
}
if(length(dy)==3){
   form <- "bi"
   ddim  <- 2
   n1 <- dy[2]
   n2 <- dy[3]
   n3 <- 1
   n <- n1*n2
   cpar$n<-n
cpar$kstar<-log(15)
hincr <- sqrt(hincr)
}
if(length(dy)==4){
   form <- "tri"
   ddim  <- 3
   n1 <- dy[2]
   n2 <- dy[3]
   n3 <- dy[4]
   n <- n1*n2*n3
   cpar$n<-n
hincr <- hincr^(1/3)
}
if(length(dy)>4)
   return("AWS for more than 3 dimensional grids is not implemented")
#
#    Initialize  list for theta
#
     tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y, fix=rep(FALSE,n))
     zobj<-list(ai=y, bi0=rep(1,n))
     bi0old<-rep(1,n)
#
#    now select the correct aws-procedure
#
#   cases:    homoskedastic or heteroskedastic depending on length(sigma2)
#
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
if(length(sigma2)==n){
# heteroskedastic
zobj <- .Fortran("vhaws",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(dd),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi2=as.double(tobj$bi2),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       double(dd),
		       double(dd),
		       as.double(vwghts),
		       PACKAGE="aws")[c("bi","bi2","bi0","ai","hakt")]
} else {
# homoskedastic
zobj <- .Fortran("vaws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(dd),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi2=as.double(tobj$bi2),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       double(dd),
		       double(dd),
		       as.double(vwghts),
		       PACKAGE="aws")[c("bi","bi2","bi0","ai","hakt")]
}
dim(zobj$ai)<-dy
dim(zobj$bi)<-dy[-1]
dim(zobj$bi2)<-dy[-1]
dim(zobj$bi0)<-dy[-1]
if(hakt>n1/2) zobj$bi0 <- hincr^ddim*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy[-1]
dim(tobj$bi2)<-dy[-1]
dim(tobj$eta)<-dy[-1]
if(graph){
if(ddim==1){ 
par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y[1,],ylim=range(y[1,],tobj$theta[1,]),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta[1,],lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
} 
if(ddim==2){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(ddim==3){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[1,,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta[1,,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
image(tobj$eta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
} 
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
###                                                                                                            
###                                                                       
if(length(sigma2)==n){
# heteroskedastic Gaussian case 
vartheta <- tobj$bi2/tobj$bi^2
} else {
vartheta <- sigma2*tobj$bi2/tobj$bi^2
}
z<-list(theta=tobj$theta,ni=tobj$bi,var=vartheta,y=y,call=args)
class(z)<-"vaws.gaussian"
z
}
vawsold <- function(y,qlambda=NULL,lkern="Triangle",sigma2=NULL,
                 hinit=NULL,hincr=NULL,hmax=NULL,heta=NULL,qtau=NULL,
                 eta0=0,u=NULL,graph=FALSE,wghts=NULL,vwghts=NULL)
{
#
#
#
#
#
  IQRdiff <- function(y) IQR(diff(y))/1.908
  updtheta<-function(zobj,tobj,cpar){
    heta<-cpar$heta
    eta0<-cpar$eta0
    tau1<-cpar$tau1
    tau2<-cpar$tau2
    kstar<-cpar$kstar
    vw<-cpar$vw
    dd<-cpar$dd
    hakt<-zobj$hakt
    tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
    bi0<-zobj$bi0
    bi<-zobj$bi
    ni2<-zobj$ni2
    n<-cpar$n
    thetanew<-zobj$ai/outer(rep(1,dd),bi)
    theta<-tobj$theta
    #thetanew[outer(rep(TRUE,dd),tobj$fix,"&")]<-theta[outer(rep(TRUE,dd),tobj$fix,"&")]

    if(hakt>heta) {
      eta<-vw%*%matrix((thetanew-theta)^2,dd,n)/tau
      dim(eta)<-dim(bi0)
      eta<-bi0*eta
      eta<-(1-eta0)*(1-pmax(0,1-eta))+eta0
    } else {
      eta <- rep(eta0,length(bi))
    }

    dim(eta)<-dim(bi)
    eta[tobj$fix]<-1
    bi <- (1-eta)*bi + eta * tobj$bi
    etadd<-outer(rep(1,dd),eta)
    theta <- (1-etadd)*thetanew + etadd * theta
    list(theta=theta,bi=bi,ni2=ni2,eta=eta,fix=(eta==1))
  }
#
#    first check arguments and initialize
#
  args <- match.call()
  spmax <- 5
  dy<-dim(y)
  if(is.null(dy)) return("dim(y) should exist for multivariate data")
  dd<-dy[1]
  if(is.null(vwghts)) vwghts<-rep(1,dd)
  if(sum(vwghts^2)!=dd) vwghts<-vwghts/sqrt(sum(vwghts^2))*dd
#  this is just to use qchisq(,dd) as an approximation for tau and lambda
  if(is.null(heta)) {
    if(length(dy)==2) heta<-2
    if(length(dy)==3) heta<-2
    if(length(dy)==4) heta<-2
  }
  if(is.null(qtau)) {
    if(length(dy)==2) qtau<-.92 
    if(length(dy)==3) qtau<-.92
    if(length(dy)==4) qtau<-.92 # not yet adjusted
  }
  if(qtau<1) tau1<-qchisq(qtau,1) else tau1<-1e50
  if(qtau<1) heta<- 1e50  
  tau2<-tau1/2
  if(is.null(hmax)){
    if(length(dy)==2) hmax<-250    # uses a maximum of about 500 points
    if(length(dy)==3) hmax<-12   # uses a maximum of about 450 points
    if(length(dy)==4) hmax<-5    # uses a maximum of about 520 points
  }
  cpar<-list(heta=heta,tau1=tau1,tau2=tau2,eta0=eta0,vw=vwghts,dd=dd)
#
#          if not provided set default value for qlambda 
#
  if(is.null(qlambda)) qlambda <- .966
# qlambda>=1  eliminates the use of the stochastic penalty
#
#  check if family is implemented and set the code for the family (used in kldist) 
#
  if(is.null(cpar$eta0)) cpar$eta0<-0
#
#    set the code for the kernel (used in lkern) and set lambda
#
  lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
#
#      get lambda as quantile of appropriate chisq,
#                rescale to be consistent with the paper in  lambda
#
  if(qlambda<1) lambda <- 2*qchisq(qlambda,dd) else lambda <- 1e50
#
#     now set hinit and hincr if not provided
#
  if(is.null(hinit)||hinit<1) hinit <- 1
  if(is.null(hincr)||hincr<=1) hincr <-1.25
#
#   estimate variance in the gaussian case if necessary
#
  if(is.null(sigma2)) sigma2 <- IQRdiff(aperm(y,c(2:length(dy),1)))^2
  if(length(sigma2)==1) sigma2<-rep(sigma2,length(y)/dd)
  if((dd*length(sigma2))!=length(y)) return("incompatible length of sigma2")
  sigma2<-1/sigma2 #  taking the invers yields simpler formulaes 
# now check which procedure is appropriate
##  this is the version on a grid
  dy <- dim(y)
  if(length(dy)==2) {
    form <- "uni"
    ddim  <- 1
    n <- dy[2]
    cpar$n<-n
    cpar$kstar<-log(100)
  }
  if(length(dy)==3){
    form <- "bi"
    ddim  <- 2
    n1 <- dy[2]
    n2 <- dy[3]
    n <- n1*n2
    cpar$n<-n
    cpar$kstar<-log(15)
    if(is.null(wghts)) wghts<-c(1,1)
    hinit<-hinit/wghts[1]
    hmax<-hmax/wghts[1]
    wghts<-(wghts[2]/wghts[1])
    hincr <- sqrt(hincr)
  }
  if(length(dy)==4){
    form <- "tri"
    ddim  <- 3
    n1 <- dy[2]
    n2 <- dy[3]
    n3 <- dy[4]
    n <- n1*n2*n3
    cpar$n<-n
    if(is.null(wghts)) wghts<-c(1,1,1)
    hinit<-hinit/wghts[1]
    hmax<-hmax/wghts[1]
    wghts<-(wghts[2:3]/wghts[1])
    hincr <- hincr^(1/3)
  }
  if(length(dy)>4)
    return("AWS for more than 3 dimensional grids is not implemented")
#
#    Initialize  list for theta
#
  tobj<-list(bi= rep(1,n), ni2= rep(1,n), theta= y, fix=rep(FALSE,n))
  zobj<-list(ai=y, bi0=rep(1,n))
  bi0old<-rep(1,n)
#
#    now select the correct aws-procedure
#
#   cases:    uni
#             bi
#             tri
#
  if(form=="uni" ){
###
###              uni
###
    hakt <- hinit
    lambda0<-lambda
    if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
    while(hakt<=hmax){
      zobj <- .Fortran("vhawsuni",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(dd),
                       as.integer(n),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
		       ni2=as.double(tobj$ni2),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       double(dd),
		       double(dd),
		       as.double(vwghts),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt","ni2")]
      if(hakt>n/2) zobj$bi0 <- hincr*biold
      dim(zobj$ai)<-c(dd,n)
      biold <- zobj$bi0
      tobj<-updtheta(zobj,tobj,cpar)
      if(graph){
        par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
        plot(y[1,],ylim=range(y,tobj$theta),col=3)
        if(!is.null(u)) lines(u,col=2)
        lines(tobj$theta[1,],lwd=2)
        title(paste("Reconstruction  h=",signif(hakt,3)))
        plot(tobj$bi/zobj$bi0,type="l",ylim=range(0,1))
        lines(tobj$eta,col=2)
        title("Sum of rel.weights and eta")
      }
      if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),
                          "   MSE: ",mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
      hakt <- hakt*hincr
      lambda0<-lambda
      gc()
    }
  }
  if(form=="bi" ){
###
###             gridded      bi
###
    hakt <- hinit
    lambda0<-lambda
    if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
    while(hakt<=hmax){
      zobj <- .Fortran("vhawsbi",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(dd),
                       as.integer(n1),
                       as.integer(n2),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
		       ni2=as.double(tobj$ni2),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       double(dd),
		       double(dd),
		       as.double(vwghts),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt","ni2")]
      dim(zobj$ai)<-c(dd,n1,n2)
      dim(zobj$bi)<-c(n1,n2)
      dim(tobj$ni2)<-c(n1,n2)
      dim(zobj$bi0)<-c(n1,n2)
      if(hakt>min(n1,n2)/2) zobj$bi0 <- hincr*hincr*biold
      biold <- zobj$bi0
      tobj<-updtheta(zobj,tobj,cpar)
      dim(tobj$theta)<-c(dd,n1,n2)
      dim(tobj$bi)<-c(n1,n2)
      dim(tobj$eta)<-c(n1,n2)
      if(graph){
        par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
        image(y[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
        title("Observed Image")
        image(tobj$theta[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=range(y))
        title(paste("Reconstruction  h=",signif(hakt,3)))
        image(tobj$bi/zobj$bi0,col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Sum of rel. weights",
                    signif(min(tobj$bi/zobj$bi0),2),"-",signif(mean(tobj$bi/zobj$bi0),2),"-",signif(max(tobj$bi/zobj$bi0),2)))
        image(tobj$eta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
        title("eta")
      }
      if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),
                          "   MSE: ",mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
      hakt <- hakt*hincr
      lambda0<-lambda
      gc()
    }
  }
  if(form=="tri" ){
###
###             gridded      tri
###
    hakt <- hinit
    lambda0<-lambda
    if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
    while(hakt<=hmax){
      zobj <- .Fortran("vhawstri",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(dd),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
		       ni2=as.double(tobj$ni2),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       double(dd),
		       double(dd),
		       as.double(vwghts),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt","ni2")]
      dim(zobj$ai)<-c(dd,n1,n2,n3)
      dim(zobj$bi)<-c(n1,n2,n3)
      dim(zobj$bi0)<-c(n1,n2,n3)
      if(hakt>min(n1,n2,n3)/2) zobj$bi0 <- hincr*hincr*hincr*biold
      biold <- zobj$bi0
      tobj<-updtheta(zobj,tobj,cpar)
      if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                          mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
      hakt <- hakt*hincr
      lambda0<-lambda
      gc()
    }
    dim(tobj$theta)<-c(dd,n1,n2,n3)
    dim(tobj$bi)<-c(n1,n2,n3)
    dim(tobj$ni2)<-c(n1,n2,n3)
    dim(tobj$eta)<-c(n1,n2,n3)
  }
###                                                                       
###            end cases                                                  
###                                                                       
  z<-list(theta=tobj$theta,ni=tobj$bi,ni2=tobj$ni2,y=y,call=args)
  class(z)<-"vaws.gaussian"
  z
}
