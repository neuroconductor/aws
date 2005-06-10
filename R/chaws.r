#
#     Tests of new adaptive control for the bivariate case and Gaussian case
#
#    R - function  caws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
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
#             uni:         heta=2   tau1=3     tau2=2
#             bi:          heta=2   tau1=1     tau2=1.5
#             tri:         heta=2   tau1=(1)   tau2=(1.5)  # not yet adjusted
#             Gaussian:      qlambda=.966
#
#
chaws <- function(y,x=NULL,qlambda=NULL,lkern="Triangle",
                 sigma2=NULL,shape=NULL,hinit=NULL,hincr=NULL,hmax=NULL,
		 heta=NULL,tau1=NULL,tau2=NULL,eta0=0,NN=FALSE,u=NULL,
                 graph=FALSE,demo=FALSE,wghts=NULL)
{
IQRdiff <- function(y) IQR(diff(y))/1.908
updtheta<-function(zobj,tobj,cpar){
heta<-cpar$heta
eta0<-cpar$eta0
tau1<-cpar$tau1
tau2<-cpar$tau2
kstar<-cpar$kstar
hakt<-zobj$hakt
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
bi<-zobj$bi
thetanew<-zobj$ai/bi
theta<-tobj$theta
thetanew[tobj$fix]<-theta[tobj$fix]
if(hakt>heta) {
eta<-zobj$bi0/tau*(thetanew-theta)^2
eta<-(1-eta0)*(1-pmax(0,1-eta))+eta0
} else {
eta <- rep(eta0,length(theta))
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
theta <- (1-eta)*thetanew + eta * theta
list(theta=theta,bi=bi,eta=eta,fix=(eta==1))
}
#
#    first check arguments and initialize
#
args <- match.call()
spmax <- 5
if(is.null(heta)) heta<-2
if(is.null(tau1)){
if(is.null(dim(y))) tau1<-3 
if(length(dim(y))==2) tau1<-1
if(length(dim(y))==3) tau1<-1 # not yet adjusted
}
if(is.null(tau2)){
if(is.null(dim(y))) tau2<-2 
if(length(dim(y))==2) tau2<-1.5
if(length(dim(y))==3) tau2<-1.5 # not yet adjusted
}
if(is.null(hmax)){
if(is.null(dim(y))) hmax<-250    # uses a maximum of about 500 points
if(length(dim(y))==2) hmax<-12   # uses a maximum of about 450 points
if(length(dim(y))==3) hmax<-5    # uses a maximum of about 520 points
}
cpar<-list(heta=heta,tau1=tau1,tau2=tau2,eta0=eta0)
#
#          if not provided set default value for qlambda 
#
if(is.null(qlambda)) qlambda <- .966
if(qlambda<.6) return("Inappropriate value of qlambda")
# qlambda>=1  eliminates the use of the stochastic penalty
#
#          check if family is implemented and set the code for the family (used in kldist) 
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
if(qlambda<1) lambda <- 2*qchisq(qlambda,1) else lambda <- 1e50
#
#   estimate variance in the gaussian case if necessary
#
    if(is.null(sigma2)) sigma2 <- IQRdiff(y)^2
    if(length(sigma2)==1) sigma2<-rep(sigma2,length(y))
    if(length(sigma2)!=length(y)) return("incompatible length of sigma2")
    sigma2<-1/sigma2 #  taking the invers yields simpler formulaes 
#
#     now set hinit and hincr if not provided
#
if(is.null(hinit)||hinit<1) hinit <- 1
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
gridded <- is.null(x)
if(gridded){
##  this is the version on a grid
dy <- dim(y)
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n <- length(y)
   cpar$kstar<-log(100)
}
if(length(dy)==2){
   form <- "bi"
   ddim  <- 2
n1 <- dy[1]
n2 <- dy[2]
n <- n1*n2
cpar$kstar<-log(15)
if(is.null(wghts)) wghts<-c(1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2]/wghts[1])
hincr <- sqrt(hincr)
}
if(length(dy)==3){
   form <- "tri"
   ddim  <- 3
n1 <- dy[1]
n2 <- dy[2]
n3 <- dy[3]
n <- n1*n2*n3
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
hincr <- hincr^(1/3)
}
if(length(dy)>3)
   return("AWS for more than 3 dimensional grids is not implemented")
} else {
# not gridded
dx <- dim(x)
ddim <- 1
if(is.null(dx)&&NN) {
#
#    order data by order of x
#
    form <- "uni"
    n <- length(x)
    cpar$kstar<-log(100)
    if(n!=length(y)) return("incompatible lengths of x and y")
    ox <- order(x)
    x <- x[ox]
    y <- y[ox]
}else {
   if(is.null(dx)){
      px <- 1
      n <- length(x)
    cpar$kstar<-log(100)
   }else{
   px <- dx[1]
   n <- dx[2]
    cpar$kstar<-log(100)
   }
   form <- "multi"
   if(n!=length(y)) return("incompatible dimensions of x and y")
   weights <- rep(1,px)
#
#  now generate matrix of nearest neighbors
#  hmax is interpreted as maximal number of neighbors
#
   if(NN){
   ihmax <- trunc(hmax)
   if(ihmax>n) ihmax <- n
   neighbors <- matrix(0,ihmax,n)
   for (i in 1:n) {
      adist <- weights%*%((x-x[,i])^2)
      neighbors[,i] <- order(adist)[1:ihmax]
      }
   } 
   }
   if(length(y)!=n) return("incompatible dimensions of x and y")
   }
#
#    Initialize  list for theta
#
     tobj<-list(bi= rep(1,n), theta= y, fix=rep(FALSE,n))
     zobj<-list(ai=y, bi0=rep(1,n))
     bi0old<-rep(1,n)
#
#    now select the correct aws-procedure
#
#   cases:    gridded      uni
#             gridded      bi
#             gridded      tri
#             !gridded     multi
#             !gridded     multi, Nearest Neighbor
#
if(gridded &&  form=="uni" ){
###
###              gridded     uni
###
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
zobj <- .Fortran("chawsuni",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt")]
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
if(graph){
par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi/zobj$bi0,type="l",ylim=range(0,1))
lines(tobj$eta,col=2)
title("Sum of rel.weights and eta")
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
}
      if(gridded &&  form=="bi" ){
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
zobj <- .Fortran("chawsbi",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt")]
gc()
dim(zobj$ai)<-c(n1,n2)
if(hakt>min(n1,n2*wghts)/2) zobj$bi0 <- hincr*hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
gc()
dim(tobj$theta)<-c(n1,n2)
dim(tobj$bi)<-c(n1,n2)
dim(tobj$eta)<-c(n1,n2)
if(graph){
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=range(y))
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi/zobj$bi0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of rel. weights",
signif(min(tobj$bi/zobj$bi0),2),"-",signif(mean(tobj$bi/zobj$bi0),2),"-",signif(max(tobj$bi/zobj$bi0),2)))
image(tobj$eta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
}
      if(gridded &&  form=="tri" ){
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
zobj <- .Fortran("chawstri",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt")]
gc()
dim(zobj$ai)<-c(n1,n2,n3)
if(hakt>min(n1,n2*wghts[1],n3*wghts[2])/2) zobj$bi0 <- hincr*hincr*hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
gc()
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
dim(tobj$theta)<-c(n1,n2,n3)
dim(tobj$bi)<-c(n1,n2,n3)
dim(tobj$eta)<-c(n1,n2,n3)
}
      if( form=="multi" ){
###
###                        multi (nongridded)    p==0 
###
if(NN){
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
# now run aws-cycle
while(hakt<=hmax){
ihakt <- min(ihmax,trunc(hakt))
zobj <- .Fortran("chawsmnn",
              as.double(y),
              as.logical(tobj$fix),
              as.double(sigma2),
              as.integer(neighbors),
              as.integer(n),
              as.integer(ihmax),
              hakt=as.integer(ihakt),
              as.double(lambda0),
              as.double(tobj$theta),
              bi=as.double(tobj$bi),
              bi0=as.double(zobj$bi0),
              ai=as.double(zobj$ai),
              as.integer(lkern),
              as.double(spmax),
              PACKAGE="aws")[c("bi","bi0","ai","hakt")]
gc()
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
gc()
if(!is.null(u)) {
cat("bandwidth: ",signif(hakt,3),"   MSE: ",
    mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
}
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
} else {
dpd <- 2
# now run aws-cycle
hakt <- hinit
while(hakt<=hmax){
ihakt <- sum(maxdist<=hakt)
zobj <- .Fortran("chawsmul",
              as.double(y),
              as.logical(tobj$fix),
              as.double(x),
              as.double(sigma2),
              as.integer(n),
	      as.integer(dx),
              hakt=as.double(hakt),
              as.double(lambda),
              as.double(tobj$theta),
              bi=as.double(tobj$bi),
              bi0=as.double(zobj$bi0),
              ai=as.double(zobj$ai),
              as.integer(lkern),
	      as.double(spmax),
	      PACKAGE="aws")[c("bi","bi0","ai","hakt")]
gc()
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
gc()
if(!is.null(u)) {
cat("bandwidth: ",signif(hakt,3),"   MSE: ",
    mean((tobj$theta-u)^2),"   MAE: ",mean(abs(tobj$theta-u)),"\n")
}
hakt <- hakt*hincr
gc()
}
}
}
###                                                                       
###            end cases                                                  
###                                                                       
z<-list(theta=tobj$theta,ni=tobj$bi,y=y,x=x,call=args)
class(z)<-"laws.gaussian"
z
}
