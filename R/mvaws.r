#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
#    emaphazises on the propagation-separation approach 
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
#             sagg:          heta=2   qtau=.95     
#             Gaussian:      qlambda=.96          
#             Bernoulli:     qlambda=.98       
#             Poisson:       qlambda=.98       
#             Exponential:   qlambda=.98       
#             Weibull:       qlambda=.98     
#             Volatility:    qlambda=.98       
#
mvaws <- function(y,qlambda=NULL,qtau=NULL,lkern="Triangle",aggkern="Uniform",
                 sigma2=NULL,df=NULL,hinit=NULL,hincr=NULL,hmax=NULL,lseq=NULL,
		 heta=NULL,eta0=NULL,u=NULL,graph=FALSE,demo=FALSE,wghts=NULL,spmax=5)
{
#
#          Auxilary functions
#
IQRdiff <- function(y) IQR(diff(y))/1.908
#
#
KLdist <- function(th1,th2,s1,s2,df){
   z<-(th2-th1)^2/s2+df*(log(s2/s1)-1+s1/s2)
   z[is.na(z)]<-0
   z
		}
#
#
updtheta<-function(zobj,tobj,cpar,aggkern){
heta<-cpar$heta
eta0<-cpar$eta0
tau1<-cpar$tau1
tau2<-cpar$tau2
shape<-cpar$shape
kstar<-cpar$kstar
hakt<-zobj$hakt
df<-cpar$df
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
hakt<-zobj$hakt
bi<-zobj$bi
bi2<-zobj$bi2
thetanew<-zobj$ai/bi
s2new<-zobj$si/bi
theta<-tobj$theta
s2<-tobj$sigma2
thetanew[tobj$fix]<-theta[tobj$fix]
s2new[tobj$fix]<-s2[tobj$fix]
if(hakt>heta) {
eta<-switch(aggkern,"Uniform"=as.numeric(zobj$bi0/tau*KLdist(thetanew,theta,s2new,s2,df)>1),
                    "Triangle"=pmin(1,zobj$bi0/tau*KLdist(thetanew,theta,s2new,s2,df)),
		    as.numeric(zobj$bi0/tau*KLdist(thetanew,theta,s2new,s2,df)>1))
} else {
eta <- rep(0,length(theta))
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
bi2 <- (1-eta)*bi2 + eta * tobj$bi2
theta <- (1-eta)*thetanew + eta * theta
s2 <- (1-eta)*s2new + eta * s2
list(theta=theta,bi=bi,bi2=bi2,eta=eta,sigma2=s2,fix=(eta==1))
}
#
#    first check arguments and initialize
#
args <- match.call()
#
#   set cut off point in K_{st}(x) = exp(-x) I_{x<spmax}
#
#spmax <- 5
#
#     set approriate defaults
#
mae<-NULL
if(is.null(heta)) heta<-max(2,hinit+.5)
if(is.null(dim(y))){
heta<-max(heta,8)
if(is.null(qlambda)) qlambda <- .975
if(is.null(lseq)) lseq<-1.3
}
if(length(dim(y))==2) {
heta<-max(heta,2)
if(is.null(qlambda)) qlambda <- .975
if(is.null(lseq)) lseq<-c(1.85,1.3,1.1,1.1)
}
if(length(dim(y))==3){
 heta<-max(heta,2)
if(is.null(qlambda)) qlambda <- .95
if(is.null(lseq)) lseq<-c(3,2,1.5,1.2,1.2,1.2)
}
if(qlambda<.9) warning("Inappropriate value of qlambda")
if(qlambda>=1){
# thats stagewise aggregation with kernel specified by aggkern
if(is.null(qtau)) qtau<-.4
if(qtau==1) tau1 <- 1e50 else tau1<-qchisq(qtau,1)
if(aggkern=="Triangle") tau1<-2.5*tau1
tau2<-tau1/2
} else {
if(is.null(qtau)) qtau<-.95
if(qtau>=1) {
tau1 <- 1e50 
} else {
tau1<-qchisq(qtau,1)
}
if(aggkern=="Triangle") tau1<-2.5*tau1
tau2<-tau1/2
}
if(is.null(hmax)){
if(is.null(dim(y))) hmax<-250    # uses a maximum of about 500 points
if(length(dim(y))==2) hmax<-12   # uses a maximum of about 450 points
if(length(dim(y))==3) hmax<-5    # uses a maximum of about 520 points
}
cpar<-list(heta=heta,tau1=tau1,tau2=tau2)
#
#          check if model is implemented and set the code for the model (used in kldist) 
#
n<-length(y)
#
#    set the code for the kernel (used in lkern) and set lambda
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,Gaussian=5,2)
#
#      get lambda as quantile of appropriate chisq,
#           rescale to be consistent with the paper in  lambda
#
if(qlambda<1) lambda <- qchisq(qlambda,2) else lambda <- 1e50
if(is.null(df)||df<0) return("Positive degrees of freedom (shape) required ")
if(length(sigma2)!=length(y)) return("Incorrect size of sigma2 ")
if(any(sigma2<=0)) return("positive elements of sigma2 required ")

cpar$df<-df
#
#     now set hinit and hincr if not provided
#
if(is.null(hinit)||hinit<1) hinit <- 1
if(lkern==5) {
    hinit <- hinit
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax*0.42445*4
    }
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
##  this is the version on a grid
dy <- dim(y)
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n1 <- n <- length(y)
   n2 <- 1
   n3 <- 1
   cpar$kstar<-log(100)
}
if(length(dy)==2){
   form <- "bi"
   ddim  <- 2
n1 <- dy[1]
n2 <- dy[2]
n3 <- 1
n <- n1*n2
cpar$kstar<-log(15)
hincr <- sqrt(hincr)
}
if(length(dy)==3){
   form <- "tri"
   ddim  <- 3
n1 <- dy[1]
n2 <- dy[2]
n3 <- dy[3]
n <- n1*n2*n3
cpar$kstar<-log(5)
hincr <- hincr^(1/3)
}
if(length(dy)>3)
   stop("AWS for more than 3 dimensional grids is not implemented")
#
#    Initialize  list for theta
#
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
     tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y, sigma2=sigma2, fix=rep(FALSE,n))
     zobj<-list(ai=y, si=sigma2, bi0= rep(1,n))
     bi0old<-rep(1,n)
###
###              gridded   ( 1D -- 3D )
###
steps<-as.integer(log(hmax/hinit)/log(hincr)+1)
if(is.null(lseq)) lseq<-1
if(length(lseq)<steps) lseq<-c(lseq,rep(1,steps-length(lseq)))
lseq<-lseq[1:steps]
k<-1
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
zobj <- .Fortran("cmvaws",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
		       as.double(df),#thats degrees of freedom in sigma^2 her
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       as.double(tobj$sigma2),
                       bi=as.double(tobj$bi),
		       bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       si=as.double(zobj$si),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       PACKAGE="aws")[c("bi","bi0","bi2","ai","si","hakt")]
gc()
dim(zobj$ai)<-dim(zobj$si)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^ddim*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar,aggkern)
gc()
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
if(graph){
if(ddim==1){ 
par(mfrow=c(1,3),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta,lwd=2)
plot(sqrt(sigma2),ylim=sqrt(range(sigma2,tobj$sigma2)),col=3)
lines(sqrt(tobj$sigma2),lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
} 
if(ddim==2){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(sqrt(tobj$sigma2),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Standard deviations  h=",signif(hakt,3)," min=",signif(min(sqrt(tobj$sigma2)),3)," max=",signif(max(sqrt(tobj$sigma2)),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
}
if(ddim==3){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(sqrt(tobj$sigma2[,,n3%/%2+1]),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Standard deviations  h=",signif(hakt,3)," min=",signif(min(sqrt(tobj$sigma2)),3)," max=",signif(max(sqrt(tobj$sigma2)),3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
} 
}
if(!is.null(u)) {
cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    signif(mean((tobj$theta-u)^2),3),"   MAE: ",signif(mean(abs(tobj$theta-u)),3)," mean(bi)=",signif(mean(tobj$bi),3),"\n")
mae<-c(mae,signif(mean(abs(tobj$theta-u)),3))
		    }
if(demo) readline("Press return")
hakt <- hakt*hincr
lambda0<-lambda*lseq[k]
k<-k+1
gc()
}
###                                                                       
###            end cases                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if eta0=0 and aggkern="Uniform", or if qtau1=1 
###   
# heteroskedastic Gaussian case 
varfactor <- tobj$bi2/tobj$bi^2
z<-list(theta=tobj$theta,ni=varfactor,var=tobj$sigma2,y=y,sigma2=sigma2,hmax=hakt/hincr,mae=mae,lseq=c(0,lseq[-steps]),call=args)
class(z)<-"aws.meanvar"
z
}
