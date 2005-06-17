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
aws <- function(y,qlambda=NULL,qtau=NULL,family="Gaussian",lkern="Triangle",aggkern="Uniform",
                 sigma2=NULL,shape=NULL,hinit=NULL,hincr=NULL,hmax=NULL,
		 heta=NULL,eta0=NULL,u=NULL,graph=FALSE,demo=FALSE,wghts=NULL)
{
#
#
IQRdiff <- function(y) IQR(diff(y))/1.908
#
#
KLdist <- function(mcode,th1,th2,bi0){
   th12<-(1-0.5/bi0)*th2+0.5/bi0*th1
   z<-switch(mcode,(th1-th2)^2,
                th1*log(th1/th12)+(1.-th1)*log((1.-th1)/(1.-th12)),
		th1*log(th1/th12)-th1+th12,
		th1/th2-1.-log(th1/th2))
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
kstar<-cpar$kstar
hakt<-zobj$hakt
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
mcode<-cpar$mcode
hakt<-zobj$hakt
bi<-zobj$bi
bi2<-zobj$bi2
thetanew<-zobj$ai/bi
theta<-tobj$theta
thetanew[tobj$fix]<-theta[tobj$fix]
if(hakt>heta) {
eta<-switch(aggkern,"Uniform"=(1-eta0)*as.numeric(zobj$bi0/tau*KLdist(mcode,thetanew,theta,max(zobj$bi0))>1)+eta0,
                    "Triangle"=(1-eta0)*pmin(1,zobj$bi0/tau*KLdist(mcode,thetanew,theta,max(zobj$bi0)))+eta0,
		    (1-eta0)*as.numeric(zobj$bi0/tau*KLdist(mcode,thetanew,theta,max(zobj$bi0))>1)+eta0)
} else {
eta <- rep(eta0,length(theta))
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
bi2 <- (1-eta)*bi2 + eta * tobj$bi2
theta <- (1-eta)*thetanew + eta * theta
list(theta=theta,bi=bi,bi2=bi2,eta=eta,fix=(eta==1))
}
#
#    first check arguments and initialize
#
args <- match.call()
#
#   set cut off point in K_{st}(x) = exp(-x) I_{x<spmax}
#
spmax <- 5
#
#          if not provided set default value for qlambda (qlambda>=1 disables the stochastic penalty)
#
if(is.null(qlambda)) qlambda <- switch(family,   Gaussian=.96,
                                                     Bernoulli=.98,
                                                     Exponential=.98,
                                                     Poisson=.98,
                                                     Weibull=.98,
                                                     Volatility=.98)
if(qlambda<.9) return("Inappropriate value of qlambda")
#
#     set approriate defaults
#
if(is.null(heta)) heta<-max(2,hinit+.5)
if(is.null(dim(y))) heta<-max(heta,switch(family,Gaussian=2,Bernoulli=2,Exponential=8,Poisson=4/min(1,mean(y)),Weibull=8,
                                                     Volatility=8))
if(length(dim(y))==2) heta<-max(heta,switch(family,Gaussian=2,Bernoulli=2,Exponential=2,Poisson=2/min(1,sqrt(mean(y))),Weibull=2,
                                                     Volatility=2))
if(length(dim(y))==3) heta<-max(heta,switch(family,Gaussian=2,Bernoulli=2,Exponential=2,Poisson=2/min(1,mean(y)^(1/3)),Weibull=2,
                                                     Volatility=2))
if(qlambda>=1){
# thats stagewise aggregation with kernel specified by aggkern
if(is.null(qtau)) qtau<-switch(family,Gaussian=.4,Bernoulli=.75,Exponential=.75,Poisson=.75,Weibull=.75,
                                                     Volatility=.75)
if(qtau==1) tau1 <- 1e50 else tau1<-qchisq(qtau,1)
if(aggkern=="Triangle") tau1<-2.5*tau1
tau2<-tau1/2
hinit<-heta
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
cpar<-list(heta=heta,tau1=tau1,tau2=tau2,eta0=eta0)
#
#          check if model is implemented and set the code for the model (used in kldist) 
#
n<-length(y)
if(family=="Weibull" && (is.null(shape) || shape<=0))
   return("Shape parameter for Weibull has to be positive")
cpar$mcode<-switch(family,Gaussian=1,Bernoulli=2,Poisson=3,Exponential=4,Weibull=4,Volatility=4,-1)
if(cpar$mcode < 0) return(paste("specified family ",family," not yet implemented"))
if(is.null(cpar$eta0)) cpar$eta0<-switch(family,Gaussian=0,Bernoulli=.25,Poisson=.5,Exponential=.25,Weibull=.25,Volatility=.25,.25)
#
#    set the code for the kernel (used in lkern) and set lambda
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
#
#      get lambda as quantile of appropriate chisq,
#           rescale to be consistent with the paper in  lambda
#
if(qlambda<1) lambda <- qchisq(qlambda,1) else lambda <- 1e50
#
#   estimate variance in the gaussian case if necessary
#
if(family=="Gaussian") {
    if(is.null(sigma2)) {
        sigma2 <- IQRdiff(y)^2
	cat("Estimated variance: ", signif(sigma2,4),"\n")
	}
    if(length(sigma2)==1){
#   homoskedastic Gaussian case
    lambda <- lambda*sigma2*2 
    cpar$tau1 <- cpar$tau1*sigma2*2 
    cpar$tau2 <- cpar$tau2*sigma2*2 
    } else {
#   heteroskedastic Gaussian case
    if(length(sigma2)!=n) {
        cat("sigma2 does not have length 1 or same length as y")
	return("sigma2 does not have length 1 or same length as y")
	}
    lambda <- lambda*2 
    cpar$tau1 <- cpar$tau1*2 
    cpar$tau2 <- cpar$tau2*2 
    sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
    }
}
#
#   specify which statistics are needed and transform data if necessary
#
weibull <- FALSE
if(family=="Weibull") {
family <- "Exponential"
y <- y^shape
weibull <- TRUE
} else {
shape <- 1
}
if(family=="Volatility"){
family <- "Exponential"
y <- y^2
lambda <- 2*lambda 
# this accounts for the additional 1/2 in Q(\hat{theta},theta)
weibull <- TRUE
shape <- 2
}
#
#     now set hinit and hincr if not provided
#
if(is.null(hinit)||hinit<1) hinit <- 1
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
   return("AWS for more than 3 dimensional grids is not implemented")
#
#    Initialize  list for theta
#
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
     tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y, fix=rep(FALSE,n))
     zobj<-list(ai=y, bi0= rep(1,n))
     bi0old<-rep(1,n)
###
###              gridded   ( 1D -- 3D )
###
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e10 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
if(family=="Gaussian"&length(sigma2)==n){
# heteroskedastic Gaussian case
zobj <- .Fortran("chaws",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt")]
} else {
# all other cases
zobj <- .Fortran("caws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt")]
}
gc()
dim(zobj$ai)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^ddim*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar,aggkern)
gc()
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
if(graph){
if(ddim==1){ 
par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
} 
if(ddim==2){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(ddim==3){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
} 
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    signif(mean((tobj$theta-u)^2),3),"   MAE: ",signif(mean(abs(tobj$theta-u)),3)," mean(bi)=",signif(mean(tobj$bi),3),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
if(weibull) tobj$theta <- tobj$theta^(1/shape)
###                                                                       
###            end cases                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if eta0=0 and aggkern="Uniform", or if qtau1=1 
###   
if( family=="Gaussian"&length(sigma2)==n){
# heteroskedastic Gaussian case 
vartheta <- tobj$bi2/tobj$bi^2
} else {
vartheta <- switch(family,Gaussian=sigma2,
                          Bernoulli=tobj$theta*(1-tobj$theta),
			  Poisson=tobj$theta,
			  Exponential=tobj$theta^2,
			  Weibull=tobj$theta^2*(gamma(2/shape+1)/gamma(1/shape+1)^2-1),
			  Volatility=2*tobj$theta,0)*tobj$bi2/tobj$bi^2
}
z<-list(theta=tobj$theta,ni=tobj$bi,var=vartheta,y=y,hmax=hakt/hincr,call=args)
class(z)<-switch(family,Gaussian="aws.gaussian",Bernoulli="aws.bernoulli",Exponential="aws.exponential",
                 Poisson="aws.poisson",Weibull="aws.weibull",Volatility="aws.vola")
z
}

