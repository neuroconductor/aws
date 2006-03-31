#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
#    emaphazises on the propagation-separation approach 
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
#     default parameters:  see function setawsdefaults
#       
aws <- function(y,hmax=NULL,qlambda=NULL,qtau=NULL,family="Gaussian",
                sigma2=NULL,scorr=0,shape=NULL,wghts=NULL,graph=FALSE,demo=FALSE,
		lkern="Triangle",skern="Triangle",aggkern="Uniform",
		spmin=0,spmax=5,lseq=NULL,u=NULL,testprop=FALSE)
{
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
#
#   set appropriate defaults
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,
	            Gaussian=5,2)
skern <- switch(skern,"Exp"=1,"Triangle"=2,2)
cpar<-setawsdefaults(dy,mean(y),family,skern,aggkern,qlambda,qtau,lseq,hmax,shape,spmax)
lambda <- cpar$lambda
hmax <- cpar$hmax
lseq <- cpar$lseq
shape <- cpar$shape
d <- cpar$d
hinit <- cpar$hinit
hincr <- cpar$hincr
spmax <- cpar$spmax
n<-length(y)
# 
#   family dependent transformations that depend on the value of family
#
zfamily <- awsfamily(family,y,sigma2,shape,scorr,lambda,cpar)
cpar <- zfamily$cpar
lambda <- zfamily$lambda
sigma2 <- zfamily$sigma2
h0 <- zfamily$h0
y <- zfamily$y
rm(zfamily)
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax*0.42445*4
    hinit <- 0.42445*4
    }
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
##  this is the version on a grid
n <- length(y)
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
#
#    Initialize  for the iteration
#
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y/shape, fix=rep(FALSE,n))
zobj<-list(ai=y, bi0= rep(1,n))
biold<-rep(1,n)
if(family=="Gaussian"&length(sigma2)==n) vred<-rep(1,n)
mae<-NULL
hakt <- hinit
hakt0 <- hinit
lambda0<-lambda
lambda0<-1e50 # that removes the stochstic term for the first step, Initialization by kernel estimates
if(testprop) {
#
#  prepare to  check for alpha in propagation condition (to adjust qlambda and lseq)
#
       if(is.null(u)) u <- 0
       cpar <- c(cpar, list(n1=n1,n2=n2,n3=n3,n=n1*n2*n3,lkern=lkern,skern=skern,wghts=wghts,spmin=spmin,spmax=spmax,family=family,u=u))
       propagation <- NULL
    } 
#
#   iteratate until maximal bandwidth is reached
#
steps <- as.integer(log(hmax/hinit)/log(hincr))
cat("Progress:")
for(k in 0:steps){
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
if(family=="Gaussian"&scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
# Correction for spatial correlation depends on h^{(k)} 
hakt0<-hakt
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
		       vred=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.integer(skern),
	               as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","bi0","bi2","vred","ai","hakt")]
vred[!tobj$fix]<-zobj$vred[!tobj$fix]
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
                       as.integer(skern),
                       as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","bi0","bi2","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^d*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
#
#  if testprop == TRUE
#  check alpha in propagation condition (to adjust qlambda and lseq)
#  
if(testprop) propagation <- awstestprop(y,family,tobj,zobj,sigma2,hakt,cpar,u,propagation)
if(graph){
#
#     Display intermediate results if graph == TRUE
#
if(d==1){ 
par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
} 
if(d==2){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(d==3){ 
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
} 
}
#
#    Calculate MAE and MSE if true parameters are given in u 
#    this is for demonstration and testing for propagation (parameter adjustments) 
#    only.
#
if(!is.null(u)) {
   cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    signif(mean((tobj$theta-u)^2),3),"   MAE: ",
		    signif(mean(abs(tobj$theta-u)),3)," mean(bi)=",
		    signif(mean(tobj$bi),3),"\n")
   mae<-c(mae,signif(mean(abs(tobj$theta-u)),3))
		    }
if(demo) readline("Press return")
#
#   Prepare for next iteration
#
hakt <- hakt*hincr
x<-1.25^k
scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
lambda0<-lambda*lseq[k+1]*scorrfactor
cat(paste(signif(sum(hincr^(2*(0:k)))/sum(hincr^(2*(0:steps)))*100,2),"% ",sep=""))
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if aggkern="Uniform", or if qtau1=1 
###   
if( family=="Gaussian"&length(sigma2)==n){
# heteroskedastic Gaussian case 
vartheta <- tobj$bi2/tobj$bi^2
} else {
vartheta <- switch(family,Gaussian=sigma2,
                          Bernoulli=tobj$theta*(1-tobj$theta),
			  Poisson=tobj$theta,
			  Exponential=tobj$theta^2,
			  Volatility=2*tobj$theta,
			  Variance=2*tobj$theta,0)*tobj$bi2/tobj$bi^2
vred<-tobj$bi2/tobj$bi^2
}
if( family=="Gaussian"){
vartheta<-vartheta/Spatialvar.gauss(hakt/0.42445/4,h0+1e-5,d)*Spatialvar.gauss(hakt/0.42445/4,1e-5,d)
}
z<-list(theta=tobj$theta,ni=tobj$bi,var=vartheta,vred=vred,y=y,
        hmax=hakt/hincr,mae=mae,lseq=c(0,lseq),call=args)
class(z)<-switch(family,Gaussian="aws.gaussian",Bernoulli="aws.bernoulli",Exponential="aws.exponential",
                 Poisson="aws.poisson",Volatility="aws.vola",Variance="aws.var")
z
}
#######################################################################################
#
#        Auxilary functions
#
#######################################################################################
#
#        Set default values
#
# default values for qlambda, lseq and qtau are chosen by propagation condition
# (strong version) with alpha=0.1
# see script aws_propagation.r
#
#######################################################################################
setawsdefaults <- function(dy,meany,family,skern,aggkern,qlambda,qtau,lseq,hmax,shape,spmax){
hinit <- 1
if(!is.null(qlambda)&&qlambda<.9){
   cat("Inappropriate value of qlambda, using defaults")
   qlambda <- NULL
}
#
#   univariate case
#
if(is.null(dy)){ 
      d<-1
      if(is.null(qlambda)) qlambda <- switch(skern,
                                             switch(family,
                                                   Gaussian=.992,
					           Bernoulli=.985,
					           Exponential=.985,
					           Poisson=.981,
                                                   Volatility=.985,
					           Variance=.985,
					           .992),
					     switch(family,
                                                    Gaussian=.98,
					            Bernoulli=.97,
					            Exponential=.965,
					            Poisson=.965,
                                                    Volatility=.965,
					            Variance=.965,
					            .97))
      if(is.null(lseq)) lseq<-switch(family,
                                     Gaussian=1.5,
				     Bernoulli=1,
				     Exponential=3.9,
				     Poisson=1,
				     Volatility=9,
				     Variance=9,
				     1.5)
     }
if(length(dy)==2) {
     d<-2
     if(is.null(qlambda)) qlambda <- switch(skern,
                                            switch(family,   
                                                   Gaussian=.992,
					           Bernoulli=.99,
					           Exponential=.99,
					           Poisson=.99,
                                                   Volatility=.99,
					           Variance=.99,
					           .992),
                                            switch(family,   
                                                   Gaussian=.97,
					           Bernoulli=.965,
					           Exponential=.965,
					           Poisson=.965,
                                                   Volatility=.965,
					           Variance=.965,
					           .97))
     if(is.null(lseq)) lseq<-switch(family,
                                    Gaussian=c(1.8,1.3,1.2,1.2,1.1,1.1,1.1),
				    Bernoulli=1,
				    Exponential=4.8,
				    Poisson=1,
				    Volatility=11,
				    Variance=11,
				    c(1.8,1.3,1.2,1.2,1.1,1.1,1.1))
	}
if(length(dy)==3){
     d<-3
     if(is.null(qlambda)) qlambda <- switch(skern,
                                            switch(family,   
                                                   Gaussian=.995,
					           Bernoulli=.995,
					           Exponential=.995,
					           Poisson=.995,
                                                   Volatility=.995,
					           Variance=.995,
					           .995),
                                            switch(family,   
                                                   Gaussian=.97,
					           Bernoulli=.97,
					           Exponential=.975,
					           Poisson=.965,
                                                   Volatility=.975,
					           Variance=.975,
					           .97))
     if(is.null(lseq)) lseq<-switch(family,
                                    Gaussian=c(1.9,1.5,1.3,1.3,1.3,1.3,rep(1.1,8)),
				    Bernoulli=1,
				    Exponential=4.7,
				    Poisson=1,
				    Volatility=11,
				    Variance=11,
				    c(1.9,1.5,1.3,1.3,1.3,1.3,rep(1.1,8)))
}
if(qlambda<1) lambda <- qchisq(qlambda,1) else lambda <- 1e50
hincr <- 1.25^(1/d)
#
#    determine heta for memory step
#
heta <- switch(family,   Gaussian=1.25,
			 Bernoulli=5/meany/(1-meany),
			 Exponential=20,
			 Poisson=max(1.25,20/meany),
                         Volatility=20,
			 Variance=20,
			 1.25)^(1/d)
kstar <- log(switch(d,100,15,5))
#
# stagewise aggregation 
#
if(qlambda>=1){
#
# force aggkern = "Triangle"  here
#
#  we need a tau for stagewise aggregation that fulfils the propagation condition
#
    aggkern <- "Triangle"
    cat("Stagewise aggregation: Triangular aggregation kernel is used\n")
#
#   this is the first bandwidth to consider for stagewise aggregation
#
    if(is.null(qtau)) qtau <- switch(d,
                                 switch(family,
                                        Gaussian=.4,
			                Bernoulli=.45,
			                Exponential=.75,
			                Poisson=.45,
                                        Volatility=.75,
			                Variance=.75,
			                .35),
			         switch(family,
                                        Gaussian=.35,
			                Bernoulli=.45,
			                Exponential=.75,
			                Poisson=.35,
                                        Volatility=.75,
			                Variance=.75,
			                .3),
                                 switch(family,
                                        Gaussian=.35,
			                Bernoulli=.45,
			                Exponential=.75,
			                Poisson=.35,
                                        Volatility=.75,
			                Variance=.75,
			                .25)) 
    if(qtau>=1){
        heta <- hmax
        cat("Neither PS nor Stagewise Aggregation is specified, compute kernel estimate with bandwidth hmax\n") 
         hinit <- heta
    }
} else {
#
#   set appropriate value for qtau (works for all families)
#
if(is.null(qtau)) qtau <- .995
if(qtau>=1) heta <- 1e10 # no memory step
}
if(qtau>=1) tau1 <- 1e50 else tau1<-qchisq(qtau,1)
#
#    adjust for different aggregation kernels
#
if(aggkern=="Triangle") tau1<-2.5*tau1
tau2<-tau1/2
#
#   set maximal bandwidth
#
if(is.null(hmax)) hmax <- switch(d,250,12,5)
# uses a maximum of about 500, 450 and 520  points, respectively.
mcode <- switch(family,
                Gaussian=1,
		Bernoulli=2,
		Poisson=3,
		Exponential=4,
		Volatility=4,
		Variance=5,
		-1)
if(mcode < 0) stop(paste("specified family ",family," not yet implemented"))
if(skern==2) {
   lambda<-lambda*1.8
   spmax <- 1
   } else {
   if(is.null(spmax)) spmax <- 5
   }
if(is.null(shape)) shape<-1
steps <- as.integer(log(hmax/hinit)/log(hincr))
lseq <- c(lseq,rep(1,steps))[1:steps]
if(qlambda<1) 
cat("Running PS with lambda=",signif(lambda,3)," hmax=",hmax," number of iterations=",steps," memory step",if(qtau>=1) "OFF" else "ON","\n")
else cat("Sequence tau_k:",signif((tau1+tau2*pmax(kstar-log(hincr^(1:steps)),0)),3),"\n")
list(heta=heta,tau1=tau1,tau2=tau2,lambda=lambda,lseq=lseq,hmax=hmax,d=d,mcode=mcode,shape=shape,aggkern=aggkern,hinit=hinit,kstar=kstar,hincr=hincr,spmax=spmax)
}
#######################################################################################
#
#    IQRdiff (for robust variance estimates
#
#######################################################################################
IQRdiff <- function(y) IQR(diff(y))/1.908
#######################################################################################
#
#    Kullback-Leibler distances
#
#######################################################################################
KLdist <- function(mcode,th1,th2,bi0){
   th12<-(1-0.5/bi0)*th2+0.5/bi0*th1
   z<-switch(mcode,(th1-th2)^2,
                th1*log(th1/th12)+(1.-th1)*log((1.-th1)/(1.-th12)),
		th1*log(th1/th12)-th1+th12,
		th1/th2-1.-log(th1/th2),
		th1/th2-1.-log(th1/th2))
   z[is.na(z)]<-0
   z
		}
####################################################################################
#
#    Memory step for local constant aws
#
####################################################################################
updtheta<-function(zobj,tobj,cpar){
heta<-cpar$heta
hakt<-zobj$hakt
bi<-zobj$bi
bi2<-zobj$bi2
thetanew<-zobj$ai/bi
if(hakt>heta) {
#
#   memory step
#
mcode<-cpar$mcode
aggkern <- cpar$aggkern
tau1<-cpar$tau1
tau2<-cpar$tau2
kstar<-cpar$kstar
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
theta<-tobj$theta
thetanew[tobj$fix]<-theta[tobj$fix]
eta<-switch(aggkern,"Uniform"=as.numeric(zobj$bi0/tau*
                     KLdist(mcode,thetanew,theta,max(zobj$bi0))>1),
                    "Triangle"=pmin(1,zobj$bi0/tau*
		     KLdist(mcode,thetanew,theta,max(zobj$bi0))),
		    as.numeric(zobj$bi0/tau*
		     KLdist(mcode,thetanew,theta,max(zobj$bi0))>1))
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
bi2 <- (1-eta)*bi2 + eta * tobj$bi2
thetanew <- (1-eta)*thetanew + eta * theta
} else {
#
#  no memory step
#
eta <- rep(0,length(thetanew))
}
list(theta=thetanew,bi=bi,bi2=bi2,eta=eta,fix=(eta==1))
}
####################################################################################
#
#    Regularize for Bernoulli and Poisson models if parameter estimates are
#    at the border of their domain
#
####################################################################################
regularize <- function(zobj,family){
if(family%in%c("Bernoulli","Poisson")){
ind <-(zobj$ai<1)
zobj$ai[ind]<-1
zobj$bi[ind]<-zobj$bi[ind]+1-zobj$ai[ind]
if(!is.null(zobj$bi0)) zobj$bi0[ind]<-zobj$bi0[ind]+1-zobj$ai[ind]
if(family=="Bernoulli") {
ind <- zobj$bi-zobj$ai < 1
zobj$bi[ind] <- 1 + zobj$ai[ind]
if(!is.null(zobj$bi0)) zobj$bi0[ind]<-zobj$bi0[ind]+1-(zobj$bi-zobj$ai)[ind]
}
}
zobj
} 
############################################################################
#
#   transformations that depend on the specified family
#
############################################################################
awsfamily <- function(family,y,sigma2,shape,scorr,lambda,cpar){
h0 <- 0
if(family=="Gaussian") {
  d <- cpar$d
  if(scorr[1]>0) {
         h0<-numeric(length(scorr))
         for(i in 1:length(h0))
         h0[i]<-geth.gauss(scorr[i])
         if(length(h0)<d) h0<-rep(h0[1],d)
         cat("Corresponding bandwiths for specified correlation:",h0,"\n")
}
    if(is.null(sigma2)) {
        sigma2 <- IQRdiff(as.vector(y))^2
        if(scorr[1]>0) sigma2<-sigma2*Varcor.gauss(h0)
	cat("Estimated variance: ", signif(sigma2,4),"\n")
	}
    if(length(sigma2)==1){
#   homoskedastic Gaussian case
    lambda <- lambda*sigma2*2 
    cpar$tau1 <- cpar$tau1*sigma2*2 
    cpar$tau2 <- cpar$tau2*sigma2*2 
    } else {
#   heteroskedastic Gaussian case
    if(length(sigma2)!=n) 
	stop("sigma2 does not have length 1 or same length as y")
    lambda <- lambda*2 
    cpar$tau1 <- cpar$tau1*2 
    cpar$tau2 <- cpar$tau2*2 
    sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
    }
}
#
#   specify which statistics are needed and transform data if necessary
#
if(family=="Volatility"){
   family <- "Exponential"
   y <- y^2
   lambda <- 2*lambda 
# this accounts for the additional 1/2 in Q(\hat{theta},theta)
}
if(family=="Variance"){
   lambda <- 2*lambda/shape 
# this accounts for the additional 1/2 in Q(\hat{theta},theta) and the degrees of freedom in chisq
   cpar$tau1 <- cpar$tau1/shape
   cpar$tau2 <- cpar$tau2/shape 
}
list(cpar=cpar,lambda=lambda,y=y,sigma2=sigma2,h0=h0)
}
##################################################################################
#
#   AWS local constant Test for propagation condition 
#
##################################################################################
awstestprop <- function(y,family,tobj,zobj,sigma2,hakt,cpar,u,propagation){
dlw<-(2*trunc(hakt/c(1,cpar$wghts))+1)[1:cpar$d]
if(family=="Gaussian"&length(sigma2)==cpar$n){
# heteroskedastic Gaussian case
pobj <- .Fortran("chaws",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(cpar$n1),
                       as.integer(cpar$n2),
                       as.integer(cpar$n3),
                       hakt=as.double(hakt),
                       as.double(1e40),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(cpar$n),
                       bi0=as.double(zobj$bi0),
		       vred=double(cpar$n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(cpar$lkern),
                       as.integer(cpar$skern),
	               as.double(cpar$spmin),
		       as.double(cpar$spmax),
		       double(prod(dlw)),
		       as.double(cpar$wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","ai","hakt")]
} else {
# all other cases
pobj <- .Fortran("caws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(cpar$n1),
                       as.integer(cpar$n2),
                       as.integer(cpar$n3),
                       hakt=as.double(hakt),
                       as.double(1e40),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(cpar$n),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(cpar$lkern),
                       as.integer(cpar$skern),
                       as.double(cpar$spmin),
		       as.double(cpar$spmax),
		       double(prod(dlw)),
		       as.double(cpar$wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) pobj<-regularize(pobj,family)
ptheta <- array(pobj$ai/pobj$bi,dim(y)) 
narisk <- sum(abs(ptheta-u))
if(narisk==0) narisk<-1e10
propagation <- c(propagation,sum(abs(tobj$theta-ptheta))/narisk)
cat("Propagation with alpha=",max(propagation),"\n")
cat("alpha values:","\n")
print(rbind(cpar$lseq[1:length(propagation[-1])],signif(propagation[-1],3)))
propagation
}
