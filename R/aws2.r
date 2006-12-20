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
aws2 <- function(y,hmax=NULL,qlambda=NULL,family="Gaussian",
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
qtau <- 1
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
    if(qlambda==1&qtau==1) hinit <- hmax
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
zobj <- .Fortran("caws2",as.double(y),
                       fix=as.logical(tobj$fix),
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
		       PACKAGE="aws",DUP=FALSE)[c("fix","bi","bi0","bi2","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^d*biold
biold <- zobj$bi0
tobj$fix <- zobj$fix
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$fix)<-dy
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
oldpar<-par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$fix*max(tobj$bi),col=2)
title("Sum of weights and eta")
} 
if(d==2){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$fix,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(d==3){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$fix[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
} 
par(oldpar)
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
