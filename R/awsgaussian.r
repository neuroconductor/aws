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
aws.gaussian <- function(y,hmax=NULL,hpre=NULL,qlambda=NULL,qtau=NULL,varmodel="Constant",
                varpar=NULL,varprop=.1,scorr=0,wghts=NULL,graph=FALSE,demo=FALSE,
		lkern="Triangle",skern="Triangle",aggkern="Uniform",
		spmin=0,spmax=5,lseq=NULL,u=NULL)
{
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
if(!(varmodel %in% c("Constant","Linear","Quadratic"))) stop("Model for variance not implemented")
#
#   set appropriate defaults
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,
	            Gaussian=5,2)
skern <- switch(skern,"Exp"=1,"Triangle"=2,2)
cpar<-setawsdefaults(dy,mean(y),"Gaussian",skern,aggkern,qlambda,qtau,lseq,hmax,1,spmax)
lambda <- cpar$lambda
hmax <- cpar$hmax
lseq <- cpar$lseq
shape <- cpar$shape
d <- cpar$d
cpar$heta <- 20^(1/d)
hinit <- cpar$hinit
hincr <- cpar$hincr
spmax <- cpar$spmax
n<-length(y)
# 
#   family dependent transformations 
#
zfamily <- awsgfamily(y,scorr,lambda,cpar)
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
vred<-rep(1,n)
mae<-NULL
hakt <- hinit*hincr
hakt0 <- hinit*hincr
lambda0<-lambda
lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
#
#   produce a presmoothed estimate to stabilze variance estimates
#
if(is.null(hpre)) hpre<-20^(1/d)
dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:d]
hobj <- .Fortran("caws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hpre),
                       as.double(1e40),
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
		       PACKAGE="aws",DUP=FALSE)[c("bi","ai")]
hobj$theta <- hobj$ai/hobj$bi
dim(hobj$theta) <- dim(hobj$bi) <- dy
#
#   iteratate until maximal bandwidth is reached
#
steps <- as.integer(log(hmax/hinit)/log(hincr))
cat("Progress:")
for(k in 1:steps){
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
if(scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
# Correction for spatial correlation depends on h^{(k)} 
hakt0<-hakt
# heteroskedastic Gaussian case
zobj <- .Fortran("cgaws",as.double(y),
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
		       gi=double(n),
		       vred=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.integer(skern),
	               as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","bi0","bi2","vred","ai","gi","hakt")]
vred[!tobj$fix]<-zobj$vred[!tobj$fix]
dim(zobj$ai)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^d*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
tobj$gi <- zobj$gi
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
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
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
} 
if(d==2){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
zlim <- quantile(tobj$theta,c(0.001,0.999))
image(array(pmax(pmin(tobj$theta,zlim[2]),zlim[1]),dy),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(d==3){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
zlim <- quantile(tobj$theta,c(0.001,0.999))
image(array(pmax(pmin(tobj$theta[,,n3%/%2+1],zlim[2]),zlim[1]),dy[-3]),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$eta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
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
#
#   Create new variance estimate
#
vobj <- awsgsigma2(y,hobj,tobj,varmodel,varpar,varprop,h0)
sigma2 <- vobj$sigma2inv
coef <- vobj$coef
rm(vobj)
hakt <- hakt*hincr
x<-1.25^(k-1)
scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
lambda0<-lambda*lseq[k]*scorrfactor
cat(paste(signif(sum(hincr^(2*(1:k)))/sum(hincr^(2*(1:steps)))*100,2),"% ",sep=""))
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if aggkern="Uniform", or if qtau1=1 
###   
if(length(sigma2)==n){
# heteroskedastic case 
vartheta <- tobj$bi2/tobj$bi^2
} else {
# homoskedastic case 
vartheta <- sigma2*tobj$bi2/tobj$bi^2
vred<-tobj$bi2/tobj$bi^2
}
vartheta<-vartheta/Spatialvar.gauss(hakt/0.42445/4,h0+1e-5,d)*Spatialvar.gauss(hakt/0.42445/4,1e-5,d)
z<-list(theta=tobj$theta,sigma2=1/sigma2,ni=tobj$bi,var=vartheta,vred=vred,y=y,
        hmax=hakt/hincr,mae=mae,lseq=c(0,lseq),call=args)
class(z)<-"aws.gaussian"
z
}
##############################################################################################
#
#   AWS - segmentation
#
###############################################################################################
aws.segment <- function(y,hmax=NULL,hpre=NULL,qlambda=NULL,sigma2=NULL,varmodel="Constant",
                varpar=NULL,varprop=.1,scorr=0,wghts=NULL,graph=FALSE,demo=FALSE,
		lkern="Triangle",skern="Triangle",aggkern="Uniform",
		spmin=0,spmax=5,lseq=NULL,u=NULL,nlevels=2,swghts=NULL,rhoseq=NULL,
		rho=.95)
{
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
if(!(varmodel %in% c("None","Constant","Linear","Quadratic"))) stop("Model for variance not implemented")
if(varmodel=="None"&is.null(sigma2)) stop("Specify either a variance model or a variance")
#
#   set appropriate defaults
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,
	            Gaussian=5,2)
skern <- switch(skern,"Exp"=1,"Triangle"=2,2)
if(is.null(swghts)) swghts<-rep(1/nlevels,nlevels)
cpar<-setawsdefaults(dy,mean(y),"Gaussian",skern,aggkern,qlambda,1,lseq,hmax,1,spmax)
lambda <- cpar$lambda
hmax <- cpar$hmax
lseq <- cpar$lseq
shape <- cpar$shape
d <- cpar$d
cpar$heta <- 20^(1/d)
#hinit <- cpar$hinit
hinit <- 1
hincr <- cpar$hincr
spmax <- cpar$spmax
n<-length(y)
# 
#   family dependent transformations 
#
zfamily <- awsgfamily(y,scorr,lambda,cpar)
cpar <- zfamily$cpar
lambda <- zfamily$lambda
if(varmodel!="None") sigma2 <- zfamily$sigma2 else sigma2 <- 1/sigma2
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
vred<-rep(1,n)
mae<-NULL
hakt <- hinit*hincr
hakt0 <- hinit*hincr
lambda0<-lambda
lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
#
#   produce a presmoothed estimate to stabilze variance estimates
#
if(is.null(hpre)) hpre<-20^(1/d)
dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:d]
hobj <- .Fortran("caws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hpre),
                       as.double(1e40),
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
		       PACKAGE="aws",DUP=FALSE)[c("bi","ai")]
hobj$theta <- hobj$ai/hobj$bi
dim(hobj$theta) <- dim(hobj$bi) <- dy
#
#   iteratate until maximal bandwidth is reached
#
steps <- as.integer(log(hmax/hinit)/log(hincr))
if(is.null(rhoseq)) rhoseq <- 1 - rho^(0:steps) 
cat("Progress:")
for(k in 1:steps){
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
if(scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
# Correction for spatial correlation depends on h^{(k)} 
hakt0<-hakt
# heteroskedastic Gaussian case
zobj <- .Fortran("cgaws",as.double(y),
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
		       gi=double(n),
		       vred=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.integer(skern),
	               as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","bi0","bi2","vred","ai","gi","hakt")]
vred[!tobj$fix]<-zobj$vred[!tobj$fix]
dim(zobj$ai)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^d*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
tobj$gi <- zobj$gi
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
#
#   Create new variance estimate
#
if(varmodel!="None") {
vobj <- awsgsigma2(y,hobj,tobj,varmodel,varpar,varprop,h0)
sigma2 <- vobj$sigma2inv
coef <- vobj$coef
rm(vobj)
}
#
#   Segmentation
#
tobj <- awssegment(tobj,nlevels,swghts,rhoseq[k])
levels <- tobj$levels
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
lines(tobj$ind*max(tobj$bi)/nlevels,col=2)
title("Sum of weights and segmentation")
} 
if(d==2){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
zlim <- quantile(tobj$theta,c(0.001,0.999))
image(array(pmax(pmin(tobj$theta,zlim[2]),zlim[1]),dy),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$ind,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Segmentation (rho=",signif(rhoseq[k],3),")"))
}
if(d==3){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
zlim <- quantile(tobj$theta,c(0.001,0.999))
image(array(pmax(pmin(tobj$theta[,,n3%/%2+1],zlim[2]),zlim[1]),dy[-3]),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$ind[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Segmentation (rho=",signif(rhoseq[k],3),")"))
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
cat("Levels and size of segments:\n")
print(rbind(levels,as.vector(table(as.vector(tobj$ind)))))
if(demo) readline("Press return")
#
#   Prepare for next iteration
#
hakt <- hakt*hincr
x<-1.25^(k-1)
scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
lambda0<-lambda*lseq[k]*scorrfactor
cat(paste(signif(sum(hincr^(2*(1:k)))/sum(hincr^(2*(1:steps)))*100,2),"% ",sep=""))
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if aggkern="Uniform", or if qtau1=1 
###   
z<-list(theta=tobj$theta,segment=tobj$ind,ni=tobj$bi,y=y,levels=tobj$levels,
        hmax=hakt/hincr,mae=mae,lseq=c(0,lseq),call=args)
class(z)<-"aws.segment"
z
}
###########################################################################
#
#   Auxialiary functions
#
############################################################################
#
#   transformations for Gaussian case with variance modelling
#
############################################################################
awsgfamily <- function(y,scorr,lambda,cpar){
d <- cpar$d
h0 <- numeric(d)
if(scorr[1]>0) {
   if(length(scorr)<d) scorr <- c(scorr,rep(0,d))[1:d]
   for(i in 1:d) h0[i]<-geth.gauss(scorr[i])
   cat("Corresponding bandwiths for specified correlation:",h0,"\n")
}
        sigma2 <- IQRdiff(as.vector(y))^2
        if(scorr[1]>0) sigma2<-sigma2*Varcor.gauss(h0)
	cat("Estimated variance: ", signif(sigma2,4),"\n")
	sigma2 <- rep(sigma2, length(y))
	dim(sigma2) <- dim(y)
    lambda <- lambda*2 
    cpar$tau1 <- cpar$tau1*2 
    cpar$tau2 <- cpar$tau2*2 
    sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
list(cpar=cpar,lambda=lambda,y=y,sigma2=sigma2,h0=h0)
}
############################################################################
#
#  estimate inverse of variances
#
############################################################################
awsgsigma2 <- function(y,hobj,tobj,varmodel,coef,varprop,h0){
if(is.null(dy <- dim(y))) dy <- length(y)
if(is.null(coef)){
if(is.null(dy)) d <- 1 else d <- length(dy)
corfactor <- 1
for( i in 1:d ) corfactor <- corfactor*sum(dnorm((-5:5),0,h0[i]*0.59+1e-10))/dnorm(0,0,h0[i]*0.59+1e-10)
ind <- tobj$gi>1
residsq <- ((y-tobj$theta)[ind]*tobj$gi[ind]/(tobj$gi[ind]-pmin(.95*tobj$gi[ind],corfactor)))^2
theta <- tobj$theta[ind]
if(varmodel=="Quadratic") theta2 <- theta^2
wght <- tobj$gi[ind]-1
coef <- switch(varmodel,
               Constant=coefficients(lm(residsq~1,weights=wght^2)),
               Linear=coefficients(lm(residsq~theta,weights=wght^2)),
	       Quadratic=coefficients(lm(residsq~theta+theta2,weights=wght^2)))
} else {
if(length(coef)<switch(varmodel,Constant=1,Linear=2,Quadratic=3)) stop("To few coefficiants for specified variance model")
} 
gamma <- pmin(tobj$gi/hobj$bi,1)
theta <- gamma*tobj$theta+(1-gamma)*hobj$theta
#
#    use smoother estimates to obtain more stable variance estimates
#
sigma2 <- switch(varmodel,
               Constant=array(coef,dy),
               Linear=coef[1]+coef[2]*theta,
	       Quadratic=coef[1]+coef[2]*theta+coef[3]*theta^2)
#varquantile <- quantile(residsq,varprop)
varquantile <- varprop*mean(sigma2)
sigma2 <- pmax(sigma2,varquantile)
cat("Estimated mean variance",signif(mean(sigma2),3)," Variance parameters:",signif(coef,3),"\n")
list(sigma2inv=1/sigma2,coef=coef)
}
############################################################################
#
#    segmentation
#
############################################################################
awssegment <- function(tobj,nlevels,swghts,rho)
{
theta <- tobj$theta
dy <- dim(theta)
dim(theta) <- NULL
squant <- cumsum(c(0,swghts))/sum(swghts)
cuts <- as.vector(quantile(theta,squant))
ind<-tapply(theta,cut(theta,cuts))
ind[is.na(ind)] <- 1 # there seems to be a problem with the minimum ... .
levels <- numeric(nlevels)
for ( i in 1:nlevels) {
if(any(ind==i)) {
   levels[i] <- mean(theta[ind==i])
   theta[ind==i] <- rho * levels[i] + (1-rho) * theta[ind==i]
   }
}
dim(theta) <- dim(ind) <- dy
tobj$theta <- theta
tobj$ind <- ind
tobj$levels <- levels
tobj
}
