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
paws <- function(y,hmax=NULL,aws=TRUE,family="Gaussian",
                lkern="Triangle",homogen=TRUE,aggkern="Uniform",
                sigma2=NULL,shape=NULL,scorr=0,spmin=0.25,
                ladjust=1,wghts=NULL,u=NULL,graph=FALSE,demo=FALSE,
                maxni=FALSE,patchsize=1)
{
#
#   patch based version (patches of patchsize neighbors in each direction)
#
#   this version uses neighborhoods with an increase in potential 
#   variance reduction by a factor of 1.25 from one iteration step 
#   to the next
#
#    wghts is interpreted as voxel extensions ..., wghts for nonexisting dimensions are are set to INFTY
#
#    first check arguments and initialize
#
args <- match.call()
memory <- FALSE
dy<-dim(y)
if(is.null(dy)) dy <- length(y)
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
#
#   set appropriate defaults
#
if(is.null(wghts)) wghts <- c(1,1,1)
wghts <- switch(length(dy),c(0,0),c(wghts[1]/wghts[2],0),wghts[1]/wghts[2:3])
if(family=="NCchi"){
varstats <- sofmchi(shape/2) # precompute table of mean, sd and var for 
#
#   NCchi for noncentral chi with shape=degrees of freedom and theta =NCP
#
}
cpar<-setawsdefaults(dy,mean(y),family,lkern,aggkern,aws,memory,ladjust,hmax,shape,wghts)
lambda <- cpar$lambda
hmax <- cpar$hmax
shape <- cpar$shape
d <- cpar$d
n<-length(y)
mc.cores <- setCores(,reprt=FALSE)

# 
#   family dependent transformations that depend on the value of family
#
zfamily <- awsfamily(family,y,sigma2,shape,scorr,lambda,cpar)
cpar <- zfamily$cpar
lambda <- zfamily$lambda
sigma2 <- zfamily$sigma2
h0 <- zfamily$h0
y <- zfamily$y
lkern <- cpar$lkern
rm(zfamily)
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
##  this is the version on a grid
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
#
#    Initialize  for the iteration
#
maxvol <- cpar$maxvol
k <- cpar$k
kstar <- cpar$kstar
tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y/shape, fix=rep(FALSE,n))
if(maxni) bi <- tobj$bi
zobj<-list(ai=y, bi0= rep(1,n))
hhom <- rep(1,n)
if(family=="Gaussian"&length(sigma2)==n) vred<-rep(1,n)
mae<-NULL
lambda0<-1e50 # that removes the stochstic term for the first step, Initialization by kernel estimates
#
#   iteratate until maximal bandwidth is reached
#
cat("Progress:")
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
while (k<=kstar) {
      hakt0 <- gethani(1,1.25*hmax,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,1.25*hmax,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"hakt",hakt,"\n")
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hakt <- hakt*0.42445*4
    }
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
if(family=="Gaussian"&scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
# Correction for spatial correlation depends on h^{(k)} 
np1 <- 2*patchsize+1
np2 <- if(n2>1) 2*patchsize+1 else 1
np3 <- if(n3>1) 2*patchsize+1 else 1
# all other cases
if(cpar$mcode!=6){
   zobj <- .Fortran("pcaws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi2=double(n),
                       bi0=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3),
                       double(np1*np2*np3*mc.cores),
                       double(np1*np2*np3*mc.cores),
                       PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt","hhom")]
   } else {
     stop("Non-central chi model not implemented")
   }                     

if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$theta)<-dy
if(maxni) bi <- tobj$bi <- pmax(bi,tobj$bi)
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
dim(tobj$fix)<-dy
if(homogen) hhom <- zobj$hhom
#
#  
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
lines(hhom/max(hhom)*max(tobj$bi),col=3)
title("Sum of weights, eta and hhom")
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
title("Estimates fixed")
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
title("Estimates fixed")
} 
par(oldpar)
}
#
#    Calculate MAE and MSE if true parameters are given in u 
#    this is for demonstration and testing for propagation (parameter adjustments) 
#    only.
#
if(!is.null(u)) {
   maxI <- max(1,diff(range(u)))
   mse <- mean((tobj$theta-u)^2)
   psnr <- 20*log(maxI,10)-10*log(mse,10)
   cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    signif(mse,3),"   MAE: ",
		    signif(mean(abs(tobj$theta-u)),3),
		    "PSNR: ",signif(psnr,3), " mean(bi)=",
		    signif(mean(tobj$bi),3),"mean hhom",signif(mean(hhom),3),"\n")
   mae<-c(mae,signif(mean(abs(tobj$theta-u)),3))
		    }
if(demo) readline("Press return")
#
#   Prepare for next iteration
#
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
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if aggkern="Uniform", or if !memory 
###   

vartheta <- switch(family,Gaussian=sigma2,
                          Bernoulli=tobj$theta*(1-tobj$theta),
                          Poisson=tobj$theta,
                          Exponential=tobj$theta^2,
                          Volatility=2*tobj$theta,
                          Variance=2*tobj$theta,0)*tobj$bi2/tobj$bi^2
vred<-tobj$bi2/tobj$bi^2
sigma2 <- switch(family,Gaussian=sigma2,
                        Bernoulli=tobj$theta*(1-tobj$theta),
                        Poisson=tobj$theta,
                        Exponential=tobj$theta^2,
                        Volatility=2*tobj$theta,
                        Variance=2*tobj$theta,0)
if( family=="Gaussian"){
vartheta<-vartheta/Spatialvar.gauss(hakt/0.42445/4,h0+1e-5,d)*Spatialvar.gauss(hakt/0.42445/4,1e-5,d)
}
awsobj(y,tobj$theta,vartheta,hakt,sigma2,lkern,lambda,ladjust,aws,memory,
   args,homogen,earlystop=FALSE,family=family,wghts=wghts,mae=mae,ni=tobj$bi)
}
