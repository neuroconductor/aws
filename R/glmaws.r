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
glmaws <- function(y,xd=NULL,p=1,family="Gaussian",qlambda=NULL,heta=NULL,tau=NULL,KL=FALSE,
                lkern="Triangle",aggkern="Uniform",sigma2=NULL,hinit=NULL,hincr=NULL,hmax=NULL,hf=2,
                lseq=NULL,iter=50,u=NULL,graph=FALSE,demo=FALSE,wghts=NULL,spmax=5,eps=1.e-8,
		showwghts=FALSE,conf=FALSE,usevar=TRUE,smoothwghts=TRUE)
{ 
#
#          Auxilary functions
#
Pardist <- function(mcode,Bi0,dtheta){
#  local polynomial uni  mcode=1
#  local polynomial bi   mcode=2
#  local linear multi    mcode=3    
   dp1 <- dim(dtheta)[1]
   dp2 <- dim(Bi0)[1]
   if(mcode==1){
      dist <- 0
      for(i in 1:dp1) for(j in 1:dp1) dist <- dist+dtheta[i,]*Bi0[i+j-1,]*dtheta[j,]
   }
   if(mcode==2){
         ind <- matrix(c(1, 2, 3, 4, 5, 6,
                         2, 4, 5, 7, 8, 9,
                         3, 5, 6, 8, 9,10,
                         4, 7, 8,11,12,13,
                         5, 8, 9,12,13,14,
                         6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
                dist <- 0
                for(i in 1:dp1) for(j in 1:dp1) dist <- dist+dtheta[i,,]*Bi0[ind[i,j],,]*dtheta[j,,]
   }
   if(mcode==3){
                for(i in 1:dp1) dist <- dist+dtheta[i,]*Bi0[i+(i-1)*i/2,]*dtheta[i,]
		for(j in 2:dp1) for(i in 1:(j-1)) dist <- dist+2*dtheta[i,]*Bi0[i+(j-1)*j/2,]*dtheta[j,]
   }
   dist
   }
updtheta <- function(zobj,tobj,cpar,aggkern){
regdiff<-function(x1,x2){
dx<-dim(x1)
x<-x1-x2
#x<-sign(x)*pmax(abs(x)-1.e-3*abs(x2),0.0)
dim(x)<-dx
x
}
heta <- cpar$heta
tau1 <- cpar$tau1
tau2 <- cpar$tau2
kstar <- cpar$kstar
hakt <- zobj$hakt
tau <- 2*(tau1+tau2*max(kstar-log(hakt),0))
mcode <- cpar$mcode
bi0 <- zobj$bi0
bi <- zobj$bi
bi2 <- zobj$bi2
theta <- tobj$theta
thetanew <- zobj$theta
dd<-dim(theta)
if(hakt>heta) {
	eta <-switch(aggkern,"Uniform"=(1-eta0)*
         as.numeric(Pardist(mcode,bi0,regdiff(thetanew,theta))>tau/2.5)+eta0,
         "Triangle"=(1-eta0)*pmin(1,Pardist(mcode,bi0,regdiff(thetanew,theta))/tau)+eta0)
	if(length(dd)>2) dim(eta)<-dd[-1]
} else {
eta <- rep(eta0,prod(dim(bi)[-1]))
dim(eta) <- dim(bi)[-1]
}
eta[tobj$fix] <- 1
dp1 <- dim(zobj$theta)[1]
dp2 <- dim(bi)[1]
etadd<-outer(rep(1,dp1),eta)
theta <- (1-etadd)*thetanew + etadd * theta
etadd<-outer(rep(1,dp2),eta)
bi <- (1-etadd)*bi + etadd * tobj$bi 
bi2 <- (1-etadd)*bi2 + etadd * tobj$bi2 
list(theta=theta,bi=bi,bi2=bi2,eta=eta,fix=(eta==1))
}
#
#          Main function body
#
#    first check arguments and initialize                                 
#
if(p==0) return("use aws for local constant models")
if(p>3)  return("no defaults for parameters available")
mfamily<-switch(family,Gaussian=1,Poisson=2,Bernoulli=3,Exponential=4,1)
args <- match.call()
if(is.null(qlambda)) {
if(is.null(dim(y))) qlambda <- switch(p,.65,.966,.966) else qlambda <- switch(p,.65,.92) 
}
if(qlambda<.6) warning("Inappropriate value of qlambda")
if(is.null(dim(y))){
if(is.null(lseq)) lseq<-switch(family,Gaussian=1.3,Bernoulli=1,Exponential=1.7,Poisson=1,1.3)
}
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
dy <- dim(y)
#  this is the version on a grid
if(is.null(hinit)||hinit<=0) hinit <- p+1
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n <- length(y)
   dp1 <- p+1
   if(is.null(xd)) xd<-1
   x<-(1:n)*xd
}
if(length(dy)==2){
   form <- "bi"
   ddim  <- 2
if(is.null(wghts)) wghts <- c(1,1)
hinit <- hinit/wghts[1]
hmax <- hmax/wghts[1]
wghts <- (wghts[2]/wghts[1])
#  only use a wght for the second component
n1 <- dy[1]
n2 <- dy[2]
n <- n1*n2
if(p>2) return("bivariate aws on a grid is not implemented for p>2")
dp1 <- switch(p+1,1,3,6)
}
if(length(dy)>2)
   return("polynomial AWS for more than 2 dimensional grids is not implemented")
#
#     now set hincr, sigma2 if not provided                               
#
mae<-NULL
if(is.null(hincr)) hincr <- 1.25^(1/ddim)  
if(is.null(sigma2)){
if(family=="Gaussian"){
sigma2 <- IQRdiff(y)^2
cat("sigma^2=",sigma2,"\n")
sigma2<-sigma2 
} else sigma2<-1
}

#
#    now generate kernel on a grid                                        
#
lkern <- switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
#
#   get lambda as quantile of appropriate chisq, rescale to be consistent 
# with the paper and multiply by 2*sigma2 to get 2*sigma2*lambda in lamakt
#
#    univariate   p=1     qlambda=.65    heta=10    tau1= 100    tau2=  500
#                 p=2     qlambda=.92    heta=25    tau1= 500    tau2= 3000
#                 p=3     qlambda=.92    heta=100   tau1=4000    tau2=80000
#
#
#    bivariate    p=1     qlambda=.65    heta=3     tau1= 4     tau2= 12
#                 p=2     qlambda=.92    heta=4     tau1=30     tau2= 50
if(is.null(dy)) {
if(is.null(heta)) heta <- switch(p,10,25,100) # 2*(p+1)^2
if(is.null(tau)) tau <- switch(p,100,500,4000)
tau2 <- (p+1)^2*tau
kstar <- log(switch(p,150,300,600))
} else {
if(is.null(heta)) heta <- switch(p,3,4)  # 
tau <- switch(p,4,12)
tau2 <- (p+1)^2*tau
kstar <- log(switch(p,15,30))
}
if(qlambda>=1) lamakt <- 1.e50 else lamakt <- 2*qchisq(qlambda,dp1)*sigma2
if(is.null(eta0)) eta0<-switch(family,Gaussian=0,Poisson=.05,Bernoulli=.05,Exponential=.05)
cpar <- list(heta=heta,eta0=eta0,tau1=tau*sigma2,tau2=tau2*sigma2,kstar=kstar)
steps<-as.integer(log(hmax/hinit)/log(hincr)+1)
if(is.null(lseq)) lseq<-1
if(length(lseq)<steps) lseq<-c(lseq,rep(1,steps-length(lseq)))
lseq<-lseq[1:steps]
k<-1
#
#    now select the correct aws-procedure                                 
#
#   cases:       uni   p>0                                     
#                bi    p=1,2                                    
#
      if( form=="uni" && p>0  ){
###
###                        uni     p>0
###
confint<-NULL
cpar$mcode <- mcode <- 1
dp1 <- p+1
dp2 <- p+dp1
dp3 <- (dp1+1)*dp1/2
dxp <- max(diff(x,p))+.1
if(is.null(hinit)||hinit<dxp) hinit <- dxp
#   generate binomial coefficients
cb <- matrix(0,dp1,dp1)
for(i in (1:dp1)) cb[i:dp1,i] <- choose((i:dp1)-1,i-1)
hakt <- hinit
tobj <- list(bi=matrix(0,dp2,n),bi2=matrix(0,dp2,n),theta=matrix(0,dp1,n),fix=rep(FALSE,n))
theta <- matrix(0,dp1,n)
tobj$theta[1,]<- theta[1,] <- switch(mfamily,Gaussian=mean(y),
                                             Poisson=log(mean(y)),
                                             Bernoulli=log(mean(y)/(1-mean(y))),
                                             Exponential=1/mean(y))
bi0old <- matrix(0,dp2,n)
bii <- matrix(0,dp3,n)
zobj <- list(bi0=bi0old,bi02=bi0old,ni=rep(hinit,n))
lamakt0 <- 1.e50
h0<-matrix(.Fortran("poisdes",
              as.double(y),
              as.integer(n),
              as.integer(hf*p),
              h0=integer(2*n),PACKAGE="aws")[["h0"]],2,n)
#
#       this defines an interval where kernel weights are used (identifiability)
#
while(hakt<=hmax){
if(showwghts){
zobj <- .Fortran("lpawsunw",
              as.integer(n),
              as.integer(dp1),
              as.integer(dp2),
              as.integer(dp3),
              as.double(y),
              as.double(xd),
	      fix=as.logical(tobj$fix),
              as.integer(mfamily),
              as.double(tobj$theta),# this dixes the weights
              theta=as.double(tobj$theta),
              as.double(bii),
              bi=as.double(tobj$bi),
              bi2=as.double(tobj$bi2),
              bi0=as.double(zobj$bi0),
              bi02=as.double(zobj$bi02),
              ai=double(dp1), 
              ni=as.double(zobj$ni),
              as.double(lamakt0),
              hakt=as.double(hakt),
              as.integer(h0),
              as.integer(lkern),
              as.double(cb),
              double(dp1*dp1),# dmat
              double(dp1), # thij, d
              double(dp2), # psix
	      wghts=double(n*n), # wghts
              double(n), # wghts0
	      as.double(spmax),
              as.integer(iter),
              as.logical(smoothwghts),
              double(n), # work
              PACKAGE="aws")[c("fix","theta","bi","bi0","bi2","bi02","ni","hakt","wghts")]
} else {
zobj <- .Fortran("lpawsuni",
              as.integer(n),
              as.integer(dp1),
              as.integer(dp2),
              as.integer(dp3),
              as.double(y),
              as.double(xd),
	      fix=as.logical(tobj$fix),
              as.integer(mfamily),
              as.double(tobj$theta),# this dixes the weights
              theta=as.double(tobj$theta),
              as.double(bii),
              bi=as.double(tobj$bi),
              bi2=as.double(tobj$bi2),
              bi0=as.double(zobj$bi0),
              bi02=as.double(zobj$bi02),
              ai=double(dp1), 
              ni=as.double(zobj$ni),
              as.double(lamakt0),
              hakt=as.double(hakt),
              as.integer(h0),
              as.integer(lkern),
              as.double(cb),
              double(dp1*dp1),# dmat
              double(dp1), # thij, d
              double(dp2), # psix
	      double(n), # wghts
              double(n), # wghts0
	      as.double(spmax),
              as.integer(iter),
              as.logical(smoothwghts),
              double(n), # work
              PACKAGE="aws")[c("fix","theta","bi","bi0","bi2","bi02","ni","hakt")]
}
dim(zobj$theta) <- c(dp1,n)
dim(zobj$bi) <- c(dp2,n)
dim(zobj$bi0) <- c(dp2,n)
dim(zobj$bi2) <- c(dp2,n)
dim(zobj$bi02) <- c(dp2,n)
failed<-(1:n)[tobj$fix!=zobj$fix]
if(length(failed)>0) cat("No convergence in", failed, "\n")
tobj$fix[failed]<-TRUE
gc()
if(hakt>n/2) {
    zobj$bi0 <- hincr*biold
    zobj$bi02 <- hincr*biold2
}
biold <- zobj$bi0
biold2 <- zobj$bi02
tobj <- updtheta(zobj,tobj,cpar,aggkern)
if(usevar) {
bii<-matrix(.Fortran("bibi2ibi",
                     as.integer(n),
                     as.integer(dp1),
                     as.integer(dp2),
                     as.integer(dp3),
                     as.double(tobj$bi),
                     as.double(tobj$bi2),
                     bii=double(dp3*n),
                     double(dp1*dp1),
                     double(dp1*dp1),
                     PACKAGE="aws")$bii,dp3,n)
} else {
bii<-tobj$bi[c(1:3,3:5,4:7)[1:dp3],]
}
gc()
if(conf){
confint<-matrix(.Fortran("confuni",
                as.integer(n),
                as.integer(dp1),
                as.integer(dp2),
                as.double(tobj$theta[1,]),
                as.double(tobj$bi),
                as.double(tobj$bi2),
                conf=double(2*n),
                double(dp1*dp1),
                PACKAGE="aws")$conf,2,n)
}
if(graph){
if(showwghts) par(mfrow=c(1,3),mar=c(3,3,3,1)) else par(mfrow=c(1,2),mar=c(3,3,3,1))
plot(x,y,ylim=range(y),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,switch(mfamily,Gaussian=tobj$theta[1,],Poisson=exp(tobj$theta[1,]),
                       Bernoulli=exp(tobj$theta[1,])/(1+exp(tobj$theta[1,])),
                       Exponential=1/tobj$theta[1,]),lwd=2)
if(conf) {
    lines(x,switch(mfamily,Gaussian=confint[1,],Poisson=exp(confint[1,]),
                       Bernoulli=exp(confint[1,])/(1+exp(confint[1,])),
                       Exponential=1/confint[1,]),lwd=1,col=4)
    lines(x,switch(mfamily,Gaussian=confint[2,],Poisson=exp(confint[2,]),
                       Bernoulli=exp(confint[2,])/(1+exp(confint[2,])),
                       Exponential=1/confint[2,]),lwd=1,col=4)
}
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(x,zobj$ni,type="l",ylim=c(0,max(zobj$ni)))
lines(x,tobj$eta*max(zobj$ni),col=2)
title("Sum of weights and eta")
if(showwghts) image(matrix(zobj$wghts,n,n),zlim=c(0,1),col=gray((0:255)/255))
}
if(!is.null(u)) {
mtheta<-switch(mfamily,Gaussian=tobj$theta[1,],Poisson=exp(tobj$theta[1,]),
                       Bernoulli=exp(tobj$theta[1,])/(1+exp(tobj$theta[1,])),
                       Exponential=1/tobj$theta[1,])
cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",mean((mtheta-u)^2),
    "   MAE: ",mean(abs(mtheta-u)),"\n")
mae<-c(mae,signif(mean(abs(mtheta-u)),3))
}
if(demo) readline("Press return")
hakt <- hakt*hincr
lamakt0<-lamakt*lseq[k]
k<-k+1
gc()
}
}
      if(  form=="bi" && p>0){
###
###               bi    p=1,2
###
cpar$mcode <- 2
dp1 <- switch(p+1,1,3,6)
dp2 <- switch(p+1,1,6,15)
dpm <- 1
ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
if(is.null(hinit)||hinit<p+.1) hinit <- p+.1
if(demo) readline("Press return")
# now run aws-cycle                                                         
hakt <- hinit
tobj <- list(bi=array(0,c(dp2,n1,n2)),theta=array(0,c(dp1,n1,n2)),fix=matrix(FALSE,n1,n2))
biold <- array(0,c(dp2,n1,n2))
zobj <- list(ai=array(0,c(dp1,n1,n2)),bi0=biold)
lamakt0 <- 1.e50
while(hakt<=hmax){
zobj <- .Fortran("cpawsbi",
              as.integer(n1),
              as.integer(n2),
              as.integer(dp1),
              as.integer(dp2),
              as.double(y),
	      as.logical(tobj$fix),
              as.double(tobj$theta),
              bi=as.double(tobj$bi),
              bi0=as.double(zobj$bi0),
              ai=as.double(zobj$ai),
              as.double(lamakt0),
              hakt=as.double(hakt),
              as.integer(lkern),
              double(dp1*dp1),
              double(dp1),
              double(dp2),
              double(dp2),
              double(dp2),
              double(dp1),
              as.integer(ind),
              as.double(wghts),
	      as.double(spmax),PACKAGE="aws")[c("ai","bi","bi0","hakt")]
gc()
dim(zobj$ai) <- c(dp1,n1,n2)
dim(zobj$bi) <- c(dp2,n1,n2)
dim(zobj$bi0) <- c(dp2,n1,n2)
if(hakt>min(n1,n2)/2) zobj$bi0 <- hincr*hincr*biold
biold <- zobj$bi0
tobj <- updtheta(zobj,tobj,cpar)
gc()
lamakt0 <- lamakt
if(graph){
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta[1,,],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
image(matrix(tobj$eta,n1,n2),col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("eta")
}
if(!is.null(u)) 
   cat("bandwidth: ",signif(hakt,3)," eta==1: ",sum(tobj$eta==1),"   MSE: ",mean((tobj$theta[1,,]-u)^2),
       "   MAE: ",mean(abs(tobj$theta[1,,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
}
###
###            end cases
###
z <- list(theta=tobj$theta,confint=confint,y=y,x=x,mae=mae,lseq=c(0,lseq[-steps]),call=args)
class(z) <- "aws"
z
}

