#
#    R - function  aws  for  Adaptive Weights Smoothing (AWS)
#    in regression models with additive sub-Gaussian errors               
#    local constant and local polynomial approach                         
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
#
#    default parameter settings
#
#    univariate   p=1     qlambda=.65    heta=8     tau1= 20    tau2= 0
#                 p=2     qlambda=.92    heta=18    tau1= 80    tau2= 40
#                 p=3     qlambda=.92    heta=32    tau1=400    tau2= 0
#
#
#    bivariate    p=1     qlambda=.65    heta=3     tau1= 4     tau2= 12
#                 p=2     qlambda=.92    heta=4     tau1=30     tau2= 50
#
cphaws <- function(y,x=NULL,p=1,sigma2=NULL,qlambda=NULL,heta=NULL,tau1=NULL,tau2=NULL,eta0=0,
                lkern="Triangle",hinit=NULL,hincr=NULL,hmax=NULL,NN=FALSE,
                u=NULL,graph=FALSE,demo=FALSE,wghts=NULL,spmin=0,spmax=5)
{ 
#
#          Auxilary functions
#
IQRdiff <- function(y) IQR(diff(y))/1.908
Pardist <- function(mcode,Bi0,dtheta){
#  local polynomial uni  mcode=1
#  local polynomial bi   MCODE=2
#  local linear multi    mcode=3    
   dp1<-dim(dtheta)[1]
   dp2<-dim(Bi0)[1]
   if(mcode==1){
      dist<-0
      for(i in 1:dp1) for(j in 1:dp1) dist<-dist+dtheta[i,]*Bi0[i+j-1,]*dtheta[j,]
   }
   if(mcode==2){
         ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
                dist<-0
                for(i in 1:dp1) for(j in 1:dp1) dist<-dist+dtheta[i,,]*Bi0[ind[i,j],,]*dtheta[j,,]
   }
   if(mcode==3){
                for(i in 1:dp1) dist<-dist+dtheta[i,]*Bi0[i+(i-1)*i/2,]*dtheta[i,]
		for(j in 2:dp1) for(i in 1:(j-1)) dist<-dist+2*dtheta[i,]*Bi0[i+(j-1)*j/2,]*dtheta[j,]
   }
   dist
   }
gettheta <- function(mcode,ai,bi){
if(mcode==1){
#  univariate 
   n<-dim(ai)[2]
   dp1<-dim(ai)[1]
   dp2<-dim(bi)[1]
   theta<-matrix(.Fortran("mpawsuni",
                  as.integer(n),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=double(dp1*n),
                  double(dp1*dp1),PACKAGE="aws")$theta,dp1,n)
}
if(mcode==2){
#  bivariate 
   n1<-dim(ai)[2]
   n2<-dim(ai)[3]
   n<-n1*n2
   dp1<-dim(ai)[1]
   dp2<-dim(bi)[1]
   ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
   theta <- array(.Fortran("mpawsbi",
                  as.integer(n),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=double(dp1*n),
                  double(dp1*dp1),
                  as.integer(ind),PACKAGE="aws")$theta,c(dp1,n1,n2))
 }
if(mcode==3){
#  multivariate
   n<-dim(ai)[2]
   dp1<-dim(ai)[1]
   dp2<-dim(bi)[1]
   theta <- matrix(.Fortran("mpawsmul",
                             as.integer(n),
                             as.integer(dp1),
                             as.integer(dp2),
                             as.double(ai),
                             as.double(bi),
                             double(dp1*n),
                             double(dp1*dp1),PACKAGE="aws")[[6]],dp1,n)
}
theta
}
updtheta<-function(zobj,tobj,cpar){
heta<-cpar$heta
eta0<-cpar$eta0
tau1<-cpar$tau1
tau2<-cpar$tau2
kstar<-cpar$kstar
hakt<-zobj$hakt
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
mcode<-cpar$mcode
bi0<-zobj$bi0
bi<-zobj$bi
thetanew<-gettheta(mcode,zobj$ai,bi)
theta<-tobj$theta
if(hakt>heta) {
	eta<-(1-eta0)*pmin(1,Pardist(mcode,bi0,thetanew-theta)/tau)+eta0
} else {
eta <- rep(eta0,prod(dim(bi)[-1]))
dim(eta)<-dim(bi)[-1]
}
eta[tobj$fix]<-1
dp1<-dim(zobj$ai)[1]
dp2<-dim(bi)[1]
etadd<-outer(rep(1,dp2),eta)
bi<-(1-etadd)*bi+etadd*tobj$bi
etadd<-outer(rep(1,dp1),eta)
theta <- (1-etadd)*thetanew + etadd * theta
#  local polynomial uni  mcode=1
#  local polynomial bi   mcode=2
#  local linear multi    mcode=3    
list(theta=theta,bi=bi,eta=eta,fix=(eta==1))
}
#
#          Main function body
#
#    first check arguments and initialize                                 
#
if(p==0) stop("use caws for local constant models")
if(p>3)  stop("no defaults for parameters available")
args <- match.call()
if(is.null(qlambda)) {
qlambda <- switch(p,.65,.92,.92)
}
if(qlambda<.6) warning("Inappropriate value of qlambda")
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
gridded <- is.null(x)
if(gridded){
#  this is the version on a grid
if(is.null(hinit)||hinit<=0) hinit <- p+.1
dy <- dim(y)
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n <- length(y)
   dp1 <- p+1
}
if(length(dy)==2){
   form <- "bi"
   ddim  <- 2
if(is.null(wghts)) wghts<-c(1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2]/wghts[1])
#  only use a wght for the second component
n1 <- dy[1]
n2 <- dy[2]
n <- n1*n2
if(p>2) stop("bivariate aws on a grid is not implemented for p>2")
dp1 <- switch(p+1,1,3,6)
}
if(length(dy)>2)
   stop("polynomial AWS for more than 2 dimensional grids is not implemented")
} else {
# not gridded
dx <- dim(x)
ddim <- 1
if(is.null(dx)) {
#
#    order data by order of x
#
    form <- "uni"
    n <- length(x)
    if(n!=length(y)) stop("incompatible lengths of x and y")
    ox <- order(x)
    x <- x[ox]
    y <- y[ox]
    dp1 <- p+1
}else {
   px <- dx[1]
   n <- dx[2]
   if(p>1) {
      p <- 1
      cat("p is set to 1, the maximal polynomial degree implemented")
      }
   form <- "multi"
   if(n!=length(y)) stop("incompatible dimensions of x and y")
   if(is.null(wghts)||length(wghts)!=px) wghts <- rep(1,px)
   dp1 <- 1+p*px
#
#  now generate matrix of nearest neighbors
#  hmax is interpreted as maximal number of neighbors
#
   if(NN){
   ihmax <- trunc(hmax)
   if(ihmax>n) ihmax <- n
   neighbors <- matrix(0,ihmax,n)
   for (i in 1:n) {
      if(px==1) adist <- (x-x[i])^2 else adist <- wghts^2%*%((x-x[,i])^2)
      neighbors[,i] <- order(adist)[1:ihmax]
      }
   } else {
   ihmax <- n
   ddim <- px
   neighbors <- distmat <- matrix(0,n,n)
   for (i in 1:n) {
      if(px==1) adist <- (x-x[i])^2 else adist <- wghts^2%*%((x-x[,i])^2)
      od <- order(adist)
      distmat[,i] <- adist[od]
      neighbors[,i] <- od
      }
#  now reduce memory used to whats needed                                 
   gc()
   distmat <- sqrt(distmat)
   maxdist <- apply(distmat,1,max)
   ihmax <- sum(maxdist<=hmax)
   distmat <- distmat[1:ihmax,]
   neighbors <- neighbors[1:ihmax,]
   maxdist <- maxdist[1:ihmax]
   gc()
   }
   }
   if(length(y)!=n) stop("incompatible dimensions of x and y")
   }
#
#   estimate variance in the gaussian case if necessary
#
    if(is.null(sigma2)) sigma2 <- IQRdiff(y)^2
    if(length(sigma2)==1) sigma2<-rep(sigma2,length(y))
    if(length(sigma2)!=length(y)) stop("incompatible length of sigma2")
    sigma2<-1/sigma2 #  taking the invers yields simpler formulaes 
#
#     now set hincr if not provided                               
#
if(is.null(hincr)) hincr <- 1.25^(1/ddim)  
#
#    now generate kernel on a grid                                        
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
#
#   get lambda as quantile of appropriate chisq, rescale to be consistent 
# with the paper and multiply by 2 to get 2*lambda in lamakt
#
#    univariate   p=1     qlambda=.65    heta=10    tau1= 20    tau2= 0
#                 p=2     qlambda=.92    heta=25    tau1= 80    tau2=40
#                 p=3     qlambda=.92    heta=100   tau1=400    tau2= 0
#
#
#    bivariate    p=1     qlambda=.65    heta=3     tau1= 4     tau2= 12
#                 p=2     qlambda=.92    heta=4     tau1=30     tau2= 50
if(is.null(dy)) {
if(is.null(heta)) heta<-switch(p,8,18,32) # 2*(p+1)^2
if(is.null(tau1)) tau1<-switch(p,20,80,400)
if(is.null(tau2)) tau2<-switch(p,0,40,0)
kstar<-log(switch(p,150,300,600))
} else {
if(is.null(heta)) heta<-switch(p,3,4)  # 
if(is.null(tau1)) tau1<-switch(p,4,12)
if(is.null(tau2)) tau2<-switch(p,30,50)
kstar<-log(switch(p,15,30))
}
if(qlambda>=1) lamakt<-1.e50 else lamakt <- 2*qchisq(qlambda,dp1)
cpar<-list(heta=heta,eta0=eta0,tau1=tau1,tau2=tau2,kstar=kstar)
#
#    now select the correct aws-procedure                                 
#
#   cases:    !gridded     uni   p>0                                     
#             gridded      bi    p=1,2                                    
#             !gridded     multi p=1                                    
#             !gridded     multi p=1  Nearest-Neighbor                  
#
      if( form=="uni" && p>0  ){
###
###                        uni     p>0
###
cpar$mcode<-1
if(gridded) x <- 1:length(y)
dp1 <- p+1
dp2 <- p+dp1
dxp <- max(diff(x,p+1))*(1+1.e-8)
if(is.null(hinit)||hinit<dxp) hinit <- dxp
#   generate binomial coefficients
cb <- matrix(0,dp1,dp1)
for(i in (1:dp1)) cb[i:dp1,i] <- choose((i:dp1)-1,i-1)
hakt <- hinit
tobj<-list(bi=matrix(0,dp2,n),theta=matrix(0,dp1,n),fix=rep(FALSE,n))
bi0old<-matrix(0,dp2,n)
zobj<-list(bi0=bi0old,ai=matrix(0,dp1,n))
lamakt0<-1.e50
while(hakt<=hmax){
zobj <- .Fortran("cphawsun",
              as.integer(n),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
	      as.logical(tobj$fix),
              as.double(sigma2),
              as.double(tobj$theta),
              bi=as.double(tobj$bi),
              as.double(zobj$bi0[1,]),
              bi0=as.double(zobj$bi0),
              ai=as.double(zobj$ai),
              as.double(lamakt0),
              hakt=as.double(hakt),
              as.integer(lkern),
              as.double(cb),
              double(dp1*dp1),
              double(dp1),
              double(dp2),
              double(dp1),
	      double(dp2),
	      as.double(spmin),
	      as.double(spmax),PACKAGE="aws")[c("ai","bi","bi0","hakt")]
dim(zobj$ai)<-c(dp1,n)
dim(zobj$bi)<-c(dp2,n)
dim(zobj$bi0)<-c(dp2,n)
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
lamakt0<-lamakt
if(graph){
par(mfrow=c(1,2),mar=c(3,3,3,1))
plot(x,y,ylim=range(y,tobj$theta[1,]),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,tobj$theta[1,],lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi[1,]/zobj$bi0[1,],type="l",ylim=range(0,1))
lines(tobj$eta,col=2)
title("Sum of rel.weights and eta")
}
if(!is.null(u)) 
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((tobj$theta[1,]-u)^2),
    "   MAE: ",mean(abs(tobj$theta[1,]-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
}
      if(gridded &&  form=="bi" && p>0){
###
###             gridded      bi    p=1,2
###
cpar$mcode<-2
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
tobj<-list(bi=array(0,c(dp2,n1,n2)),theta=array(0,c(dp1,n1,n2)),fix=matrix(FALSE,n1,n2))
biold<-array(0,c(dp2,n1,n2))
zobj<-list(ai=array(0,c(dp1,n1,n2)),bi0=biold)
lamakt0<-1.e50
while(hakt<=hmax){
zobj <- .Fortran("cphawsbi",
              as.integer(n1),
              as.integer(n2),
              as.integer(dp1),
              as.integer(dp2),
              as.double(y),
	      as.logical(tobj$fix),
              as.double(sigma2),
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
	      as.double(spmin),
	      as.double(spmax),PACKAGE="aws")[c("ai","bi","bi0","hakt")]
dim(zobj$ai)<-c(dp1,n1,n2)
dim(zobj$bi)<-c(dp2,n1,n2)
dim(zobj$bi0)<-c(dp2,n1,n2)
if(hakt>min(n1,n2)/2) zobj$bi0 <- hincr*hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
lamakt0<-lamakt
if(graph){
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta[1,,],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(tobj$bi[1,,]/zobj$bi0[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of rel. weights",
signif(min(tobj$bi[1,,]/zobj$bi0[1,,]),2),"-",
       signif(mean(tobj$bi[1,,]/zobj$bi0[1,,]),2),"-",
       signif(max(tobj$bi[1,,]/zobj$bi0[1,,]),2)))
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
      if( form=="multi" ){
###
###                        multi (nongridded)    p==0 or p==1
###
cpar$mcode<-3
dp1 <- 1+p*px
dp2 <- dp1*(dp1+1)/2
bi <- matrix(0,dp2,n)
theta <- ai <- matrix(0,dp1,n)
if(NN){
if(is.null(hinit)||hinit<(p+1)) hinit <- p+1
info <- 1
while(info>0){
ihinit <- trunc(hinit)
z <- .Fortran("ipawsmnn",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihinit,]),
              as.integer(ihinit),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              theta=as.double(theta),
              as.integer(lkern),
              double(dp1*dp1),
              double(dp1),
              info=as.integer(info),PACKAGE="aws")[c("ai","bi","theta","info")]
info <- z$info
hinit <- hinit+1
}
theta <- matrix(z$theta,dp1,n)
bi <- bi0 <- matrix(z$bi,dp2,n)
ai <- z$ai
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,]-u)^2),
                "   MAE: ",mean(abs(theta[1,]-u)),"\n")
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
ihakt <- min(ihmax,trunc(hakt))
z <- .Fortran("cpawsmnn",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihakt,]),
              as.integer(ihakt),
              as.double(theta),
              as.double(bi),
              bi=as.double(bi),
              bi0=as.double(bi0),
              as.double(ai),
              ai=as.double(ai),
              as.double(lamakt),
              as.double(hakt),
              as.integer(lkern),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1),
              double(dp1),
              double(dp1),
	      as.double(spmax),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- matrix((1-eta)*z$ai + eta * ai,dp1,n)
    bi <- matrix((1-eta)*z$bi + eta * bi,dp2,n)
    bi0 <- matrix((1-eta)*z$bi0 + eta * bi0,dp2,n)
    theta<-gettheta(mcode,ai,bi)
if(!is.null(u))
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,]-u)^2),
                 "   MAE: ",mean(abs(theta[1,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
} else {
dpd <- dp1+1
if(is.null(hinit)) hinit <- maxdist[dpd]
info <- 1
while(info>0){
if(hinit<=maxdist[dpd]) hinit <- maxdist[dpd]
ihinit <- sum(maxdist<=hinit)
z <- .Fortran("ipawsmul",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihinit,]),
              as.double(distmat[1:ihinit,]),
              as.integer(ihinit),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              theta=as.double(theta),
              as.integer(lkern),
              double(dp1*dp1),
              double(dp1),
              info=integer(1),PACKAGE="aws")[c("ai","bi","theta","info")]
info <- z$info
dpd <- dpd+1
}
theta <- matrix(z$theta,dp1,n)
bi <- bi0 <- matrix(z$bi,dp2,n)
ai <- z$ai
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,]-u)^2),
                "   MAE: ",mean(abs(theta[1,]-u)),"\n")
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
ihakt <- sum(maxdist<=hakt)
z <- .Fortran("cpawsmul",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihakt,]),
              as.double(distmat[1:ihakt,]),
              as.integer(ihakt),
              as.double(theta),
              as.double(bi),
              bi=as.double(bi),
              bi0=as.double(bi0),
              as.double(ai),
              ai=as.double(ai),
              as.double(lamakt),
              as.double(hakt),
              as.integer(lkern),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1),
              double(dp1),
              double(dp1),
	      as.double(spmax),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- matrix((1-eta)*z$ai + eta * ai,dp1,n)
    bi <- matrix((1-eta)*z$bi + eta * bi,dp2,n)
    bi0 <- matrix((1-eta)*z$bi0 + eta * bi0,dp2,n)
    theta <- gettheta(mcode,ai,bi)
if(!is.null(u))
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,]-u)^2),
                  "   MAE: ",mean(abs(theta[1,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
}
}
###
###            end cases
###
z<-list(theta=tobj$theta,y=y,x=x,call=args)
class(z)<-"aws"
z
}
