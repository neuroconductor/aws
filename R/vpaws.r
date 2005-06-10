#
#    R - function  vpaws  for  Adaptive Weights Smoothing (AWS)
#    in multivariate regression models with additive sub-Gaussian errors           
#    vectorized local polynomial approach                         
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
vpaws <- function(y,p=1,sigma2=NULL,qlambda=NULL,heta=NULL,qtau=NULL,eta0=0,
                lkern="Triangle",hinit=NULL,hincr=NULL,hmax=NULL,
                u=NULL,graph=FALSE,wghts=NULL,vwghts=NULL)
{ 
#
#          Auxilary functions
#
IQRdiff <- function(y) IQR(diff(y))/1.908
Pardist <- function(mcode,Bi0,dtheta,vwghts){
#  local polynomial uni  mcode=1
#  local polynomial bi   MCODE=2
#  local linear multi    mcode=3    
   dd<-dim(dtheta)[2]
   dp1<-dim(dtheta)[1]
   dp2<-dim(Bi0)[1]
   if(mcode==1){
      dist<-0
      for(m in 1:dd) {
      dist0<-0
      for(i in 1:dp1) for(j in 1:dp1)
                dist0<-dist0+dtheta[i,m,]*Bi0[i+j-1,]*dtheta[j,m,]
      dist<-dist+dist0*vwghts[m]
      }
   }
   if(mcode==2){
         ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
         dist<-0
         for(m in 1:dd) {
            dist0<-0
            for(i in 1:dp1) for(j in 1:dp1)
		dist0<-dist0+dtheta[i,m,,]*Bi0[ind[i,j],,]*dtheta[j,m,,]
            dist<-dist+dist0*vwghts[m]
         }		
   }
   dist
   }
gettheta <- function(mcode,ai,bi){
if(mcode==1){
#  univariate 
   n<-dim(ai)[2]
   dp1<-dim(ai)[1]
   dp2<-dim(bi)[1]
   dd<-dim(ai)[2]
   theta<-ai
   theta<-array(.Fortran("vmpawsun",
                as.integer(n),
                as.integer(dd),
                as.integer(dp1),
                as.integer(dp2),
                as.double(ai),
                as.double(bi),
                theta=double(dp1*dd*n),
                double(dp1*dp1),PACKAGE="aws")$theta,c(dp1,dd,n))
}
if(mcode==2){
#  bivariate 
   n1<-dim(ai)[3]
   n2<-dim(ai)[4]
   n<-n1*n2
   dd<-dim(ai)[2]
   dp1<-dim(ai)[1]
   dp2<-dim(bi)[1]
   ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
   theta <- array(.Fortran("vmpawsbi",
                  as.integer(n),
		  as.integer(dd),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=double(dp1*dd*n),
                  double(dp1*dp1),
                  as.integer(ind),PACKAGE="aws")$theta,c(dp1,dd,n1,n2))
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
vw<-cpar$vw
dd<-cpar$dd
bi<-zobj$bi
thetanew<-gettheta(mcode,zobj$ai,bi)
theta<-tobj$theta
if(hakt>heta) {
	eta<-(1-eta0)*pmin(1,Pardist(mcode,zobj$bi0,thetanew-theta,vw)/tau)+eta0
} else {
eta <- rep(eta0,prod(dim(bi)[-1]))
}
dim(eta)<-dim(bi)[-1]
eta[tobj$fix]<-1
dp1<-dim(zobj$ai)[1]
dp2<-dim(bi)[1]
etadd<-outer(rep(1,dp2),eta)
bi<-(1-etadd)*bi+etadd*tobj$bi
etadd<-outer(matrix(1,dp1,dd),eta)
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
args <- match.call()
spmax <- 5
if(p==0) return("use vaws for local constant models")
if(p>3)  return("no defaults for parameters available")
if(is.null(qlambda)) {
qlambda <- switch(p,.65,.92,.92)
}
if(qlambda<.6) return("Inappropriate value of qlambda")
# now check which procedure is appropriate
#  this is the version on a grid
if(is.null(hinit)||hinit<=0) hinit <- p+1
dy <- dim(y)
if(length(dy)==2) {
   form <- "uni"
   ddim  <- 1
   dd<- dy[1]
   n <- dy[2]
   dp1 <- p+1
   dp2 <- p+dp1
   kstar<-log(switch(p,150,300,600))
}
if(length(dy)==3){
   form <- "bi"
   ddim  <- 2
if(is.null(wghts)) wghts<-c(1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2]/wghts[1])
#  only use a wght for the second component
dd<- dy[1]
n1 <- dy[2]
n2 <- dy[3]
n <- n1*n2
kstar<-log(switch(p,15,40))
if(p>2) return("bivariate aws on a grid is not implemented for p>2")
dp1 <- switch(p+1,1,3,6)
dp2 <- switch(p+1,1,6,15)
}
if(length(dy)>3)
   return("polynomial AWS for more than 2 dimensional grids is not implemented")
if(is.null(vwghts)) vwghts<-rep(1,dd)
if(sum(vwghts^2)!=dd) vwghts<-vwghts/sqrt(sum(vwghts^2))*dd
#
#   estimate variance in the gaussian case if necessary
#
    if(is.null(sigma2)) sigma2 <- IQRdiff(y)^2
    if(length(sigma2)==1) sigma2<-rep(sigma2,n)
    if(length(sigma2)!=n) return("incompatible length of sigma2")
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
if(is.null(heta)) {
if(length(dy)==2) heta<-switch(p,8,18,32)
if(length(dy)==3) heta<-switch(p,3,4)
}
if(is.null(qtau)) {
if(length(dy)==2) qtau<-.92 
if(length(dy)==3) qtau<-.92
}
if(qtau<1) tau1<-qchisq(qtau,dp1*dp1) else tau1<-1e50
if(qtau==1) heta<- 1e50  
tau2<-2*tau1
if(is.null(hmax)){
if(length(dy)==2) hmax<-250    # uses a maximum of about 500 points
if(length(dy)==3) hmax<-12   # uses a maximum of about 450 points
}
if(qlambda>=1) lamakt<-1.e50 else lamakt <- 2*qchisq(qlambda,dp1*dd)
cpar<-list(heta=heta,eta0=eta0,tau1=tau1,tau2=tau2,kstar=kstar,dd=dd,vw=vwghts)
#
#    now select the correct aws-procedure                                 
#
#   cases:        uni   p>0                                     
#                 bi    p=1,2                                    
      if( form=="uni" && p>0  ){
###
###                        uni     p>0
###
cpar$mcode<-1
x <- 1:length(y)
dxp <- max(diff(x,p+1))*(1+1.e-8)
if(is.null(hinit)||hinit<dxp) hinit <- dxp
#   generate binomial coefficients
cb <- matrix(0,dp1,dp1)
for(i in (1:dp1)) cb[i:dp1,i] <- choose((i:dp1)-1,i-1)
hakt <- hinit
tobj<-list(bi=matrix(0,dp2,n),theta=array(0,c(dp1,dd,n)),fix=rep(FALSE,n))
bi0old<-matrix(0,dp2,n)
zobj<-list(bi0=bi0old,ai=array(0,c(dp1,dd,n)))
lamakt0<-1.e50
while(hakt<=hmax){
zobj <- .Fortran("vphawsun",
              as.integer(n),
	      as.integer(dd),
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
              double(d*dp1),
              double(dp2),
              double(d*dp1),
	      double(dp2),
	      as.double(spmax),
	      as.double(vwghts),PACKAGE="aws")[c("ai","bi","bi0","hakt")]
dim(zobj$ai)<-c(dp1,dd,n)
dim(zobj$bi)<-c(dp2,n)
dim(zobj$bi0)<-c(dp2,n)
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
lamakt0<-lamakt
if(graph){
par(mfrow=c(1,2),mar=c(3,3,3,1))
plot(x,y[1,],ylim=range(y[1,],tobj$theta[1,1,]),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,tobj$theta[1,1,],lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi[1,]/zobj$bi0[1,],type="l",ylim=range(0,1))
lines(tobj$eta,col=2)
title("Sum of rel.weights and eta")
}
if(!is.null(u)) 
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((tobj$theta[1,,]-u)^2),
    "   MAE: ",mean(abs(tobj$theta[1,,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
}
      if(form=="bi" && p>0){
###
###                   bi    p=1,2
###
cpar$mcode<-2
ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
if(is.null(hinit)||hinit<p+1) hinit <- p+1
# now run aws-cycle                                                         
hakt <- hinit
tobj<-list(bi=array(0,c(dp2,n1,n2)),theta=array(0,c(dp1,dd,n1,n2)),fix=matrix(FALSE,n1,n2))
biold<-array(0,c(dp2,n1,n2))
zobj<-list(ai=array(0,c(dp1,dd,n1,n2)),bi0=biold)
lamakt0<-1.e50
while(hakt<=hmax){
zobj <- .Fortran("vphawsbi",
              as.integer(n1),
              as.integer(n2),
	      as.integer(dd),
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
              double(dd*dp1),
              as.integer(ind),
              as.double(wghts),
	      as.double(spmax),
	      as.double(vwghts),PACKAGE="aws")[c("ai","bi","bi0","hakt")]
dim(zobj$ai)<-c(dp1,dd,n1,n2)
dim(zobj$bi)<-c(dp2,n1,n2)
dim(zobj$bi0)<-c(dp2,n1,n2)
if(hakt>min(n1,n2)/2) zobj$bi0 <- hincr*hincr*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
lamakt0<-lamakt
if(graph){
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$theta[1,1,,],col=gray((0:255)/255),zlim=range(y[1,,]),xaxt="n",yaxt="n")
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
   cat("bandwidth: ",signif(hakt,3)," eta==1: ",sum(tobj$eta==1),"   MSE: ",mean((tobj$theta[1,,,]-u)^2),
       "   MAE: ",mean(abs(tobj$theta[1,,,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
}
###
###            end cases
###
z<-list(theta=tobj$theta,y=y,call=args)
class(z)<-"vpaws.gaussian"
z
}
