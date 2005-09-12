#
#
#       General Gaussian Model - Heteroskedastic regression (unknown mu and sigma2)
#
#
awsnorm<-function(y,qlambda=NULL,qtau=NULL,lkern="Triangle",
        aggkern="Uniform",hinit=NULL,hincr=NULL,hmax=NULL,heta=NULL,
	eta0=NULL,u=NULL,graph=FALSE,wghts=NULL){
#KLnorm <- function(mu1,mu2,sig1,sig2){
#    log(sig2/sig1)-1+(sig1+(mu1-mu2)^2)/sig2
KLnorm <- function(mu1,mu2,sig1,sig2){
       log(sig1/sig2)-1+(sig2+(mu1-mu2)^2)/sig1
}
updtheta<-function(zobj,tobj,cpar,aggkern){
heta<-cpar$heta
eta0<-cpar$eta0
tau1<-cpar$tau1
kstar<-cpar$kstar
hakt<-zobj$hakt
tau<-tau1*(2+max(kstar-log(hakt),0))
hakt<-zobj$hakt
bi0<-zobj$bi0
bi<-zobj$bi
n<-length(bi)
munew<-zobj$ami/bi
sigmanew<-zobj$asi/bi-munew*munew
sigmanew<-sigmanew*bi/(bi-1)+cpar$eps # this corrects for discretization error
sigma<-tobj$sigma
sigmanew[is.na(sigmanew)]<-sigma[is.na(sigmanew)]
mu<-tobj$mu
sigmanew[tobj$fix]<-sigma[tobj$fix]
munew[tobj$fix]<-mu[tobj$fix]
if(hakt>heta) {
eta<-switch(aggkern,"Uniform"=(1-eta0)*as.numeric(bi0/tau*
   KLnorm((1-eta0)*munew+eta0*mu,mu,
          (1-eta0)*sigmanew+eta0*sigma,sigma)>1)+eta0,
        "Triangle"=(1-eta0)*pmin(1,bi0/tau*
   KLnorm((1-eta0)*munew+eta0*mu,mu,
          (1-eta0)*sigmanew+eta0*sigma,sigma))+eta0)
} else {
eta <- rep(eta0,n)
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
mu <- (1-eta)*munew + eta * mu
sigma <- (1-eta)*sigmanew + eta * sigma
list(mu=mu,sigma=sigma,bi=bi,eta=eta,fix=(eta==1))
}
args <- match.call()
spmax <- 5
#
#    first construct sufficient statistics  Y*Y^T
#
dy<-dim(y)
yyt<-y*y
#
#   Initialize parameters
#
if(is.null(qlambda)) qlambda<-.96
if(is.null(qtau)) if(qlambda==1) {
if(is.null(dy)) qtau<-.8 else qtau<-.25 
} else qtau<-.92
d1<-2
if(qtau<1) tau1<-qchisq(qtau,d1) else tau1<-1e50
if(aggkern=="Triangle") tau1<-2.5*tau1
if(is.null(eta0)) eta0<-.0
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
if(qlambda<1) lambda <- 2*qchisq(qlambda,2) else lambda <- 1e50
#
#   now run aws
#
if(is.null(dy)) {
form="uni"
n1<-n<-length(y)
n2<-n3<-1
ddim<-1
kstar<-log(100)
}
if(length(dy)==2) {
form="bi"
n1<-dy[1]
n2<-dy[2]
n3<-1
n<-n1*n2
ddim<-2
kstar<-log(15)
hincr <- sqrt(hincr)
}
if(length(dy)==3) {
form="tri"
n1<-dy[1]
n2<-dy[2]
n3<-dy[3]
n<-n1*n2*n3
n1 <- dy[1]
n2 <- dy[2]
n3 <- dy[3]
n <- n1*n2*n3
ddim<-3
kstar<-log(5)
hincr <- hincr^(1/3)
}
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
if(is.null(wghts)) wghts<-c(1,1,1)
if(is.null(hinit)||hinit<1.2^(1/ddim)) hinit <- 2^(1/ddim)
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(is.null(heta)) heta<-max(4,hinit+1)
if(is.null(hmax)){
if(is.null(dim(y))) hmax<-250    # uses a maximum of about 500 points
if(length(dim(y))==2) hmax<-12   # uses a maximum of about 450 points
if(length(dim(y))==3) hmax<-5    # uses a maximum of about 520 points
}
cpar<-list(heta=heta,tau1=tau1,eta0=eta0,eps=1.e-6,kstar=kstar)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
tobj<-list(bi= rep(1,n), mu=y, sigma= yyt, fix=rep(FALSE,n))
zobj<-list(asi=yyt, ami=y, bi0=rep(1,n))
bi0old<-rep(1,n)
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e50 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
zobj <- .Fortran("cawsnorm",
                       as.double(y),
                       as.double(yyt),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$mu),
                       as.double(tobj$sigma),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
                       ami=as.double(zobj$ami),
                       asi=as.double(zobj$asi),
                       as.integer(lkern),
		       as.double(spmax),
		       as.double(wghts),
		       PACKAGE="aws")[c("bi","bi0","ami","asi","hakt")]
if(hakt>n1/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
dim(zobj$ami)<-dy
dim(zobj$asi)<-dy
dim(zobj$bi)<-dy
dim(zobj$bi0)<-dy
#cat("Now update for h=",hakt,"\n")
tobj<-updtheta(zobj,tobj,cpar,aggkern)
dim(tobj$mu)<-dy
dim(tobj$sigma)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
if(graph){
if(ddim==1){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
plot(y)
lines(tobj$mu,col=2)
title(paste("Observed data and estimated mean (h=",signif(hakt,3),")"))
plot(sqrt(tobj$sigma),col=2,type="l")
title(paste("Estimated standard deviation(h=",signif(hakt,3),")"))
plot(tobj$bi,type="l")
title(paste("Sum of weights  min=",signif(min(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
}
if(ddim==2){
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$mu,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(sqrt(tobj$sigma),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Estimated SD  h=",signif(hakt,3)," min=",signif(sqrt(min(tobj$sigma)),3)," max=",signif(sqrt(max(tobj$sigma)),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights  min=",signif(min(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
}
if(ddim==3){
n3h<-n3%/%2
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3h],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$mu[,,n3h],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(sqrt(tobj$sigma[,,n3h]),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Estimated SD  h=",signif(hakt,3)," min=",signif(sqrt(min(tobj$sigma[,,n3h])),3)," max=",signif(sqrt(max(tobj$sigma[,,n3h])),3)))
image(tobj$bi[,,n3h],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights  min=",signif(min(tobj$bi[,,n3h]),3)," max=",signif(max(tobj$bi[,,n3h]),3)))
}
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
    mean((tobj$mu-u)^2),"   MAE: ",mean(abs(tobj$mu-u)),"\n")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
list(mu=tobj$mu,sigma=tobj$sigma,bi=tobj$bi,args=args)
}

