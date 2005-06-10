#
#
#       General Gaussian Model - Heteroskedastic regression (unknown mu and sigma2)
#
#
cnorm<-function(y,hmax=100,model="full",qlambda=.966,qtau1=.92,lkern="Triangle",
        heta=NULL,eta0=0,hinit=NULL,hincr=NULL,graph=FALSE,u=NULL,wghts=NULL,eps=5.1e-6){
#KLnorm <- function(mu1,mu2,sig1,sig2){
#    log(sig2/sig1)-1+(sig1+(mu1-mu2)^2)/sig2
KLnorm <- function(mu1,mu2,sig1,sig2){
       log(sig1/sig2)-1+(sig2+(mu1-mu2)^2)/sig1
}
updtheta<-function(zobj,tobj,cpar){
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
eta<-(1-eta0)*pmin(1,bi0/tau*
   KLnorm((1-eta0)*munew+eta0*mu,mu,
          (1-eta0)*sigmanew+eta0*sigma,sigma))+eta0
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
if(is.null(qtau1)) qtau1<-.92
if(qtau1<1) tau1<-qchisq(qtau1,2) else tau1<-1e50
if(is.null(eta0)) eta0<-.25
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
if(qlambda<1) lambda <- 2*qchisq(qlambda,2) else lambda <- 1e50
if(is.null(hinit)||hinit<1.25) hinit <- 1.25
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(is.null(heta)) heta<-max(4,hinit+1)
cpar<-list(heta=heta,tau1=tau1,eta0=eta0,model=model,eps=eps)
#
#   now run aws
#
if(is.null(dy)) {
form="uni"
n<-length(y)
cpar$kstar<-log(100)
}
if(length(dy)==2) {
form="bi"
n1<-dy[1]
n2<-dy[2]
n<-n1*n2
cpar$kstar<-log(15)
if(is.null(wghts)) wghts<-c(1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2]/wghts[1])
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
if(is.null(wghts)) wghts<-c(1,1,1)
cpar$kstar<-log(5)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
hincr <- hincr^(1/3)
}
if(length(dy)>3) return("AWS for more than 3 dimensional grids is not implemented")
tobj<-list(bi= rep(1,n), mu=y, sigma= yyt, fix=rep(FALSE,n))
zobj<-list(asi=yyt, ami=y, bi0=rep(1,n))
bi0old<-rep(1,n)
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e50 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
if(form=="uni"){
while(hakt<=hmax){
zobj <- .Fortran("cawsnuni",
                       as.double(y),
                       as.double(yyt),
                       as.logical(tobj$fix),
                       as.integer(n),
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
		       PACKAGE="aws")[c("bi","bi0","ami","asi","hakt")]
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
#cat("Now update for h=",hakt,"\n")
tobj<-updtheta(zobj,tobj,cpar)
if(graph){
par(mfrow=c(1,3),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$mu),col=3)
lines(tobj$mu,lwd=2)
title(paste("Mean  h=",signif(hakt,3)))
plot(yyt-tobj$mu^2,ylim=range(yyt-tobj$mu^2,tobj$sigma),col=3)
lines(tobj$sigma,lwd=2)
title(paste("Variance  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    mean((tobj$mu-u)^2),"   MAE: ",mean(abs(tobj$mu-u)),"\n")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
}
if(form=="bi"){
while(hakt<=hmax){
zobj <- .Fortran("cawsnbi",
                       as.double(y),
                       as.double(yyt),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
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
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
dim(zobj$ami)<-c(n1,n2)
dim(zobj$asi)<-c(n1,n2)
dim(zobj$bi)<-c(n1,n2)
dim(zobj$bi0)<-c(n1,n2)
#cat("Now update for h=",hakt,"\n")
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$mu)<-c(n1,n2)
dim(tobj$sigma)<-c(n1,n2)
dim(tobj$bi)<-c(n1,n2)
dim(tobj$eta)<-c(n1,n2)
if(graph){
par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(tobj$mu,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(sqrt(tobj$sigma),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Estimated SD  h=",signif(hakt,3)," min=",signif(sqrt(min(tobj$sigma)),3)," max=",signif(sqrt(max(tobj$sigma)),3)))
image(sqrt(tobj$bi),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights  min=",signif(min(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
title("Sum of weights")
}
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    mean((tobj$mu-u)^2),"   MAE: ",mean(abs(tobj$mu-u)),"\n")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
}
if(form=="tri"){
while(hakt<=hmax){
zobj <- .Fortran("cawsntri",
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
if(hakt>n/2) zobj$bi0 <- hincr*biold
biold <- zobj$bi0
dim(zobj$ami)<-c(n1,n2,n3)
dim(zobj$asi)<-c(n1,n2,n3)
dim(zobj$bi)<-c(n1,n2,n3)
dim(zobj$bi0)<-c(n1,n2,n3)
#cat("Now update for h=",hakt,"\n")
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$mu)<-c(n1,n2,n3)
dim(tobj$sigma)<-c(n1,n2,n3)
dim(tobj$bi)<-c(n1,n2,n3)
dim(tobj$eta)<-c(n1,n2,n3)
if(graph){
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
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
    mean((tobj$mu-u)^2),"   MAE: ",mean(abs(tobj$mu-u)),"\n")
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
}
list(mu=tobj$mu,sigma=tobj$sigma,bi=tobj$bi,args=args)
}

