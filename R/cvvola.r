#
#
#       Multivariate Volatility estimation
#
#
cvvola<-function(y,hmax=100,model="full",qlambda=.966,qtau1=.92,lkern="Triangle",
             heta=NULL,eta0=0,hinit=dim(y)[1],hincr=NULL,graph=FALSE){
KLvola <- function(th1,th2){
    ds<-dim(th1)[1]
    n<-dim(th1)[2]
    d<-(sqrt(1+8*ds)-1)/2         
    .Fortran("klvolan",
             as.double(th1),
	     as.double(th2),
	     as.integer(d),
             as.integer(ds),
	     as.integer(n),
	     double(d*d),
	     double(d*d),
	     double(d*d),
	     kld=double(n),
	     PACKAGE="aws")$kld
}
updtheta<-function(zobj,tobj,cpar){
heta<-cpar$heta
eta0<-cpar$eta0
tau1<-cpar$tau1
kstar<-cpar$kstar
d<-cpar$d
hakt<-zobj$hakt
tau<-tau1*(2+max(kstar-log(hakt),0))
hakt<-zobj$hakt
bi0<-zobj$bi0
bi<-zobj$bi
yyt<-tobj$yyt
if(cpar$model=="full"){
thetanew<-t(t(zobj$ai)/bi)
} else {
thetanew<-t(t(zobj$ai)/bi)
ind<-(1:d)*((1:d)+1)/2
m<-1
for(i in (1:d)) for(j in 1:i){
if(i!=j) {
thetanew<-mean(yyt[m]/sqrt(thetanew[i,]*thetanew[j,]))*sqrt(thetanew[i,]*thetanew[j,])
}
}
}
theta<-tobj$theta
n<-dim(theta)[2]
thetanew[,tobj$fix]<-theta[,tobj$fix]
if(hakt>heta) {
eta<-(1-eta0)*pmin(1,bi0/tau*KLvola((1-eta0)*thetanew+eta0*theta,theta))+eta0
} else {
eta <- rep(eta0,n)
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
for( i in 1:dim(theta)[1]){
theta[i,] <- (1-eta)*thetanew[i,] + eta * theta[i,]
}
list(theta=theta,bi=bi,eta=eta,fix=(eta==1),yyt=yyt)
}
args <- match.call()
spmax <- 5
#
#    first construct sufficient statistics  Y*Y^T
#
d<-dim(y)[1]
n<-dim(y)[2]
ds<-d*(d+1)/2
yyt<-matrix(0,ds,n)
m<-1
for(i in 1:d) for(j in 1:i) {
yyt[m,]<-y[i,]*y[j,]
m<-m+1
}
#
#   Initialize parameters
#
if(is.null(qtau1)) qtau1<-.92
if(qtau1<1) tau1<-qchisq(qtau1,ds) else tau1<-1e50
if(is.null(eta0)) eta0<-.25
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
if(qlambda<1) lambda <- 2*qchisq(qlambda,ds) else lambda <- 1e50
if(is.null(hinit)||hinit<d) hinit <- d
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(is.null(heta)) heta<-max(2*ds,hinit+1)
cpar<-list(heta=heta,tau1=tau1,eta0=eta0,model=model,kstar=log(d*100),d=d)
#
#   now run aws
#
tobj<-list(yyt=yyt,bi= rep(1,n), theta= yyt, fix=rep(FALSE,n))
zobj<-list(ai=yyt, bi0=rep(1,n))
bi0old<-rep(1,n)
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e50 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
zobj <- .Fortran("cvawsvol",as.double(yyt),
                       as.logical(tobj$fix),
                       as.integer(n),
                       as.integer(d),
                       as.integer(ds),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(lkern),
		       as.double(spmax),
		       double(ds),
		       double(d*d),
		       double(d*d),
		       double(d*d),
		       PACKAGE="aws")[c("bi","bi0","ai","hakt")]
if(hakt>n/2) zobj$bi0 <- hincr*biold
dim(zobj$ai)<-c(ds,n)
biold <- zobj$bi0
cat("Now update for h=",hakt,"\n")
tobj<-updtheta(zobj,tobj,cpar)
if(graph){
par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(yyt[1,],ylim=range(yyt[1,],tobj$theta[1,]),col=3)
lines(tobj$theta[1,],lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
}
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
list(theta=tobj$theta,bi=tobj$bi,args=args)
}

