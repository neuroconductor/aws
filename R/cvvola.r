#
#
#       Multivariate Volatility estimation
#
#
vawsvola<-function(y,qlambda=NULL,qtau=NULL,model="full",lkern="Triangle",
             aggkern="Uniform",hinit=dim(y)[1],hincr=NULL,hmax=NULL,heta=NULL,
	     eta0=NULL,graph=FALSE){
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
updtheta<-function(zobj,tobj,cpar,aggkern){
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
thetanew<-t(t(zobj$ai)/bi)
ind<-(1:d)*((1:d)+1)/2
for(i in ind) {
msi<-mean(thetanew[i,])
thetanew[i,thetanew[i,]<1.e-5*msi]<-1.e-5*msi
}
if(cpar$model!="full"){
# just to avoid zero variances
sn<-sqrt(thetanew[ind,])
m<-1
## estimate and use global correlation
for(i in (1:d)) for(j in 1:i){
if(i!=j) {
thetanew[m,]<-mean(yyt[m,]/(sn[i,]*sn[j,]))*sn[i,]*sn[j,]
}
m<-m+1
}
}
theta<-tobj$theta
n<-dim(theta)[2]
thetanew[,tobj$fix]<-theta[,tobj$fix]
if(hakt>heta) {
eta<-switch(aggkern,"Uniform"=(1-eta0)*as.numeric(bi0/tau*
             KLvola((1-eta0)*thetanew+eta0*theta,theta)>1)+eta0,
            "Triangle"=(1-eta0)*pmin(1,bi0/tau*
	     KLvola((1-eta0)*thetanew+eta0*theta,theta))+eta0)
} else {
eta <- rep(0,n)
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
if(is.null(qlambda)) qlambda<-.98
if(is.null(qtau)) if(qlambda==1) qtau<-.85 else qtau<-.96
if(model=="full") d1<-(d+1)/2 else d1<-1
if(qtau<1) tau1<-qchisq(qtau,2) else tau1<-1e50
if(aggkern=="Triangle") tau1<-2.5*tau1
if(is.null(eta0)) eta0<-.0
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
if(qlambda<1) lambda <- 2*qchisq(qlambda,d1*d) else lambda <- 1e50
if(is.null(hinit)||hinit< max(10,d)) hinit <- max(10,d)
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(is.null(hmax)) hmax <- 100*d
if(is.null(heta)) heta<-max(2*(d+d1),hinit+1)
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
tobj<-updtheta(zobj,tobj,cpar,aggkern)
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

