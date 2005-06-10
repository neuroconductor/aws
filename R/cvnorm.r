#
#
#       Multivariate Volatility estimation
#
#
cvnorm<-function(y,hmax=100,model="full",qlambda=.966,qtau1=.92,lkern="Triangle",
             heta=NULL,eta0=0,hinit=dim(y)[1],hincr=NULL,graph=FALSE){
KLnorm <- function(mu1,mu2,sig1,sig2){
    ds<-dim(sig1)[1]
    n<-dim(sig1)[2]
    d<-dim(mu1)[1]
    .Fortran("klmnormn",
             as.double(mu1),
             as.double(mu2),
             as.double(sig1),
             as.double(sig2),
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
n<-length(bi)
yyt<-tobj$yyt
sigma<-tobj$sigma
mu<-tobj$mu
munew<-t(t(zobj$ami)/bi)
mumut<-matrix(0,ds,n)
m<-1
for(i in 1:d)
for(j in 1:i){
mumut[m,]<-munew[i,]*munew[j,]
m<-m+1
}
sigmanew<-t(t(zobj$asi)/bi)-mumut
sigmanew<-t(t(sigmanew)*bi/(bi-d))
m<-1
for(i in 1:d) for(j in 1:i){
sigmanew[m,bi<(d+1)]<-sigma[m,bi<(d+1)]
m<-m+1
}
if(cpar$model!="full"){
sn<-sqrt(sigmanew[(1:d)*((1:d)+1)/2,])
m<-1
for(i in 1:d) for(j in 1:i){
if(i!=j) sigmanew[m,]<-mean((yyt[m,]-mumut[m,])/(sn[i,]*sn[j,]))*sn[i,]*sn[j,]
m<-m+1
}
}
sigmanew[,tobj$fix]<-sigma[,tobj$fix]
munew[,tobj$fix]<-mu[,tobj$fix]
if(hakt>heta) {
eta<-(1-eta0)*pmin(1,bi0/tau*
   KLnorm((1-eta0)*munew+eta0*mu,mu,
          (1-eta0)*sigmanew+eta0*sigma,sigma))+eta0
} else {
eta <- rep(eta0,n)
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
for( i in 1:d){
mu[i,] <- (1-eta)*munew[i,] + eta * mu[i,]
}
for( i in 1:dim(sigma)[1]){
sigma[i,] <- (1-eta)*sigmanew[i,] + eta * sigma[i,]
}
list(mu=mu,sigma=sigma,bi=bi,eta=eta,fix=(eta==1),yyt=yyt)
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
if(qtau1<1) tau1<-qchisq(qtau1,d+1) else tau1<-1e50
if(is.null(eta0)) eta0<-.25
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,2)
if(model=="full") {
if(qlambda<1) lambda <- 2*qchisq(qlambda,ds+d) else lambda <- 1e50
} else {
if(qlambda<1) lambda <- 2*qchisq(qlambda,2*d) else lambda <- 1e50
}
if(is.null(hinit)||hinit<d+1.5) hinit <- d+1.5
if(is.null(hincr)||hincr<=1) hincr <-1.25
if(is.null(heta)) heta<-max(2*(ds+d),hinit+1)
cpar<-list(heta=heta,tau1=tau1,eta0=eta0,model=model,kstar=log(d*100),d=d)
#
#   now run aws
#
tobj<-list(yyt=yyt,bi= rep(1,n), mu=y, sigma= yyt, fix=rep(FALSE,n))
zobj<-list(asi=yyt, ami=y, bi0=rep(1,n))
bi0old<-rep(1,n)
hakt <- hinit
lambda0<-lambda
if(hinit>1) lambda0<-1e50 # that removes the stochstic term for the first step
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
zobj <- .Fortran("cvawsnun",
                       as.double(y),
                       as.double(yyt),
                       as.logical(tobj$fix),
                       as.integer(n),
                       as.integer(d),
                       as.integer(ds),
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
                       double(ds),
                       double(d),
                       double(d*d),
                       double(d*d),
                       double(d*d),
                       PACKAGE="aws")[c("bi","bi0","ami","asi","hakt")]
if(hakt>n/2) zobj$bi0 <- hincr*biold
dim(zobj$ami)<-c(d,n)
dim(zobj$asi)<-c(ds,n)
biold <- zobj$bi0
cat("Now update for h=",hakt,"\n")
tobj<-updtheta(zobj,tobj,cpar)
if(graph){
par(mfrow=c(1,3),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y[1,],ylim=range(y[1,],tobj$mu[1,]),col=3)
lines(tobj$mu[1,],lwd=2)
title(paste("Mean  h=",signif(hakt,3)))
plot(yyt[1,]-tobj$mu[1,]^2,ylim=range(yyt[1,]-tobj$mu[1,]^2,tobj$sigma[1,]),col=3)
lines(tobj$sigma[1,],lwd=2)
title(paste("Variance  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
title("Sum of weights and eta")
}
hakt <- hakt*hincr
lambda0<-lambda
gc()
}
list(mu=tobj$mu,sigma=tobj$sigma,bi=tobj$bi,args=args)
}

