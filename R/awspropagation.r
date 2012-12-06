awstestprop <- function(dy,hmax,theta=1,family="Gaussian",
                 lkern="Triangle",aws=TRUE,memory=FALSE,shape=2,
                 homogeneous=TRUE,varadapt=FALSE,ladjust=1,seed=1,
                 minlevel=1e-6,maxz=25,diffz=.5){
if(length(dy)>3) return("maximum array dimension is 3")
nnn <- prod(dy)
if(minlevel < 5/nnn) {
minlevel <- 5/nnn
cat("minlevel reset to",minlevel,"due to insufficient size of test sample\n")
}
set.seed(seed)
par(mfrow=c(1,1),mar=c(3,3,3,3),mgp=c(2,1,0))
y <- array(switch(family,"Gaussian"=rnorm(nnn),
                         "Poisson"=rpois(nnn,theta),
                         "Exponential"=rexp(nnn,1),
                         "Bernoulli"=rbinom(nnn,1,theta),
                         "Volatility"=rnorm(nnn),
                         "Variance"=rchisq(nnn,shape)/shape),dy)
z <- seq(0,maxz,diffz)
elevel <- trunc(log(1e-6,10))
levels <- as.vector(outer(c(.5,.2,.1),10^(-0:elevel)))
levels <- levels[levels>=minlevel]
wghts <- switch(length(dy),c(0,0),c(1,0),c(1,1))
cpar<-setawsdefaults(dy,mean(y),family,lkern,"Uniform",aws,memory,ladjust,hmax,shape,wghts)
lambda <- cpar$lambda
hmax <- cpar$hmax
shape <- cpar$shape
d <- cpar$d
n<-length(y)
if(!homogeneous&family=="Gaussian"){
sigma2 <- array(rchisq(prod(dy),shape)/shape,dy)
} else sigma2 <- 1
zfamily <- awsfamily(family,y,sigma2,shape,0,lambda,cpar)
cpar <- zfamily$cpar
lambda <- zfamily$lambda
sigma2 <- zfamily$sigma2
h0 <- zfamily$h0
y <- zfamily$y
lkern <- cpar$lkern
rm(zfamily)
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
maxvol <- cpar$maxvol
# k <- cpar$k
k <- 1  
kstar <- cpar$kstar
h <- numeric(kstar)
if(k>1) h[1:(k-1)] <- 1+(0:(k-2))*.001
exceedence  <- exceedencena  <- matrix(0,length(z),kstar) # this is used to store exceedence probabilities for adaptive and nonadaptive estimates
zobj<-zobj0<-list(ai=y, bi=rep(1,n))
yhat <- y/shape
hhom <- rep(1,n)
lambda0<-1e50
cat("Progress:")
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
#
#  get initial conditions for a comparison
#
exceedence0 <- numeric(length(z))
if(family=="Bernoulli") y0 <- (10*y+1)/12
if(family=="Poisson") y0 <- y+.1
#
#  this corresponds to the regularization used to avoid Inf distances
#
KLdist0 <- switch(family,"Gaussian"=y^2/2,
                         "Poisson"=(theta-y0+y0*(log(y0)-log(theta))),
                         "Exponential"=(log(y)-1+1/y),
                         "Bernoulli"=(y0*log(y0/theta)+
                                     (1-y0)*log((1-y0)/(1-theta))),
                         "Volatility"=(log(y)-1+1/y)/2,
                         "Variance"=shape/2*(log(y)-1+1/y))
#KLdist0[is.na(KLdist0)] <- 1e8
for(i in 1:length(z)) exceedence0[i] <- mean(KLdist0>z[i])
#
#  now iterate
#
while (k<=kstar) {
      hakt0 <- gethani(1,1.25*hmax,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,1.25*hmax,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"hakt",hakt,"\n")
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hakt <- hakt*0.42445*4
    }
h[k] <- hakt
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
#
#   get nonadaptive estimate
#
if(!homogeneous&family=="Gaussian"){
zobj0 <- .Fortran("chaws1",as.double(y),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       bi=double(n),
                       double(n),
                       double(n),
                       double(n),#vred
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","ai","hakt")]
} else {
zobj0 <- .Fortran("caws1",as.double(y),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       bi=double(n),
                       double(n),
                       double(n),
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj0<-regularize(zobj0,family)
yhat0 <- zobj0$ai/zobj0$bi
dim(yhat0) <- dy
ih <- as.integer(hakt)
ind1 <- (ih+1):(dy[1]-ih)
if(length(dy)>1) ind2 <- (ih+1):(dy[2]-ih)
if(length(dy)>2) ind3 <- (ih+1):(dy[3]-ih)
yhat0 <- switch(length(dy),yhat0[ind1],yhat0[ind1,ind2],yhat0[ind1,ind2,ind3])
ni <- max(zobj0$bi)
KLdist0 <- switch(family,"Gaussian"=yhat0^2/2,
                         "Poisson"=(theta-yhat0+yhat0*(log(yhat0)-log(theta))),
                         "Exponential"=(log(yhat0)-1+1/yhat0),
                         "Bernoulli"=(yhat0*log(yhat0/theta)+
                                     (1-yhat0)*log((1-yhat0)/(1-theta))),
                         "Volatility"=(log(yhat0)-1+1/yhat0)/2,
                         "Variance"=shape/2*(log(yhat0)-1+1/yhat0))
for(i in 1:length(z)) exceedencena[i,k] <- mean(KLdist0>z[i]/ni)
#
#   get adaptive estimate
#
if(!homogeneous&family=="Gaussian"){
zobj <- .Fortran("chaws",as.double(y),
                       as.logical(rep(FALSE,n)),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(yhat),
                       bi=as.double(zobj$bi),
                       bi2=double(n),
                       double(n),
                       double(n),#vred
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi2","ai","hakt")]
} else {
zobj <- .Fortran("caws",as.double(y),
                       as.logical(rep(FALSE,n)),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(yhat),
                       bi=as.double(zobj$bi),
                       bi2=double(n),
                       double(n),
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi2","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
bi <- zobj$bi
yhat <-zobj$ai/bi
dim(yhat)<-dy
if(varadapt) bi <- bi^2/zobj$bi2
lambda0 <- lambda
yhat0 <- switch(length(dy),yhat[ind1],yhat[ind1,ind2],yhat[ind1,ind2,ind3])
KLdist1 <- switch(family,"Gaussian"=yhat0^2/2,
                         "Poisson"=(theta-yhat0+yhat0*(log(yhat0)-log(theta))),
                         "Exponential"=(log(yhat0)-1+1/yhat0),
                         "Bernoulli"=(yhat0*log(yhat0/theta)+
                                     (1-yhat0)*log((1-yhat0)/(1-theta))),
                         "Volatility"=(log(yhat0)-1+1/yhat0)/2,
                         "Variance"=shape/2*(log(yhat0)-1+1/yhat0))
for(i in 1:length(z)) exceedence[i,k] <- mean(KLdist1>z[i]/ni)
contour(z,0:k,cbind(exceedence0,exceedence[,1:k]),levels=levels,ylab="step",xlab="z",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," Exceed. Prob."))
contour(z,0:k,cbind(exceedence0,exceedencena[,1:k]),levels=levels,ylab="step",xlab="z",
       add=TRUE,col=2,lty=3)
       yaxp <- par("yaxp")
       at <- unique(as.integer(seq(yaxp[1],yaxp[2],length=yaxp[3]+1)))
       at <- at[at>0&at<=k]
       axis(4,at=at,labels=as.character(signif(h[at],3)))
       mtext("bandwidth",4,1.8)
if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
     }
tpar <- if(family%in%c("Bernoulli","Poisson"))  paste("theta=",theta) else  ""
cat(family,"(dim:",length(dy),tpar,") ni=",ni,"\n")
#cat(family,"(dim:",length(dy),tpar,") ni=",ni,"e-prob:",signif(exceedence[4*(1:15)+1,k],3),"\n")
cat("Quantile KLdist0 (.5,.75,.9,.95,.99,.995,.999,1)",signif(quantile(KLdist0,c(.5,.75,.9,.95,.99,.995,.999,1)),3),"\n")
cat("Quantile KLdist1 (.5,.75,.9,.95,.99,.995,.999,1)",signif(quantile(KLdist1,c(.5,.75,.9,.95,.99,.995,.999,1)),3),"\n")
k <- k+1
gc()
}
list(h=h,z=z,prob=exceedence,probna=exceedencena,levels=levels)
}
