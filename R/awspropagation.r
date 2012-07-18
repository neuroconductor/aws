awstestprop <- function(dy,hmax,theta=1,family="Gaussian",
                 lkern="Triangle",aws=TRUE,memory=FALSE,shape=2,
                 homogeneous=TRUE,ladjust=1,seed=1){
if(length(dy)>3) return("maximum array dimension is 3")
nnn <- prod(dy)
set.seed(seed)
par(mfrow=c(1,2),mar=c(3,3,3,1),mgp=c(2,1,0))
y <- array(switch(family,"Gaussian"=rnorm(nnn),
                         "Poisson"=rpois(nnn,theta),
                         "Exponential"=rexp(nnn,1),
                         "Bernoulli"=rbinom(nnn,1,theta),
                         "Volatility"=rnorm(nnn),
                         "Variance"=rchisq(nnn,shape)/shape),dy)
z <- seq(0,30,.5)
alpha <- exp(seq(-10,3.5,.25))
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
k <- cpar$k
kstar <- cpar$kstar
h <- numeric(kstar)
if(k>1) h[1:(k-1)] <- 1+(0:(k-2))*.001
fix <- rep(FALSE,n)
exceedence  <- exceedencena  <- matrix(0,61,kstar) # this is used to store exceedence probabilities for adaptive and nonadaptive estimates
pofalpha <- matrix(0,55,kstar) # this is used to store exceedence probabilities 
zobj<-zobj0<-list(ai=y, bi0= rep(1,n), bi=rep(1,n),theta= y/shape)
hhom <- rep(1,n)
lambda0<-1e50
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
                       bi=as.double(zobj0$bi),
                       bi2=double(n),
                       bi0=as.double(zobj0$bi0),
                       double(n),#vred
                       ai=as.double(zobj0$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt")]
} else {
zobj0 <- .Fortran("caws1",as.double(y),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       bi=as.double(zobj0$bi),
		                 bi2=double(n),
                       bi0=as.double(zobj0$bi0),
                       ai=as.double(zobj0$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
		                 double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj0<-regularize(zobj0,family)
yhat0 <- zobj0$ai/zobj0$bi
dim(yhat0) <- dy
ih <- as.integer(hakt)
ind1 <- (ih+1):(dy[1]-ih)
if(length(dy)>1) ind2 <- (ih+1):(dy[2]-ih)
if(length(dy)>2) ind3 <- (ih+1):(dy[3]-ih)
yhat0 <- switch(length(dy),yhat0[ind1],yhat0[ind1,ind2],yhat0[ind1,ind2,ind3])
ni <- max(zobj0$bi0)
KLdist0 <- switch(family,"Gaussian"=yhat0^2/2,
                         "Poisson"=(theta-yhat0+yhat0*(log(yhat0)-log(theta))),
                         "Exponential"=(log(yhat0)-1+1/yhat0),
                         "Bernoulli"=(yhat0*log(yhat0/theta)+
                                     (1-yhat0)*log((1-yhat0)/(1-theta))),
                         "Volatility"=(log(yhat0)-1+1/yhat0)/2,
                         "Variance"=shape/2*(log(yhat0)-1+1/yhat0))
for(i in 1:61) exceedencena[i,k] <- mean(ni*KLdist0>z[i])
#
#   get adaptive estimate
#
if(!homogeneous&family=="Gaussian"){
zobj <- .Fortran("chaws",as.double(y),
                       as.logical(fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi=as.double(zobj$bi),
                       bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       double(n),#vred
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt")]
} else {
zobj <- .Fortran("caws",as.double(y),
                       as.logical(fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi=as.double(zobj$bi),
		                 bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
		                 double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt")]
}
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
biold <- zobj$bi0
zobj$theta <-zobj$ai/zobj$bi
dim(zobj$theta)<-dy
dim(zobj$bi)<-dy
lambda0 <- lambda
yhat <- zobj$theta
bi <- zobj$bi
yhat <- switch(length(dy),yhat[ind1],yhat[ind1,ind2],yhat[ind1,ind2,ind3])
KLdist1 <- switch(family,"Gaussian"=yhat^2/2,
                         "Poisson"=(theta-yhat+yhat*(log(yhat)-log(theta))),
                         "Exponential"=(log(yhat)-1+1/yhat),
                         "Bernoulli"=(yhat*log(yhat/theta)+
                                     (1-yhat)*log((1-yhat)/(1-theta))),
                         "Volatility"=(log(yhat)-1+1/yhat)/2,
                         "Variance"=shape/2*(log(yhat)-1+1/yhat))
bi <- switch(length(dy),bi[ind1],bi[ind1,ind2],bi[ind1,ind2,ind3])
for(i in 1:55) pofalpha[i,k] <- mean(KLdist1 > (1+alpha[i])*KLdist0)
if(k>1) contour(log(alpha),h[1:k],pofalpha[,1:k],levels=c(.5,.2,.1,.05,.02,.01,.005,
                .002,.001,.0005,.0002,.0001,.00005,.00002,.00001,.000005,.000002,.000001),
                       ylab="h",xlab="ln(alpha)",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," p"))
for(i in 1:61) exceedence[i,k] <- mean(ni*KLdist1>z[i])
if(k>1){
contour(z,h[1:k],exceedence[,1:k],levels=c(.5,.2,.1,.05,.02,.01,.005,
                .002,.001,.0005,.0002,.0001,.00005,.00002,.00001,.000005,.000002,.000001,.0000005,.0000002,.0000001),ylab="h",xlab="z",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," Exceed. Prob."))
contour(z,h[1:k],exceedencena[,1:k],levels=c(.5,.2,.1,.05,.02,.01,.005,
                .002,.001,.0005,.0002,.0001,.00005,.00002,.00001,.000005,.000002,.000001,.0000005,.0000002,.0000001),ylab="h",xlab="z",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," Exceed. Prob."),add=TRUE,col=2,lty=3)
       }
if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
     }
tpar <- if(family%in%c("Bernoulli","Poisson"))  paste("theta=",theta) else  ""
cat(family,"(dim:",length(dy),tpar,") ni=",ni,"e-prob:",signif(exceedence[4*(1:15)+1,k],3),"\n")
cat("p",signif(pofalpha[5*(1:11),k],3),"\n")
cat("Quantile KLdist0 (.5,.75,.9,.95,.99,.995,.999,1)",signif(quantile(KLdist0,c(.5,.75,.9,.95,.99,.995,.999,1)),3),"\n")
cat("Quantile KLdist1 (.5,.75,.9,.95,.99,.995,.999,1)",signif(quantile(KLdist1,c(.5,.75,.9,.95,.99,.995,.999,1)),3),"\n")
k <- k+1
gc()
}
list(h=h,z=z,prob=exceedence,probna=exceedencena,pofalpha=pofalpha)
}
