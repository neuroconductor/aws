awstestprop <- function(dy,hmax,theta=1,family="Gaussian",
                 lkern="Triangle",aws=TRUE,memory=FALSE,shape=2,ladjust=1,seed=1){
if(length(dy)>3) return("maximum array dimension is 3")
nnn <- prod(dy)
set.seed(seed)
par(mfrow=c(1,4),mar=c(3,3,3,1),mgp=c(2,1,0))
y <- array(switch(family,"Gaussian"=rnorm(nnn),
                         "Poisson"=rpois(nnn,theta),
                         "Exponential"=rexp(nnn,1),
                         "Bernoulli"=rbinom(nnn,1,theta),
                         "Volatility"=rnorm(nnn),
                         "Variance"=rchisq(nnn,shape)/shape),dy)
z <- 1:30
alpha <- exp(seq(-10,3.5,.25))
wghts <- switch(length(dy),c(0,0),c(1,0),c(1,1))
cpar<-setawsdefaults(dy,mean(y),family,lkern,"Uniform",aws,memory,ladjust,hmax,shape,wghts)
lambda <- cpar$lambda
hmax <- cpar$hmax
shape <- cpar$shape
d <- cpar$d
n<-length(y)
sigma2 <- 1
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
exceedence  <- matrix(0,30,kstar) # this is used to store exceedence probabilities 
pofalpha <- matrix(0,55,kstar) # this is used to store exceedence probabilities 
tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y/shape, fix=fix)
zobj<-zobj0<-list(ai=y, bi0= rep(1,n), bi=rep(1,n))
hhom <- rep(1,n)
biold<-rep(1,n)
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
zobj0 <- .Fortran("caws",as.double(y),
                       as.logical(fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(1e50),
                       as.double(tobj$theta),
                       bi=as.double(zobj0$bi),
		       bi2=double(n),
                       bi0=as.double(zobj0$bi0),
                       ai=as.double(zobj0$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
		       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt")]
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
plot(density(sqrt(KLdist0)))
#
#   get adaptive estimate
#
zobj <- .Fortran("caws",as.double(y),
                       as.logical(fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
		       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt")]
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
dim(tobj$fix)<-dy
lambda0 <- lambda
cat("lambda0",lambda0,"\n")
yhat <- tobj$theta
bi <- tobj$bi
yhat <- switch(length(dy),yhat[ind1],yhat[ind1,ind2],yhat[ind1,ind2,ind3])
KLdist1 <- switch(family,"Gaussian"=yhat^2/2,
                         "Poisson"=(theta-yhat+yhat*(log(yhat)-log(theta))),
                         "Exponential"=(log(yhat)-1+1/yhat),
                         "Bernoulli"=(yhat*log(yhat/theta)+
                                     (1-yhat)*log((1-yhat)/(1-theta))),
                         "Volatility"=(log(yhat)-1+1/yhat)/2,
                         "Variance"=shape/2*(log(yhat)-1+1/yhat))
lines(density(sqrt(KLdist1)),col=2)
bi <- switch(length(dy),bi[ind1],bi[ind1,ind2],bi[ind1,ind2,ind3])
plot(density(bi))
#lines(density(zobj$bi),col=2)
cat("quantiles diff bi",quantile(zobj0$bi-zobj$bi),"\n")
#plot(density(KLdist1/KLdist0,na.rm=TRUE))
for(i in 1:55) pofalpha[i,k] <- mean(KLdist1 > (1+alpha[i])*KLdist0)
if(k>1) contour(log(alpha),h[1:k],pofalpha[,1:k]^.2,n=30,ylab="h",xlab="ln(alpha)",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," p^(.2)"))
#KLdist <- switch(family,"Gaussian"=ni*yhat^2/2,
#                         "Poisson"=ni*(theta-yhat+yhat*(log(yhat)-log(theta))),
#                         "Exponential"=ni*(log(yhat)-1+1/yhat),
#                         "Bernoulli"=ni*(yhat*log(yhat/theta)+
#                                     (1-yhat)*log((1-yhat)/(1-theta))),
#                         "Volatility"=ni*(log(yhat)-1+1/yhat)/2,
#                         "Variance"=ni*shape/2*(log(yhat)-1+1/yhat))
for(i in 1:30) exceedence[i,k] <- mean(ni*KLdist1>z[i])
if(k>1) contour(z,h[1:k],exceedence[,1:k]^.2,n=30,ylab="h",xlab="z",
       main=paste(family,length(dy),"-dim. ladj=",ladjust," Exceed. Prob.^(.2)"))
if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
     }
tpar <- if(family%in%c("Bernoulli","Poisson"))  paste("theta=",theta) else  ""
cat(family,"(dim:",length(dy),tpar,")",signif(exceedence[,k],3),"\n")
cat("p",signif(pofalpha[,k],3),"\n")
cat("Quantile KLdist0",signif(quantile(KLdist0,c(0,.25,.5,.75,.9,.95,.99,.995,.999)),3),"\n")
cat("Quantile KLdist1",signif(quantile(KLdist1,c(0,.25,.5,.75,.9,.95,.99,.995,.999)),3),"\n")
k <- k+1
gc()
}
list(h=h,z=z,prob=exceedence,pofalpha=pofalpha)
}
