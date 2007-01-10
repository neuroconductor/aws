require(aws)
x<-seq(-1,1,length=256)
fbi<-function(x,y,k,r){
z1<-sqrt(x^2+y^2)
theta<-asin(x/z1)
z<-sin(k*theta)
z[z1<r]<-sin(pi*z1/r)[z1<r]
z<-sign((x+y)*(x-y))*z
(z-min(z))/(max(z)-min(z))
}
k <- as.integer(readline("Select number of waves: Enter for 7, provide positive integer otherwise "))
if(is.na(as.numeric(k))) k <- 7 else k <- max(1,as.integer(k))
im0<-outer(x,x,"fbi",k,.5)
image(im0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image (k=",k,")"))
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=.25, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- .25 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- .25
y <- im0+rnorm(im0,0,sigma)
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
degree <- readline("Degree of polynomial model:\n Press 'Enter' for degree =2, otherwise provide degree:")
if(is.na(as.numeric(degree))) degree <- 2 else degree <- as.numeric(degree)
if(!(degree %in% 0:2)) degree <- 2
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=15, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 15 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 15
qtau <- readline("Memory control(N/Y) :")
if(qtau %in% c("n","N")) qtau <- 1 else qtau <- NULL
risk <- readline("Report risks (N/Y):")
if(risk %in% c("y","Y")) u <-im0 else u <- NULL
cat("Run aws \n")
if(k<25) {
yhat <- lpaws(y,degree=degree,hmax=hmax,graph=TRUE,qtau=qtau,u=u)
} else {
yhat <- lpaws(y,degree=degree,hmax=hmax,graph=TRUE,qtau=qtau,sigma2=sigma^2,u=u)
}
readline("Press ENTER to show results")
oldpar <- par(mfrow=c(2,2),mar=c(1,1,2,.25),mgp=c(2,1,0))
image(im0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image (k=",k,")"))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
image(yhat$theta[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction (degree=",degree," hmax=",signif(yhat$hmax,3),")"))
image(yhat$ni,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights (min:",signif(min(yhat$ni),3)," max:",signif(max(yhat$ni),3),")"))
par(oldpar)
rm(fbi,x,im0,y,sigma,hmax,yhat,degree,u,risk,qtau)
