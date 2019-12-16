require(aws)
if(exists("X11")) X11(,12,4.5)
par(mfrow=c(1,3),mar=c(3,3,3,.25),mgp=c(2,1,0))
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
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=15, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 15 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 15
psize <- readline("Patchsize:\n Press 'Enter' for psize=2, otherwise provide value of psize:")
if(is.na(as.numeric(psize))) psize <- 2 else psize <- as.integer(pmin(4,pmax(0,psize)))
risk <- readline("Report risks (N/Y):")
if(risk %in% c("n","N")) u <- NULL else u <- im0
cat("Run aws \n")
if(k<25) {
yhat <- paws(y,hmax=hmax,graph=TRUE,u=u,patchsize=2)
} else {
yhat <- paws(y,hmax=hmax,graph=TRUE,sigma2=sigma^2,u=u,patchsize=psize)
}
readline("Press ENTER to show results")
oldpar <- par(mfrow=c(1,4),mar=c(3,3,3,.25),mgp=c(2,1,0))
image(im0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image (k=",k,")"))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
image(awsdata(yhat,"est"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction (patchsize=",2," hmax=",signif(yhat@hmax,3),")"))
readline("Enter for next plot:")
image(awsdata(yhat,"sd"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Standard deviation of estimates (min:",signif(min(awsdata(yhat,"sd")),3)," max:",signif(max(awsdata(yhat,"sd")),3),")"))
par(oldpar)
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){
rm(fbi,x,im0,y,sigma,hmax,yhat,u,risk)
dev.off()
}
