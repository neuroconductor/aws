aws.irreg <- function(y,x,d=2,hmax=NULL,hpre=NULL,qlambda=NULL,qtau=NULL,varmodel="Constant",
                sigma2=NULL,varprop=.1,graph=FALSE,
		lkern="Triangle",skern="Triangle",aggkern="Uniform",
		spmin=0,spmax=5,lseq=NULL,nbins=100,henv=NULL)
{
#
#    first check arguments and initialize
#
args <- match.call()
n<-length(y)
dx <- dim(x)
if(!(d %in% 1:2)) stop("this version is for 1D and 2D only")
if((d==1 && length(x)!=length(y))||(d==2 && (is.null(dx)||dx[1]!=n||dx[2]!=2))) stop("incorrect size of x")
if(!(varmodel %in% c("Constant","Linear","Quadratic"))) stop("Model for variance not implemented")
#
#   2D binning
#
require(sm)
zbins<-binning(x,y,nbins=n^(1/d)/2)
given.var<-!is.null(sigma2)
if(!given.var) {
sigma20 <- mean(zbins$devs[zbins$x.freq>1]/(zbins$x.freq[zbins$x.freq>1]-1))
cat("Preliminary variance estimate:",sigma20,"\n")
} else {
coef <- sigma2[1]
}
zbins<-binning(x,y,nbins=nbins)
ni <- t(zbins$table.freq)
mask <- ni>0 
if(!is.null(henv)) mask <- .Fortran("mask",
                                    as.logical(mask),
				    mask=as.logical(mask),
				    as.integer(nbins),
				    as.integer(switch(d,1,nbins)),
				    as.integer(max(0,henv)),
				    PACKAGE="aws",DUP=FALSE)$mask
yy <- rep(mean(y),length(mask))
dim(yy)<-dim(mask) <- dim(ni)
yy[ni>0] <- zbins$means
nn <- length(yy)
if(given.var) {
if(length(sigma2)!=nn) sigma2<-rep(sigma2[1],nn)
sigma2 <- 1/sigma2
} else {
sigma2 <- 1/rep(sigma20,nn)
}
if(d==2) wghts<-diff(range(x[,1]))/diff(range(x[,2])) else wghts<-1
if(d==2) dy<-dim(yy)<-dim(sigma2)<-c(nbins,nbins)
#
#   set appropriate defaults
#
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,
	            Gaussian=5,2)
skern <- switch(skern,"Exp"=1,"Triangle"=2,2)
cpar<-setawsdefaults(dim(yy),mean(y),"Gaussian",skern,aggkern,qlambda,qtau,lseq,hmax,1,spmax)
lambda <- 2*cpar$lambda # Gaussian case
cpar$tau1 <- cpar$tau1*2 
cpar$tau2 <- cpar$tau2*2 
hmax <- cpar$hmax
lseq <- cpar$lseq
shape <- cpar$shape
cpar$heta <- 20^(1/d)
hinit <- cpar$hinit
hincr <- cpar$hincr
spmax <- cpar$spmax
cpar$heta <- 1e10
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax*0.42445*4
    hinit <- 0.42445*4
    }
# now check which procedure is appropriate
##  this is the version on a grid
n1 <- nbins
n2 <- switch(d,1,nbins)
dy <- switch(d,NULL,c(n1,n2))
#
#    Initialize  for the iteration
#  
tobj<-list(bi= ni, bi2= ni^2, theta= yy/shape, fix=rep(FALSE,nn))
zobj<-list(ai=yy, bi0= rep(1,nn))
biold<-ni
vred<-ni
hakt <- hinit*hincr
hakt0 <- hinit*hincr
lambda0<-lambda
lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
#
#   produce a presmoothed estimate to stabilze variance estimates
#
if(is.null(hpre)) hpre<-(20*nn/n)^(1/d)
dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:d]
hobj <- .Fortran("cawsmask",as.double(yy),
                       as.logical(ni>0),# bins where we need estimates 
                       as.integer(ni),# contains number of points in bin
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       hakt=as.double(hpre),
                       as.double(1e40),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(nn),
                       bi0=as.double(zobj$bi0),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.integer(skern),
                       as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","ai")]
hobj$theta <- hobj$ai/hobj$bi
hobj$theta[ni==0]<-mean(hobj$theta[ni>0])
dim(hobj$theta) <- dim(hobj$bi) <- dy
#
#   iteratate until maximal bandwidth is reached
#
steps <- as.integer(log(hmax/hinit)/log(hincr))
cat("Progress:")
for(k in 1:steps){
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
# Correction for spatial correlation depends on h^{(k)} 
hakt0<-hakt
# heteroskedastic Gaussian case
zobj <- .Fortran("cgawsmas",as.double(yy),
                       as.logical(mask),# bins where we need estimates 
                       as.integer(ni),# contains number of points in bin
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
		       bi2=double(nn),
                       bi0=as.double(zobj$bi0),
		       vred=double(nn),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.integer(skern),
	               as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="aws",DUP=FALSE)[c("bi","bi0","bi2","vred","ai","hakt")]
vred[!tobj$fix]<-zobj$vred[!tobj$fix]
dim(zobj$ai)<-dy
if(hakt>n1/2) zobj$bi0 <- hincr^d*biold
biold <- zobj$bi0
tobj<-updtheta(zobj,tobj,cpar)
tobj$vred <- vred
tobj$theta[tobj$bi==0]<-mean(tobj$theta[ni>0])
dim(tobj$vred)<-dy
dim(tobj$theta)<-dy
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
if(graph){
#
#     Display intermediate results if graph == TRUE
#
if(d==1){ 
oldpar<-par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot((1:nn)[ni>0],yy[ni>0],ylim=range(yy,tobj$theta[mask]),col=3)
points((1:nn)[mask&ni==0],yy[mask&ni==0],col=4)
lines((1:nn)[mask],tobj$theta[mask],lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot((1:nn),tobj$bi,type="l",ylim=range(0,tobj$bi))
points((1:nn)[ni>0],max(tobj$bi)/max(ni)*ni[ni>0],col=3)
points((1:nn)[mask],rep(0,sum(mask)),col=4)
title("Sum of weights, ni and mask")
} 
if(d==2){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(yy,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(yy[mask]),3)," max=",signif(max(yy[mask]),3)))
zlim <- quantile(tobj$theta,c(0.001,0.999))
image(array(pmax(pmin(tobj$theta,zlim[2]),zlim[1]),dy),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta[mask]),3)," max=",signif(max(tobj$theta[mask]),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi[mask]),3)," mean=",signif(mean(tobj$bi[mask]),3)," max=",signif(max(tobj$bi),3)))
image(mask,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("mask")
}
par(oldpar)
}
#
#   Prepare for next iteration
#
#
#   Create new variance estimate
#
if(!given.var){
sigma2 <- awsisigma2(yy,hobj,tobj,ni,sigma20,varmodel,varprop)
}
hakt <- hakt*hincr
x<-1.25^(k-1)
lambda0<-lambda*lseq[k]
cat(paste(signif(sum(hincr^(2*(1:k)))/sum(hincr^(2*(1:steps)))*100,2),"% ",sep=""))
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if aggkern="Uniform", or if qtau1=1 
###   
if(length(sigma2)==nn){
# heteroskedastic case 
vartheta <- tobj$bi2/tobj$bi^2
} else {
# homoskedastic case 
vartheta <- sigma2*tobj$bi2/tobj$bi^2
vred<-tobj$bi2/tobj$bi^2
}
z<-list(theta=tobj$theta,sigma2=1/sigma2,bi=tobj$bi,var=vartheta,vred=vred,y=yy,ni=ni,varcoef=coef,
        hmax=hakt/hincr,lseq=c(0,lseq),call=args,zbins=zbins,x=x)
class(z)<-"aws.gaussian"
z
}
############################################################################
#
#  estimate inverse of variances
#
############################################################################
awsisigma2 <- function(y,hobj,tobj,ni,sigma20,varmodel,varprop){
if(is.null(dy <- dim(y))) dy <- length(y)
vredinv <- 1/tobj$vred
vredinv[is.na(vredinv)]<-0
vredinv[vredinv>1e10]<-0
ind <- vredinv>ni&ni>0
residsq <- pmax(1,ni[ind]-1)*((y-tobj$theta)[ind]*vredinv[ind]/(vredinv[ind]-ni[ind]))^2
theta <- tobj$theta[ind]
if(varmodel=="Quadratic") theta2 <- theta^2
wght <- (vredinv-ni)[ind]
coef <- switch(varmodel,
               Constant=coefficients(lm(residsq~1,weights=wght)),
               Linear=coefficients(lm(residsq~theta,weights=wght)),
	       Quadratic=coefficients(lm(residsq~theta+theta2,weights=wght)))
gamma <- pmin(vredinv/hobj$bi,1)
gamma[is.na(gamma)]<-0
theta <- gamma*tobj$theta+(1-gamma)*hobj$theta
#
#    use smoother estimates to obtain more stable variance estimates
#
eta <- sum(vredinv[ind]/ni[ind]-1)/(5*sum(ni[ni>0]-1)+sum(vredinv[ind]/ni[ind]-1))
#cat("eta",sum(vredinv[ind]-1),sum(ni[ni>0]-1),eta,"\n")
sigma2 <- switch(varmodel,
               Constant=array(coef,dy),
               Linear=coef[1]+coef[2]*theta,
	       Quadratic=coef[1]+coef[2]*theta+coef[3]*theta^2)
varquantile <- quantile(residsq,varprop)
sigma2 <- eta*pmax(sigma2,varquantile)+(1-eta)*sigma20
cat("Estimated mean variance",signif(mean(sigma2[ni>0]),3)," Variance parameters:",signif(coef,3),"\n")
1/sigma2
}
