#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2006 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#
#     default parameters:  see function setawsdefaults
#       
awscorr <- function(y,hmax=NULL,aws=TRUE,family="Gaussian",
                lkern="Triangle",homogen=TRUE,
                sigma2=NULL,shape=NULL,scorr=0,spmin=0.25,
                ladjust=1,wghts=NULL,u=NULL,graph=FALSE,demo=FALSE,
                testprop=FALSE,maxni=FALSE)
{
  #
  #   this version uses neighborhoods with an increase in potential 
  #   variance reduction by a factor of 1.25 from one iteration step 
  #   to the next
  #
  #    wghts is interpreted as voxel extensions ..., wghts for nonexisting dimensions are are set to INFTY
  #
  #    first check arguments and initialize
  #
  if(all(scorr==0)) stop("no spatial correlation specified, please use aws insted")
  if(length(scorr)<3) scorr <- c(scorr,rep(1e-5,3-length(scorr)))
  cwidth <- as.integer(log(5e-2)/log(scorr))
  acor <- corarray(scorr,cwidth)
  args <- match.call()
  dy<-dim(y)
  if(is.null(dy)) dy <- length(y)
  if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
  #
  #   set appropriate defaults
  #
  if(is.null(wghts)) wghts <- c(1,1,1)
  wghts <- switch(length(dy),c(0,0),c(wghts[1]/wghts[2],0),wghts[1]/wghts[2:3])
  if(family=="NCchi"){
    varstats <- sofmchi(shape/2) # precompute table of mean, sd and var for 
    #
    #   NCchi for noncentral chi with shape=degrees of freedom and theta =NCP
    #
  }
  cpar<-setawsdefaults(dy,mean(y),family,lkern,"Uniform",aws,FALSE,ladjust,hmax,shape,wghts)
  lambda <- cpar$lambda
  hmax <- cpar$hmax
  shape <- cpar$shape
  d <- cpar$d
  n<-length(y)
  # 
  #   family dependent transformations that depend on the value of family
  #
  zfamily <- awsfamily(family,y,sigma2,shape,scorr,lambda,cpar)
  cpar <- zfamily$cpar
  lambda <- zfamily$lambda
  sigma2 <- zfamily$sigma2
  h0 <- zfamily$h0
  y <- zfamily$y
  lkern <- cpar$lkern
  rm(zfamily)
  if(demo&& !graph) graph <- TRUE
  # now check which procedure is appropriate
  ##  this is the version on a grid
  n1 <- switch(d,n,dy[1],dy[1])
  n2 <- switch(d,1,dy[2],dy[2])
  n3 <- switch(d,1,1,dy[3])
  #
  #    Initialize  for the iteration
  #
  maxvol <- cpar$maxvol
  k <- cpar$k
  kstar <- cpar$kstar
  tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y/shape, fix=rep(FALSE,n))
  if(maxni) bi <- tobj$bi
  zobj<-list(ai=y, bi0= rep(1,n))
  hhom <- rep(1,n)
  if(family=="Gaussian"&length(sigma2)==n) vred<-rep(1,n)
  mae<-NULL
  lambda0<-1e50 # that removes the stochstic term for the first step, Initialization by kernel estimates
  if(testprop) {
    #
    #  prepare to  check for alpha in propagation condition (to adjust value of lambda using parameter ladjust)
    #
    if(is.null(u)) u <- 0
    cpar <- c(cpar, list(n1=n1,n2=n2,n3=n3,n=n1*n2*n3,family=family,u=u))
    propagation <- NULL
  } 
  #
  #   iteratate until maximal bandwidth is reached
  #
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
    dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
#aws::    if(family=="Gaussian"&scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
    # Correction for spatial correlation depends on h^{(k)} 
    if(family=="Gaussian"&length(sigma2)==n){
      # heteroskedastic Gaussian case
      zobj <- .Fortran("corchaws",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(acor),# local correlation array
                       as.integer(cwidth[1]+1),# dimensionality of 1 sector of local correlation array
                       as.integer(cwidth[2]+1),# dimensionality of 1 sector of local correlation array
                       as.integer(cwidth[3]+1),# dimensionality of 1 sector of local correlation array
                       integer(3*prod(dlw)),#array for relative indices of active points
                       double(prod(dlw)), # array for weights for active points
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi2=double(n),
                       bi0=double(n),
                       vred=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws")[c("bi","bi0","bi2","vred","ai","hakt")]
      vred[!tobj$fix]<-zobj$vred[!tobj$fix]
    } else {
      # all other cases
      if(cpar$mcode!=6){
        zobj <- .Fortran("corcaws",as.double(y),
                         as.logical(tobj$fix),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         hakt=as.double(hakt),
                         as.double(acor),
                         as.integer(cwidth[1]+1),
                         as.integer(cwidth[2]+1),# dimensionality of 1 sector of local correlation array
                         as.integer(cwidth[3]+1),# dimensionality of 1 sector of local correlation array
                         integer(3*prod(dlw)),#array for relative indices of active points
                         double(prod(dlw)), # array for weights for active points
                         hhom=as.double(hhom),
                         as.double(lambda0),
                         as.double(tobj$theta),
                         bi=as.double(tobj$bi2),
                         bi2=double(n),
                         bi0=double(n),
                         ai=as.double(zobj$ai),
                         as.integer(cpar$mcode),
                         as.integer(lkern),
                         as.double(spmin),
                         double(prod(dlw)),
                         as.double(wghts),
                         PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt","hhom")]
      } else {
        zobj <- .Fortran("caws6",as.double(y),
                         as.logical(tobj$fix),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         hakt=as.double(hakt),
                         hhom=as.double(hhom),
                         as.double(lambda0),
                         as.double(tobj$theta),
                         as.double(fncchiv(tobj$theta,varstats)/2),
                         bi=as.double(tobj$bi2),
                         bi2=double(n),
                         bi0=double(n),
                         ai=as.double(zobj$ai),
                         as.integer(cpar$mcode),
                         as.integer(lkern),
                         as.double(spmin),
                         double(prod(dlw)),
                         as.double(wghts),
                         PACKAGE="aws")[c("bi","bi0","bi2","ai","hakt","hhom")]
      }                     
    }
    if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
    dim(zobj$ai)<-dy
    tobj<-updtheta(zobj,tobj,cpar)
    dim(tobj$theta)<-dy
    if(maxni) bi <- tobj$bi2 <- pmax(bi,tobj$bi2)
    dim(tobj$bi)<-dy
    dim(tobj$eta)<-dy
    dim(tobj$fix)<-dy
    if(homogen) hhom <- zobj$hhom
    #
    #  if testprop == TRUE
    #  check alpha in propagation condition (to adjust value of lambda)
    #  
    if(testprop) propagation <- awstestprop(y,family,tobj,zobj,sigma2,hakt,cpar,u,propagation)
    if(graph){
      #
      #     Display intermediate results if graph == TRUE
      #
      if(d==1){ 
        oldpar<-par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
        plot(y,ylim=range(y,tobj$theta),col=3)
        if(!is.null(u)) lines(u,col=2)
        lines(tobj$theta,lwd=2)
        title(paste("Reconstruction  h=",signif(hakt,3)))
        plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
        lines(tobj$eta*max(tobj$bi),col=2)
        lines(hhom/max(hhom)*max(tobj$bi),col=3)
        title("Sum of weights, eta and hhom")
      } 
      if(d==2){ 
        oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
        image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
        image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
        image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
        image(tobj$fix,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
        title("Estimates fixed")
      }
      if(d==3){ 
        oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
        image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
        image(tobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
        image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
        title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
        image(tobj$fix[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
        title("Estimates fixed")
      } 
      par(oldpar)
    }
    #
    #    Calculate MAE and MSE if true parameters are given in u 
    #    this is for demonstration and testing for propagation (parameter adjustments) 
    #    only.
    #
    if(!is.null(u)) {
      cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
          signif(mean((tobj$theta-u)^2),3),"   MAE: ",
          signif(mean(abs(tobj$theta-u)),3)," mean(bi)=",
          signif(mean(tobj$bi),3),"mean hhom",signif(mean(hhom),3),"\n")
      mae<-c(mae,signif(mean(abs(tobj$theta-u)),3))
    }
    if(demo) readline("Press return")
    #
    #   Prepare for next iteration
    #
#    x<-1.25^k
#    scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
    lambda0<-lambda#*scorrfactor
    if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
    }
    k <- k+1
    gc()
  }
  cat("\n")
  ###                                                                       
  ###            end iterations now prepare results                                                  
  ###                                 
  if( family=="Gaussian"&length(sigma2)==n){
    # heteroskedastic Gaussian case 
    vartheta <- tobj$bi2/tobj$bi^2
    #  pointwise variances are reflected in weights
  } else {
    vartheta <- switch(family,Gaussian=sigma2,
                       Bernoulli=tobj$theta*(1-tobj$theta),
                       Poisson=tobj$theta,
                       Exponential=tobj$theta^2,
                       Volatility=2*tobj$theta,
                       Variance=2*tobj$theta,0)*tobj$bi2/tobj$bi^2
    vred<-tobj$bi2/tobj$bi^2
  }
  sigma2 <- switch(family,Gaussian=sigma2,
                   Bernoulli=tobj$theta*(1-tobj$theta),
                   Poisson=tobj$theta,
                   Exponential=tobj$theta^2,
                   Volatility=2*tobj$theta,
                   Variance=2*tobj$theta,0)
  if( family=="Gaussian"){
    vartheta<-vartheta/Spatialvar.gauss(hakt/0.42445/4,h0+1e-5,d)*Spatialvar.gauss(hakt/0.42445/4,1e-5,d)
  }
  awsobj(y,tobj$theta,vartheta,hakt,sigma2,lkern,lambda,ladjust,aws,FALSE,
         args,homogen,earlystop=FALSE,family=family,wghts=wghts,mae=mae,ni=tobj$bi)
}


corarray <- function(scorr,width=rep(3,3)){
  ## assumes autocorrelated (order 1) errors in all directions
  cfield <- array(1,width+1)
  cfield <- sweep(cfield,1,scorr[1]^(0:width[1]),"*")
  cfield <- sweep(cfield,2,scorr[2]^(0:width[2]),"*")
  sweep(cfield,3,scorr[3]^(0:width[3]),"*")
}