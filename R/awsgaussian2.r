aws.gaussian2 <- function (y, hmax=4, aws=TRUE, varmodel="Constant",
                      ladjust=1.0 , scorr=TRUE, lkern="Triangle", 
                      demo=FALSE, graph=FALSE, zext=1) {

  #
  #          Auxilary functions
  #
  IQRdiff <- function(data) IQR(diff(data))/1.908

  #####################################################################################
  ###    
  ###    function body
  ###    
  ###    first check arguments and initialize
  ###
  #####################################################################################    
  args <- match.call()
  bcf1 <- c(-.2408,1.9618,.5945,-1.356,-0.4299,0)
  bcf2 <- c(-.4724,1.1969,.22167,.48691,.13731,-1.0648)
  bcf3 <- c(-.4719,0.9642,.2291,0,.0,0)
  bcf <- cbind(bcf1,bcf2,bcf3)
# coefficients for  bias correction for spatial correlation
  if(!(toupper(varmodel) %in% c("CONSTANT","LINEAR"))) stop("varmodel incorrect")
  nvarpar <- switch(varmodel,Constant=1,Linear=2,1)
  #
  #   Check image type
  #
  dimg <- dimg0 <- dim(y)
  d <- max(1,length(dimg))
  n1 <- switch(d,length(y),dimg[1],dimg[1])
  n2 <- switch(d,1,dimg[2],dimg[2])
  n3 <- switch(d,1,1,dimg[3])
  n <- n1*n2*n3
  spcorr <- numeric(d)
  h0 <- rep(0,d)
  hpre <- switch(d,5,2.,2.)
  dlw<-(2*trunc(hpre)+1)
  #
  #     set approriate defaults
  #
  if (aws) qlambda <- .95 else qlambda <- 1
  ladjust <- max(1,ladjust)
  lseq <- pmax(1,c(1,1,1,1,1,1,1,1.6,1.5,1.4,1.3,1.2,1.1)/ladjust)
  lkern <- switch(lkern,
                  Triangle=2,
                  Quadratic=3,
                  Cubic=4,
                  Uniform=1,
                  2)
  if (is.null(hmax)) hmax <- 4
  if (qlambda<1) lambda <- ladjust*2*qchisq(qlambda,1) else lambda <- 1e50
  #
  #      in case of colored noise get the corresponding bandwidth (for Gaussian kernel)
  #
  sigma2 <- IQRdiff(y)^2
  cat("Estimated variance (assuming independence): ", signif(sigma2,4),"\n")
    lambda <- 2*lambda
  #     now set hinit and hincr if not provided
  if(aws) hinit <- 1 else {
    cat("No adaptation method specified. Calculate kernel estimate with bandwidth hmax.\n")
    hinit <- hmax
  }
  hincr <- 1.25^(1/d)
  if (demo && !graph) graph <- TRUE
  if(graph){
    oldpar <- if(d==1) par(mfrow=c(1,1),mar=c(3,3,3,1),mgp=c(2,1,0))
              else par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
  }
  # now check which procedure is appropriate
  #
  #    Initialize  list for theta
  #
  hmax <- hmax
  bi <- rep(1,n)
  theta <- y 
  bi0 <- 1
  #
  #  if varmodel specified prepare for initial variance estimation
  #
  coef <- matrix(0,nvarpar)
  coef[1,] <- sigma2
  vobj <- list(coef=coef,meanvar=sigma2)
  imgq995 <- quantile(y,.995)
  if(scorr){
    twohp1 <- 2*trunc(hpre)+1
    twohp3 <- if(n3==1) 1 else 2*trunc(hpre)+1
    pretheta <- .Fortran("gaws0",
                         as.double(y),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         as.double(hpre),
                         bi=as.double(bi),
                         as.double(bi0),
                         theta=double(n),
                         as.integer(lkern),
                         double(twohp1*twohp1*twohp3),# array for location weights
                         as.double(zext),
                         DUP=FALSE,
                         PACKAGE="aws")$theta
    dim(pretheta) <- dimg
    spchcorr <- .Fortran(switch(d,"estcorr1","estcorr2","estcorr3"),
                         as.double(y-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         scorr=double(d),
                         DUP=FALSE,
                         PACKAGE="aws")["scorr"]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hpre) 
    spcorr <- pmin(.9,spcorr+
                          bcf[1,d]/srh+bcf[2,d]/hpre+
                          bcf[3,d]*spcorr/srh+bcf[4,d]*spcorr/hpre+
                          bcf[5,d]*spcorr^2/srh+bcf[6,d]*spcorr^2/hpre)
    #  bias correction for spatial correlation
    cat("Estimated spatial correlation:",signif(spcorr,2),"\n")
  } else {
    spcorr <- rep(0,d)
  }
  #
  #         fix values of the image in inactiv pixel
  #
  ###
  ###              gridded   1D -- 3D
  ###
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  hakt0 <- hakt <- hinit
  if(aws) hakt <- hakt*hincr
  step <- 1
  lambda0 <- 1e50
  progress <- 0
  total <- (hincr^(2*ceiling(log(hmax/hinit)/log(hincr)))-1)/(hincr^2-1)
  if (total == 0) total <- hincr^(2*step) # for (hmax == hinit)
  #
  #   run single steps to display intermediate results
  #
  while (hakt<=hmax) {
    twohp1 <- 2*trunc(hakt)+1
    twohp3 <- if(n3==1) 1 else 2*trunc(hakt/zext)+1
    if(any(spcorr)>0) {
      h0<-numeric(length(spcorr))
      for(i in 1:length(h0))
        h0[i]<-geth.gauss(spcorr[i])
      if(length(h0)<d) h0<-rep(h0[1],d)
      # cat("Corresponding bandwiths for specified correlation:",h0,"\n")
    }
    if (any(spcorr>0)) {
      lcorr <- Spatialvar.gauss(hakt0,h0,d)/
        Spatialvar.gauss(h0,1e-5,d)/Spatialvar.gauss(hakt0,1e-5,d)
      # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)} 
      lambda0 <-lambda0*lcorr
      if(varmodel=="None") 
        lambda0 <- lambda0*Varcor.gauss(h0)
    } 
    hakt0 <- hakt
    zobj <- .Fortran("gvaws",
                      as.double(y),
                      as.integer(n1),
                      as.integer(n2),
                      as.integer(n3),
                      as.double(vobj$coef),
                      as.integer(nvarpar),
                      as.double(vobj$meanvar),
                      hakt=as.double(hakt),
                      as.double(lambda0),
                      as.double(theta),
                      bi=as.double(bi),
                      bi0=as.double(bi0),# just take a scalar here
                      theta=double(n),
                      as.integer(lkern),
                      as.double(spmin),		       
                      as.double(zext),
                      double(twohp1*twohp1*twohp3),# array for location weights
                      DUP=FALSE,
                      PACKAGE="aws")[c("bi","bi0","theta","hakt")]
    theta <- zobj$theta
    bi <- zobj$bi
    bi0 <- zobj$bi0
    rm(zobj)
    gc()
    dim(bi) <- dim(theta) <- dimg
    if (graph) {
       if(d==1){
          plot(1:n1,y,xlab="index",ylab="y")
          lines(1:n,theta,col=2)
          lines(1:n,theta+1.96/bi,lty=2,col=2)
          lines(1:n,theta-1.96/bi,lty=2,col=2)
          title("Data and fitted curve")
       }
       if(d==2) {
          image(y,col=grey((0:255)/255))
          title("Data")
          image(theta,col=grey((0:255)/255))
          title("Reconstruction")
          image(bi,col=grey((0:255)/255))
          title("ni")
       }
       if(d==3) {
          image(y[,,n3/2],col=grey((0:255)/255))
          title("Data")
          image(theta[,,n3/2],col=grey((0:255)/255))
          title("Reconstruction")
          image(bi[,,n3/2],col=grey((0:255)/255))
          title("ni")
       }
    }
    progress <- progress + hincr^(2*step)
    cat("Bandwidth",signif(hakt,3)," Progress",signif(progress/total,2)*100,"% \n")
    if (scorr) {
      #
      #   Estimate Correlations  (keep old estimates until hmax > hpre)
      #
      if(hakt > hpre){
      spchcorr <- .Fortran(switch(d,"estcorr1","estcorr2","estcorr3"),
                         as.double(y-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         scorr=double(d),,
                         DUP=FALSE,
                         PACKAGE="aws")["scorr"]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hakt) 
    spcorr <- pmin(.9,spcorr+
                         bcf[1,d]/srh+bcf[2,d]/hakt+
                         bcf[3,d]*spcorr/srh+bcf[4,d]*spcorr/hakt+
                         bcf[5,d]*spcorr^2/srh+bcf[6,d]*spcorr^2/hakt)
      #  bias correction for spatial correlation
      cat("Estimated spatial correlation:",signif(spcorr,2),"\n")
      } 
    } else {
      spcorr <- rep(0,d)
    }
      #
      #   Create new variance estimate
      #
      vobj <- .Fortran(switch(varmodel,Constant="esigmac",Linear="esigmal"),
                       as.double(y),
                       as.integer(n1*n2*n3),
                       as.double(theta),
                       as.double(bi),
                       as.double(imgq995),
                       coef=double(nvarpar),
                       meanvar=double(1),
                       DUP=FALSE,
                       PACKAGE="aws")[c("coef","meanvar")]
      cat("Estimated mean variance",signif(vobj$meanvar,3),"\n")
    step <- step + 1
    if (demo) readline("Press return")
    hakt <- hakt*hincr
    lambda0 <- lambda*lseq[step]
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 .....................................
  if(graph) par(oldpar)
  list(theta=theta,ni=bi,vcoef=vobj$coef,spcorr=spcorr,args=args)
}
