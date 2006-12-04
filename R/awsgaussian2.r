aws.gaussian2 <- function (y, hmax=4, aws=TRUE, varmodel="Constant",
                      ladjust=1.0 , mask = NULL, scorr=TRUE, lkern="Triangle", 
                      skern="Triangle", demo=FALSE, graph=FALSE) {

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
  if(!is.null(mask)) {
    varmodel <- "None"
    scorr <- FALSE
    clip <- FALSE
  }
  bcf <- c(-.4724,1.1969,.22167,.48691,.13731,-1.0648)
# coefficients for  bias correction for spatial correlation
  estvar <- toupper(varmodel) %in% c("CONSTANT","LINEAR")
  hpre <- 2.
  if(estvar) {
    dlw<-(2*trunc(hpre)+1)
    nvarpar <- switch(varmodel,Constant=1,Linear=2,1)
  }
  #
  #   Check image type
  #
  dimg <- dimg0 <- dim(y)
  d <- length(dimg)
  n1 <- dimg[1]
  n2 <- switch(d,1,dimg[2],dimg[2])
  n3 <- switch(d,1,1,dimg[3])
  spcorr <- numeric(d)
  h0 <- rep(0,d)
  #
  #     set approriate defaults
  #
  if (aws) qlambda <- .9 else qlambda <- 1
  ladjust <- max(1,ladjust)
  lseq <- pmax(1,c(1,1,1,1,1,1,1,1.6,1.5,1.4,1.3,1.2,1.1)/ladjust)
  lkern <- switch(lkern,
                  Triangle=2,
                  Quadratic=3,
                  Cubic=4,
                  Uniform=1,
                  2)
  skern <- switch(skern,
                  Triangle=1,
                  Exp=2,
                  1)
  if (is.null(hmax)) hmax <- 4
  wghts <- wghts/sum(wghts)
  dgf <- sum(wghts)^2/sum(wghts^2)
  if (qlambda<1) lambda <- ladjust*2*qchisq(qlambda,dgf)/dgf else lambda <- 1e50
  #
  #      in case of colored noise get the corresponding bandwidth (for Gaussian kernel)
  #
  sigma2 <- IQRdiff(y)^2)
  cat("Estimated variance (assuming independence): ", signif(sigma2,4),"\n")
  if (!estvar) {
    if (length(sigma2)==1) {
      #   homoskedastic Gaussian case
      lambda <- lambda*sigma2 
    } 
  }
  if (skern==1) {
    # set the support of the statistical kernel to (0,2) for skern==1, set spmin and spmax
    lambda <- 2*lambda
    spmin <- 0
    spmax <- 1
  } else {
    #   set cut off point in K_{st}(x) = exp(-x) I_{x<spmax}
    spmin <- 0
    spmax <- 3  
  }
  #     now set hinit and hincr if not provided
  if(aws) hinit <- 1 else {
    cat("No adaptation method specified. Calculate kernel estimate with bandwidth hmax.\n")
    hinit <- hmax
  }
  hincr <- sqrt(1.25)
  if (demo && !graph) graph <- TRUE
  if(graph){
    oldpar <- par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
    on.exit(par(oldpar))
    graphobj0 <- object[-(1:length(object))[names(object)=="img"]]
    graphobj0$dim <- c(n1,n2)
  }
  # now check which procedure is appropriate
  #
  #    Initialize  list for theta
  #
  hmax <- hmax
  bi <- rep(1,n)
  theta <- switch(imgtype,greyscale=object$img[xind,yind],rgb=object$img[xind,yind,])
  bi0 <- 1
  #
  #  if varmodel specified prepare for initial variance estimation
  #
  if(estvar){
    coef <- matrix(0,nvarpar)
    coef[1,] <- sigma2
    vobj <- list(coef=coef,meanvar=sigma2)
    imgq995 <- quantile(y,.995)
  }
  if(scorr){
    twohp1 <- 2*trunc(hpre)+1
    pretheta <- .Fortran("gaws0",
                         as.double(y),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         as.double(hpre),
                         as.double(theta),
                         bi=as.double(bi),
                         as.double(bi0),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         double(twohp1*twohp1),# array for location weights
                         as.double(wghts),DUP=FALSE,
                         PACKAGE="aws")$theta
    dim(pretheta) <- dimg
    spchcorr <- .Fortran("estcorr",
                         as.double(switch(imgtype,
                                          greyscale=object$img[xind,yind],
                                          rgb=object$img[xind,yind,])-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         scorr=double(2*dv),
                         chcorr=double(max(1,dv*(dv-1)/2)),
                         as.double(hpre),
                         DUP=FALSE,
                         PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
#    spcorr <- matrix(pmin(.9,0.8817*spcorr+0.231/hpre+6.018*spcorr/hpre^2+
#                            1.753*spcorr^2/hpre-10.622*spcorr^2/hpre^2),2,dv)
    srh <- sqrt(hpre) 
    spcorr <- matrix(pmin(.9,spcorr+
                          bcf[1]/srh+bcf[2]/hpre+
                          bcf[3]*spcorr/srh+bcf[4]*spcorr/hpre+
                          bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hpre),2,dv)
    #  bias correction for spatial correlation
    chcorr <- spchcorr$chcorr
    for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
    if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
        " (1,3):",signif(chcorr[2],2),
        " (2,3):",signif(chcorr[3],2),"\n")
  } else {
    spcorr <- matrix(0,2,dv)
    chcorr <- numeric(max(1,dv*(dv-1)/2))
  }
  #
  #         fix values of the image in inactiv pixel
  #
  if(!is.null(mask)) fix <- !mask[xind,yind]
  ###
  ###              gridded   2D
  ###
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  hakt0 <- hakt <- hinit
  if(estvar) { 
    if(aws) hakt <- hakt*hincr
    step <- 1
  } else step <- 0 
  lambda0 <- 1e50
  progress <- 0
  total <- (hincr^(2*ceiling(log(hmax/hinit)/log(hincr)))-1)/(hincr^2-1)
  if (total == 0) total <- hincr^(2*step) # for (hmax == hinit)
  #
  #   run single steps to display intermediate results
  #
  while (hakt<=hmax) {
    twohp1 <- 2*trunc(hakt)+1
    spcorr <- pmax(apply(spcorr,1,mean),0)
    if(any(spcorr)>0) {
      h0<-numeric(length(spcorr))
      for(i in 1:length(h0))
        h0[i]<-geth.gauss(spcorr[i])
      if(length(h0)<2) h0<-rep(h0[1],2)
      # cat("Corresponding bandwiths for specified correlation:",h0,"\n")
    }
    if (any(spcorr>0)) {
      lcorr <- Spatialvar.gauss(hakt0,h0)/
        Spatialvar.gauss(h0,1e-5)/Spatialvar.gauss(hakt0,1e-5)
      # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)} 
      lambda0 <-lambda0*lcorr
      if(varmodel=="None") 
        lambda0 <- lambda0*Varcor.gauss(h0)
    } 
    hakt0 <- hakt
    if(is.null(mask)){
      if(estvar){
        varcase <- 1
        zobj <- .Fortran("awsvimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(vobj$coef),
                         as.integer(nvarpar),
                         as.double(vobj$meanvar),
                         as.double(chcorr),
                         hakt=as.double(hakt),
                         as.double(lambda0),
                         as.integer(theta),
                         bi=as.double(bi),
                         bi0=as.double(bi0),# just take a scalar here
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(skern),
                         as.double(spmin),		       
                         as.double(spmax),
                         as.double(sqrt(wghts)),
                         double(twohp1*twohp1),# array for location weights
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","bi0","theta","hakt")]
      } else {
        # all other cases
        varcase <- 3
        zobj <- .Fortran("awsimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         hakt=as.double(hakt),
                         as.double(lambda0),
                         as.integer(theta),
                         bi=as.double(bi),
                         bi0=as.double(bi0),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(skern),
                         as.double(spmin),
                         as.double(spmax),
                         double(twohp1*twohp1),# array for location weights
                         as.double(wghts),
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","bi0","theta","hakt")]
      }
    } else {
      # all other cases
      varcase <- 3
      zobj <- .Fortran("mawsimg",
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.logical(fix),
		       as.logical(mask[xind,yind]),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(dv),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.integer(theta),
                       bi=as.double(bi),
                       bi0=as.double(bi0),
                       theta=integer(prod(dimg)),
                       as.integer(lkern),
		       as.integer(skern),
                       as.double(spmin),
                       as.double(spmax),
                       double(twohp1*twohp1),# array for location weights
                       as.double(wghts),
                       double(dv),DUP=FALSE,
                       PACKAGE="adimpro")[c("bi","bi0","theta","hakt")]
    }
    theta <- zobj$theta
    bi <- zobj$bi
    bi0 <- zobj$bi0
    rm(zobj)
    gc()
    dim(bi) <- dimg[1:2]
    if (graph) {
      graphobj <- graphobj0
      class(graphobj) <- "adimpro"
      graphobj$img <- switch(imgtype,
                             greyscale=object$img[xind,yind],
                             rgb=object$img[xind,yind,])
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title("Observed Image")
      graphobj$img <- array(as.integer(theta),dimg)
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt,3)))
      graphobj$img <- matrix(as.integer(65534*bi/bi0),n1,n2)
      graphobj$type <- "greyscale"
      graphobj$gamma <- FALSE
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Adaptation (rel. weights):",signif(mean(bi)/bi0,3)))
      rm(graphobj)
      gc()
    }
    progress <- progress + hincr^(2*step)
    cat("Bandwidth",signif(hakt,3)," Progress",signif(progress/total,2)*100,"% \n")
    if (scorr) {
      #
      #   Estimate Correlations  (keep old estimates until hmax > hpre)
      #
      if(hakt > hpre){
      spchcorr <- .Fortran("estcorr",
                           as.double(switch(imgtype,
                                            greyscale=object$img[xind,yind],
                                            rgb=object$img[xind,yind,]) - theta),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(dv),
                           scorr=double(2*dv),
                           chcorr=double(max(1,dv*(dv-1)/2)),
                           DUP=FALSE,
                           PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
#    spcorr <- matrix(pmin(.9,0.8817*spcorr+0.231/hakt+6.018*spcorr/hakt^2+
#                            1.753*spcorr^2/hakt-10.622*spcorr^2/hakt^2),2,dv)
    srh <- sqrt(hakt) 
    spcorr <- matrix(pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hakt+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hakt+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hakt),2,dv)
      #  bias correction for spatial correlation
      chcorr <- spchcorr$chcorr
      for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
      if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
                   " (1,3):",signif(chcorr[2],2),
                   " (2,3):",signif(chcorr[3],2),"\n")
      } else {
         spcorr <- matrix(pmin(.9,.06+0.9*spchcorr$scorr+1.108/hpre),2,dv)
      #  bias correction for spatial correlation
         chcorr <- spchcorr$chcorr
      }
    } else {
      spcorr <- matrix(0,2,dv)
      chcorr <- numeric(max(1,dv*(dv-1)/2))
    }
    if (estvar) {
      #
      #   Create new variance estimate
      #
      vobj <- .Fortran(switch(varmodel,Constant="esigmac",Linear="esigmal"),
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.integer(n1*n2),
                       as.integer(dv),
                       as.integer(theta),
                       as.double(bi),
                       as.integer(imgq995),
                       coef=double(nvarpar*dv),
                       meanvar=double(dv),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("coef","meanvar")]
      dim(vobj$coef) <- c(nvarpar,dv)
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
    }
    step <- step + 1
    if (demo) readline("Press return")
    hakt <- hakt*hincr
    lambda0 <- lambda*lseq[step]
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 .....................................
  if(graph) par(oldpar)
  if(imgtype=="rgb") {
    object$img[xind,yind,] <- theta
  } else if(imgtype=="greyscale") {
    object$img[xind,yind] <- theta
  }
  #  if(dv==1) dim(img) <- dim(img)[1:2]
  ni <- array(1,dimg0[1:2])
  ni[xind,yind] <- bi
  object$ni <- ni
  object$ni0 <- bi0
  object$hmax <- hakt/hincr
  object$call <- args
  if(estvar) {
    if(varmodel=="Linear") {
       vobj$coef[1,] <- vobj$coef[1,]/65535^2
       vobj$coef[2,] <- vobj$coef[2,]/65535
    } else vobj$coef <- vobj$coef/65535^2
    object$varcoef <- vobj$coef
    object$wghts <- vobj$wghts
  }
  if(scorr) {
    object$scorr <- spcorr
    object$chcorr <- chcorr
  }
  invisible(if(compress) compress.image(object) else object)
}
