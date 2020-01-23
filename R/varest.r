estimateSigmaCompl <- function(magnitude,phase,mask,kstar=20,kmin=8,hsig=5,lambda=12,verbose=TRUE){
  ## kmin = 10 corresponds to an initial bandwidth of 1.47 giving positive weight to direct neighbors and
  ## 2D diagonal neigbors
  args <- sys.call(-1)
  sdim <- dim(mask)
#  if(!is.numeric(magnitude)){
#    if (verbose) cat("reading Magnitude file ... ")
#    R <- readNIfTI(magnitude, reorient = FALSE)
#  } else {
    R <- magnitude
#  }
#  if(!is.numeric(phase)){
#    if (verbose) cat("reading Phase file ... ")
#    Ph <- readNIfTI(phase, reorient = FALSE)
#  } else {
    Ph <- phase
#  }
  ComplImg <- array(0,c(2,sdim))
  ComplImg[1,,,] <- R*cos(Ph)
  ComplImg[2,,,] <- R*sin(Ph)
  ## find the number of usable cores
  mc.cores <- setCores(, reprt = FALSE)
  ##
  ##  start smoothing and variance estimation
  ##
  n <- prod(sdim)
  lambda0 <- 1e40
  sigma2 <- array(1e10,sdim)
  # just inilitialize with something large, first step is nonadaptive due to lambda0
  k <- kmin
  hmax <- 1.25^(kstar/3)
  ## preparations for median smoothing
  parammd <- getparam3d(hsig,c(1,1))
  nwmd <- length(parammd$w)

  if (verbose) pb <- txtProgressBar(min = 0, max = kstar-kmin+1, style = 3)
  bi <- array(1,sdim)
  zobj <- list(theta=ComplImg, bi=bi)
  if (verbose) {
    mae <- NULL
    protocol <- matrix("", kstar-kmin+1, 1, dimnames = list(paste("step", kmin:kstar), "protocol"))
  }
  while (k <= kstar) {
    ## determine the actual bandwidth for this step
    hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, c(1,1), 1e-4)

    ## we need the (approx.) size of the weigthing scheme array
    dlw <- (2*trunc(hakt/c(1, 1, 1))+1)[1:3]

    ## perform the actual adaptive smoothing
    zobj <- .Fortran(C_cplxawss,
                     as.double(ComplImg),
                     as.logical(mask),
                     as.integer(2),
                     as.integer(sdim[1]),
                     as.integer(sdim[2]),
                     as.integer(sdim[3]),
                     hakt = as.double(hakt),
                     as.double(lambda0),
                     as.double(zobj$theta),
                     as.double(sigma2),
                     bi = as.double(zobj$bi),
                     theta = double(2*n),
                     sigma2 = double(n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(c(1,1)),
                     double(2 * mc.cores))[c("bi", "theta", "hakt","sigma2")]
    ##
    ##  now get local median variance estimates
    ##
    dim(zobj$sigma2) <- sdim
    sigma2 <- .Fortran(C_mediansm,
                       as.double(zobj$sigma2),
                       as.logical(mask),
                       as.integer(sdim[1]),
                       as.integer(sdim[2]),
                       as.integer(sdim[3]),
                       as.integer(parammd$ind),
                       as.integer(nwmd),
                       double(nwmd*mc.cores), # work(nw,nthreds)
                       as.integer(mc.cores),
                       sigma2n = double(n))$sigma2n/0.6931
    # sigma2n containes sum of 2 independent squared residuals
    # 0.6931 approximates  median \chi_2 /2
    # needed to get correct results
    ## use maximum ni
    bi <- zobj$bi <- pmax(bi, zobj$bi)

    ## some verbose stuff
    if (verbose) {
      protocol[k-kmin+1,1] <- paste("bandwidth: ", signif(hakt, 3),
                                    "sigma: mean: ", signif(sqrt(mean(sigma2[mask])),3),
                                    "median: ", signif(sqrt(median(sigma2[mask])),3),
                                    "sd: ", signif(sd(sqrt(sigma2[mask])),3),
                                    "median(bi):", signif(median(zobj$bi[mask]),3),
                                    "max(bi):", signif(max(zobj$bi[mask]),3))
      setTxtProgressBar(pb, k-kmin+1)
    }

    ## go for next iteration
    k <- k+1
    lambda0 <- lambda
    gc()
  }
  dim(zobj$theta) <- c(2,sdim)
  # return estimated parameters of rician distribution
  list(sigma=array(sqrt(sigma2),sdim),
       theta=array(sqrt(zobj$theta[1,,,]^2+zobj$theta[2,,,]^2),sdim),
       sigmal=array(sqrt(zobj$sigma2),sdim),mask=mask,
       protocol=protocol,args=args)
}

awsLocalSigma <- function(y, steps, mask, ncoils, vext=c(1,1),
      lambda = 5, minni=2, hsig=5, sigma=NULL, family=c("NCchi","Gauss"),
      verbose = FALSE, trace = FALSE, u=NULL){
if(trace) tergs <- array(0,c(steps,4,sum(mask))) else tergs <- NULL
family <- match.arg(family[1],c("NCchi","Gauss","Gaussian"))

if(family == "NCchi"){
  varstats <- sofmchi(ncoils)
  th <- seq(0,30,.01)
  z <- fncchiv(th,varstats)
  minz <- min(z)
  th <- th[z>min(z)]
  minth <- min(th)
  z <- z[z>min(z)]
  nz <- nls(z~(a*th+b*th*th+c*th*th*th+d)/(a*th+b*th*th+c*th*th*th+1+d),data=list(th=th,z=z),start=list(a=1,b=1,c=1,d=1))
  vpar <- c(minth,minz,coef(nz))
  ##  this provides an excellent approximation for the variance reduction in case of low ncp
}

if(length(vext)==3) vext <- vext[2:3]/vext[1]
## dimension and size of cubus
ddim <- dim(y)
n <- prod(ddim)

## test dimension
if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")

## check mask
if (is.null(mask)) mask <- array(TRUE, ddim)
if(length(mask) != n) stop("dimensions of data array and mask should coincide")

## initial value for sigma_0
if(is.null(sigma)){
  # sigma <- sqrt( mean( y[mask]^2) / 2 / ncoils)
  sigma <- IQQdiff( y, mask, .25, verbose=verbose)
  if("NCchi" == family){
## iterative improvement for NC chi-distribution
    sigma <- iniSigmaNCchi( y, mask, .25, ncoils, sigma)
    sigma <- iniSigmaNCchi( y, mask, .25, ncoils, sigma)
    sigma <- iniSigmaNCchi( y, mask, .25, ncoils, sigma)
  }
}
##
##   Prepare for diagnostics plots
##
if(verbose){
  mslice <-  (ddim[3]+1)/2
  ymslice <- y[,,mslice]
  ymslice[!mask[,,mslice]] <- 0
  if(!is.null(u)&&"NCchi" == family){
    par(mfrow=c(2,4),mar=c(3,3,3,1),mgp=c(2,1,0))
  } else {
    par(mfrow=c(2,3),mar=c(3,3,3,1),mgp=c(2,1,0))
  }
} else {
  cat("step")
}
## define initial arrays for parameter estimates and sum of weights (see PS)
th <- y
ksi <- array( y^2, ddim)
ni <- array( 1, ddim)
sigma <- array(sigma, ddim)
# initialize array for local sigma by global estimate
mc.cores <- setCores(,reprt=FALSE)
## preparations for median smoothing
parammd <- getparam3d(hsig, vext)
  ## iterate PS starting with bandwidth h0
if(family=="NCchi"){
      #  precompute values of log(besselI) for interpolation
  nfb <- 200 * ncoils  # for arguments > nfb  asymp. approx can be used
  x <- 1:nfb
  flb <- x + log(besselI(x, ncoils-1, TRUE))
}
for (i in 1:steps) {

  h <- 1.25^((i-1)/3)
  param <- getparam3d(h, vext)
  nw <- length(param$w)
        ## perform one step PS with bandwidth h
  if(family=="NCchi"){
    z <- .Fortran(C_awslchi2,
                  as.double(y),# data
                  as.double(ksi),# \sum w_ij S_j^2
                  ni = as.double(ni),
                  as.double(sigma),
                  as.double(vpar),# parameters for var. reduction
                  as.double(ncoils),
                  as.integer(mask),
                  as.integer(ddim[1]),
                  as.integer(ddim[2]),
                  as.integer(ddim[3]),
                  as.integer(param$ind),
                  as.double(param$w),
                  as.integer(nw),
                  as.double(minni),
                  double(nw*mc.cores), # wad(nw,nthreds)
                  double(nw*mc.cores), # sad(nw,nthreds)
                  as.double(lambda),
                  as.integer(mc.cores),
                  as.integer(floor(ncoils)),
                  double(floor(ncoils)*mc.cores), # work(L,nthreds)
                  th = double(n),
                  sigman = double(n),
                  ksi = double(n),
                  as.double(flb),
                  as.integer(nfb))[c("ni","ksi","th","sigman")]
    thchi <- z$th
    ksi <- z$ksi
    thchi[!mask] <- 0
  } else{
    ## assume Gaussian data
    z <- .Fortran(C_awslgaus,
                  as.double(y),# data
                  as.double(th),# \sum w_ij S_j
                  ni = as.double(ni),
                  as.double(sigma),
                  as.integer(mask),
                  as.integer(ddim[1]),
                  as.integer(ddim[2]),
                  as.integer(ddim[3]),
                  as.integer(param$ind),
                  as.double(param$w),
                  as.integer(nw),
                  as.double(minni),
                  as.double(lambda),
                  th = double(n),
                  sigman = double(n))[c("ni","th","sigman")]

  }
    ## extract sum of weigths (see PS) and consider only voxels with ni larger then mean
  th <- array(z$th,ddim)
  ni <- array(z$ni,ddim)
  z$sigman[z$sigman==0] <- median(z$sigman[z$sigman>0])
  nmask <- sum(mask)
  if(verbose) cat("local estimation in step ",i," completed",format(Sys.time()),"\n")
  ##
  ##  nonadaptive smoothing of estimated standard deviations
  ##
  if(any(ni[mask]>minni)){
     z$sigman[mask] <- z$sigman[mask]*(1+runif(nmask,-.0001,.0001))
  ##
  ##  avoid ties in local neighborhood cubes of z$sigman
  ##
     nwmd <- length(parammd$w)
     sigma <- .Fortran(C_mediansm,
                    as.double(z$sigman),
                    as.integer(mask),
                    as.integer(ddim[1]),
                    as.integer(ddim[2]),
                    as.integer(ddim[3]),
                    as.integer(parammd$ind),
                    as.integer(nwmd),
                    double(nwmd*mc.cores), # work(nw,nthreds)
                    as.integer(mc.cores),
                    sigman = double(n))$sigman
  }
  dim(sigma) <- ddim
  mask[sigma==0] <- FALSE
  if(verbose) cat("local median smoother in step ",i," completed",format(Sys.time()),"\n")
  ##
  ##  diagnostics
  ##
  if(verbose){
    meds <- median(sigma[mask])
    means <- mean(sigma[mask])
    image(ymslice,col=grey(0:255/255))
    title(paste("S  max=",signif(max(y[mask]),3)," median=",signif(median(y[mask]),3)))
    image(th[,,mslice],col=grey(0:255/255))
    title(paste("E(S)  max=",signif(max(th[mask]),3)," median=",signif(median(th[mask]),3)))
    image(sigma[,,mslice],col=grey(0:255/255),zlim=c(0,max(sigma[mask])))
    title(paste("sigma max=",signif(max(sigma[mask]),3)," median=",signif(meds,3)))
    image(ni[,,mslice],col=grey(0:255/255))
    title(paste("Ni    max=",signif(max(ni[mask]),3)," median=",signif(median(ni[mask]),3)))
    plot(density(sigma[mask]),main="density of sigma")
    plot(density(ni[mask]),main="density of Ni")
    cat("mean sigma",means,"median sigma",meds,"sd sigma",sd(sigma[mask]),"\n")
    if(!is.null(u)&&"NCchi" == family){
      thchims <- fncchir(th/sigma,varstats)*sigma
      thchims[!mask] <- u[!mask]
      image(abs(thchims-u)[,,mslice],col=grey(0:255/255))
      title("abs Error in thchims")
      plot(density((thchims-u)[mask]),main="density of thchims-u")
      cat("MAE(th)",mean(abs(thchi-u)[mask]),"RMSE(th)",sqrt(mean((thchi-u)[mask]^2)),"MAE(thms)",mean(abs(thchims-u)[mask]),"RMSE(thms)",sqrt(mean((thchims-u)[mask]^2)),"\n")
    }
  } else {
    cat(" ",i)
  }
  if(trace) {
      tergs[i,1,] <- ni[mask]
      tergs[i,2,] <- th[mask]
      tergs[i,3,] <- z$sigma[mask]
      tergs[i,4,] <- sigma[mask]
  }
}
## END PS iteration
if(!verbose) cat("\n")
if(family=="NCchi"){
   thchi <- fncchir(th/sigma,varstats)*sigma
   thchi[!mask] <- 0
} else {
  thchi <- NULL
}
## this is the result (th is expectation, not the non-centrality parameter !!!)
invisible(list(sigma = sigma,
               sigmal = array(z$sigman,ddim),
               theta = th,
               thchi = thchi,
               ni  = ni,
               tergs = tergs,
               mask = mask))
      }

IQQ <- function (x, q = .25, na.rm = FALSE, type = 7){
  cqz <- qnorm(.05)/qnorm(q)
  x <- as.numeric(x)
  z <- diff(quantile(x, c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))
  z0 <- 0
  while(abs(z-z0)>1e-5*z0){
    # outlier removal
    z0 <- z
    #       cat(sum(x>(z0*cqz))," ")
    x <- x[x<(z0*cqz)]
    z <- diff(quantile(x, c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))
  }
  z
}

IQQdiff <- function(y, mask, q = .25, verbose = FALSE) {
  cq <- qnorm(1-q)*sqrt(2)*2
  sx <- IQQ( diff(y[mask]), q)/cq
  sy <- IQQ( diff(aperm(y,c(2,1,3))[aperm(mask,c(2,1,3))]), q)/cq
  sz <- IQQ( diff(aperm(y,c(3,1,2))[aperm(mask,c(3,1,2))]), q)/cq
  if(verbose) cat( "Pilot estimates of sigma", sx, sy, sz, "\n")
  min( sx, sy, sz)
}

iniSigmaNCchi <- function(y, mask, q, L, sigma, verbose=FALSE){
#
#  robust improvement for sigma using moment estimates and IQQdiff
#
  meany <- mean(y[mask]^2)
  eta <- sqrt(max( 0, meany/sigma^2-2*L))
  m <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)*hg1f1(-1/2,L,-eta^2/2)
  v <- max( .01, 2*L+eta^2-m^2)
  if(verbose) cat( eta, m, v, "\n")
  IQQdiff( y, mask, q)/sqrt(v)
}

AFLocalSigma <- function(y,ncoils,level=NULL,mask=NULL,h=2,hadj=1,vext = c( 1, 1)){
  ##
  ##   estimate effective sigma and effective ncoils (L) according to Aja-Fernandez MRI 2013
  ##
  #    sh2B stands for \hat{\sigma}^2_{n,B}
  #    sh2L stands for \hat{\sigma}^2_{nL}
  #    sh2eB stands for \hat{\sigma}^2_{eff,B}
  #    sh2eS stands for \hat{\sigma}^2_{eff,S}
  #    LeB stands for \hat{L}_{eff,B}
  #    vmlb  -  Var(M_L(x)|x \in B)
  #    vmlbx -  Var(M_L(x)|x \in B)_x
  #    phix  -  \hat{\Phi}_n(x)
  #    seff - 3D array of local variance estimates
  #    Leff - 3D array of effective L
  ddim <- dim(y)
  n <- prod(ddim)
  if(is.null(level)&is.null(mask)){
    warning("no background definition, need either level or mask")
    return(invisible(NULL))
  }
  if(is.null(level)) indB <- !mask else indB <- y<level
  sh2Lsimple <- mean(y[indB]^2)/2
  mask1 <- array(TRUE,ddim)
  ##
  ##  local variance estimates in vx
  ##
  vx <- .Fortran(C_afmodevn,
                 as.double(y),
                 as.integer(ddim[1]),
                 as.integer(ddim[2]),
                 as.integer(ddim[3]),
                 as.integer(mask1),
                 as.double(h),
                 as.double(vext),
                 sigma = double(n))$sigma
  dim(vx) <- ddim
  vxb <- vx[indB]
  vxs <- vx[!indB]
  dv2b <- density( vxb[vxb>0], n = 4092, adjust = hadj, to = min( max(vxb[vxb>0]), median(vxb[vxb>0])*5) )
  sh2B <- dv2b$x[dv2b$y == max(dv2b$y)][1]
  dv2s <- density( vxs[vxs>0], n = 4092, adjust = hadj, to = min( max(vxs[vxs>0]), median(vxs[vxs>0])*5) )
  sh2eS <- dv2s$x[dv2s$y == max(dv2s$y)][1]
  m2 <- .Fortran(C_afmodem2,
                 as.double(y),
                 as.integer(ddim[1]),
                 as.integer(ddim[2]),
                 as.integer(ddim[3]),
                 as.integer(mask1),
                 as.double(h),
                 as.double(vext),
                 sm = double(n))$sm
  dim(m2) <- ddim
  m2b <- m2[indB]
  dm2b <- density( m2b[m2b>0], n = 4092, adjust = hadj, to = min( max(m2b[m2b>0]), median(m2b[m2b>0])*5) )
  sh2L <- dm2b$x[dm2b$y == max(dm2b$y)][1]/2
  cat("sh2B",sh2B,"sh2Lsimple",sh2Lsimple,"sh2L",sh2L,"\n")
  ##
  ##  now get sh2eB and LeB
  ##
  LeBn <- ncoils
  sh2eBn <- sh2B/2/(LeBn-gamma(LeBn+.5)^2/gamma(LeBn)^2)
  LeB <- 1
  sh2eB <- 0
  while(abs(sh2eB-sh2eBn)>1e-5||abs(LeB-LeBn)>1e-5){
    LeB <- LeBn
    sh2eB <- sh2eBn
    LeBn <- sh2L/sh2eB
    sh2eBn <- sh2B/2/(LeBn-gamma(LeBn+.5)^2/gamma(LeBn)^2)
    cat("LeB",LeB,"LeBn",LeBn,"sh2eB",sh2eB,"sh2eBn",sh2eBn,"\n")
  }
  sh2eB <- sh2eBn
  LeB <- sh2L/sh2eB
  phix <- pmax(0,pmin(1,sh2L/(m2-sh2L)))
  dim(phix) <- dim(m2)
  seff <- sqrt((1-phix)*sh2eS+phix*sh2eB)
  Leff <- sh2L/seff^2
  list(sigma=seff,Leff=Leff,sh2L=sh2L,sh2B=sh2B,sh2eS=sh2eS,sh2eB=sh2eB,LeB=LeB,m2=m2,vx=vx,indB=indB)
}

estGlobalSigma <- function(y, mask=NULL, ncoils=1, steps=16, vext=c(1,1),
                  lambda=20, hinit=2, hadj=1, q=.25, level=NULL,
                  sequence=FALSE,method=c("awsVar","awsMAD",
                  "AFmodevn","AFmodem1chi","AFbkm2chi","AFbkm1chi")){
##
##  estimate global scale parameter sigma for NCchi distributed data
##
  method <- match.arg(method)
    qni <- .8
## methods using the propagation-separation (PS) approach
  if(method %in% c("awsVar","awsMAD")){
    varstats <- sofmchi(ncoils)
    if(length(vext)==3) vext <- vext[2:3]/vext[1]
    ## dimension and size of cubus
    ddim <- dim(y)
    n <- prod(ddim)

    ## test dimension
    if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")

    ## check mask
    if (is.null(mask)) mask <- array(TRUE, ddim)
    if(length(mask) != n) stop("dimensions of data array and mask should coincide")

    ## initial value for sigma_0
    # sigma <- sqrt( mean( y[mask]^2) / 2 / ncoils)
    sigma <- IQQdiff( y, mask, q)
    #  cat( "sigmahat1", sigma, "\n")
    sigma <- iniSigmaNCchi( y, mask, q, ncoils, sigma)
    #  cat( "sigmahat2", sigma, "\n")
    sigma <- iniSigmaNCchi( y, mask, q, ncoils, sigma)
    #  cat( "sigmahat3", sigma, "\n")
    sigma <- iniSigmaNCchi( y, mask, q, ncoils, sigma)
    #  cat( "sigmahat4", sigma,"\n")

    ## define initial arrays for parameter estimates and sum of weights (see PS)
    th <- array( 1, ddim)
    ni <- array( 1, ddim)
    #  y <- y/sigma # rescale to avoid passing sigma to awsvchi
    minlev <- sqrt(2)*gamma(ncoils+.5)/gamma(ncoils)
    if (sequence) sigmas <- lind <- minni <- numeric(steps)
    mc.cores <- setCores(,reprt=FALSE)
    ## iterate PS starting with bandwidth hinit
    for (i in 1:steps) {

      h <- hinit * 1.25^((i-1)/3)
      param <- getparam3d(h,vext)
      nw <- length(param$w)
      fncchi <- fncchiv(th/sigma,varstats)
      ## correction factor for variance of NC Chi distribution
      ## perform one step PS with bandwidth h
      if(method=="awsVAR"){
        z <- .Fortran(C_awsvchi,
                      as.double(y),        # data
                      as.double(th),       # previous estimates
                      ni = as.double(ni),
                      as.double(fncchi/2),
                      as.integer(mask),
                      as.integer(ddim[1]),
                      as.integer(ddim[2]),
                      as.integer(ddim[3]),
                      as.integer(param$ind),
                      as.double(param$w),
                      as.integer(nw),
                      as.double(lambda),
                      as.double(sigma),
                      th = double(n),
                      sy = double(n))[c("ni","th","sy")]
      } else { # method=="awsMAD"
        z <- .Fortran(C_awsadchi,
                      as.double(y),        # y(n1,n2,n3)
                      as.double(th),       # th(n1,n2,n3)
                      ni = as.double(ni),  # ni(n1,n2,n3)
                      as.double(fncchi/2), # fns(n1,n2,n3)
                      as.integer(mask),    # mask(n1,n2,n3)
                      as.integer(ddim[1]), # n1
                      as.integer(ddim[2]), # n2
                      as.integer(ddim[3]), # n3
                      as.integer(param$ind), # ind(3,nw)
                      as.double(param$w), # w(nw)
                      as.integer(nw), # nw
                      as.double(lambda), # lambda
                      as.double(sigma), # sigma
                      double(nw*mc.cores), # wad(nw,nthreds)
                      as.integer(mc.cores), # nthreds
                      th = double(n), # thn(n1*n2*n3)
                      sy = double(n))[c("ni","th","sy")]
      }
      ## extract sum of weigths (see PS) and consider only voxels with ni larger then mean
      th <- z$th
      ni <- z$ni
      ni[!mask]<-1
      ind <- (ni > .9999*quantile(ni[ni>1],qni))#&(z$th>sigma*minlev)
      ## use correction factor for sd of NC Chi distribution
      sy1 <- z$sy[ind]
      th1 <- th[ind]
      sy1 <- sy1/fncchis(th1/sigma,varstats)
      ## use the maximal mode of estimated local sd parameters, exclude largest values for better precision
      dsigma <- density( sy1[sy1>0], n = 4092, adjust = hadj, to = min( max(sy1[sy1>0]), median(sy1[sy1>0])*5) )
      sigma <- dsigma$x[dsigma$y == max(dsigma$y)][1]

      if (sequence) {
        sigmas[i] <- sigma
        lind[i] <- sum(ind)
        minni[i] <- min(ni[ind])
      }

    }
    ## END PS iteration
    ## this is the result (th is expectation, not the non-centrality parameter !!!)
    result <- list(sigma = if(sequence) sigmas else sigma,
                   theta = th,
                   lind  = if(sequence) lind else sum(ind),
                   minni  = if(sequence) minni else min(ni[ind]))
  }
##
##  estimation of sigma from background or brain area
##  various methods from table 2 in Aja-Fernandez 2009 following Aja-Fernandez
##
  if(method%in%c("AFmodevn","AFmodem1chi","AFbkm2chi","AFbkm1chi")){
    if(method!="AFmodevn"&is.null(level)&is.null(mask)) {
      stop("need information on background using either level or mask")
    }
    ## dimension and size of cubus
    ddim <- dim(y)
    n <- prod(ddim)

    ## test dimension
    if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")

    ## check mask
    if (is.null(mask)) mask <- array(TRUE, ddim)
    if(!is.null(level)){
      if (method=="AFmodevn"){
        mask[y<level] <- FALSE
      } else {
        mask[y>level] <- FALSE
      }
    }
    if(length(mask) != n) stop("dimensions of data array and mask should coincide")

    ## let FORTRAN do the calculation
    if(method%in%c("AFmodevn","AFmodem1chi")){
      if(method=="AFmodevn"){
        sigma <- .Fortran(C_afmodevn,
                          as.double(y),
                          as.integer(ddim[1]),
                          as.integer(ddim[2]),
                          as.integer(ddim[3]),
                          as.integer(mask),
                          as.double(hinit),
                          as.double(vext),
                          sigma = double(n))$sigma
        sigma <- sigma/2/(ncoils-gamma(ncoils+.5)^2/gamma(ncoils)^2)
        sigma <- array( sqrt(sigma), ddim)
      } else {
        afactor <- sqrt(1/2)*gamma(ncoils)/gamma(ncoils+.5)
        sigma <- .Fortran(C_afmodem1,
                          as.double(y),
                          as.integer(ddim[1]),
                          as.integer(ddim[2]),
                          as.integer(ddim[3]),
                          as.integer(mask),
                          as.double(hinit),
                          as.double(vext),
                          sigma = double(n))$sigma
        sigma <- array( afactor*sigma, ddim)
      }
      ##  use the maximal mode of estimated local variance parameters, exclude largest values for better precision
      dsigma <- density( sigma[sigma>0], n = 4092, adjust = hadj, to = min( max(sigma[sigma>0]), median(sigma[sigma>0])*5) )
      sigmag <- dsigma$x[dsigma$y == max(dsigma$y)][1]
    } else {
      if(method=="AFbkm2chi"){
        sigmag <- sqrt(mean(y[mask]^2)/2/ncoils)
      } else {
        sigmag <- mean(y[mask])*sqrt(ncoils/2)*gamma(ncoils)/gamma(ncoils+.5)/sqrt(ncoils)
      }
    }
    result <- list(sigma = sigmag)
  }
  class(result) <- "NCchiGlobalSigma"
  result
}

awslinsd <- function(y,hmax=NULL,hpre=NULL,h0=NULL,mask=NULL,
                     ladjust=1,wghts=NULL,varprop=.1,A0,A1)
{
  #
  #    first check arguments and initialize
  #
  homogen <- TRUE
  wghts <- NULL

  args <- match.call()
  dy<-dim(y)
  if(length(dy)!=3) stop("Image should be 3D")
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  #
  #   set appropriate defaults
  #
  if(is.null(wghts)) wghts <- c(1,1,1)
  wghts <- wghts[1]/wghts[2:3]
  cpar<-setawsdefaults(mean(y),ladjust,hmax,wghts)
  qlambda <- .98
  if(is.null(hmax)) hmax <- 5
  lambda <- ladjust*qchisq(qlambda,1)*2
  maxvol <- getvofh(hmax,c(1,0,0,1,0,1),c(1,wghts))
  kstar <- as.integer(log(maxvol)/log(1.25))
  k <- 6
  cat("Estimating variance model using PS with with lambda=",signif(lambda,3)," hmax=",hmax,"number of iterations:",kstar-k+1,"\n")
  if(is.null(mask)) {
    if(length(dy)==0) mask <- rep(TRUE,length(y)) else mask <- array(TRUE,dy)
  }
  n<-length(y)
  dmask <- dim(mask)
  nvoxel <- sum(mask)
  position <- array(0,dmask)
  position[mask] <- 1:nvoxel
  y <- y[mask]
  #
  #   family dependent transformations
  #
  sigma2 <- max(1,IQRdiff(as.vector(y))^2)
  if(any(h0)>0) sigma2<-sigma2*Varcor.gauss(h0)
  #  cat("Estimated variance: ", signif(sigma2,4),"\n")
  sigma2 <- rep(sigma2, nvoxel)
  sigma2 <- 1/sigma2
  #  taking the invers yields simpler formulaes
  # now check which procedure is appropriate
  ##  this is the version on a grid
  #
  #    Initialize  for the iteration
  #
  zobj<-list(ai=y, bi= rep(1,nvoxel), theta= y)
  mae<-NULL
  lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
  #
  #   produce a presmoothed estimate to stabilze variance estimates
  #
  if(is.null(hpre)) hpre<-20^(1/3)
  dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:3]
  hobj <- .Fortran(C_smooth3d,
                     as.double(y[mask]),
                     as.double(rep(1,nvoxel)),
                     as.integer(position),
                     as.integer(0L),
                     as.integer(nvoxel),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(1L),
                     as.double(hpre),
                     theta=as.double(zobj$theta),
                     bi=as.double(zobj$bi),
                     as.integer(2L), # lkern
                     double(prod(dlw)),
                     as.double(wghts),
                     double(1))[c("bi","theta")]
  theta <- bi <- array(0,dy)
  theta[mask] <- hobj$theta
  bi[mask] <- hobj$bi
  hobj <- list(bi=bi, theta=theta)
    #
  #   iteratate until maximal bandwidth is reached
  #
  #  cat("Progress:")
  #  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  pb <- txtProgressBar(0, kstar, style = 3)
  while (k<=kstar) {
    hakt0 <- gethani(1,10,1.25^(k-1),c(1,0,0,1,0,1),c(1,wghts),1e-4)
    hakt <- gethani(1,10,1.25^k,c(1,0,0,1,0,1),c(1,wghts),1e-4)
    dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:3]
    if(any(h0>0)) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,3)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
    # Correction for spatial correlation depends on h^{(k)}
    hakt0<-hakt
    # heteroskedastic Gaussian case
    zobj <- .Fortran(C_cgaws,
                       as.double(y),
                       as.integer(position),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi = as.double(zobj$bi),
                       bi2 = double(n),
                       bi0 = as.double(zobj$bi0),
                       gi = double(n),
                       gi2 = double(n),
                       ai = as.double(zobj$ai),
                       as.integer(2L),
                       as.double(0.25),
                       double(prod(dlw)),
                       as.double(wghts)
                     )[c("bi", "bi0", "bi2", "gi2", "ai", "gi", "hakt")]
    zobj$theta <- (zobj$ai/zobj$bi)
    #
    #    Calculate MAE and MSE if true parameters are given in u
    #    this is for demonstration and testing for propagation (parameter adjustments)
    #    only.
    #
    #   Prepare for next iteration
    #
    #
    #   Create new variance estimate
    #
    vobj <- awsgsigmasd(y,hobj,zobj,varprop,h0,A0,A1)
    sigma2 <- vobj$sigma2inv
    coef <- vobj$coef
    rm(vobj)
    lambda0<-lambda
    setTxtProgressBar(pb, k)
    #    if (max(total) >0) {
    #      cat(signif(total[k],2)*100,"% . ",sep="")
    #    }
    k <- k+1
    gc()
  }
  close(pb)
  # cat("\n")
  theta <-array(0, dmask)
  theta[mask] <- zobj$theta
  ###
  ###            end iterations now prepare results
  ###
  list(theta = theta, vcoef=coef, mask=mask)
}

############################################################################
#
#  estimate inverse of variances, uses nonadaptive hobj to stabilize,
#    based on absolute residuals, linear
#  expects only voxel within mask
#
############################################################################
awsgsigmasd <- function(y,hobj,tobj,varprop,h0,thmin,thmax){
  ## specify sd to be linear to mean
  thrange <- range(y)
  #    thmax <- thrange[2]-.2*diff(thrange)
  #    thmin <- thrange[1]+.1*diff(thrange)
  ind <- tobj$gi>1.5&tobj$theta>thmin&tobj$theta<thmax
  absresid <- abs((y-tobj$theta)[ind]*tobj$gi[ind]/sqrt(tobj$gi[ind]^2-tobj$gi2[ind]))/.8
  #     absresid <- abs((y-tobj$theta)[ind])/.8
  theta <- tobj$theta[ind]
  wght <- (tobj$gi[ind]^2-tobj$gi2[ind])/tobj$gi[ind]^2
  coef <- coefficients(lm(absresid~theta,weights=wght^2))
  # force positive variance for positive mean by increasing variance estimate
  if(coef[2] < 0){
    coef[1] <- coefficients(lm(absresid~1,weights=wght^2))
    coef[2] <- 0
  }
  if(coef[1] <0.5){
    coef[2] <- coefficients(lm(absresid~theta-1,weights=wght^2))
    coef[1] <- 0
  }
  gamma <- pmin(tobj$gi/hobj$bi,1)
  theta <- gamma*tobj$theta+(1-gamma)*hobj$theta
  sigma2 <- (coef[1]+coef[2]*theta)^2
  varquantile <- varprop*mean(sigma2)
  sigma2 <- pmax(sigma2,varquantile)
  #  cat("Estimated mean variance",signif(mean(sigma2),3)," Variance parameters:",signif(coef,3),"\n")
  list(sigma2inv=1/sigma2,coef=coef)
}
