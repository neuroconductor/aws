nlmeans <- function(x,lambda,sigma,patchhw=1,searchhw=7,pd=NULL){
#
#  Non-Local-Means Algorithm, see
#  Buades, B. Coll, J.f Morel, A review of image denoising algorithms,
#  with a new one, SIAM Multiscale Modeling and Simulation, Vol 4 (2),
#  pp: 490-530, 2005.
#
   dimx <- dim(x)
   if(is.null(dimx)) d <- 1 else d <- length(dimx)
   if(d>3) stop("no implementation for more than 3 dimensions")
   if(d==1) dimx <- length(x)
   if(d<3) dimx <- c(dimx,rep(1,3-d))
   n <- prod(dimx)
   psize <- (2*patchhw+1)^d
   if(is.null(pd)) pd <- psize
#
#  Create matrix of patches
#
   patchmat <- switch(d,.Fortran(C_fillpat1,
                                as.double(x),
                                as.integer(dimx),
                                as.integer(patchhw),
                                as.integer(psize),
                                pmat = double(psize*n))$pmat,
                       .Fortran(C_fillpat2,
                                as.double(x),
                                as.integer(dimx[1]),
                                as.integer(dimx[2]),
                                as.integer(patchhw),
                                as.integer(psize),
                                pmat = double(psize*n))$pmat,
                       .Fortran(C_fillpat3,
                                as.double(x),
                                as.integer(dimx[1]),
                                as.integer(dimx[2]),
                                as.integer(dimx[3]),
                                as.integer(patchhw),
                                as.integer(psize),
                                pmat = double(psize*n))$pmat)
  if(pd<psize){
     dim(patchmat) <- c(n,psize)
     patchmat <- prcomp(patchmat, rank. = pd)$x
  }
  dim(patchmat) <- c(n,pd)
  xhat <- .Fortran(C_nlmeans,
                   as.double(x),
                   as.integer(dimx[1]),
                   as.integer(dimx[2]),
                   as.integer(dimx[3]),
                   as.double(t(patchmat)),
                   as.integer(pd),
                   as.integer(searchhw),
                   as.double(lambda*sigma),
                   xhat=double(n))$xhat
   dim(xhat) <- dim(x)
   z <- list(theta=xhat, lambda=lambda, sigma=sigma, patchhw=patchhw, pd=psize, searchhw=searchhw)
   class(z) <- "nlmeans"
   z
}
