qmeasures <- function(img,ref,
     which=c("PSNR","MAE","MSE","RMSE","SSIM","MAGE","RMSGE"),mask=FALSE){
     results <- list(NULL)
     if("PSNR" %in% which) results$PSNR <- getPSNR(img,ref,mask)
     if("MAE" %in% which) results$MAE <- mean(abs(img-ref)[mask])
     if("MSE" %in% which) results$MSE <- mean(((img-ref)^2)[mask])
     if("RMSE" %in% which) results$RMSE <- sqrt(mean(((img-ref)^2)[mask]))
     if("SSIM" %in% which) results$SSIM <- getSSIM(img,ref,mask)
     if("MAGE" %in% which) results$MAGE <- getMAGE(img,ref)
     if("RMSGE" %in% which) results$RMSGE <- getRMSGE(img,ref)
     results
   }

getPSNR <- function(img,ref,mask){
  if(!is.null(mask)){
     if(length(mask) == length(img)){
     img <- img[mask]
     ref <- ref(mask)
   }
    }
  drref <- diff(range(ref))
  20*log(drref,10)-10*log(mean((img-ref)^2),10)
}

getSSIM <- function(img,ref,alpha=1,beta=1,gamma=1,mask){

  if(!is.null(mask)){
     if(length(mask) != length(img)){
     img <- img[mask]
     ref <- ref(mask)
   }
    }
  m1 <- mean(img)
  m2 <- mean(ref)
  v1 <- var(as.vector(img))
  v2 <- var(as.vector(ref))
  c12 <- cov(as.vector(img),as.vector(ref))
  ll <- (2*m1*m2+.01^2)/(m1^2+m2^2+.01^2)
  cc <- (2*sqrt(v1)*sqrt(v2)+.03^2)/(v1+v2+.03^2)
  ss <- (c12+.03^2/2)/(sqrt(v1*v2)+.03^2/2)
  ll^alpha*cc^beta*ss^gamma
}

edges3d <- function(img){
  dimg <- dim(img)
  dx <- (img[-c(1:2),,]-img[-c(dimg[1]-1,dimg[1]),,])[,-c(1,dimg[2]),-c(1,dimg[3])]/2
  dy <- (img[,-c(1:2),]-img[,-c(dimg[2]-1,dimg[2]),])[-c(1,dimg[1]),,-c(1,dimg[3])]/2
  dz <- (img[,,-c(1:2)]-img[,,-c(dimg[3]-1,dimg[3])])[-c(1,dimg[1]),-c(1,dimg[2]),]/2
  array(sqrt(dx*dx+dy*dy+dz*dz),dimg-2)
}

edges2d <- function(img){
  dimg <- dim(img)
  dx <- (img[-c(1:2),]-img[-c(dimg[1]-1,dimg[1]),])[,-c(1,dimg[2])]/2
  dy <- (img[,-c(1:2)]-img[,-c(dimg[2]-1,dimg[2])])[-c(1,dimg[1]),]/2
  array(sqrt(dx*dx+dy*dy),dimg-2)
}

getMAGE <- function(img,refimg,tau=0){
   dimg <- dim(img)
   ret <- NULL
   if(length(dimg==2)){
     gimg <- edges2d(img)
     grefimg <- edges2d(refimg)
   } else if(length(dimg==3)) {
     if(dimg[3]==3) ret <- "Not implemented for color images"
     gimg <- edges3d(img)
     grefimg <- edges3d(refimg)
   } else {
     ret <- "illegal image dimensions"
   }
   if(is.null(ret)) ret <- mean(abs(gimg-grefimg)[gimg>tau])
   ret
}

getRMSGE <- function(img,refimg,tau=0){
   dimg <- dim(img)
   ret <- NULL
   if(length(dimg==2)){
     gimg <- edges2d(img)
     grefimg <- edges2d(refimg)
   } else if(length(dimg==3)) {
     if(dimg[3]==3) ret <- "Not implemented for color images"
     gimg <- edges3d(img)
     grefimg <- edges3d(refimg)
   } else {
     ret <- "illegal image dimensions"
   }
   if(is.null(ret)) ret <- sqrt(mean((gimg-grefimg)[gimg>tau]^2))
   ret
}

edgediff <- function(img,refimg,tau=0){
   gimg <- edges3d(img)
   grefimg <- edges3d(refimg)
   paste0("MAGE=",signif(mean(abs(gimg-grefimg)[gimg>tau]),3)," RMSGE=",signif(sqrt(mean((gimg-grefimg)[gimg>tau]^2)),3))
}
