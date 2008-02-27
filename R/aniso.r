#
#    R - function  aws  for estimating the flow direction in an image
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2008 Weierstrass-Institut fuer
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
IQRdiff <- function(y) IQR(diff(y))/1.908
awsflow2 <- function(y,kstar,h0=1.12,h1=1.12,ch=sqrt(1.25),ch1=sqrt(1.25),cg=1,lambda1=10,lambda2=10,graph=TRUE,adapt=TRUE,
           u=NULL,ugamma=NULL,demo=FALSE){
args <- match.call()
dy<-dim(y)
n <- prod(dy)
# Initialization
z1 <- .Fortran("initflob",
                as.double(y),
                as.integer(dy[1]),
                as.integer(dy[2]),
                as.double(h0),
                ahat=double(n),
                chat=double(n),
                gamma=double(2*n),
                Ni=double(n),
                Bmat=double(3*n),
                bvec=double(2*n),
                Svec=double(2*n),
                Si=double(n),
                Pvec=double(2*n),
                Qmat=double(3*n),
                PACKAGE="aws")[c("ahat","chat","gamma","Ni","Bmat",
                         "bvec","Svec","Si","Pvec","Qmat")]
ahat <- matrix(z1$ahat,dy[1],dy[2])
chat <- matrix(z1$chat,dy[1],dy[2])
Ni <- matrix(z1$Ni,dy[1],dy[2])
Si <- matrix(z1$Si,dy[1],dy[2])
Bmat <- array(z1$Bmat,c(3,dy))
bvec <- array(z1$bvec,c(2,dy))
Svec <- array(z1$Svec,c(2,dy))
Qmat <- array(z1$Qmat,c(3,dy))
Pvec <- array(z1$Pvec,c(2,dy))
gamma <- array(z1$gamma,c(2,dy))
#gamma <- .Fortran("ingamma2",
#                  as.double(y),
#                  as.integer(dy[1]),
#                  as.integer(dy[2]),
#                  gamma=double(2*n),
#                  PACKAGE="aws")$gamma
#dim(gamma) <- c(2,dy)
ngamma <- sqrt(gamma[1,,]^2+gamma[2,,]^2)
if(graph) {
par(mfrow=c(2,2),mar=c(3,3,3,.1),mgp=c(2,1,0))
image(atan(gamma[1,,]/gamma[2,,])+pi*(gamma[2,,]<0),col=rainbow(512),zlim=c(-pi/2,3/2*pi))
title(paste("Orientation h=",signif(h0,3)))
image(Ni, col=grey(0:255/255))
title(paste("Ni",signif(min(Ni),3),signif(mean(Ni),3),signif(max(Ni),3)))
#image(matrix(z1$ahat,dy[1],dy[2]), col=grey(0:255/255),zlim=quantile(z1$ahat,c(.01,.99)))
#title(paste("ahat",signif(quantile(z1$ahat,.01),3),signif(quantile(z1$ahat,.99),3)))
#image(matrix(z1$chat,dy[1],dy[2]), col=grey(0:255/255),zlim=quantile(z1$chat,c(.01,.99)))
#title(paste("chat",signif(quantile(z1$chat,.01),3),signif(quantile(z1$chat,.99),3)))
image(matrix(z1$ahat,dy[1],dy[2]), col=grey(0:255/255))
title(paste("ahat",signif(min(z1$ahat),3),signif(max(z1$ahat),3)))
image(matrix(z1$chat,dy[1],dy[2]), col=grey(0:255/255))
title(paste("chat",signif(max(z1$chat),3)))
}
if(!is.null(u)) cat("h0=",signif(h0,2),"MAE",signif(mean(abs(z1$ahat-u)),3),"MSE",signif(mean((z1$ahat-u)^2),3),"\n")
if(!is.null(ugamma)) cat("h0=",signif(h0,2),"MSE(gamma)",signif(mean(1-abs(gamma[1,,]*ugamma[1,,]+gamma[2,,]*ugamma[2,,])/ngamma),3),"\n")
if(demo) tmp <- readline("Press enter for next step")
#Pvec <- double(2*n)
#Qmat <- double(3*n)
#Ni <- rep(1,n)
sigma2 <- min(IQRdiff(as.vector(y))^2,IQRdiff(as.vector(t(y)))^2)
cat("estimated sigma",signif(sqrt(sigma2),2),"\n")
# Iterate
ll <- length(lambda1)
if(ll<kstar+1) lambda1 <- c(lambda1,rep(lambda1[ll],kstar-ll+1))
ll <- length(lambda2)
if(ll<kstar+1) lambda2 <- c(lambda2,rep(lambda2[ll],kstar-ll+1))
step <- 0
hakt <- h0*ch^step
hakt1 <- h1*ch1^step
gakt <- h0*cg^step
etaaktinv <- sqrt(hakt^2-gakt^2)/(hakt*gakt)
haktinv <- 1/hakt
while(step<kstar+1){
#  step 1
z1 <- .Fortran("itstep1b",
              as.double(y),
              as.integer(dy[1]),
              as.integer(dy[2]),
              as.double(gamma),
              as.double(haktinv),
              as.double(etaaktinv),
              as.double(sigma2),
              as.double(lambda1[step+1]),
              Bmat=as.double(Bmat),
              bvec=as.double(bvec),
              Svec=as.double(Svec),
              Ni=as.double(Ni),
              Si=as.double(Si),
              chat=as.double(chat),
              ahat=as.double(ahat),
              as.double(ngamma),
              as.double(h0*h0),# do not use adapive weights if distance is less than h0
              PACKAGE="aws")[c("Bmat","bvec","Svec","Ni","Si","chat","ahat")]
#  step 2 
z1$ahat[z1$Ni==0] <- y[z1$Ni==0]
z1$Ni[z1$Ni==0] <- 1
ahat <- matrix(z1$ahat,dy[1],dy[2])
chat <- matrix(z1$chat,dy[1],dy[2])
Si <- matrix(z1$Si,dy[1],dy[2])
Bmat <- array(z1$Bmat,c(3,dy))
bvec <- array(z1$bvec,c(2,dy))
Svec <- array(z1$Svec,c(2,dy))
z2 <- .Fortran("itfstep2",
              as.double(gamma),
              as.integer(dy[1]),
              as.integer(dy[2]),
              as.double(hakt1),
              as.double(sigma2),
              as.double(lambda2[step+1]),
              Pvec=as.double(Pvec),
              Qmat=as.double(Qmat),
              as.double(Ni),
              as.double(chat),
              as.double(ahat),
              as.double(Bmat),
              as.double(bvec),
              as.double(Svec),
              gamma=double(2*n),
              ngamma=as.double(ngamma),
              as.double(h0*h0),# do not use adapive weights if distance is less than h0
              PACKAGE="aws")[c("Qmat","Pvec","gamma","ngamma")]
#  organize everything for the next step
Ni <- matrix(z1$Ni,dy[1],dy[2])
Qmat <- array(z2$Qmat,c(3,dy))
Pvec <- array(z2$Pvec,c(2,dy))
gamma <- array(z2$gamma,c(2,dy))
ngamma <- matrix(z2$ngamma,dy[1],dy[2])
if(graph) {
par(mfrow=c(2,2),mar=c(3,3,3,.1),mgp=c(2,1,0))
image(atan(gamma[1,,]/gamma[2,,])+pi*(gamma[2,,]<0),col=rainbow(512),zlim=c(-pi/2,3/2*pi))
title(paste("Orientation h=",signif(hakt,3),"g=",signif(gakt,3)))
image(matrix(Ni,dy[1],dy[2]), col=grey(0:255/255))
title(paste("Ni",signif(min(Ni),3),signif(mean(Ni),3),signif(max(Ni),3)))
image(ahat, col=grey(0:255/255))
title(paste("ahat",signif(min(ahat),3),signif(max(ahat),3)))
image(chat, col=grey(0:255/255))
title(paste("chat",signif(max(chat),3)))
}
if(!is.null(u)) cat("h=",signif(hakt,2),"g=",signif(gakt,2),"MAE",signif(mean(abs(z1$ahat-u)),3),"MSE",signif(mean((z1$ahat-u)^2),3),"\n")
if(!is.null(ugamma)) cat("h1=",signif(hakt1,2),"MSE(gamma)",signif(mean(1-abs(gamma[1,,]*ugamma[1,,]+gamma[2,,]*ugamma[2,,])/ngamma),3),"\n")
cat("step=",step+1,"h=",signif(hakt,2),"lambda1=",signif(lambda1[step+1],3),"lambda2=",signif(lambda2[step+1],3),"\n")
step<-step+1
hakt <- h0*ch^step
hakt1 <- h1*ch1^step
gakt <- h0*cg^step
etaaktinv <- sqrt(hakt^2-gakt^2)/(hakt*gakt)
haktinv <- 1/hakt
if(demo) tmp <- readline("Press enter for next step")
}
gamma <- sweep(gamma,2:3,ngamma,"/")
list(gamma=gamma,ahat=matrix(z1$ahat,dy[1],dy[2]),chat=matrix(z1$chat,dy[1],dy[2]),Ni=Ni)
}
#
#    R - function  aws  for estimating the flow direction in an image
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2008 Weierstrass-Institut fuer
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
IQRdiff <- function(y) IQR(diff(y))/1.908
awsflow <- function(y,kstar,h0=1.12,ch=sqrt(1.25),cg=1,lambda=10,graph=TRUE,adapt=TRUE,
           u=NULL,ugamma=NULL,demo=FALSE,sel1=1,sel2=1){
args <- match.call()
dy<-dim(y)
n <- prod(dy)
# Initialization
z1 <- .Fortran("initflow",
                as.double(y),
                as.integer(dy[1]),
                as.integer(dy[2]),
                as.double(h0),
                ahat=double(n),
                chat=double(n),
                gamma=double(2*n),
                Ni=double(n),
                Pvec=double(2*n),
                Qmat=double(3*n),
                PACKAGE="aws")[c("ahat","chat","gamma","Ni","Pvec","Qmat")]
Ni <- matrix(z1$Ni,dy[1],dy[2])
Qmat <- array(z1$Qmat,c(3,dy))
Pvec <- array(z1$Pvec,c(2,dy))
gamma <- array(z1$gamma,c(2,dy))
#gamma <- .Fortran("ingamma2",
#                  as.double(y),
#                  as.integer(dy[1]),
#                  as.integer(dy[2]),
#                  gamma=double(2*n),
#                  PACKAGE="aws")$gamma
#dim(gamma) <- c(2,dy)
ngamma <- sqrt(gamma[1,,]^2+gamma[2,,]^2)
if(graph) {
par(mfrow=c(2,2),mar=c(3,3,3,.1),mgp=c(2,1,0))
image(atan(gamma[1,,]/gamma[2,,])+pi*(gamma[2,,]<0),col=rainbow(512))
title(paste("Orientation h=",signif(h0,3)))
image(Ni, col=grey(0:255/255))
title(paste("Ni",signif(min(Ni),3),signif(mean(Ni),3),signif(max(Ni),3)))
#image(matrix(z1$ahat,dy[1],dy[2]), col=grey(0:255/255),zlim=quantile(z1$ahat,c(.01,.99)))
#title(paste("ahat",signif(quantile(z1$ahat,.01),3),signif(quantile(z1$ahat,.99),3)))
#image(matrix(z1$chat,dy[1],dy[2]), col=grey(0:255/255),zlim=quantile(z1$chat,c(.01,.99)))
#title(paste("chat",signif(quantile(z1$chat,.01),3),signif(quantile(z1$chat,.99),3)))
image(matrix(z1$ahat,dy[1],dy[2]), col=grey(0:255/255))
title(paste("ahat",signif(min(z1$ahat),3),signif(max(z1$ahat),3)))
image(matrix(z1$chat,dy[1],dy[2]), col=grey(0:255/255))
title(paste("chat",signif(max(z1$chat),3)))
}
if(!is.null(u)) cat("h0=",signif(h0,2),"MAE",signif(mean(abs(z1$ahat-u)),3),"MSE",signif(mean((z1$ahat-u)^2),3),"\n")
if(!is.null(ugamma)) cat("h0=",signif(h0,2),"MSE(gamma)",signif(mean(1-abs(gamma[1,,]*ugamma[1,,]+gamma[2,,]*ugamma[2,,])/ngamma),3),"\n")
if(demo) tmp <- readline("Press enter for next step")
#Pvec <- double(2*n)
#Qmat <- double(3*n)
#Ni <- rep(1,n)
sigma2 <- min(IQRdiff(as.vector(y))^2,IQRdiff(as.vector(t(y)))^2)
cat("estimated sigma",signif(sqrt(sigma2),2),"\n")
# Iterate
ll <- length(lambda)
if(ll<kstar+1) lambda <- c(lambda,rep(lambda[ll],kstar-ll+1))
step <- 0
hakt <- h0*ch^step
etaakt <- h0^2*ch^step/(h0^2*(ch*cg)^step-1)
while(step<kstar+1){
#  step 1
z1 <- .Fortran("itfstep1",
              as.double(y),
              as.integer(dy[1]),
              as.integer(dy[2]),
              as.double(gamma),
              as.double(Ni),
              as.double(Qmat),
              as.double(Pvec),
              as.double(hakt),
              as.double(etaakt),
              as.double(sigma2),
              as.double(lambda[step+1]),
              Bmat=double(3*n),
              bvec=double(2*n),
              Svec=double(2*n),
              Ni=double(n),
              chat=double(n),
              ahat=double(n),
              as.double(ngamma),
              as.double(h0*h0),
              as.integer(sel1),# do not use adapive weights if distance is less than h0
              PACKAGE="aws")[c("Bmat","bvec","Svec","Ni","chat","ahat")]
#  step 2 
z1$ahat[z1$Ni==0] <- y[z1$Ni==0]
z1$Ni[z1$Ni==0] <- 1
z2 <- .Fortran("itfstep2",
              as.double(gamma),
              as.integer(dy[1]),
              as.integer(dy[2]),
              as.double(hakt),
              as.double(sigma2),
              as.double(lambda[step+1]),
              Pvec=as.double(Pvec),
              Qmat=as.double(Qmat),
              as.double(Ni),
              as.double(z1$chat),
              as.double(z1$ahat),
              as.double(z1$Bmat),
              as.double(z1$bvec),
              as.double(z1$Svec),
              gamma=double(2*n),
              ngamma=as.double(ngamma),
              as.double(h0*h0),# do not use adapive weights if distance is less than h0
              PACKAGE="aws")[c("Qmat","Pvec","gamma","ngamma")]
#  organize everything for the next step
Ni <- matrix(z1$Ni,dy[1],dy[2])
Qmat <- array(z2$Qmat,c(3,dy))
Pvec <- array(z2$Pvec,c(2,dy))
gamma <- array(z2$gamma,c(2,dy))
ngamma <- matrix(z2$ngamma,dy[1],dy[2])
if(graph) {
par(mfrow=c(2,2),mar=c(3,3,3,.1),mgp=c(2,1,0))
image(atan(gamma[1,,]/gamma[2,,])+pi*(gamma[2,,]<0),col=rainbow(512))
title(paste("Orientation h=",signif(hakt,3),"eta=",signif(etaakt,3)))
image(matrix(Ni,dy[1],dy[2]), col=grey(0:255/255))
title(paste("Ni",signif(min(Ni),3),signif(mean(Ni),3),signif(max(Ni),3)))
#image(matrix(z1$ahat,dy[1],dy[2]), col=grey(0:255/255),zlim=quantile(z1$ahat,c(.01,.99)))
#title(paste("ahat",signif(quantile(z1$ahat,.01),3),signif(quantile(z1$ahat,.99),3)))
#image(matrix(z1$chat,dy[1],dy[2]), col=grey(0:255/255),zlim=quantile(z1$chat,c(.01,.99)))
#title(paste("chat",signif(quantile(z1$chat,.01),3),signif(quantile(z1$chat,.99),3)))
image(matrix(z1$ahat,dy[1],dy[2]), col=grey(0:255/255))
title(paste("ahat",signif(min(z1$ahat),3),signif(max(z1$ahat),3)))
image(matrix(z1$chat,dy[1],dy[2]), col=grey(0:255/255))
title(paste("chat",signif(max(z1$chat),3)))
}
if(!is.null(u)) cat("h=",signif(hakt,2),"g=",signif(h0*cg^step,2),"MAE",signif(mean(abs(z1$ahat-u)),3),"MSE",signif(mean((z1$ahat-u)^2),3),"\n")
if(!is.null(ugamma)) cat("h=",signif(hakt,2),"g=",signif(h0*cg^step,2),"MSE(gamma)",signif(mean(1-abs(gamma[1,,]*ugamma[1,,]+gamma[2,,]*ugamma[2,,])/ngamma),3),"\n")
cat("step=",step+1,"h=",signif(hakt,2),"lambda=",signif(lambda[step+1],3),"\n")
step<-step+1
hakt <- h0*ch^step
etaakt <- h0^2*ch^step/(h0^2*(ch*cg)^step-1)
if(demo) tmp <- readline("Press enter for next step")
}
gamma <- sweep(gamma,2:3,ngamma,"/")
list(gamma=gamma,ahat=matrix(z1$ahat,dy[1],dy[2]),chat=matrix(z1$chat,dy[1],dy[2]),Ni=Ni)
}