#
#    R - function  sagg  
#    
#    this function provides a unified iterface to the 
#    stagewise aggregation for a variety of models
#    and problems
#    depending on the parameters specified the following functions are 
#    called internally
#
#    model=="Regression", homoskedastic=TRUE,    p = 0                         psawslike
#    model=="Regression", family=="Gaussian", homoskedastic=FALSE,   p = 0     psawshlike
#    model=="Regression", family=="Gaussian", homoskedastic=TRUE,    p > 0     psawspoly
#    model=="Regression", family=="Gaussian", homoskedastic=FALSE,   p > 0     psawshpoly
#
#    emaphazises on the propagation-separation approach.
#
#    Parameter defaults that are difficult to access, like 
#
#           qlambda, tau1, tau2, heta, eta0, hinit, hincr
# 
#    can only be specified when the individual functions are called. 
#    
#
#    Copyright (C) 2002 Weierstrass-Institut fuer
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
sagg<-function(y,x=NULL,p=0,hmax=NULL,family="Gaussian",homoskedastic=TRUE,sigma2=NULL,shape=NULL,
                lkern="Triangle",NN=FALSE,u=NULL,graph=FALSE,demo=FALSE,control=TRUE,
                wghts=NULL,heta=1,tau1=1,tau2=1){
ergs<-"not implemented"
if(control) heta<-NULL else heta<-1e10 
   if(p==0){
         if(is.null(dim(y))) {
         heta<-switch(family,"Gaussian"=2,heta)
         tau1<-switch(family,"Gaussian"=.8,.8)
         } else {
         heta<-switch(family,"Gaussian"=1.5,"Poisson"=2,heta)
         tau1<-switch(family,"Gaussian"=.6,"Poisson"=2.4,.6)
       }
      if(homoskedastic) {
         ergs<-caws(y,x,lkern=lkern,family=family,sigma2=sigma2,qlambda=1,
                    shape=shape,hmax=hmax,heta=heta,tau1=tau1,tau2=tau1/2,NN=NN,u=u,graph=graph,
                    demo=demo,wghts=wghts)
         } else {
         ergs<-chaws(y,x,lkern=lkern,sigma2=sigma2,qlambda=1,
                     shape=shape,hmax=hmax,heta=heta,tau1=tau1,tau2=tau1/2,NN=NN,u=u,graph=graph,
                     demo=demo,wghts=wghts)	 
         }
      } else {
         if(is.null(dim(y))) {
            heta<-p+1.5
          tau1<-switch(p,10,50,250)
         } else {
         heta<-p+1.5
         tau1<-switch(p,3.2,10)
         }
      if(homoskedastic) {
         ergs<-cpaws(y,x,p=p,lkern=lkern,sigma2=sigma2,qlambda=1,
                     hmax=hmax,heta=heta,tau1=tau1,tau2=tau2,NN=NN,u=u,graph=graph,
                     demo=demo,wghts=wghts)
            } else {
         ergs<-cphaws(y,x,p=p,lkern=lkern,sigma2=sigma2,qlambda=1,
                      hmax=hmax,heta=heta,tau1=tau1,tau2=tau2,NN=NN,u=u,graph=graph,
                      demo=demo,wghts=wghts)	 
            }
      } 
ergs
}
