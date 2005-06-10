#
#    R - function  psaws  for regression
#    
#    this function provides a unified iterface to the 
#    propagation-separation approach to AWS for a variety of models
#    and problems
#    depending on the parameters specified the following functions are 
#    called internally
#
#    homoskedastic=TRUE,    p = 0                         psawslike
#    family=="Gaussian", homoskedastic=FALSE,   p = 0     psawshlike
#    family=="Gaussian", homoskedastic=TRUE,    p > 0     psawspoly
#    family=="Gaussian", homoskedastic=FALSE,   p > 0     psawshpoly
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
psaws<-function(y,x=NULL,p=0,hmax=NULL,model="Regression",family="Gaussian",
                homoskedastic=TRUE,sigma2=NULL,shape=NULL,lkern="Triangle",NN=FALSE,
		u=NULL,graph=FALSE,demo=FALSE,aws=TRUE,control=TRUE,wghts=NULL){
ergs<-"not implemented"
if(aws) qlambda<-NULL else qlambda<-1
if(control) heta<-NULL else heta<-1e10 
   if(p==0){
      if(homoskedastic) {
         ergs<-caws(y,x,lkern=lkern,family=family,sigma2=sigma2,qlambda=qlambda,
                    heta=heta,shape=shape,hmax=hmax,NN=NN,u=u,graph=graph,
	            demo=demo,wghts=wghts,demo=demo)
	 } else {
         ergs<-chaws(y,x,lkern=lkern,sigma2=sigma2,qlambda=qlambda,
                     heta=heta,shape=shape,hmax=hmax,NN=NN,u=u,graph=graph,
	             demo=demo,wghts=wghts,demo=demo)	 
	 }
      } else {
      if(homoskedastic) {
         ergs<-cpaws(y,x,p=p,lkern=lkern,sigma2=sigma2,qlambda=qlambda,
                     heta=heta,hmax=hmax,NN=NN,u=u,graph=graph,
	             demo=demo,wghts=wghts,demo=demo)
	 } else {
         ergs<-cphaws(y,x,p=p,lkern=lkern,sigma2=sigma2,qlambda=qlambda,
                      heta=heta,hmax=hmax,NN=NN,u=u,graph=graph,
	              demo=demo,wghts=wghts,demo=demo)	 
	 }
      }
ergs
}
