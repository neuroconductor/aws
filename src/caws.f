C
C    Copyright (C) 2004 Weierstrass-Institut fuer 
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines are part of the aws package and contain  
C  FORTRAN 77 code needed in R functions aws, vaws, 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    location penalty for multivariate non-gridded aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lmkern(kern,dx,xi,xj,h2)
      implicit logical (a-z)
      external lkern
      integer kern,dx,i
      real*8 xi(dx),xj(dx),h2,z,zd,lkern
      z=0.d0
      do 1 i=1,dx
         zd=xi(i)-xj(i)
         z=z+zd*zd
	 if(z.gt.h2) goto 2
1     continue
      lmkern=lkern(kern,z/h2)
      goto 999
2     lmkern=0.d0
999   return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws on a grid      
C
C   this is a reimplementation of the original aws procedure
C
C   should be slightly slower for non-Gaussian models (see function kldist)
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Model=1    Gaussian   
C          Model=2    Bernoulli   
C          Model=3    Poisson   
C          Model=4    Exponential   
C
C     computing dlog(theta) and dlog(1.d0-theta) outside the AWS-loops 
C     will reduces computational costs at the price of readability
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldist(model,thi,thj,bi0)
      implicit logical (a-z)
      integer model
      real*8 thi,thj,z,thij,bi0,eta,tthi
      IF (model.eq.1) THEN
C        Gaussian
         z=thi-thj
         kldist=z*z
      ELSE IF (model.eq.2) THEN
C        Bernoulli
         kldist=0.d0
         eta=0.5d0/bi0
         thij=(1.d0-eta)*thj+eta*thi
         tthi=(1.d0-thi)
         IF (thi.gt.1.d-10) kldist=kldist+thi*dlog(thi/thij)
         IF (tthi.gt.1.d-10) kldist=kldist+tthi*dlog(tthi/(1.d0-thij))
      ELSE IF (model.eq.3) THEN
C        Poisson
         kldist=0.d0
         eta=0.5d0/bi0
         thij=(1.d0-eta)*thj+eta*thi
         IF (thi.gt.1.d-10) kldist=thi*dlog(thi/thij)-thi+thij
      ELSE IF (model.eq.4) THEN
C        Exponential
         kldist=thi/thj-1.d0-dlog(thi/thj)
      ELSE IF (model.eq.5) THEN
C        Exponential
         kldist=thi/thj-1.d0-dlog(thi/thj)
      ELSE
C        use Gaussian
         z=thi-thj
         kldist=z*z
      ENDIF
      RETURN
      END
      real*8 function kldist0(model,thi,thj)
      implicit logical (a-z)
      integer model
      real*8 thi,thj,z
      IF (model.eq.1) THEN
C        Gaussian
         z=thi-thj
         kldist0=z*z
      ELSE IF (model.eq.2) THEN
C        Bernoulli
         kldist0=thi*dlog(thi/thj)+(1.d0-thi)*dlog((1.d0-thi)/
     1           (1.d0-thj))
      ELSE IF (model.eq.3) THEN
C        Poisson
         kldist0=thi*dlog(thi/thj)-thi+thj
      ELSE IF (model.eq.4) THEN
C        Exponential
         kldist0=thi/thj-1.d0-dlog(thi/thj)
      ELSE
C        use Gaussian
         z=thi-thj
         kldist0=z*z
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         lkern=1.d0
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=dexp(-xsq*8.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        Compute truncated Exponential Kernel 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function skern(x,xmax)
      implicit logical (a-z)
      real*8 x,xmax
      IF (x.gt.xmax) THEN
         skern=0.d0
      ELSE
         skern=dexp(-x)
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws0(y,fix,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                   bi0,ai,model,kern,spmax,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,
     1        iind,jind,jind3,jind2
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               ja3=max0(1,i3-ih3)
               je3=min0(n3,i3+ih3)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               swjy=0.d0
               DO j3=ja3,je3
                  z3=(i3-j3)*wght(2)
                  z3=z3*z3
                  ih2=dsqrt(hakt2-z3)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
		  jind3=(j3-1)*n1*n2
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     ja1=max0(1,i1-ih1)
                     je1=min0(n1,i1+ih1)
		     jind2=jind3+(j2-1)*n1
                     DO j1=ja1,je1
C  first stochastic term
                        jind=j1+jind2
                        z1=(i1-j1)
                        wj=lkern(kern,(z1*z1+z2)/hakt2)
                        swj0=swj0+wj
                        IF (aws) THEN
                  sij=bii*kldist(model,thetai,theta(jind),bii0)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj
                        swjy=swjy+wj*y(jind)
                     END DO
                  END DO
               END DO
               ai(iind)=swjy
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws(y,fix,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                   bi0,ai,model,kern,spmax,lwght,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt,lwght(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      DO j3=1,dlw3
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=dsqrt(hakt2-z3)/wght(1)
         jind3=(j3-1)*dlw1*dlw2
         DO j2=clw2-ih2,clw2+ih2
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=dsqrt(hakt2-z2)
            jind2=jind3+(j2-1)*dlw1
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               swjy=0.d0
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
		  jwind3=(jw3-1)*dlw1*dlw2
	          jind3=(j3-1)*n1*n2
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=dsqrt(hakt2-z3)/wght(1)
                  DO jw2=clw2-ih2,clw2+ih2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
	             jind2=(j2-1)*n1+jind3
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jind=j1+jind2
                        wj=lwght(jw1+jwind2)
                        swj0=swj0+wj
                        IF (aws) THEN
                  sij=bii*kldist(model,thetai,theta(jind),bii0)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj
                        swjy=swjy+wj*y(jind)
                     END DO
                  END DO
               END DO
               ai(iind)=swjy
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws(y,fix,si2,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                   bi0,ai,model,kern,spmax,lwght,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt,lwght(1),si2(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      DO j3=1,dlw3
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         if(n2.gt.1) ih2=dsqrt(hakt2-z3)/wght(1)
         jind3=(j3-1)*dlw1*dlw2
         DO j2=clw2-ih2,clw2+ih2
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=dsqrt(hakt2-z2)
            jind2=jind3+(j2-1)*dlw1
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               swjy=0.d0
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
	          jind3=(j3-1)*n1*n2
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=dsqrt(hakt2-z3)/wght(1)
		  jwind3=(jw3-1)*dlw1*dlw2
                  DO jw2=clw2-ih2,clw2+ih2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
	             jind2=(j2-1)*n1+jind3
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
		     jwind2=jwind3+(jw2-1)*dlw1
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jind=j1+jind2
                        wj=lwght(jw1+jwind2)*si2(jind)
                        swj0=swj0+wj
                        IF (aws) THEN
                  sij=bii*kldist(model,thetai,theta(jind),bii0)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj
                        swjy=swjy+wj*y(jind)
                     END DO
                  END DO
               END DO
               ai(iind)=swjy
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws0(y,fix,si2,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                 bi0,ai,model,kern,spmax,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   bi2      \sum  Wi^2*sigma2   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt,si2(1),hakt2
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,
     1        iind,jind,jind3,jind2
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,wj0,bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               ja3=max0(1,i3-ih3)
               je3=min0(n3,i3+ih3)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               swjy=0.d0
               DO j3=ja3,je3
                  z3=(i3-j3)*wght(2)
                  z3=z3*z3
                  ih2=dsqrt(hakt2-z3)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
		  jind3=(j3-1)*n1*n2
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     ja1=max0(1,i1-ih1)
                     je1=min0(n1,i1+ih1)
		     jind2=jind3+(j2-1)*n1
                     DO j1=ja1,je1
C  first stochastic term
                        jind=j1+jind2
                        z1=(i1-j1)
                        wj0=lkern(kern,(z1*z1+z2)/hakt2)
			wj=wj0*si2(jind)
                        swj0=swj0+wj
                        IF (aws) THEN
                   sij=bii*kldist(model,thetai,theta(jind),bii0)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj0
                        swjy=swjy+wj*y(jind)
                     END DO
                  END DO
               END DO
               ai(iind)=swjy
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance (for vaws)
C
C          Model=1    Gaussian   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function vkldist(model,d,thi,thj,bi0,vw)
      implicit logical (a-z)
      integer model,d,i
      real*8 thi(d),thj(d),z,z0,bi0,vw(d)
C        Gaussian case only in the moment
         z=0.d0
	 DO i=1,d
            z0=thi(i)-thj(i)
	    z=z+z0*z0*vw(i)
	 END DO
         vkldist=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vaws(y,fix,d,n1,n2,n3,ns,hakt,lambda,theta,bi,bi2,
     1                 bi0,ai,kern,spmax,wght,thetai,swjy,vw)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   ns       length of sigma
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   bi2      \sum  Wi^2*sigma2   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external vkldist,lkern
      real*8 vkldist,lkern
      integer d,n1,n2,n3,model,kern,ns
      logical aws,fix(1),hetero
      real*8 y(d,1),theta(d,1),bi(1),bi0(1),ai(d,1),lambda(1),spmax,
     1       wght(2),bi2(1),hakt,hakt2,thetai(d),swjy(d),vw(d)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,
     1        iind,jind,jind3,jind2,k,indis
      real*8 bii,sij,swj,swj2,swj0,z1,z2,z3,wj,bii0
      model=1
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda(1).lt.1d40
      hetero=ns.gt.1
      indis=1
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
	       if(hetero) indis=iind
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               bii=bi(iind)/lambda(indis)
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               ja3=max0(1,i3-ih3)
               je3=min0(n3,i3+ih3)
               DO k=1,d
                  thetai(k)=theta(k,iind)
		  swjy(k)=0.d0
	       END DO
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               DO j3=ja3,je3
                  z3=(i3-j3)*wght(2)
                  z3=z3*z3
                  ih2=dsqrt(hakt2-z3)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
		  jind3=(j3-1)*n1*n2
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     ja1=max0(1,i1-ih1)
                     je1=min0(n1,i1+ih1)
		     jind2=jind3+(j2-1)*n1
                     DO j1=ja1,je1
C  first stochastic term
                        jind=j1+jind2
                        z1=(i1-j1)
                        wj=lkern(kern,(z1*z1+z2)/hakt2)
                        swj0=swj0+wj
                        IF (aws) THEN
                        sij=bii*vkldist(model,d,thetai,theta(1,jind),
     1                                  bii0,vw)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj
			DO k=1,d
                           swjy(k)=swjy(k)+wj*y(k,jind)
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,d
                  ai(k,iind)=swjy(k)
	       END DO
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vhaws(y,fix,si2,d,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                 bi0,ai,kern,spmax,wght,thetai,swjy,vw)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   bi2      \sum  Wi^2*sigma2   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external vkldist,lkern
      real*8 vkldist,lkern
      integer d,n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(d,1),theta(d,1),bi(1),bi0(1),ai(d,1),lambda,spmax,
     1       wght(2),bi2(1),hakt,si2(1),hakt2,thetai(d),swjy(d),vw(d)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,
     1        iind,jind,jind3,jind2,k
      real*8 bii,sij,swj,swj2,swj0,z1,z2,z3,wj,wj0,bii0,esij
      model=1
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               ja3=max0(1,i3-ih3)
               je3=min0(n3,i3+ih3)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               DO k=1,d
                  thetai(k)=theta(k,iind)
		  swjy(k)=0.d0
	       END DO
               DO j3=ja3,je3
                  z3=(i3-j3)*wght(2)
                  z3=z3*z3
                  ih2=dsqrt(hakt2-z3)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
		  jind3=(j3-1)*n1*n2
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     ja1=max0(1,i1-ih1)
                     je1=min0(n1,i1+ih1)
		     jind2=jind3+(j2-1)*n1
                     DO j1=ja1,je1
C  first stochastic term
                        jind=j1+jind2
                        z1=(i1-j1)
                        wj0=lkern(kern,(z1*z1+z2)/hakt2)
			wj=wj0*si2(jind)
                        swj0=swj0+wj
                        IF (aws) THEN
                        sij=bii*vkldist(model,d,thetai,theta(1,jind),
     1                                  bii0,vw)
                           IF (sij.gt.spmax) CYCLE
			   esij=exp(-sij)
                           wj=wj*esij
                           wj0=wj0*esij
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj0
			DO k=1,d
                           swjy(k)=swjy(k)+wj*y(k,jind)
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,d
                  ai(k,iind)=swjy(k)
	       END DO
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cpawsuni(n,dp1,dp2,x,y,fix,theta,bi,bi0,
     1        ai,lam,h,kern,cb,dmat,thij,psix,psiy,bii,spmax)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy                working memory 
C     
      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,l,je,ja,nwij,kern 
      logical aws,fix(n)
      real*8 lkern,x(n),y(n),psix(dp2),psiy(dp1),theta(dp1,n),
     1 bi(dp2,n),ai(dp1,n),lam,spmax,dmat(dp1,dp1),bii(dp2),
     2 thij(dp1),cb(dp1,dp1),bi0(dp2,n),lambda,h,
     3 sij,wij,z,xij,ha,ha2,eps,xi,thijl
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=1.e-10
      aws=lam.lt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ha=h
         xi=x(i)
C     first search for ja and je
1099     DO j=i,1,-1
            IF (x(j).le.xi-ha) EXIT
            ja=j
         END DO
         DO j=i,n
            IF (x(j).gt.xi+ha) EXIT
            je=j
         END DO
         IF (je-ja-1.le.dp1) THEN
            ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
            goto 1099
         END IF
         ha2=ha*ha
8999     nwij=0
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,dp1
                  dmat(j,k)=bi(j+k-1,i)/lambda
	       END DO
            END DO
         END IF
         DO l=1,dp2
            bii(l)=0.d0
            bi0(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l,i)=0.d0
         END DO
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased 
C      (to reduce panalization by stochastic and influence term)
C
C             loop over local neighborhood
C
         DO j=ja,je
C
C              get location weights
C
            xij=(x(j)-x(i))
            wij=lkern(kern,xij*xij/ha2)
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               IF (k.le.dp1) psiy(k)=z*y(j)
               z=z*xij
            END DO
            DO k=1,dp2
               bi0(k,i)=bi0(k,i)+wij*psix(k)
	    END DO
            IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
               DO l=1,dp1
                 thij(l)=theta(l,j)
               END DO
               IF (dp1.gt.1) THEN
                  z=1.d0
                  DO l=2,dp1
                     z=-z*xij
                     DO k=1,dp1-l+1
                        thij(k)=thij(k)+cb(k+l-1,k)*z*theta(k+l-1,j)
                     END DO
                  END DO
               END IF
C     thij contains theta_{j} in the model centered at xi
C
C                 now get sij
C
	       DO l=1,dp1
                  thij(l)=theta(l,i)-thij(l)
               END DO
C
C   thats the difference between thetai and thetaij
C
               sij=0.0d0
               DO l=1,dp1
	          thijl=thij(l)
		  sij=sij+dmat(l,l)*thijl*thijl
		  IF (l.eq.dp1) CYCLE
                  DO k=l+1,dp1
                     sij=sij+2.d0*dmat(k,l)*thijl*thij(k)
                  END DO
               END DO
C
C     now we have everything to compute  w_{ij}
C       
               IF (sij.gt.spmax) CYCLE
               wij=wij*dexp(-sij)
            END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),bi0(i),ai(i)  
            DO k=1,dp2
               bii(k)=bii(k)+wij*psix(k)
               IF (k.le.dp1) ai(k,i)=ai(k,i)+wij*psiy(k)
            END DO    
         END DO
C    this was the j - loop
         IF (nwij.lt.dp1.and.aws) THEN
            lambda=lambda*1.25
C    increase lambda to weaken stochastic penalization
            goto 8999
	 END IF
	 DO k=1,dp2 
	    bi(k,i)=bii(k)
         END DO	 
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        search for sufficiently long interval to fit a GLM for Poisson
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine poisdes(y,n,p,h)
      integer n,p,h(2,n)
      real*8 y(n)
      integer i,j,k
      real*8 z
      DO i=1,n
         k=0
         h(1,i)=i
         h(2,i)=i
         IF(y(i).gt.0.d0) k=k+1
         DO j=1,n
            IF(k.gt.p) CYCLE
            z=0.d0
            IF(i.gt.j) z=y(i-j)
            IF(z.gt.0.d0) THEN
               h(1,i)=i-j
               k=k+1
            ENDIF
            z=0.d0
            IF(i+j.le.n) z=y(i+j)
            IF(z.gt.0.d0) THEN
               h(2,i)=i+j
               k=k+1
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpawsuni(n,dp1,dp2,dp3,y,xd,fix,mfamily,theta0,theta,
     1                    bii,bi,bi2,bi0,bi02,ai,ni,lam,h,h0,kern,cb,
     2                    dmat,thij,psix,wghts,wghts0,spmax,iter,smw,
     3                    work)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     y          observed values at design points
C     xd         gridlength
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy                working memory 
C     
      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,l,je,ja,kern,mfamily,h0(2,n),ih,iter,
     1        dp3,info
      logical aws,fix(n),smw
      real*8 y(n),xd,psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),lkern,thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n)
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ih=h
C         ja=min0(h0(1,i),max0(1,i-ih))
C         je=max0(h0(2,i),min0(n,i+ih))
         ja=max0(1,h0(1,i)-ih)
         je=min0(n,h0(2,i)+ih)
         ha=max0(i-ja,je-i)
         ha=dmax1(h,ha+h-ih)*xd
         ha2=ha*ha
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,j
                  d=bii(l,i)/lambda
                  dmat(k,j)=d
                  dmat(j,k)=d
                  l=l+1
	       END DO
            END DO
         END IF
C
C          Now get weights in wghts and wghts0
C
         call lpuniwgt(n,i,ja,je,dp1,dp2,xd,ha2,h0(1,i),aws,kern,
     1                theta0,spmax,psix,thij,dmat,cb,wghts,wghts0)
         IF (smw) THEN
            call smwghtun(n,i,ja,je,wghts,work,h0(1,i),kern)
C            call smwghtun(n,i,ja,je,wghts0,work,h0(1,i),kern)
         END IF
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         call lpuniit(mfamily,n,dp1,dp2,i,ja,je,y,xd,theta(1,i),psix,
     1                dmat,wghts,wghts0,thij,bi(1,i),ai,bi0(1,i),
     2                bi2(1,i),bi02(1,i),iter,info)
         if(info.gt.0) fix(i)=.TRUE.
C   singularity or failed convergence keep old estimate in this case
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpawsunw(n,dp1,dp2,dp3,y,xd,fix,mfamily,theta0,theta,
     1                    bii,bi,bi2,bi0,bi02,ai,ni,lam,h,h0,kern,cb,
     2                    dmat,thij,psix,wghts,wghts0,spmax,iter,smw,
     3                    work)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     y          observed values at design points
C     xd         gridlength
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy                working memory 
C     
      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,l,je,ja,kern,mfamily,h0(2,n),ih,iter,
     1        dp3,info
      logical aws,fix(n),smw
      real*8 y(n),xd,psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),lkern,thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n,n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n)
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ih=h
         ja=min0(h0(1,i),max0(1,i-ih))
         je=max0(h0(2,i),min0(n,i+ih))
         ha=max0(i-ja,je-i)
         ha=dmax1(h,ha)*xd
         ha2=ha*ha
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,j
                  d=bii(l,i)/lambda
                  dmat(k,j)=d
                  dmat(j,k)=d
                  l=l+1
	       END DO
            END DO
         END IF
C
C          Now get weights in wghts and wghts0
C
         call lpuniwgt(n,i,ja,je,dp1,dp2,xd,ha2,h0(1,i),aws,kern,
     1              theta0,spmax,psix,thij,dmat,cb,wghts(1,i),wghts0)
         IF (smw) THEN
            call smwghtun(n,i,ja,je,wghts(1,i),work,h0(1,i),kern)
C            call smwghtun(n,i,ja,je,wghts0,work,h0(1,i),kern)
         END IF
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j,i)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         call lpuniit(mfamily,n,dp1,dp2,i,ja,je,y,xd,theta(1,i),psix,
     1                dmat,wghts(1,i),wghts0,thij,bi(1,i),ai,bi0(1,i),
     2                bi2(1,i),bi02(1,i),iter,info)
         if(info.gt.0) fix(i)=.TRUE.
C   singularity or failed convergence keep old estimate in this case
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    smooth univariate weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghtu0(n,i,ja,je,wghts,work,h,kern)
      integer n,kern,i,ja,je,h(2)
      real*8 wghts(n),work(n)
      integer j,k,ka,ke,ih,kj
      real*8 ha,maxwght,lkern,x,d
      external lkern
      ih=max0(i-h(1),h(2)-i)
      ha=ih+0.5d0
      ha=ha*ha
      DO j=ja,je
         work(j)=0.d0
      END DO
      DO k=-ih,ih
         x=lkern(kern,k*k/ha)
         DO j=ja,je
            kj=k+j
            if(kj.lt.ja.or.kj.gt.je) CYCLE
            work(j)=work(j)+x*wghts(kj)
         END DO
      END DO
      maxwght=0.d0
      DO j=ja,je
         if(work(j).gt.maxwght) maxwght=work(j)
      END DO
      DO j=ja,je
         wghts(j)=work(j)/maxwght
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    smooth univariate weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghtu1(n,i,ja,je,wghts,work,h,kern)
      integer n,kern,i,ja,je,h(2)
      real*8 wghts(n),work(n)
      integer j,k,ka,ke,ih
      real*8 ha,maxwght,lkern,x
      external lkern
      ih=max0(i-h(1),h(2)-i)
      DO j=ja,je
         ka=max0(ja,j-ih)
         ke=min0(je,j+ih)
         work(j)=0.d0
         DO k=ka,ke
            work(j)=work(j)+wghts(k)
         END DO
      END DO
      maxwght=work(i)
C      DO j=ja,je
C         if(work(j).gt.maxwght) maxwght=work(j)
C      END DO
      DO j=ja,je
         wghts(j)=dmax1(1.d0,work(j)/maxwght)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    smooth univariate weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghtun(n,i,ja,je,wghts,work,h,kern)
      integer n,kern,i,ja,je,h(2)
      real*8 wghts(n),work(n)
      integer j,k,ka,ke,ih
      real*8 ha,maxwght,lkern,d
      external lkern
      ih=max0(i-h(1),h(2)-i)
      d=0.d0
      DO k=ja,ja+ih
         if(k.le.je) d=d+wghts(k)
      END DO
      DO j=ja,je
         work(j)=d
         ka=j-ih
         ke=j+ih
         if(ke.le.je) d=d+wghts(ke)
         if(ka.ge.ja) d=d-wghts(ka)
      END DO
      maxwght=work(i)
C      DO j=ja,je
C         if(work(j).gt.maxwght) maxwght=work(j)
C      END DO
      DO j=ja,je
         wghts(j)=dmin1(1.d0,work(j)/maxwght)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    one iteration for local GLM estimates
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpuniit(mfamily,n,dp1,dp2,i,ja,je,y,xd,theta,psix,
     1           dmat,wghts,wghts0,d,bi,ai,bi0,bi2,bi02,iter,info)
      implicit logical (a-z)
      integer n,dp1,dp2,mfamily,i,ja,je,iter,info
      real*8 y(n),theta(dp1),psix(dp2),bi(dp2),ai(dp1),bi0(dp2),
     1       bi2(dp2),bi02(dp2),d(dp1),wghts(n),wghts0(n),xd,
     2       dmat(dp1,dp1)
      integer it,j,k,l
      real*8 z,dist,xij
      DO it=1,iter
         DO l=1,dp2
            bi(l)=0.d0
            bi0(l)=0.d0
            bi2(l)=0.d0
            bi02(l)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l)=0.d0
         END DO
         DO j=ja,je
            xij=(j-i)*xd
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               z=z*xij
            END DO
C        now compute contributions to bi(i),bi0(i),ai(i)  
            call guniaibi(mfamily,dp1,dp2,psix,theta,y(j),bi,
     1           ai,bi0,bi2,bi02,wghts(j),wghts0(j))
         END DO
         call mpawsun0(dp1,dp2,ai,bi,dmat,d,info)
         IF(info.gt.0) goto 999
         dist=0.d0
         DO k=1,dp1
            DO l=1,dp1
               dist=dist+d(k)*bi(k+l-1)*d(l)
            END DO
         END DO
         DO k=1,dp1
            theta(k)=theta(k)+d(k)
         END DO
         IF(dist.lt.1.d-12) goto 999
         info=10
      END DO
999   RETURN
      END
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    compute weights for local GLM estimates
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpuniwgt(n,i,ja,je,dp1,dp2,xd,ha2,h0,aws,kern,theta,
     1                    spmax,psix,thij,dmat,cb,wghts,wghts0)
      implicit logical (a-z)
      integer n,i,ja,je,dp2,dp1,kern,h0(2)
      real*8 xd,ha2,wghts(n),wghts0(n),psix(dp2),theta(dp1,n),
     1       thij(dp1),dmat(dp1,dp1),cb(dp1,dp1),spmax
      logical aws
      integer j,k,l
      real*8 xij,z,sij,thijl,wij,wij0,lkern
      external lkern
      DO j=ja,je
         xij=(j-i)*xd
         wij0=lkern(kern,xij*xij/ha2)
         wghts0(j)=wij0
         z=1.d0
         DO k=1,dp2
            psix(k)=z
            z=z*xij
         END DO
C         IF ((j.lt.h0(1).or.j.gt.h0(2)).and.aws) THEN
         IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            DO l=1,dp1
               thij(l)=theta(l,j)
            END DO
            IF (dp1.gt.1) THEN
               z=1.d0
               DO l=2,dp1
                  z=-z*xij
                  DO k=1,dp1-l+1
                     thij(k)=thij(k)+cb(k+l-1,k)*z*theta(k+l-1,j)
                  END DO
               END DO
            END IF
C     thij contains theta_{j} in the model centered at xi
C
C                 now get sij
C
	    DO l=1,dp1
               thij(l)=theta(l,i)-thij(l)
            END DO
C
C   thats the difference between thetai and thetaij
C
            sij=0.0d0
            DO l=1,dp1
	       thijl=thij(l)
	       sij=sij+dmat(l,l)*thijl*thijl
	       IF (l.eq.dp1) CYCLE
               DO k=l+1,dp1
                  sij=sij+2.d0*dmat(k,l)*thijl*thij(k)
               END DO
            END DO
C
C     now we have everything to compute  w_{ij}
C       
            IF (sij.le.spmax) THEN
               wij=wij0*dexp(-sij)
               wghts(j)=wij
            ELSE
               wghts(j)=0.d0
            ENDIF
         ELSE
            wghts(j)=wij0
         END IF
      END DO
      RETURN
      END     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate contribution of Y_j to ai and bi in univariate local polynomial aws (GLM)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine guniaibi(mfamily,dp1,dp2,psix,theta,y,bii,ai,
     1                    bi0,bi2,bi02,wij,wij0)
      implicit logical (a-z)
      integer mfamily,dp1,dp2
      real*8 psix(dp1),theta(dp1),bii(dp2),ai(dp1),bi0(dp2),wij,wij0,
     1       y,ymz1,z,bi2(dp2),bi02(dp2),wij2,wij02
      integer i
      real*8 psith,z1,z2
      psith=0.0d0
      wij2=wij*wij
      wij02=wij0*wij0
      DO i=1,dp1
         psith=psith+psix(i)*theta(i)
      ENDDO
      IF (mfamily.eq.1) THEN
C  Gaussian case
         DO i=1,dp2
            bii(i)=bii(i)+psix(i)*wij
            bi0(i)=bi0(i)+psix(i)*wij0
            bi2(i)=bi2(i)+psix(i)*wij2
            bi02(i)=bi02(i)+psix(i)*wij02
            IF (i.le.dp1) ai(i)=ai(i)+(y-psith)*psix(i)*wij
         ENDDO 
      ELSE IF (mfamily.eq.2) THEN
C  Poisson case
         IF (psith.gt.5.d0) THEN 
             z1=psith-5.d0
C             z1=dexp(5.d0)*(1.d0+z1*(1.d0+.5d0*z1))
             z1=148.4132d0*(1.d0+z1*(1.d0+.5d0*z1*(1.d0+z1/6.d0)))
             z2=z1*dexp(-psith)
C             wij=0.d0
             IF(wij.gt.0.d0) ymz1=z2*y-z1
C             call dblepr("psith",5,psith,1)
         ELSE 
             z1=dexp(psith)
             IF(wij.gt.0.d0) ymz1=y-z1
         ENDIF
         DO i=1,dp2
            z=psix(i)*z1
            bi0(i)=bi0(i)+z*wij0
            bi02(i)=bi02(i)+z*wij02
            if(wij.le.0.d0) CYCLE
            bii(i)=bii(i)+z*wij
            bi2(i)=bi2(i)+z*wij2
            IF (i.le.dp1) ai(i)=ai(i)+ymz1*psix(i)*wij
         ENDDO 
      ELSE IF (mfamily.eq.3) THEN
C  Bernoulli case
         z2=dexp(psith)
         z1=z2/(1+z2)
         z2=z1/(1+z2)
         DO i=1,dp2
            bii(i)=bii(i)+z1*psix(i)*wij
            bi0(i)=bi0(i)+z1*psix(i)*wij0
            IF (i.le.dp1) ai(i)=ai(i)+(y-z2)*psix(i)*wij
         ENDDO 
      ELSE IF (mfamily.eq.4) THEN
C  Exponential case
         z1=-1/psith
         z2=-z1/psith
         DO i=1,dp2
            bii(i)=bii(i)+z1*psix(i)*wij
            bi0(i)=bi0(i)+z1*psix(i)*wij0
            IF (i.le.dp1) ai(i)=ai(i)+(y-z2)*psix(i)*wij
         ENDDO 
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cphawsun(n,dp1,dp2,x,y,fix,si2,theta,bi,si0,bi0,
     1   ai,lam,h,kern,cb,dmat,thij,psix,psiy,bii,spmax)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     bin        \sum \Psi^T Wi^k \Psi        (output)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat,        working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy                working memory 
C     
C      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,l,je,ja,nwij,kern 
      logical aws,fix(n)
      real*8 lkern,skern,x(n),y(n),psix(dp2),psiy(dp1),theta(dp1,n),
     1 bi(dp2,n),bii(dp2),ai(dp1,n),lam,spmax,si2(n),dmat(dp1,dp1),
     2 thij(dp1),cb(dp1,dp1),bi0(dp2,n),lambda,h,si0(n),
     4 sij,wij,z,xij,ha,ha2,eps,xi,wijs,thijl
      external lkern,skern
      eps=1.e-10
      aws=lam.lt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ha=h
         xi=x(i)
C     first search for ja and je
1099     DO j=i,1,-1
            IF (x(j).le.xi-ha) EXIT
            ja=j
         END DO
         DO j=i,n
            IF (x(j).gt.xi+ha) EXIT
            je=j
         END DO
         IF (je-ja-1.le.dp1) THEN
            ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
            goto 1099
         END IF
         ha2=ha*ha
8999     nwij=0
          IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat 
C    expand Bi in dmat
           l=1
            DO j=1,dp1
               DO k=1,dp1
                  dmat(j,k)=bi(j+k-1,i)/lambda
	       END DO
            END DO
         END IF
         DO l=1,dp2
            bii(l)=0.d0
            bi0(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l,i)=0.d0
         END DO
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased 
C      (to reduce panalization by stochastic and influence term)
C
C             loop over local neighborhood
C
         DO j=ja,je
C
C              get location weights
C
            xij=(x(j)-x(i))
            wij=lkern(kern,xij*xij/ha2)*si2(j)
	    wijs=1.d0
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               IF (k.le.dp1) psiy(k)=z*y(j)
               z=z*xij
            END DO
            DO k=1,dp2
               bi0(k,i)=bi0(k,i)+wij*psix(k)
            END DO    
            IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
               DO l=1,dp1
                 thij(l)=theta(l,j)
               END DO
               IF (dp1.gt.1) THEN
                  z=1.d0
                  DO l=2,dp1
                     z=-z*xij
                     DO k=1,dp1-l+1
                        thij(k)=thij(k)+cb(k+l-1,k)*z*theta(k+l-1,j)
                     END DO
                  END DO
               END IF
C     thij contains theta_{j} in the model centered at xi
C
C                 now get sij
C
	       DO l=1,dp1
                  thij(l)=theta(l,i)-thij(l)
               END DO
C
C   thats the difference between thetai and thetaij
C
               sij=0.0d0
               DO l=1,dp1
	          thijl=thij(l)
		  sij=sij+dmat(l,l)*thijl*thijl
		  IF (l.eq.dp1) CYCLE
                  DO k=l+1,dp1
                     sij=sij+2.d0*dmat(k,l)*thijl*thij(k)
                  END DO
               END DO
C
C     now we have everything to compute  w_{ij}
C       
               IF (sij.gt.spmax) CYCLE
               wij=wij*skern(sij,spmax)
            END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),bi0(i),ai(i)  
            DO k=1,dp2
               bii(k)=bii(k)+wij*psix(k)
               IF (k.le.dp1) ai(k,i)=ai(k,i)+wij*psiy(k)
            END DO    
         END DO    
C    this was the j - loop
         IF (nwij.lt.dp1.and.aws) THEN
            lambda=lambda*1.25
C    increase lambda to weaken stochastic penalization
            goto 8999
	 END IF
	 DO k=1,dp2 
	    bi(k,i)=bii(k)
         END DO	 
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in bivariate local polynomial aws (gridded) 
C
C   p > 0   only  !!!! 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cpawsbi(n1,n2,dp1,dp2,y,fix,theta,bi,bi0,ai,
     1    lam,h,kern,dmat,thij,psix,si,si0,siy,ind,wght,spmax)
C
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     theta      estimates from step (k-1)   (input)
C     bi         \sum \Psi^T Wi \Psi  from step (k-1)
C     bin        \sum \Psi^T Wi \Psi  (output)
C     bi0        \sum \Psi^T Wi0 \Psi  from step (k-1)
C     ai        \sum \Psi^T Wi Y     (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat       working arrays  dp1 times dp1
C     thij               vector for parameter differences
C     psix                     working memory for Psij
C     si         working array  for \sum \Psi^T Wi \Psi
C     si0        working array  for \sum \Psi^T Wi0 \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C
      implicit logical (a-z)
      integer n1,n2,dp1,dp2,i1,i2,j1,j2,k,l,info,je1,ja1,je2,ja2,
     1        j,ih1,ih2,ind(dp1,dp1),m,kern
      logical aws,fix(n1,n2)
      real*8 bi(dp2,n1,n2),bi0(dp2,n1,n2),lkern,
     1       ai(dp1,n1,n2),wght,theta(dp1,n1,n2),siy(dp1),
     2       dmat(dp1,dp1),thij(dp1),spmax,y(n1,n2),h,lam,
     3       psix(dp2),si(dp2),si0(dp2),sii,wijy,d,z11,z12,z22,ha2,
     4       lambda,wij,s0i,z1,z2,thijl
      external lkern
C     
C     in case of dp1==1  lawsbi should be called  (p=0)
C
      ha2=h*h
      ih1=h
      aws=lam.lt.1.d20
      DO i1=1,n1
         DO i2=1,n2
	    IF (fix(i1,i2)) CYCLE
C    nothing to do, estimate in (i1,i2) is already fixed by control
            lambda=lam
C  first fill si and siy with 0's
C  fields are used to sum components of ai, bin and bi0
8999        DO j=1,dp2
               si(j)=0.d0
               si0(j)=0.d0
               IF (j.le.dp1) siy(j)=0.d0
            END DO
            IF (aws) THEN
               DO j=1,dp1
                  DO k=1,dp1
                     m=ind(j,k)
                     dmat(j,k)=bi(m,i1,i2)/lambda
                  END DO
               END DO
               s0i=bi0(1,i1,i2)
               sii=bi(1,i1,i2)
            END IF
C
C     Prepare for loop over j's
C
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            DO j1=ja1,je1
               z1=(i1-j1)
               z11=z1*z1
               ih2=dsqrt(ha2-z11)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
                  z2=(i2-j2)*wght
                  z22=z2*z2
                  z12=z1*z2
C          first compute location part
                  wij=lkern(kern,(z11+z22)/ha2)
                  si0(1)=si0(1)+wij
                  si0(2)=si0(2)-z1*wij
                  si0(3)=si0(3)-z2*wij
                  si0(4)=si0(4)+z11*wij
                  si0(5)=si0(5)+z12*wij
                  si0(6)=si0(6)+z22*wij
                  IF (dp1.gt.3) THEN
                     si0(7)=si0(7)-z11*z1*wij
                     si0(8)=si0(8)-z11*z2*wij
                     si0(9)=si0(9)-z1*z22*wij
                     si0(10)=si0(10)-z2*z22*wij
                     si0(11)=si0(11)+z11*z11*wij
                     si0(12)=si0(12)+z11*z12*wij
                     si0(13)=si0(13)+z11*z22*wij
                     si0(14)=si0(14)+z12*z22*wij
                     si0(15)=si0(15)+z22*z22*wij
		  END IF
C          this is the location penalty, now fill si0
                  IF (aws) THEN 
C          now fill psix 
C           now translate thetaj into model centered in xi
C
                     thij(1)=theta(1,j1,j2)
                     thij(2)=theta(2,j1,j2)
                     thij(3)=theta(3,j1,j2)
                     thij(1)=thij(1)+thij(2)*z1+thij(3)*z2
                     IF (dp1.gt.3) THEN
                        thij(4)=theta(4,j1,j2)
                        thij(5)=theta(5,j1,j2)
                        thij(6)=theta(6,j1,j2)
                        thij(1)=thij(1)+thij(4)*z11+
     +                           thij(5)*z12+thij(6)*z22
                        thij(2)=thij(2)+thij(5)*z2+2.d0*thij(4)*z1
                        thij(3)=thij(3)+thij(5)*z1+2.d0*thij(6)*z2
		     END IF
C  
C           get difference of thetas
C
                     DO l=1,dp1
                        thij(l)=theta(l,i1,i2)-thij(l)
                     END DO
C
C           get stochastic penalty
C
                     d=0.d0
                     DO l=1,dp1
		        thijl=thij(l)
		        d=d+dmat(l,l)*thijl*thijl
			IF(l.eq.dp1) CYCLE
                        DO k=l+1,dp1
                           d=d+2.d0*dmat(k,l)*thijl*thij(k)
                        END DO
                     END DO
C
C           now compute weights 
C           stochastic and extension penalty can be added because kernel is exp
C
                     IF (d.gt.spmax) CYCLE
                     wij=wij*dexp(-d)
	          END IF
C           now compute contributions to bi(i),ai(i)  
                  wijy=y(j1,j2)*wij
                  si(1)=si(1)+wij
                  si(2)=si(2)-z1*wij
                  si(3)=si(3)-z2*wij
                  si(4)=si(4)+z11*wij
                  si(5)=si(5)+z12*wij
                  si(6)=si(6)+z22*wij
                  siy(1)=siy(1)+wijy
                  siy(2)=siy(2)-z1*wijy
                  siy(3)=siy(3)-z2*wijy
        	  IF (dp1.le.3) CYCLE
                  si(7)=si(7)-z11*z1*wij
                  si(8)=si(8)-z11*z2*wij
                  si(9)=si(9)-z1*z22*wij
                  si(10)=si(10)-z2*z22*wij
                  si(11)=si(11)+z11*z11*wij
                  si(12)=si(12)+z11*z12*wij
                  si(13)=si(13)+z11*z22*wij
                  si(14)=si(14)+z12*z22*wij
                  si(15)=si(15)+z22*z22*wij
                  siy(4)=siy(4)+z11*wijy
                  siy(5)=siy(5)+z12*wijy
		  siy(6)=siy(6)+z22*wijy
               END DO
            END DO
C        prepare matrix to test for singularity of Bi 
C        this should be changed to SVD at some point
            DO k=1,dp1
               DO j=1,dp1
                  IF (j.gt.k) THEN
                     dmat(j,k)=0.d0
                  ELSE
                     dmat(j,k)=si(ind(j,k))
                  ENDIF
               END DO
            END DO
C
C     compute choleski decomposition
C
            call invers(dmat,dp1,info)
            IF (info.gt.0) THEN
C
C          if singular relax stochastic and extension penalty 
C
               lambda=1.5*lambda
               goto 8999
            END IF
C     
C     now fill ai, bi and bi0
C
            DO j=1,dp1
               ai(j,i1,i2)=siy(j)
            END DO
            DO j=1,dp2
               bi(j,i1,i2)=si(j)
               bi0(j,i1,i2)=si0(j)
            END DO
         END DO     
      END DO     
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in bivariate local polynomial aws (gridded) 
C
C   p > 0   only  !!!! 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cphawsbi(n1,n2,dp1,dp2,y,fix,si2,theta,bi,bi0,ai,
     1    lam,h,kern,dmat,thij,psix,si,si0,siy,ind,wght,spmax)
C
C
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     theta      estimates from step (k-1)   (input)
C     bi         \sum \Psi^T Wi \Psi  from step (k-1)
C     bi0        \sum \Psi^T Wi0 \Psi  from step (k-1)
C     ai        \sum \Psi^T Wi Y     (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat       working arrays  dp1 times dp1
C     thij               vector for parameter differences
C     psix                     working memory for Psij
C     si         working array  for \sum \Psi^T Wi \Psi
C     si0        working array  for \sum \Psi^T Wi0 \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C
      implicit logical (a-z)
      integer n1,n2,dp1,dp2,i1,i2,j1,j2,k,l,info,je1,ja1,je2,ja2,
     1        j,ih1,ih2,ind(dp1,dp1),m,kern
      logical aws,fix(n1,n2)
      real*8 bi(dp2,n1,n2),bi0(dp2,n1,n2),lkern,skern,ai(dp1,n1,n2),
     1       wght,si2(n1,n2),theta(dp1,n1,n2),dmat(dp1,dp1),siy(dp1),
     2       thij(dp1),spmax,y(n1,n2),lam,psix(dp2),si(dp2),si0(dp2),
     3       sii,wijy,d,z11,z12,z22,ha2,h,lambda,wij,s0i,z1,z2,thijl
      external lkern,skern
C     
C     in case of dp1==1  lawsbi should be called  (p=0)
C
      ha2=h*h
      ih1=h
      aws=lam.lt.1.d20
      DO i1=1,n1
         DO i2=1,n2
	    IF (fix(i1,i2)) CYCLE
C    nothing to do, estimate in (i1,i2) is already fixed by control
            lambda=lam
C  first fill si and siy with 0's
C  fields are used to sum components of ai, bi and bi0
8999        DO j=1,dp2
               si(j)=0.d0
               si0(j)=0.d0
               IF (j.le.dp1) siy(j)=0.d0
            END DO
            IF (aws) THEN
               DO j=1,dp1
                  DO k=1,dp1
                     m=ind(j,k)
                     dmat(j,k)=bi(m,i1,i2)
                  END DO
               END DO
               s0i=bi0(1,i1,i2)
               sii=bi(1,i1,i2)
            END IF
C
C     Prepare for loop over j's
C
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            DO j1=ja1,je1
               z1=(i1-j1)
               z11=z1*z1
               ih2=dsqrt(ha2-z11)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
                  z2=(i2-j2)*wght
                  z22=z2*z2
                  z12=z1*z2
C          first compute location part
                  wij=lkern(kern,(z11+z22)/ha2)*si2(j1,j2)
                  si0(1)=si0(1)+wij
                  si0(2)=si0(2)-z1*wij
                  si0(3)=si0(3)-z2*wij
                  si0(4)=si0(4)+z11*wij
                  si0(5)=si0(5)+z12*wij
                  si0(6)=si0(6)+z22*wij
                  IF (dp1.gt.3) THEN
                     si0(7)=si0(7)-z11*z1*wij
                     si0(8)=si0(8)-z11*z2*wij
                     si0(9)=si0(9)-z1*z22*wij
                     si0(10)=si0(10)-z2*z22*wij
                     si0(11)=si0(11)+z11*z11*wij
                     si0(12)=si0(12)+z11*z12*wij
                     si0(13)=si0(13)+z11*z22*wij
                     si0(14)=si0(14)+z12*z22*wij
                     si0(15)=si0(15)+z22*z22*wij
		  END IF
C          this is the location penalty, now fill si0
                  IF (aws) THEN 
C          now fill psix 
C           now translate thetaj into model centered in xi
C
                     thij(1)=theta(1,j1,j2)
                     thij(2)=theta(2,j1,j2)
                     thij(3)=theta(3,j1,j2)
                     thij(1)=thij(1)+thij(2)*z1+thij(3)*z2
                     IF (dp1.gt.3) THEN
                        thij(4)=theta(4,j1,j2)
                        thij(5)=theta(5,j1,j2)
                        thij(6)=theta(6,j1,j2)
                        thij(1)=thij(1)+thij(4)*z11+
     +                           thij(5)*z12+thij(6)*z22
                        thij(2)=thij(2)+thij(5)*z2+2.d0*thij(4)*z1
                        thij(3)=thij(3)+thij(5)*z1+2.d0*thij(6)*z2
		     END IF
C  
C           get difference of thetas
C
                     DO l=1,dp1
                        thij(l)=theta(l,i1,i2)-thij(l)
                     END DO
C
C           get stochastic penalty
C
                     d=0.d0
                     DO l=1,dp1
		        thijl=thij(l)
		        d=d+dmat(l,l)*thijl*thijl
			IF(l.eq.dp1) CYCLE
                        DO k=l+1,dp1
                           d=d+2.d0*dmat(k,l)*thijl*thij(k)
                        END DO
                     END DO
C
C           now compute weights 
C           stochastic and extension penalty can be added because kernel is exp
C
                     wij=wij*skern(d/lambda,spmax)
C           now compute contributions to bi(i),ai(i)  
                     wijy=y(j1,j2)*wij
                     si(1)=si(1)+wij
                     si(2)=si(2)-z1*wij
                     si(3)=si(3)-z2*wij
                     si(4)=si(4)+z11*wij
                     si(5)=si(5)+z12*wij
                     si(6)=si(6)+z22*wij
                     siy(1)=siy(1)+wijy
                     siy(2)=siy(2)-z1*wijy
                     siy(3)=siy(3)-z2*wijy
                     IF (dp1.le.3) CYCLE
                     si(7)=si(7)-z11*z1*wij
                     si(8)=si(8)-z11*z2*wij
                     si(9)=si(9)-z1*z22*wij
                     si(10)=si(10)-z2*z22*wij
                     si(11)=si(11)+z11*z11*wij
                     si(12)=si(12)+z11*z12*wij
                     si(13)=si(13)+z11*z22*wij
                     si(14)=si(14)+z12*z22*wij
                     si(15)=si(15)+z22*z22*wij
                     siy(4)=siy(4)+z11*wijy
                     siy(5)=siy(5)+z12*wijy
		     siy(6)=siy(6)+z22*wijy
	          END IF
               END DO
            END DO
C        prepare matrix to test for singularity of Bi 
C        this should be changed to SVD at some point
            DO k=1,dp1
               DO j=1,dp1
                  IF (j.gt.k) THEN
                     dmat(j,k)=0.d0
                  ELSE
                     dmat(j,k)=si(ind(j,k))
                  ENDIF
               END DO
            END DO
C
C     compute choleski decomposition
C
            call invers(dmat,dp1,info)
            IF (info.gt.0) THEN
C
C          if singular relax stochastic and extension penalty 
C
               lambda=1.5*lambda
               goto 8999
            END IF
C     
C     now fill ai, bin and bi0
C
            DO j=1,dp1
               ai(j,i1,i2)=siy(j)
            END DO
            DO j=1,dp2
               bi(j,i1,i2)=si(j)
               bi0(j,i1,i2)=si0(j)
            END DO
         END DO     
      END DO     
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsun0(dp1,dp2,ai,bi,dmat,d,info)
C    
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer dp1,dp2,info
      real*8 ai(dp1),bi(dp2),dmat(dp1,dp1),d(dp1)
      integer j,k
      DO j=1,dp1
         DO k=1,dp1
            IF (j.gt.k) then 
               dmat(j,k)=0.0d0
            ELSE
               dmat(j,k)=bi(j+k-1)
            END IF
         END DO
      END DO
      call invers(dmat,dp1,info)
C      now dmat contains inverse of B_i 
C      now calculate theta as B_i^{-1} A_i
      DO j=1,dp1
         d(j)=0.0d0
         DO k=1,dp1
            d(j)=d(j)+dmat(j,k)*ai(k)  
         END DO
      END DO
C     just keep the old estimate if info > 0
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsuni(n,dp1,dp2,ai,bi,theta,dmat)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer n,dp1,dp2
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1)
      integer i,j,k,info
      real*8 d
      DO i=1,n
         DO j=1,dp1
            DO k=1,dp1
               IF (j.gt.k) then 
                  dmat(j,k)=0.0d0
               ELSE
                  dmat(j,k)=bi(j+k-1,i)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.ne.0) CYCLE 
C      now dmat contains inverse of B_i 
C      now calculate theta as B_i^{-1} A_i
         DO j=1,dp1
            d=0.0d0
            DO k=1,dp1
               d=d+dmat(j,k)*ai(k,i)  
            END DO
            theta(j,i)=d
         END DO
C     just keep the old estimate if info > 0
      END DO
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine confuni(n,dp1,dp2,theta,bi,bi2,conf,dmat)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer n,dp1,dp2
      real*8 bi(dp2,n),bi2(dp2,n),theta(n),dmat(dp1,dp1),conf(2,n)
      integer i,j,k,info
      real*8 d
      DO i=1,n
         DO j=1,dp1
            DO k=1,dp1
               IF (j.gt.k) then 
                  dmat(j,k)=0.0d0
               ELSE
                  dmat(j,k)=bi(j+k-1,i)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.ne.0) CYCLE 
         d=0.d0
         DO j=1,dp1
            DO k=1,dp1
               d=d+dmat(1,k)*dmat(1,j)*bi2(j+k-1,i)
            END DO
         END DO
         d=dsqrt(2*d)*1.96
         conf(1,i)=theta(i)-d
         conf(2,i)=theta(i)+d
      END DO
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bibi2ibi(n,dp1,dp2,dp3,bi,bi2,erg,dmat,dmat2)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
      implicit logical(a-z)
      integer n,dp1,dp2,dp3
      real*8 bi(dp2,n),bi2(dp2,n),erg(dp3,n),dmat(dp1,dp1),
     1       dmat2(dp1,dp1)
      integer ii,i,j,k,l,m,info
      real*8 d
      DO ii=1,n
         DO j=1,dp1
            DO k=1,dp1
               IF (j.gt.k) then 
                  dmat(j,k)=0.0d0
               ELSE
                  dmat(j,k)=bi2(j+k-1,ii)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.ne.0) CYCLE 
         m=1
         DO i=1,dp1
            DO j=1,i
               d=0.d0
               DO k=1,dp1
                  DO l=1,dp1
                     d=d+dmat(k,l)*bi(i+k-1,ii)*bi(j+l-1,ii)
                  END DO
               END DO
               erg(m,ii)=d
               m=m+1
            END DO
         END DO
      END DO
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (bivariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsbi(n,dp1,dp2,ai,bi,theta,dmat,ind)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working arrays
C
C      implicit logical (a-z)
      integer n,dp1,dp2
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1)
      integer i,j,k,info,ind(dp1,dp1)
      real*8 d
      DO i=1,n
         DO k=1,dp1
            DO j=1,dp1
               IF (j.gt.k) then
                  dmat(j,k)=0.d0
               ELSE
                  dmat(j,k)=bi(ind(j,k),i)
               END IF
            END DO
C            dmat(k,k)=dmat(k,k)*1.001
	 END DO
         call invers(dmat,dp1,info)
         IF (info.gt.0) CYCLE  
C     just keep the old estimate
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} A_i
         DO j=1,dp1
            d=0.0d0
            DO k=1,dp1
               d=d+dmat(j,k)*ai(k,i)  
            END DO
            theta(j,i)=d
         END DO
      END DO
      RETURN
      END      
c
c     invers computes the inverse of a certain
c     double precision symmetric positive definite matrix (see below)
c     
c     this code is based on choleski decomposition and
c     integrates code from linpack (dpofa.f and dpodi.f, version of 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.)
c     and BLAS
C
c     on entry
c
c        a       double precision(n, n)
c
c        n       integer
c                the order of the matrix  a .
c
c
c     on return
c
c        a       invers produces the upper half of inverse(a) .
c
c     subroutines and functions
c
c     fortran dsqrt
c
      subroutine invers(a, n, info)
      integer n,info
      double precision a(n,n)
c
c     internal variables
c
      double precision ddot,t
      double precision s
      integer j,jm1,k,l,i,kp1,im1
c     begin block with ...exits to 940 if singular (info.ne.0)
c
c
      DO j = 1, n
         info = j
         s = 0.0d0
         jm1 = j - 1
         IF (jm1 .ge. 1) THEN
            DO k = 1, jm1
               ddot=0.0d0
               DO l=1,k-1
                  ddot=ddot+a(l,k)*a(l,j)
               END DO           
               t = a(k,j) - ddot
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
            END DO
	 END IF
         s = a(j,j) - s
c     .....exit
         IF (s .le. 1.d-100) RETURN
         a(j,j) = dsqrt(s)
      END DO
      info = 0
c
c     now we have the choeski decomposition in a     
c
c     next code from dpodi to compute inverse
c
      DO k = 1, n
         a(k,k) = 1.0d0/a(k,k)
         t = -a(k,k)
         DO l=1,k-1
            a(l,k)=t*a(l,k)
         END DO          
         kp1 = k + 1
         IF (n .ge. kp1) THEN
            DO j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               DO l=1,k
                  a(l,j)=a(l,j)+t*a(l,k)
               END DO       
            END DO
         END IF
      END DO
c
c        form  inverse(r) * trans(inverse(r))
c
      DO j = 1, n
         jm1 = j - 1
         IF (jm1 .ge. 1) THEN
            DO k = 1, jm1
               t = a(k,j)
               DO l=1,k
                  a(l,k)=a(l,k)+t*a(l,j)
               END DO       
            END DO
	 END IF
         t = a(j,j)
         DO l=1,j
            a(l,j)=t*a(l,j)
         END DO          
      END DO
c
c     now fill lower triangle       
c  
      DO i=1,n
         im1 = i-1
         DO j=1,im1
            a(i,j) = a(j,i)
         END DO
      END DO
      RETURN
      END

